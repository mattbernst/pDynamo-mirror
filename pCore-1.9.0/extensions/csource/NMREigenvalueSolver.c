/*------------------------------------------------------------------------------
! . File      : NMREigenvalueSolver.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module implements the Nino, Munoz-Caro, Reyes eigenvalue solver for symmetric matrices.
!=================================================================================================================================*/

/* . Use JDEigenvalueSolver as model eventually. */

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "ExecutionEnvironment.h"
# include "Memory.h"
# include "NMREigenvalueSolver.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters and constants.
!---------------------------------------------------------------------------------------------------------------------------------*/
const static Real epsM          = 1.0e-17 ; /* . Machine epsilon for a real in double precision. */
const static Real eigenLimit    = 1.0e-10 ; /* . Number of significant decimal figures in the eigenvalues. */
const static Real degLimit      = 1.0e-07 ; /* . This defines the limit for two eigenvalues to be considered degenerate (equivalent to eigenLimit * 1.0e+3). */
const static Real eVectorsLimit = 1.0e-05 ; /* . Minimum difference between the components of two eigenvectors to be considered equal. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void NMREigenvalueSolver_Eigenvalues   ( const Integer n, const Integer nEigen, const Real infNorm, Real *alpha, Real *beta, Real *eValue ) ;
static void NMREigenvalueSolver_Eigenvectors  ( const Integer n, const Integer nEigen, Real *a, Real *alpha, Real *beta, Real *eValue, Real *eVectors ) ;
static Real NMREigenvalueSolver_Householder   ( const Boolean computeEigenvectors, const Integer n, Real *a, Real *alpha, Real *beta, Real *u ) ;
static void NMREigenvalueSolver_ReorderMatrix ( const Integer n, SymmetricMatrix *self, Real *a ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the solution.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NMREigenvalueSolver_CheckSolution (       SymmetricMatrix *self                 ,
                                         const Integer          numberOfEigenvalues  ,
                                               Real1DArray     *eigenvalues          ,
                                               Real1DArray     *referenceEigenvalues ,
                                               Real2DArray     *eigenvectors         ,
                                               Real            *eigenvalueError      ,
                                               Real            *eigenvectorError     ,
                                               Real            *normalizationError   )
{
    if ( eigenvalueError    != NULL ) (*eigenvalueError   ) = 0.0e+00 ;
    if ( eigenvectorError   != NULL ) (*eigenvectorError  ) = 0.0e+00 ;
    if ( normalizationError != NULL ) (*normalizationError) = 0.0e+00 ;
    if ( ( self != NULL ) && ( eigenvalues != NULL ) && ( numberOfEigenvalues > 0 ) )
    {
        if ( eigenvectors != NULL )
        {
            Integer     i, n ;
            Real        deviation1 = 0.0e+00, deviation2 = 0.0e+00 ;
            Real1DArray column, *temporary = NULL ;

            /* . Initialization. */
            n         = SymmetricMatrix_Dimension ( self ) ;
            temporary = Real1DArray_Allocate ( n, NULL ) ;
            if ( temporary != NULL )
            {
                /* . Loop over eigenvectors. */
                for ( i = 0 ; i < numberOfEigenvalues ; i++ )
                {
                    /* . Get A * X - mu * X. */
                    Real2DArray_ColumnSlice ( eigenvectors, i, &column, NULL ) ;
                    SymmetricMatrix_VectorMultiply ( self, &column, temporary, NULL ) ;
                    Real1DArray_AddScaledArray ( temporary, - Real1DArray_Item ( eigenvalues, i ), &column, NULL ) ;
                    deviation1 = Maximum ( deviation1, Real1DArray_AbsoluteMaximum ( temporary ) ) ;
                    deviation2 = Maximum ( deviation2, Real1DArray_Norm2 ( &column ) - 1.0e+00 ) ;
                }
                /* . Finish up. */
                if ( eigenvectorError   != NULL ) (*eigenvectorError  ) = deviation1 ;
                if ( normalizationError != NULL ) (*normalizationError) = deviation2 ;
                Real1DArray_Deallocate ( &temporary ) ;
            }
        }
        if ( ( eigenvalueError != NULL ) && ( referenceEigenvalues != NULL ) )
        {
            auto Integer i ;
            auto Real    deviation3 = 0.0e+00 ;
            for ( i = 0 ; i < numberOfEigenvalues ; i++ )
            {
                deviation3 = Maximum ( deviation3, fabs ( Real1DArray_Item ( eigenvalues, i ) - Real1DArray_Item ( referenceEigenvalues, i ) ) ) ;
            }
            (*eigenvalueError) = deviation3 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Do a double diagonalization to check the method.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define PRINTTOLERANCE 1.0e-08
void NMREigenvalueSolver_FullCheck ( SymmetricMatrix *self, Real1DArray *eigenvalues, Real2DArray *eigenvectors, Status *status )
{
    Boolean          doPrinting ;
    Integer          n ;
    Real             eigenvalueError, eigenvectorError, eigenvectorError0, normalizationError, normalizationError0 ;
    Real1DArray     *nmrEigenvalues  ;
    Real2DArray     *nmrEigenvectors ;
    SymmetricMatrix *copy ;

    /* . Clone all data. */
    n               = SymmetricMatrix_Dimension ( self ) ;
    copy            = SymmetricMatrix_Clone ( self ) ;
    nmrEigenvalues  = Real1DArray_Clone ( eigenvalues , status ) ;
    nmrEigenvectors = Real2DArray_Clone ( eigenvectors, status ) ;

    /* . NMR solver. */
    NMREigenvalueSolver_Solve ( copy, nmrEigenvalues, nmrEigenvectors, status ) ;

    /* . Reference solver. */
    SymmetricMatrix_Diagonalize ( self, eigenvalues, eigenvectors, status ) ;

    /* . Check the solutions. */
    NMREigenvalueSolver_CheckSolution ( copy, n, eigenvalues   , NULL       , eigenvectors   , NULL            , &eigenvectorError0, &normalizationError0 ) ;
    NMREigenvalueSolver_CheckSolution ( copy, n, nmrEigenvalues, eigenvalues, nmrEigenvectors, &eigenvalueError, &eigenvectorError , &normalizationError  ) ;

    /* . Finish up. */
    Real1DArray_Deallocate     ( &nmrEigenvalues  ) ;
    Real2DArray_Deallocate     ( &nmrEigenvectors ) ;
    SymmetricMatrix_Deallocate ( &copy            ) ;

    /* . Print out. */
    doPrinting = ( ( eigenvalueError >= PRINTTOLERANCE ) || ( eigenvectorError >= PRINTTOLERANCE ) || ( normalizationError >= PRINTTOLERANCE ) ) ;
    if ( doPrinting )
    {
        printf ( "\nNMR Diagonalization Errors (with references):\n" ) ;
        if ( eigenvalueError    >= PRINTTOLERANCE ) printf ( "Eigenvalue    error = %20.10g\n", eigenvalueError ) ;
        if ( eigenvectorError   >= PRINTTOLERANCE ) printf ( "Eigenvector   error = %20.10g (%20.10g)\n", eigenvectorError  , eigenvectorError0   ) ;
        if ( normalizationError >= PRINTTOLERANCE ) printf ( "Normalization error = %20.10g (%20.10g)\n", normalizationError, normalizationError0 ) ;
    }
}
# undef PRINTTOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Solve for the eigenvalues and eigenvectors.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NMREigenvalueSolver_Solve ( SymmetricMatrix *self, Real1DArray *eigenvalues, Real2DArray *eigenvectors, Status *status )
{
    if ( ( self != NULL ) && ( eigenvalues != NULL ) )
    {
        /* . Various checks - self and eigenvalues need to be compact as does dimension 1 of eigenvectors. */
             if ( ! ( SymmetricMatrix_IsCompact ( self ) && Real1DArray_IsCompact ( eigenvalues ) ) ) Status_Set ( status, Status_NonCompactArray     ) ;
        else if ( ( eigenvectors != NULL ) && ( ! Real2DArray_IsCompact ( eigenvectors, 1, NULL ) ) ) Status_Set ( status, Status_NonCompactDimension ) ;
        else
        {
            Boolean computeEigenvectors ;
            Integer n, size ;
            Real   *a, *alpha, *beta, fNorm, infNorm, *u ;

            /* . Initialization. */
            computeEigenvectors = ( eigenvectors != NULL ) ;
            n                   = SymmetricMatrix_Dimension ( self ) ;
            size                = SymmetricMatrix_Size      ( self ) ;

            /* . Allocate space. */
            a     = Memory_Allocate_Array_Real ( size ) ; /* . Reordered matrix. */
            alpha = Memory_Allocate_Array_Real ( n    ) ; /* . Diagonal elements of the tridiagonal matrix. */
            beta  = Memory_Allocate_Array_Real ( n    ) ; /* . Off-diagonal elements of the tridiagonal matrix and auxiliary workspace. */
            u     = Memory_Allocate_Array_Real ( n    ) ; /* . Auxiliary workspace. */
            if ( ( a != NULL ) && ( alpha != NULL ) && ( beta != NULL ) && ( u != NULL ) )
            {
                auto Integer i ;
                auto Real    scale ;

                /* . Compute the Frobenius (Euclidean) norm of the matrix. */
                fNorm = sqrt ( SymmetricMatrix_Trace2 ( self, self, NULL ) ) ;

                /* . Reorder the matrix. */
                NMREigenvalueSolver_ReorderMatrix ( n, self, a ) ;

                /* . Scale the matrix. */
                scale = 1.0e+00 / fNorm ;
                for ( i = 0 ; i < size ; i++ ) a[i] *= scale ;

                /* . Determine the tridiagonal matrix. */
                infNorm = NMREigenvalueSolver_Householder ( computeEigenvectors, n, a, alpha, beta, u );

                /* . Compute the eigenvalues from the tridiagonal matrix. */
                NMREigenvalueSolver_Eigenvalues ( n, n, infNorm, alpha, beta, Real1DArray_Data ( eigenvalues ) ) ;

                /* . Computing eigenvectors of the original matrix if required. */
                if ( computeEigenvectors )
                {
                    NMREigenvalueSolver_Eigenvectors (  n, n, a, alpha, beta, Real1DArray_Data ( eigenvalues ), Real2DArray_Data ( eigenvectors ) ) ;
                    Real2DArray_Transpose ( eigenvectors, status ) ;
                }

                /* . Scale back the eigenvalues. */
                Real1DArray_Scale ( eigenvalues, fNorm ) ;
            }
            else Status_Set ( status, Status_MemoryAllocationFailure ) ;

            /* . Finish up. */
            Memory_Deallocate_Real ( &alpha ) ;
            Memory_Deallocate_Real ( &beta  ) ;
            Memory_Deallocate_Real ( &u     ) ;
        }
    }
}

/* . Check on errors? */
/* if ( ifail != 0 ) { Status_Set ( status, Status_DiagonalizationFailure ) ; printf ( "\nDiagonalization Error = %d\n", ifail ) ; } */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Compute the eigenvalues.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void NMREigenvalueSolver_Eigenvalues ( const Integer n, const Integer nEigen, const Real infNorm, Real *alpha, Real *beta, Real *eValue )
{
    Integer i, k, s;
    Real    aux, two = 2.0E+00, x, xi, y, yi, z, zero = 0.0E+00 ;
# ifdef USEOPENMP
    Integer numberOfThreads ;
    Real    q ;
# endif

    /*
    ! . Starting computation of selected eigenvalues using the Sturm sequence.
    ! . Select the whole interval within which all the eigenvalues lie using the Gershgorin theorem.
    ! . At the end, all the eigenvalues lie in the x-y interval.
    */
    if (fabs(beta[0])<epsM) beta[0]=beta[0]>0 ? infNorm*epsM: -infNorm*epsM;
    xi=alpha[0]-fabs(beta[0]);
    yi=alpha[0]+fabs(beta[0]);

    for (i=1; i<n-1; ++i)
    {
        if (fabs(beta[i])<epsM) beta[i]=beta[i]>0 ? infNorm*epsM: -infNorm*epsM;
        aux=alpha[i]-fabs(beta[i])-fabs(beta[i-1]);
        if (aux<xi) xi=aux;
        aux=alpha[i]+fabs(beta[i])+fabs(beta[i-1]);
        if (aux>yi) yi=aux;
    }
    aux=alpha[n-1]-fabs(beta[n-2]);
    if (aux<xi) xi=aux;
    aux=alpha[n-1]+fabs(beta[n-2]);
    if (aux>yi) yi=aux;

    /* . Use the bisection method and a Sturm sequence to locate an interval that brackets the desired eigenvalue. */
    y = yi;
# ifdef USEOPENMP
    numberOfThreads = omp_get_num_threads ( ) ;
    q               = ( ( Real ) nEigen ) / ( Real ) numberOfThreads ;
    #pragma omp parallel private(x, z, s, aux, k, i) firstprivate(y)
# endif
    {
        /* #pragma omp parallel for private(x, z, s, aux, k, i) firstprivate(y) schedule(dynamic) */
# ifdef USEOPENMP
        auto Integer start, end, thread ;
        thread = omp_get_thread_num ( ) ;
        start  = nEigen - 1 - ( Integer ) (   thread       * q ) ;
        end    = nEigen     - ( Integer ) ( ( thread + 1 ) * q ) ;
        for ( k = start ; k >= end ; --k )
        {
# else
        for ( k = ( nEigen - 1 ); k >= 0 ; --k )
        {
# endif
            x = xi;
            while ((y-x)>eigenLimit)
            {
                z=(x+y)/two;
                s=0;
                aux=alpha[0]-z;
                if (aux<zero) s++;
                for (i=1; i<n; ++i)
                {
                    aux=alpha[i]-z-beta[i-1]*beta[i-1]/aux;
                    if (aux<zero) s++;
                }
                if (s<=k) { x=z; }
                else      { y=z; }
            }
            aux=(x+y)/two;
            eValue[k]=aux;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Compute the eigenvectors.
! . Use inverse iteration and rotation of the eigenvectors to the original matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void NMREigenvalueSolver_Eigenvectors ( const Integer n, const Integer nEigen, Real *a, Real *alpha, Real *beta, Real *eValue, Real *eVectors )
{
    Integer i, j, k, index, index0, size, nuroll, ipre;
    Real   aux, aux2, aux3, lambda = 0.0e+00, one = 1.0E+00, two = 2.0E+00, zero = 0.0E+00 ;

    unsigned long g, m, random ; /* . Variables for the MINSTD random number generator. */

    /* . Initialization. */
    nuroll = 5 ;
    size   = ( n*(n+1) )/2 ;

    /*
    ! . Defining parameters for the Park-Miller MINSTD random number generator
    ! . S.K. Park and K.W. Miller (1988). "Random Number Generators: Good Ones
    ! . Are Hard To Find". Communications of the ACM, 31 (10): 1192-1201.
    ! . Values are 7^5, 2^31-1 and the generator seed.
    */
    g      =      16807 ;
    m      = 2147483647 ;
    random =          1 ;

    /* . Compute the eigenvectors of the tridiagonal matrix using inverse iteration. */
# ifdef USEOPENMP
    #pragma omp parallel private (i, j, aux, aux3, k, lambda, ipre, aux2, index, index0) firstprivate(random, zero, one, two)
# endif
    {
        auto Real *alf, *c, *gamma, *vAux ;
        alf   = Memory_Allocate_Array_Real ( n ) ;
        c     = Memory_Allocate_Array_Real ( n ) ;
        gamma = Memory_Allocate_Array_Real ( n - 1 ) ;
        vAux  = Memory_Allocate_Array_Real ( n ) ;

# ifdef USEOPENMP
        #pragma omp for schedule(guided)
# endif
        /* . #pragma omp for private (i, j, aux, aux3, k, lambda, ipre, aux2, index, index0) schedule(dynamic) . */
        for (i=0; i<nEigen; ++i)
        {
            auto Integer iRow ;
            iRow = i * n ;
            /* . Random initialization of the eigenvector using the MINSTD linear congruential generator for random numbers in the (-1, 1) interval. */
            aux3=zero ;
            for (j=0; j<n; ++j)
            {
                vAux[j]=zero;
                random=(random*g)%m;
                aux=one-two*((Real)random)/(m-1);
                eVectors[iRow+j]=aux;
                aux3+=aux*aux;
            }
            aux3=one/sqrt(aux3);
            for (j=0; j<n; ++j) { eVectors[iRow+j]*=aux3; }

            /* . Use LU decomposition to solve the linear system of equations of form Ax=b where A is tridiagonal matrix (a-lambda*I). */
            if (i>0 && (eValue[i]-eValue[i-1])<degLimit) { lambda += two * degLimit ; }
            else                                         { lambda  = eValue[i]      ; }
            alf[0]=alpha[0]-lambda;
            gamma[0]=beta[0]/alf[0];

            ipre=(n-2)%nuroll+1;
            for (j=1; j<ipre; j++)
            {
                alf[j]=alpha[j]-lambda-beta[j-1]*gamma[j-1];
                gamma[j]=beta[j]/alf[j];
            }
            for (j=ipre; j<n-1; j+=nuroll)
            {
                alf[j]=alpha[j]-lambda-beta[j-1]*gamma[j-1];
                gamma[j]=beta[j]/alf[j];
                alf[j+1]=alpha[j+1]-lambda-beta[j]*gamma[j];
                gamma[j+1]=beta[j+1]/alf[j+1];
                alf[j+2]=alpha[j+2]-lambda-beta[j+1]*gamma[j+1];
                gamma[j+2]=beta[j+2]/alf[j+2];
                alf[j+3]=alpha[j+3]-lambda-beta[j+2]*gamma[j+2];
                gamma[j+3]=beta[j+3]/alf[j+3];
                alf[j+4]=alpha[j+4]-lambda-beta[j+3]*gamma[j+3];
                gamma[j+4]=beta[j+4]/alf[j+4];
            }
            alf[n-1]=alpha[n-1]-lambda-beta[n-2]*gamma[n-2];

            /* . Start the inverse iterations using the previous LU decomposition. */
            k=0;
            do
            {
                c[0]=eVectors[iRow]/alf[0];

                ipre=(n-1)%nuroll+1;
                for (j=1; j<ipre; j++)
                {
                    c[j]=(eVectors[iRow+j]-beta[j-1]*c[j-1])/alf[j];
                }
                for (j=ipre; j<n; j+=nuroll)
                {
                    c[j]  =(eVectors[iRow+j]-beta[j-1]*c[j-1])/alf[j];
                    c[j+1]=(eVectors[iRow+j+1]-beta[j]*c[j])/alf[j+1];
                    c[j+2]=(eVectors[iRow+j+2]-beta[j+1]*c[j+1])/alf[j+2];
                    c[j+3]=(eVectors[iRow+j+3]-beta[j+2]*c[j+2])/alf[j+3];
                    c[j+4]=(eVectors[iRow+j+4]-beta[j+3]*c[j+3])/alf[j+4];
                }

                eVectors[iRow+n-1]=c[n-1];
                aux2=c[n-1]*c[n-1];

                ipre=n-2-(n-1)%nuroll;
                for (j=n-2; j>ipre; --j)
                {
                    aux=c[j]-gamma[j]*eVectors[iRow+j+1];
                    eVectors[iRow+j]=aux;
                    aux2+=aux*aux;
                }
                for (j=ipre; j>=0; j-=nuroll)
                {
                    aux=c[j]-gamma[j]*eVectors[iRow+j+1];
                    eVectors[iRow+j]=aux;
                    aux2+=aux*aux;
                    aux=c[j-1]-gamma[j-1]*eVectors[iRow+j];
                    eVectors[iRow+j-1]=aux;
                    aux2+=aux*aux;
                    aux=c[j-2]-gamma[j-2]*eVectors[iRow+j-1];
                    eVectors[iRow+j-2]=aux;
                    aux2+=aux*aux;
                    aux=c[j-3]-gamma[j-3]*eVectors[iRow+j-2];
                    eVectors[iRow+j-3]=aux;
                    aux2+=aux*aux;
                    aux=c[j-4]-gamma[j-4]*eVectors[iRow+j-3];
                    eVectors[iRow+j-4]=aux;
                    aux2+=aux*aux;
                }
                aux2=one/sqrt(aux2);
                aux3=zero;

                ipre=n%nuroll;
                for (j=0; j<ipre; j++)
                {
                    eVectors[iRow+j]*=aux2;
                    aux3+=fabs(eVectors[iRow+j])-fabs(vAux[j]);
                    vAux[j]=eVectors[iRow+j];
                }
                for (j=ipre; j<n; j+=nuroll)
                {
                    eVectors[iRow+j]*=aux2;
                    aux3+=fabs(eVectors[iRow+j])-fabs(vAux[j]);
                    vAux[j]=eVectors[iRow+j];
                    eVectors[iRow+j+1]*=aux2;
                    aux3+=fabs(eVectors[iRow+j+1])-fabs(vAux[j+1]);
                    vAux[j+1]=eVectors[iRow+j+1];
                    eVectors[iRow+j+2]*=aux2;
                    aux3+=fabs(eVectors[iRow+j+2])-fabs(vAux[j+2]);
                    vAux[j+2]=eVectors[iRow+j+2];
                    eVectors[iRow+j+3]*=aux2;
                    aux3+=fabs(eVectors[iRow+j+3])-fabs(vAux[j+3]);
                    vAux[j+3]=eVectors[iRow+j+3];
                    eVectors[iRow+j+4]*=aux2;
                    aux3+=fabs(eVectors[iRow+j+4])-fabs(vAux[j+4]);
                    vAux[j+4]=eVectors[iRow+j+4];
                }
                ++k;
            } while(fabs(aux3)>eVectorsLimit && k<10); /* . Maximum iterations? */


            /* . Rotate the eigenvector to the original matrix. */
            index0=size-5; /* . Index of first u matrix element in this row. */
            for (j=n-3; j>=0; --j)
            {
                index=index0;
                aux=zero;
                ipre=(n-j-1)%nuroll+j+1;
                for (k=j+1; k<ipre; ++k)
                {
                    aux+=a[index]*eVectors[iRow+k];
                    ++index;
                }
                for (k=ipre; k<n; k+=nuroll)
                {
                    aux+=a[index]*eVectors[iRow+k];
                    aux+=a[++index]*eVectors[iRow+k+1];
                    aux+=a[++index]*eVectors[iRow+k+2];
                    aux+=a[++index]*eVectors[iRow+k+3];
                    aux+=a[++index]*eVectors[iRow+k+4];
                    ++index;
                }
                aux*=a[index0-1];
                index=index0;
                for (k=j+1; k<ipre; ++k)
                {
                    eVectors[iRow+k]-=aux*a[index];
                    ++index;
                }
                for (k=ipre; k<n; k+=nuroll)
                {
                    eVectors[iRow+k]  -=aux*a[index];
                    eVectors[iRow+k+1]-=aux*a[++index];
                    eVectors[iRow+k+2]-=aux*a[++index];
                    eVectors[iRow+k+3]-=aux*a[++index];
                    eVectors[iRow+k+4]-=aux*a[++index];
                    ++index;
                }
                index0-=n-j+1;
            }
        }

        /* . Finish up. */
        Memory_Deallocate_Real ( &alf   ) ;
        Memory_Deallocate_Real ( &c     ) ;
        Memory_Deallocate_Real ( &gamma ) ;
        Memory_Deallocate_Real ( &vAux  ) ;
    }

    /* . Orthonormalize degenerate eigenvectors using a numerically stable modified Gram-Schmidt. */
    for (i=0; i<nEigen; ++i)
    {
        if (i>0 && (eValue[i]-eValue[i-1])<degLimit)
        {
            auto Integer iRow, jRow ;
            iRow = i * n ;
            j=i-1;
            while (j>=0 && (eValue[i]-eValue[j])<degLimit)
            {
                jRow = j * n ;
                aux=zero;
# ifdef USEOPENMP
                #pragma omp parallel for reduction(+:aux) schedule(dynamic)
# endif
                for (k=0; k<n; ++k) { aux+=eVectors[iRow+k]*eVectors[jRow+k]; }

# ifdef USEOPENMP
                #pragma omp parallel for schedule(dynamic)
# endif
                for (k=0; k<n; ++k) { eVectors[iRow+k]+=-aux*eVectors[jRow+k]; }
                --j;
            }

            /* . Normalize. */
            aux=zero;
# ifdef USEOPENMP
            #pragma omp parallel for reduction(+:aux) schedule(dynamic)
# endif
            for (j=0; j<n; ++j) { aux+=eVectors[iRow+j]*eVectors[iRow+j]; }
            aux=one/sqrt(aux);

# ifdef USEOPENMP
            #pragma omp parallel for schedule(dynamic)
# endif
            for (j=0; j<n; ++j) { eVectors[iRow+j]*=aux; }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Tridiagonalize a real symmetrical matrix using Householder reflections.
! . a is the matrix on input and overwritten on output.
! . alpha is the diagonal of the tridiagonal matrix.
! . beta is the secondary diagonal of the tridiagonal matrix.
! . u is workspace.
! . infNorm is the infinite norm of the tridiagonal matrix (the greatest sum of row elements).
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real NMREigenvalueSolver_Householder ( const Boolean computeEigenvectors, const Integer n, Real *a, Real *alpha, Real *beta, Real *u )
{
    Integer i, index, index0, ipre, j, k, nuroll, piv, size ;
    Real    aux, aux2, aux3, aux4, aux5, infNorm, one = 1.0e+00, two = 2.0e+00, zero = 0.0e+00 ;

    /* . Initialization. */
    infNorm = -1.0e+00    ;
    nuroll  = 5           ; /* . Unroll loops 5 times. */
    size    = (n*(n+1))/2 ;

    /* . Loop over the columns of the matrix. */
    for (i=0; i<n-2; ++i)
    {
        piv     = i+1;
        index0  = i*n-i*(i-1)/2;
        aux     = zero;
        beta[i] = zero;
        u[i]    = zero;

# ifdef USEOPENMP
        #pragma omp parallel private (ipre, j, k, aux4, aux5, index) firstprivate(nuroll, i, piv, size, index0)
# endif
        {
            auto Real *p;

            p = Memory_Allocate_Array_Real ( n ) ;
            index=i*n-i*(i+1)/2;

# ifdef USEOPENMP
            #pragma omp for reduction(+:aux) schedule (static)
# endif
            /* . Initial vector u = [0,a(i,1),...,a(i,n-1)] and its norm aux. */
            for (j=piv; j<n; ++j)
            {
                beta[j]=zero;
                aux5=a[index+j];
                u[j]=aux5;
                aux+=aux5*aux5;
            }

# ifdef USEOPENMP
            #pragma omp single
# endif
            {
                aux3=u[piv]>=zero ? sqrt(aux) : -sqrt(aux);
                u[piv]+=aux3;
                aux=aux+aux3*a[index0+1];
                aux2=zero;
                aux3=one/aux;
            }

            for (j=i; j<n; ++j) { p[j]=zero; }

# ifdef USEOPENMP
            #pragma omp for schedule (guided)
# endif
            for (j=i; j<n; ++j)
            {
                index=j*(n+n-j+1)/2;
                p[j]+=a[index]*u[j];
                ++index;

                ipre=(n-j-1)%nuroll+j+1;
                for (k=j+1; k<ipre; ++k)
                {
                    p[j]+=a[index]*u[k];
                    p[k]+=a[index]*u[j];
                    ++index;
                }
                for (k=ipre; k<n; k+=nuroll)
                {
                    p[j]+=a[index]*u[k];
                    p[k]+=a[index]*u[j];
                    ++index;
                    p[j]+=a[index]*u[k+1];
                    p[k+1]+=a[index]*u[j];
                    ++index;
                    p[j]+=a[index]*u[k+2];
                    p[k+2]+=a[index]*u[j];
                    ++index;
                    p[j]+=a[index]*u[k+3];
                    p[k+3]+=a[index]*u[j];
                    ++index;
                    p[j]+=a[index]*u[k+4];
                    p[k+4]+=a[index]*u[j];
                    ++index;
                }
            }

# ifdef USEOPENMP
            #pragma omp critical
# endif
            {
                for (j=i; j<n; ++j) beta[j]+=p[j]*aux3;
            }

# ifdef USEOPENMP
            #pragma omp barrier
# endif

# ifdef USEOPENMP
            #pragma omp single
# endif
            {
                for (j=i; j<n; ++j) aux2+=u[j]*beta[j];
            }

# ifdef USEOPENMP
            #pragma omp single
# endif
            {
                aux2*=(aux3/two);
            }

# ifdef USEOPENMP
            #pragma omp for schedule(guided)
# endif
            for (j=n-1; j>=i; --j) { beta[j]-=aux2*u[j]; }

# ifdef USEOPENMP
            #pragma omp for schedule(guided)
# endif
            for (j=n-1; j>=i; --j)
            {
                index=j*n-j*(j+1)/2+n-1;
                aux5=u[j];
                aux4=beta[j];

                ipre=n-1-(n-j)%nuroll;
                for (k=n-1; k>ipre; --k)
                {
                    a[index]+=-aux5*beta[k]-aux4*u[k];
                    --index;
                }
                for (k=ipre; k>=j; k-=nuroll)
                {
                    a[index]+=-aux5*beta[k]-aux4*u[k];
                    --index;
                    a[index]+=-aux5*beta[k-1]-aux4*u[k-1];
                    --index;
                    a[index]+=-aux5*beta[k-2]-aux4*u[k-2];
                    --index;
                    a[index]+=-aux5*beta[k-3]-aux4*u[k-3];
                    --index;
                    a[index]+=-aux5*beta[k-4]-aux4*u[k-4];
                    --index;
                }
            }

# ifdef USEOPENMP
            #pragma omp single
# endif
            {
                alpha[i]=a[index0];
                beta[i]=a[index0+1];
                aux=fabs(alpha[i])+fabs(beta[i]);
                if (i>0) aux+=fabs(beta[i-1]);
                if(aux>infNorm) infNorm=aux;
                if (computeEigenvectors) a[index0]=aux3;
            }

            if ( computeEigenvectors )
            {
                index=size-(n-i)*(n-i+1)/2-i;
# ifdef USEOPENMP
                #pragma omp for schedule(guided)
# endif
                for (j=i+1; j<n; ++j) { a[index+j]=u[j]; }
            }

            /* . Clear up. */
            Memory_Deallocate_Real ( &p ) ;
        }
    }

    /* . Finish up. */
    index=(n-2)*(n+3)/2 ;   /* . Before last diagonal element. */
    alpha[n-2]=a[index] ;   /* . Last diagonal elements. */
    beta [n-2]=a[++index] ; /* . Last off-diagonal element. */
    alpha[n-1]=a[++index] ;

    aux=fabs(alpha[n-2])+fabs(beta[n-2])+fabs(beta[n-3]);
    if(aux>infNorm) infNorm=aux;
    aux=fabs(alpha[n-1])+fabs(beta[n-2]);
    if(aux>infNorm) infNorm=aux;

    /* . Finish up. */
    return infNorm ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Reorder the input matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void NMREigenvalueSolver_ReorderMatrix ( const Integer n, SymmetricMatrix *self, Real *a )
{
    auto Integer i, index, j ;
    index = 0 ;
    for ( i = 0 ; i < n ; i++ )
    {
        for ( j = i ; j < n ; j++ ) { a[index++] = SymmetricMatrix_Item ( self, j, i ) ; }
    }
}
