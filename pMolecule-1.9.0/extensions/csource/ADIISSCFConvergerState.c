/*------------------------------------------------------------------------------
! . File      : ADIISSCFConvergerState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . The state for ADIIS SCF convergers.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "ADIISSCFConvergerState.h"
# include "Memory.h"
# include "Status.h"

/*# define DIISDEBUG*/

/*
! . A flexible approach is needed for EDIIS and ODA with UHF/UKS.
! . At the moment the expansion coefficients are constrained to
! . be identical for the alpha and beta density matrices.
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void ADIISSCFConvergerOnePDM_Finalize   ( ADIISSCFConvergerOnePDM *self ) ;
static void ADIISSCFConvergerOnePDM_Initialize ( ADIISSCFConvergerOnePDM *self, const Integer sizeSM, const Integer sizeASM, const Boolean doError, Status *status ) ;

static void ADIISSCFConvergerState_ADIISSetUpDerivatives ( ADIISSCFConvergerState *self, const Integer newest, const Integer increment, const Real2DArray *adiisTraces, const Boolean useEDIIS ) ;
static Real ADIISSCFConvergerState_DIISSaveOne           ( ADIISSCFConvergerState *self, const Integer newest, const QCOnePDM *density, ADIISSCFConvergerOnePDM *densities, Real2DArray *adiisTraces ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . ADIIS setup.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConvergerState_ADIISSetUp ( ADIISSCFConvergerState *self, const Boolean useEDIIS )
{
    if ( ( self != NULL ) && ( self->history > 1 ) )
    {
        auto Integer dimension, newest ;

        /* . Allocate space. */
        dimension = self->history ;
        if ( ( ! useEDIIS ) && ( self->densitiesQ != NULL ) ) dimension *= 2 ;
        self->adiisAlphas    = Real1DArray_Allocate      ( dimension, NULL ) ; Real1DArray_Set     ( self->adiisAlphas   , 0.0e+00 ) ;
        self->adiisGradients = Real1DArray_Allocate      ( dimension, NULL ) ; Real1DArray_Set     ( self->adiisGradients, 0.0e+00 ) ;
        self->adiisHessian   = SymmetricMatrix_AllocateN ( dimension, NULL ) ; SymmetricMatrix_Set ( self->adiisHessian  , 0.0e+00 ) ;

        /* . Evaluate the gradient and Hessian (with scaling). */
        newest = Integer1DArray_Item ( self->historyIndices, 0 ) ;
        ADIISSCFConvergerState_ADIISSetUpDerivatives ( self, newest, 0            , self->adiisTracesP, useEDIIS ) ;
        ADIISSCFConvergerState_ADIISSetUpDerivatives ( self, newest, self->history, self->adiisTracesQ, useEDIIS ) ;
        SymmetricMatrix_Scale ( self->adiisHessian, 0.5e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . ADIIS density and Fock matrix update.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConvergerState_ADIISUpdate ( ADIISSCFConvergerState *self )
{
    if ( ( self != NULL ) && ( self->history > 1 ) )
    {
        auto Integer   m, n ;
        auto QCOnePDM *densityP, *densityQ ;
        auto Real      f ;

        /* . Aliases. */
        densityP = self->densityP ;
        densityQ = self->densityQ ;

# ifdef DIISDEBUG
printf ( "\nADIIS Updating (%d):\n", self->history ) ;
# endif

        /* . Calculate the new density and Fock matrices. */
        SymmetricMatrix_Set ( densityP->density, 0.0e+00 ) ;
        SymmetricMatrix_Set ( densityP->fock   , 0.0e+00 ) ;
        for ( m = 0 ; m < self->history ; m++ )
        {
            f = Real1DArray_Item ( self->adiisAlphas, m ) ;
            n = Integer1DArray_Item ( self->historyIndices, m ) ;
# ifdef DIISDEBUG
printf ( "M, N, F, E: %5d %5d %25.15f %25.15f\n", m, n, f, Real1DArray_Item ( self->energies, n ) ) ;
# endif
            SymmetricMatrix_AddScaledMatrix ( densityP->density, f, self->densitiesP[n].density ) ;
            SymmetricMatrix_AddScaledMatrix ( densityP->fock   , f, self->densitiesP[n].fock    ) ;
        }
        if ( densityQ != NULL )
        {
            auto Integer increment ;
            if ( self->adiisAlphas->length > self->history ) increment = self->history ;
            else                                             increment = 0 ;
            SymmetricMatrix_Set ( densityQ->density, 0.0e+00 ) ;
            SymmetricMatrix_Set ( densityQ->fock   , 0.0e+00 ) ;
            for ( m = 0 ; m < self->history ; m++ )
            {
                f = Real1DArray_Item ( self->adiisAlphas, m+increment ) ;
                n = Integer1DArray_Item ( self->historyIndices, m ) ;
                SymmetricMatrix_AddScaledMatrix ( densityQ->density, f, self->densitiesQ[n].density ) ;
                SymmetricMatrix_AddScaledMatrix ( densityQ->fock   , f, self->densitiesQ[n].fock    ) ;
            }
        }

        /* . Free space. */
        Real1DArray_Deallocate     ( &(self->adiisAlphas   ) ) ;
        Real1DArray_Deallocate     ( &(self->adiisGradients) ) ;
        SymmetricMatrix_Deallocate ( &(self->adiisHessian  ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ADIISSCFConvergerState *ADIISSCFConvergerState_Allocate ( void )
{
    ADIISSCFConvergerState *self = NULL ;
    self = ( ADIISSCFConvergerState * ) Memory_Allocate ( sizeof ( ADIISSCFConvergerState ) ) ;
    if ( self != NULL )
    {
        /* . Scalars. */
        self->iterationType  = ADIISSCFConvergerIterationType_Null ;
        self->isConverged    = False ;
        self->diisActive     = 0 ;
        self->history        = 0 ;
        self->iteration      = 0 ;
        self->maximumHistory = 0 ;
        self->diisError      = 0.0e+00 ;
        self->energyChange   = 0.0e+00 ;
        self->odaOldEnergy   = 0.0e+00 ;
        self->oldEnergy      = 0.0e+00 ;
        self->odaMu          = 1.0e+00 ; /* . Default first value. */
        self->rmsDifference  = 0.0e+00 ;
        /* . Aliases. */
        self->densityP       = NULL ;
        self->densityQ       = NULL ;
        self->orthogonalizer = NULL ;
        self->overlap        = NULL ;
        /* . Allocated arrays. */
        self->densitiesP     = NULL ;
        self->densitiesQ     = NULL ;
        self->odaDensityP    = NULL ;
        self->odaDensityQ    = NULL ;
        self->asmWork        = NULL ;
        self->historyIndices = NULL ;
        self->adiisAlphas    = NULL ;
        self->adiisGradients = NULL ;
        self->energies       = NULL ;
        self->adiisTracesP   = NULL ;
        self->adiisTracesQ   = NULL ;
        self->diisTraces     = NULL ;
        self->adiisHessian   = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for convergence.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ADIISSCFConvergerState_Converged ( const ADIISSCFConvergerState *self )
{
    Boolean isConverged = False ;
    if ( self != NULL ) isConverged = self->isConverged ;
    return isConverged ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConvergerState_Deallocate ( ADIISSCFConvergerState **self )
{
    if ( (*self) != NULL )
    {
        auto Integer i ;
        /* . Basic arrays. */
        AntisymmetricMatrix_Deallocate ( &((*self)->asmWork       ) ) ;
        Integer1DArray_Deallocate      ( &((*self)->historyIndices) ) ;
        Real1DArray_Deallocate         ( &((*self)->adiisAlphas   ) ) ;
        Real1DArray_Deallocate         ( &((*self)->adiisGradients) ) ;
        Real1DArray_Deallocate         ( &((*self)->energies      ) ) ;
        Real2DArray_Deallocate         ( &((*self)->adiisTracesP  ) ) ;
        Real2DArray_Deallocate         ( &((*self)->adiisTracesQ  ) ) ;
        Real2DArray_Deallocate         ( &((*self)->diisTraces    ) ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->adiisHessian  ) ) ;
        /* . ADIIS history. */
        if ( (*self)->densitiesP != NULL )
        {
            for ( i = 0 ; i < (*self)->maximumHistory ; i++ ) ADIISSCFConvergerOnePDM_Finalize ( &((*self)->densitiesP[i]) ) ;
            Memory_Deallocate ( (*self)->densitiesP ) ;
        }
        if ( (*self)->densitiesQ != NULL )
        {
            for ( i = 0 ; i < (*self)->maximumHistory ; i++ ) ADIISSCFConvergerOnePDM_Finalize ( &((*self)->densitiesQ[i]) ) ;
            Memory_Deallocate ( (*self)->densitiesQ ) ;
        }
        /* . ODA history. */
        if ( (*self)->odaDensityP != NULL )
        {
            ADIISSCFConvergerOnePDM_Finalize ( (*self)->odaDensityP ) ;
            Memory_Deallocate ( (*self)->odaDensityP ) ;
        }
        if ( (*self)->odaDensityQ != NULL )
        {
            ADIISSCFConvergerOnePDM_Finalize ( (*self)->odaDensityQ ) ;
            Memory_Deallocate ( (*self)->odaDensityQ ) ;
        }
        /* . finish up. */
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply the DIIS procedure.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DEVIATION 1.0e-6
void ADIISSCFConvergerState_DIISIterate ( ADIISSCFConvergerState *self )
{
    if ( ( self != NULL ) && ( self->history > 1 ) )
    {
        auto Integer          i, j, m, n          ;
        auto QCOnePDM        *densityP, *densityQ ;
        auto Real             deviation, f        ;
        auto Real1DArray     *b                   ;
        auto Status           localStatus         ;
        auto SymmetricMatrix *a                   ;

        /* . Aliases. */
        densityP = self->densityP ;
        densityQ = self->densityQ ;

# ifdef DIISDEBUG
printf ( "\nHistory Indices (%d %d):\n", self->history, self->maximumHistory ) ;
Integer1DArray_Print ( self->historyIndices ) ;
printf ( "\nDIIS Traces:\n" ) ;
Real2DArray_Print ( self->diisTraces ) ;
# endif

        /* . Set up and solve the linear equations. */
        /* . Loop until there is success. */
        self->diisActive = self->history ;
        while ( True )
        {
            /* . Allocate the matrix a and the rhs b. */
            Status_Set ( &localStatus, Status_Continue ) ;
            a = SymmetricMatrix_AllocateN ( self->diisActive + 1, &localStatus ) ; SymmetricMatrix_Set ( a, 0.0e+00 ) ;
            b = Real1DArray_Allocate      ( self->diisActive + 1, &localStatus ) ; Real1DArray_Set     ( b, 0.0e+00 ) ;

            /* . Fill the matrix a. */
            for ( i = 0 ; i < self->diisActive ; i++ )
            {
                m = Integer1DArray_Item ( self->historyIndices, i ) ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    n = Integer1DArray_Item ( self->historyIndices, j ) ;
                    SymmetricMatrix_Set_Component ( a, i, j, Real2DArray_Item ( self->diisTraces, m, n ) ) ;
                }
            }

            /* . Fill the remaining elements of a and b. */
            for ( i = 0 ; i < self->diisActive ; i++ )
            {
                SymmetricMatrix_Set_Component ( a, self->diisActive, i, -1.0e+00 ) ;
            }
            Real1DArray_Item ( b, self->diisActive ) = -1.0e+00 ;

# ifdef DIISDEBUG
printf ( "\nRight Hand Side (%d):\n", self->diisActive ) ;
Real1DArray_Print ( b ) ;
printf ( "\nLeft Hand Side:\n" ) ;
SymmetricMatrix_Print ( a ) ;
# endif

            /* . Solve the matrix equation. */
            SymmetricMatrix_SolveLinearEquations ( a, b, &localStatus ) ;

            /* . Determine the sum. */
            Real1DArray_Item ( b, self->diisActive ) = -1.0e+00 ;
            deviation = Real1DArray_Sum ( b ) ;

# ifdef DIISDEBUG
printf ( "\nSolution:\n" ) ;
Real1DArray_Print ( b ) ;
# endif

            /* . Success. */
            if ( Status_OK ( &localStatus ) && ( fabs ( deviation ) < DEVIATION ) ) break ;
            /* . An ill-conditioned solution. */
            else
            {
               /* . Free space. */
                Real1DArray_Deallocate     ( &b ) ;
                SymmetricMatrix_Deallocate ( &a ) ;

                /* . Remove the oldest matrix. */
                self->diisActive-- ;

                /* . Leave if there are not enough matrices. */
                if ( self->diisActive <= 1 ) return ;
            }
        }

        /* . Calculate the new density and Fock matrices. */
        SymmetricMatrix_Set ( densityP->density, 0.0e+00 ) ;
        SymmetricMatrix_Set ( densityP->fock   , 0.0e+00 ) ;
# ifdef DIISDEBUG
printf ( "\nSummation:\n" ) ;
# endif
        for ( m = 0 ; m < self->diisActive ; m++ )
        {
            f = Real1DArray_Item ( b, m ) ;
            n = Integer1DArray_Item ( self->historyIndices, m ) ;
# ifdef DIISDEBUG
printf ( "m, n, f = %d %d %25.15f\n", m, n, f ) ;
# endif
            SymmetricMatrix_AddScaledMatrix ( densityP->density, f, self->densitiesP[n].density ) ;
            SymmetricMatrix_AddScaledMatrix ( densityP->fock   , f, self->densitiesP[n].fock    ) ;
        }
        if ( densityQ != NULL )
        {
            SymmetricMatrix_Set ( densityQ->density, 0.0e+00 ) ;
            SymmetricMatrix_Set ( densityQ->fock   , 0.0e+00 ) ;
            for ( m = 0 ; m < self->diisActive ; m++ )
            {
                f = Real1DArray_Item ( b, m ) ;
                n = Integer1DArray_Item ( self->historyIndices, m ) ;
                SymmetricMatrix_AddScaledMatrix ( densityQ->density, f, self->densitiesQ[n].density ) ;
                SymmetricMatrix_AddScaledMatrix ( densityQ->fock   , f, self->densitiesQ[n].fock    ) ;
            }
        }

        /* . Free space. */
        Real1DArray_Deallocate     ( &b ) ;
        SymmetricMatrix_Deallocate ( &a ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Save the data necessary for the ADIIS and SCF procedures.
! . Return the maximum absolute value of an element of the error vectors.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real ADIISSCFConvergerState_DIISSave ( ADIISSCFConvergerState *self, const Real energy )
{
    Real maximumError = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer i, n, newest  ;
        auto Real    maximumErrorQ ;

        /* . Update the indexing variables. */
        /* . The newest index is always first in historyIndices and the oldest last. */
        if ( self->history < self->maximumHistory ) self->history++ ;
        else                                        self->history = self->maximumHistory ;
        Integer1DArray_RightCircularShift ( self->historyIndices ) ;
        newest = Integer1DArray_Item ( self->historyIndices, 0 ) ;

        /* . Energies. */
        Real1DArray_Item ( self->energies, newest ) = energy ;

        /* . Densities and traces. */
        /* . For the unrestricted case scale the error by 2 to get rough equivalence with the restricted case result. */
        /* . Zero out DIIS elements. */
        for ( i = 0 ; i < self->history ; i++ )
        {
           n = Integer1DArray_Item ( self->historyIndices, i ) ;
           Real2DArray_Item ( self->diisTraces, n, newest ) = 0.0e+00 ;
           if ( n != newest ) Real2DArray_Item ( self->diisTraces, newest, n ) = 0.0e+00 ;
        }
        maximumError = ADIISSCFConvergerState_DIISSaveOne ( self, newest, self->densityP, self->densitiesP, self->adiisTracesP ) ;
        if ( self->densityQ != NULL )
        {
            /* . Do everything separately here as possible to have multiple calls due to macro. */
            maximumErrorQ = ADIISSCFConvergerState_DIISSaveOne ( self, newest, self->densityQ, self->densitiesQ, self->adiisTracesQ ) ;
            maximumError  = Maximum ( maximumError, maximumErrorQ ) ;
            maximumError *= 2.0e+00 ;
        }
    }
    return maximumError ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Perform the ODA procedure.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConvergerState_ODAIterate ( ADIISSCFConvergerState *self, const Real energy, const Real minimumMu )
{
    if ( self != NULL )
    {
        auto QCOnePDM *densityP, *densityQ ;
        auto Real      a, b, c, d, fac, fmin, mu, mubmin, mut, mut1, mut2, nu ;

        /* . Aliases. */
        densityP = self->densityP ;
        densityQ = self->densityQ ;

        /* . Calculate some traces. */
        a = SymmetricMatrix_Trace2 ( densityP->fock         , densityP->density, NULL ) - SymmetricMatrix_Trace2 ( densityP->fock         , self->odaDensityP->density, NULL ) ;
        c = SymmetricMatrix_Trace2 ( self->odaDensityP->fock, densityP->density, NULL ) - SymmetricMatrix_Trace2 ( self->odaDensityP->fock, self->odaDensityP->density, NULL ) ;
        if ( densityQ != NULL )
        {
            a += SymmetricMatrix_Trace2 ( densityQ->fock         , densityQ->density, NULL ) - SymmetricMatrix_Trace2 ( densityQ->fock         , self->odaDensityQ->density, NULL ) ;
            c += SymmetricMatrix_Trace2 ( self->odaDensityQ->fock, densityQ->density, NULL ) - SymmetricMatrix_Trace2 ( self->odaDensityQ->fock, self->odaDensityQ->density, NULL ) ;
        }

        /* . Calculate the remaining polynomial coefficients. */
        d  = self->odaOldEnergy ;
        a += c + 2.0e+00 * ( d - energy ) ;
        b  = energy - a - c - d ;

        /* . Find mu by minimizing the cubic polynomial in the range [0,1]. */
        /* . Find the mu at the boundary with the smallest energy (either 0 or 1). */
        fmin   = energy  ;
        mubmin = 1.0e+00 ;
        if ( self->odaOldEnergy < energy ) { fmin = self->odaOldEnergy ; mubmin = 0.0e+00 ; }

        /* . Solve the polynomial. */
        if ( a == 0.0e+00 )
        {
            /* . Line. */
            if ( b == 0.0e+00 ) mu = mubmin ;
            /* . Quadratic. */
            else
            {
                mut = - c / ( 2.0e+00 * b ) ;
                if ( ( b > 0.0e+00 ) && ( mut > 0.0e+00 ) && ( mut < 1.0e+00 ) ) mu = mut ;
                else mu = mubmin ;
            }
        }
        /* . Cubic. */
        else
        {
            fac = b * b - 3.0e+00 * a * c ;
            /* . Real turning points. */
            if ( fac >= 0.0e+00 )
            {
                mut1 = ( - b + sqrt ( fac ) ) / ( 3.0e+00 * a ) ;
                mut2 = ( - b - sqrt ( fac ) ) / ( 3.0e+00 * a ) ;
                if ( ( mut1 > 0.0e+00 ) && ( mut1 < 1.0e+00 ) )
                {
                    fac = ( ( a * mut1 + b ) * mut1 + c ) * mut1 + d ;
                    if ( fac < fmin ) { fmin = fac ; mubmin = mut1 ; }
                }
                if ( ( mut2 > 0.0e+00 ) && ( mut2 < 1.0e+00 ) )
                {
                    fac = ( ( a * mut2 + b ) * mut2 + c ) * mut2 + d ;
                    if ( fac < fmin ) { fmin = fac ; mubmin = mut2 ; }
                }
                mu = mubmin ;
            }
            /* . No real turning points. */
            else mu = mubmin ;
        }

        /* . Check mu. */
        mu = Maximum ( mu, minimumMu ) ;

        /* . Form the final matrices. */
        nu = 1.0e+00 - mu ;
        SymmetricMatrix_Scale           ( densityP->density, mu ) ;
        SymmetricMatrix_AddScaledMatrix ( densityP->density, nu, self->odaDensityP->density ) ;
        SymmetricMatrix_Scale           ( densityP->fock   , mu ) ;
        SymmetricMatrix_AddScaledMatrix ( densityP->fock   , nu, self->odaDensityP->fock    ) ;
        if ( densityQ != NULL )
        {
            SymmetricMatrix_Scale           ( densityQ->density, mu ) ;
            SymmetricMatrix_AddScaledMatrix ( densityQ->density, nu, self->odaDensityQ->density ) ;
            SymmetricMatrix_Scale           ( densityQ->fock   , mu ) ;
            SymmetricMatrix_AddScaledMatrix ( densityQ->fock   , nu, self->odaDensityQ->fock    ) ;
        }

        /* . Save the mu parameter. */
        self->odaMu = mu ;

        /* . Store the interpolated energy (assuming the interpolation is a reasonable approximation). */
        self->odaOldEnergy = ( ( a * mu + b ) * mu + c ) * mu + d ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Save data for the ODA procedure.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConvergerState_ODASave ( ADIISSCFConvergerState *self )
{
    if ( self != NULL )
    {
        if ( ( self->densityP != NULL ) && ( self->odaDensityP != NULL ) )
        {
            SymmetricMatrix_CopyTo ( self->densityP->density, self->odaDensityP->density ) ;
            SymmetricMatrix_CopyTo ( self->densityP->fock   , self->odaDensityP->fock    ) ;
        }
        if ( ( self->densityQ != NULL ) && ( self->odaDensityQ != NULL ) )
        {
            SymmetricMatrix_CopyTo ( self->densityQ->density, self->odaDensityQ->density ) ;
            SymmetricMatrix_CopyTo ( self->densityQ->fock   , self->odaDensityQ->fock    ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set up a state.
!---------------------------------------------------------------------------------------------------------------------------------*/
ADIISSCFConvergerState *ADIISSCFConvergerState_SetUp ( const Boolean useODA, const Integer maximumHistory, QCOnePDM *densityP, QCOnePDM *densityQ,
                                                                       SymmetricMatrix *overlap, Real2DArray *orthogonalizer, Status *status )
{
    ADIISSCFConvergerState *self = NULL ;
    if ( densityP != NULL )
    {
        auto Integer i, m, n ;
        auto Status  localStatus ;

        /* . Initialization. */
        Status_Set ( &localStatus, Status_Continue ) ;

        /* . Allocate the structure. */
        self = ADIISSCFConvergerState_Allocate ( ) ;

        /* . Options and aliases. */
        self->maximumHistory = maximumHistory ;
        self->densityP       = densityP       ;
        self->densityQ       = densityQ       ;
        self->overlap        = overlap        ;
        self->orthogonalizer = orthogonalizer ;

        /* . Dimensions. */
        m = densityP->density->dimension ;
        /* . Orthogonalizer specific items. */
        if ( ( overlap != NULL ) && ( orthogonalizer != NULL ) )
        {
            n = orthogonalizer->length1 ;
            self->asmWork = AntisymmetricMatrix_Allocate ( densityP->density->dimension, &localStatus ) ;
        }
        else n = m ;

        /* . The DIIS history indices. */
        self->historyIndices = Integer1DArray_Allocate ( self->maximumHistory, &localStatus ) ;
        for ( i = 0 ; i < self->maximumHistory ; i++ ) Integer1DArray_Item ( self->historyIndices, i ) = i ;
        Integer1DArray_LeftCircularShift ( self->historyIndices ) ;

        /* . The DIIS traces. */
        self->diisTraces = Real2DArray_Allocate ( self->maximumHistory, self->maximumHistory, &localStatus ) ;
        Real2DArray_Set ( self->diisTraces, 0.0e+00 ) ;

        /* . The histories. */
        self->energies   = Real1DArray_Allocate ( self->maximumHistory, &localStatus ) ;
        self->densitiesP = ( ADIISSCFConvergerOnePDM * ) Memory_Allocate_Array ( self->maximumHistory, sizeof ( ADIISSCFConvergerOnePDM ) ) ;
        for ( i = 0 ; i < self->maximumHistory ; i++ ) ADIISSCFConvergerOnePDM_Initialize ( &(self->densitiesP[i]), m, n, True, &localStatus ) ;
        if ( densityQ != NULL )
        {
            self->densitiesQ = ( ADIISSCFConvergerOnePDM * ) Memory_Allocate_Array ( self->maximumHistory, sizeof ( ADIISSCFConvergerOnePDM ) ) ;
            for ( i = 0 ; i < self->maximumHistory ; i++ ) ADIISSCFConvergerOnePDM_Initialize ( &(self->densitiesQ[i]), m, n, True, &localStatus ) ;
        }

        /* . ODA. */
        if ( useODA )
        {
            self->odaDensityP = ( ADIISSCFConvergerOnePDM * ) Memory_Allocate ( sizeof ( ADIISSCFConvergerOnePDM ) ) ;
            ADIISSCFConvergerOnePDM_Initialize ( self->odaDensityP, m, n, False, &localStatus ) ;
            if ( densityQ != NULL )
            {
                self->odaDensityQ = ( ADIISSCFConvergerOnePDM * ) Memory_Allocate ( sizeof ( ADIISSCFConvergerOnePDM ) ) ;
                ADIISSCFConvergerOnePDM_Initialize ( self->odaDensityQ, m, n, False, &localStatus ) ;
            }
        }
        /* . ADIIS. */
        else
        {
            self->adiisTracesP = Real2DArray_Allocate ( self->maximumHistory, self->maximumHistory, &localStatus ) ;
            Real2DArray_Set ( self->adiisTracesP, 0.0e+00 ) ;
            if ( densityQ != NULL )
            {
                self->adiisTracesQ = Real2DArray_Allocate ( self->maximumHistory, self->maximumHistory, &localStatus ) ;
                Real2DArray_Set ( self->adiisTracesQ, 0.0e+00 ) ;
            }
        }

        /* . Finish up. */
        Status_Set ( status, localStatus ) ;
        if ( ! Status_OK ( &localStatus ) ) ADIISSCFConvergerState_Deallocate ( &self ) ;
    }
    return self ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Finalize.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConvergerOnePDM_Finalize ( ADIISSCFConvergerOnePDM *self )
{
    if ( self != NULL )
    {
        AntisymmetricMatrix_Deallocate ( &(self->error  ) ) ;
        SymmetricMatrix_Deallocate     ( &(self->density) ) ;
        SymmetricMatrix_Deallocate     ( &(self->fock   ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConvergerOnePDM_Initialize ( ADIISSCFConvergerOnePDM *self, const Integer sizeSM, const Integer sizeASM, const Boolean doError, Status *status )
{
    if ( self != NULL )
    {
        self->error   = NULL ;
        self->density = NULL ;
        self->fock    = NULL ;
        if ( sizeSM > 0 )
        {
            self->density = SymmetricMatrix_AllocateN ( sizeSM, status ) ;
            self->fock    = SymmetricMatrix_AllocateN ( sizeSM, status ) ;
        }
        if ( ( sizeASM > 0 ) && doError ) self->error = AntisymmetricMatrix_Allocate ( sizeASM, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . ADIIS setup for derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void ADIISSCFConvergerState_ADIISSetUpDerivatives ( ADIISSCFConvergerState *self, const Integer newest, const Integer increment, const Real2DArray *adiisTraces, const Boolean useEDIIS )
{
    if ( ( self != NULL ) && ( adiisTraces != NULL ) )
    {
        auto Integer i, j, m, n ;
        auto Real    dnfn ;

        /* . EDIIS. */
        if ( useEDIIS )
        {
            /* . Initial guess. */
            Real1DArray_Item ( self->adiisAlphas, 0 ) = 1.0e+00 ;

            /* . Diagonal Hessian terms zero. */
            for ( m = 0 ; m < self->history ; m++ )
            {
                i = Integer1DArray_Item ( self->historyIndices, m ) ;
                Real1DArray_Item ( self->adiisGradients, m ) = Real1DArray_Item ( self->energies, i ) ;
                for ( n = 0 ; n < m ; n++ )
                {
                    j = Integer1DArray_Item ( self->historyIndices, n ) ;
                    SymmetricMatrix_Item ( self->adiisHessian, m, n ) += ( Real2DArray_Item ( adiisTraces, i, j ) + Real2DArray_Item ( adiisTraces, j, i ) -
                                                                           Real2DArray_Item ( adiisTraces, i, i ) - Real2DArray_Item ( adiisTraces, j, j ) ) ;
                }
            }
        }
        /* . ADIIS. */
        else
        {
            /* . Initial guess. */
            Real1DArray_Item ( self->adiisAlphas, increment ) = 1.0e+00 ;

            /* . All terms involving newest are zero. */
            dnfn = Real2DArray_Item ( adiisTraces, newest, newest ) ;
            for ( m = 1 ; m < self->history ; m++ )
            {
                i = Integer1DArray_Item ( self->historyIndices, m ) ;
                Real1DArray_Item ( self->adiisGradients, m+increment ) = ( Real2DArray_Item ( adiisTraces, i, newest ) - dnfn ) ;
                for ( n = 1 ; n <= m ; n++ )
                {
                    j = Integer1DArray_Item ( self->historyIndices, n ) ;
                    SymmetricMatrix_Item ( self->adiisHessian, m+increment, n+increment ) = ( Real2DArray_Item ( adiisTraces, i, j      ) + Real2DArray_Item ( adiisTraces, j     , i ) -
                                                                                              Real2DArray_Item ( adiisTraces, i, newest ) - Real2DArray_Item ( adiisTraces, newest, i ) -
                                                                                              Real2DArray_Item ( adiisTraces, j, newest ) - Real2DArray_Item ( adiisTraces, newest, j ) +
                                                                                              2.0e+00 * dnfn ) ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Save the data necessary for the ADIIS and SCF procedures.
! . Return the maximum absolute value of elements of the error vectors.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real ADIISSCFConvergerState_DIISSaveOne ( ADIISSCFConvergerState *self, const Integer newest, const QCOnePDM *density, ADIISSCFConvergerOnePDM *densities, Real2DArray *adiisTraces )
{
    Real maximumError = 0.0e+00 ;
    if ( ( self != NULL ) && ( density != NULL ) && ( densities != NULL ) && ( newest >= 0 ) )
    {
        auto Integer m, n ;
        auto Real    f ;

        /* . Save the data. */
        SymmetricMatrix_CopyTo ( density->density, densities[newest].density ) ;
        SymmetricMatrix_CopyTo ( density->fock   , densities[newest].fock    ) ;

        /* . Calculate ( f*p*s - s*p*f ) and transform the error vector to the orthogonal basis. */
        if ( ( self->overlap == NULL ) || ( self->orthogonalizer == NULL ) )
        {
            AntisymmetricMatrix_CommutatorSS_Reference ( densities[newest].error, density->fock, density->density, NULL ) ;
        }
        else
        {
            AntisymmetricMatrix_CommutatorSSS ( self->asmWork, density->fock, density->density, self->overlap, NULL ) ;
            AntisymmetricMatrix_Transform     ( self->asmWork, self->orthogonalizer, False, densities[newest].error, NULL ) ;
        }

        /* . Find the biggest absolute element of the error vector. */
        maximumError = AntisymmetricMatrix_AbsoluteMaximum ( densities[newest].error ) ;

# ifdef DIISDEBUG
printf ( "\nMaximum Error = %25.15f\n", maximumError ) ;
# endif

        /* . Evaluate the traces required for ADIIS and DIIS. */
        for ( m = 0 ; m < self->history ; m++ )
        {
            n = Integer1DArray_Item ( self->historyIndices, m ) ;
            /* . ADIIS. */
            if ( adiisTraces != NULL )
            {
                f = SymmetricMatrix_Trace2 ( densities[n].density, densities[newest].fock, NULL ) ;
                Real2DArray_Item ( adiisTraces, n, newest ) = f ;
                if ( n != newest )
                {
                    f = SymmetricMatrix_Trace2 ( densities[n].fock, densities[newest].density, NULL ) ;
                    Real2DArray_Item ( adiisTraces, newest, n ) = f ;
                }
            }
            /* . DIIS. */
            f = - AntisymmetricMatrix_Trace2 ( densities[n].error, densities[newest].error, NULL ) ;

# ifdef DIISDEBUG
printf ( "Increment%d  %d %d %25.15f %25.15f\n", m, n, newest, f, Real2DArray_Item ( self->diisTraces, n, newest ) ) ;
# endif

            Real2DArray_Item ( self->diisTraces, n, newest ) += f ;
            if ( n != newest ) Real2DArray_Item ( self->diisTraces, newest, n ) += f ;
        }
    }
    return maximumError ;
}
