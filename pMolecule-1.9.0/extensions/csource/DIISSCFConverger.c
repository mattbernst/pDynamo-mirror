/*------------------------------------------------------------------------------
! . File      : DIISSCFConverger.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . An scf converger combining damping or the ODA/RCA method with DIIS.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "DIISSCFConverger.h"
# include "Memory.h"

/* . To avoid problems it may be better to have valid densities on entry. */

/* . The DIIS equations very easily become ill-conditioned. Great care should be taken! */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Defaults for some variables. */
# define DEFAULT_DENSITYTOLERANCE 1.0e-8
# define DEFAULT_DIISONSET        0.2e+00
# define DEFAULT_ENERGYTOLERANCE  2.0e-4
# define DEFAULT_MAXIMUMSCFCYCLES 100
# define DEFAULT_NDIISMATRICES    10
# define DEFAULT_QUSERCA          True

/* . Other options. */
# define DAMP_TOLERANCE           1.0e+00
# define ENERGY_TOLERANCE         2.0e-4
# define EXTRAPOLATION_FREQUENCY  15
# define MINIMUM_MU               0.01e+00
# define RCA_ONSET                0.8e+00

/* . Compile options. */
/* # define NOUPDATEDENSITY */
/*# define DIISDEBUG*/
/*# define ODADEBUG*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Davidson_Damping_Factor    ( const Real de, const Real dep, Real deavg, Real *damp ) ;
static void Davidson_Damping_Fock_Save ( SymmetricMatrix *d0, SymmetricMatrix *d1, SymmetricMatrix *d2 ) ;
static void Davidson_Damping_Iterate   ( const Integer iteration, const Real damp, SymmetricMatrix *d0, SymmetricMatrix *d1, SymmetricMatrix *d2 ) ;

static Real DIIS_Fock_Save ( DIISSCFConvergerState *self, QCOnePDM *densityp, QCOnePDM *densityq, SymmetricMatrix *overlap, Real2DArray *orthogonalizer ) ;
static void DIIS_Iterate   ( DIISSCFConvergerState *self, QCOnePDM *densityp, QCOnePDM *densityq ) ;

static void RCA_Iterate      ( DIISSCFConvergerState *self, QCOnePDM *densityp, QCOnePDM *densityq, Real *energy ) ;
static void RCA_Save         ( QCOnePDM *df, SymmetricMatrix *density, SymmetricMatrix *fock ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DIISSCFConverger *DIISSCFConverger_Allocate ( void )
{
    DIISSCFConverger *self = NULL ;
    self = ( DIISSCFConverger * ) Memory_Allocate ( sizeof ( DIISSCFConverger ) ) ;
    if ( self != NULL )
    {
        self->densityTolerance = DEFAULT_DENSITYTOLERANCE ;
        self->diisonset        = DEFAULT_DIISONSET        ;
        self->energytolerance  = DEFAULT_ENERGYTOLERANCE  ;
        self->maximumSCFCycles = DEFAULT_MAXIMUMSCFCYCLES ;
        self->ndiis            = DEFAULT_NDIISMATRICES    ;
        self->QUSERCA          = DEFAULT_QUSERCA          ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
DIISSCFConverger *DIISSCFConverger_Clone ( const DIISSCFConverger *self )
{
    DIISSCFConverger *new = NULL ;
    if ( self != NULL )
    {
        new = DIISSCFConverger_Allocate ( ) ;
        new->densityTolerance = self->densityTolerance ;
        new->diisonset        = self->diisonset        ;
        new->energytolerance  = self->energytolerance  ;
        new->maximumSCFCycles = self->maximumSCFCycles ;
        new->ndiis            = self->ndiis            ;
        new->QUSERCA          = self->QUSERCA          ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for continuation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean DIISSCFConverger_Continue ( const DIISSCFConverger *self, const DIISSCFConvergerState *convergerstate )
{
    Boolean QCONTINUE = False ;
    if ( ( self != NULL ) && ( convergerstate != NULL ) )
    {
        QCONTINUE = ( convergerstate->iteration < self->maximumSCFCycles ) && ( ! convergerstate->isConverged ) ;
    }
# ifdef DIISDEBUG
    printf ( "DIIS Continue = %d\n", QCONTINUE ) ;
# endif
    return QCONTINUE ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for convergence.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean DIISSCFConverger_Converged ( const DIISSCFConverger *self, DIISSCFConvergerState *convergerstate )
{
    Boolean isConverged = False ;
    if ( ( self != NULL ) && ( convergerstate != NULL ) )
    {
        isConverged = ( convergerstate->iteration > 0 ) && ( convergerstate->rmsdifference <= self->densityTolerance ) && ( convergerstate->deltaeold <= self->energytolerance ) ;
        convergerstate->isConverged = isConverged ;
# ifdef DIISDEBUG
printf ( "Iteration, RMS Tolerance, DeltaE: %d %25.15f %25.15f %25.15f %25.15f\n", convergerstate->iteration, convergerstate->rmsdifference, self->densityTolerance, convergerstate->deltaeold, self->energytolerance ) ;
# endif
    }
# ifdef DIISDEBUG
else printf ( "Self or State is NULL\n" ) ;
# endif
    return isConverged ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DIISSCFConverger_Deallocate ( DIISSCFConverger **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Iterate with no density construction.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DIISSCFConverger_IterateWithoutDensities ( const DIISSCFConverger *self, DIISSCFConvergerState *convergerstate, const Real energy )
{
    if ( ( self != NULL ) && ( convergerstate != NULL ) )
    {
        auto Real deltae, etemp, maximum_error = 0.0e+00 ;

        /* . Energies. */
        deltae = energy - convergerstate->eold ;
        etemp  = energy ;

        /* . Check for a valid density. */
        if ( convergerstate->densityp->isValid )
        {
            /* . Determine the average energy and damping factor and update the old deltae. */
            if ( ! self->QUSERCA )
            {
                if ( convergerstate->iteration > 0 )
                {
                   if      ( convergerstate->iteration == 1 ) convergerstate->deavg =   fabs ( deltae ) ;
                   else if ( convergerstate->iteration >  1 ) convergerstate->deavg = ( fabs ( deltae ) + fabs ( convergerstate->deltaeold ) + ( 0.2e+00 * convergerstate->deavg ) ) / 2.2e+00 ;
                   Davidson_Damping_Factor ( deltae, convergerstate->deltaeold, convergerstate->deavg, &(convergerstate->damp) ) ;
                }
            }

            /* . Save the input data for the DIIS procedure. */
            maximum_error = DIIS_Fock_Save ( convergerstate, convergerstate->densityp, convergerstate->densityq, convergerstate->overlap, convergerstate->orthogonalizer ) ;

            /* . Start the convergence procedure. */
            if ( convergerstate->iteration > 0 )
            {

                /* . Check to see which procedure to do. */
                if ( self->QUSERCA )
                {
                    if ( convergerstate->QDIIS && ( maximum_error > RCA_ONSET ) )
                    {
                        convergerstate->QDIIS = False ;
                        convergerstate->QRCA  = True  ;
                    }
                    else if ( ( convergerstate->iteration > 1 ) && convergerstate->QRCA && ( maximum_error < self->diisonset ) )
                    {
                        convergerstate->QDIIS = True  ;
                        convergerstate->QRCA  = False ;
                    }
                }
                else
                {
                    /* . Check to see if extrapolation is to be used. */
                    if ( convergerstate->QDAMP )
                    {
                        convergerstate->QDIIS = ( fabs ( deltae ) < ENERGY_TOLERANCE ) && ( convergerstate->damp < DAMP_TOLERANCE ) ;
                        if ( convergerstate->QDIIS ) convergerstate->QDAMP = False ;
                    }
                    /* . Check to see if damping is to be switched back on. */
                    else
                    {
                        convergerstate->QDAMP = ( fabs ( deltae ) >= ENERGY_TOLERANCE ) ;
                        if ( convergerstate->QDAMP )
                        {
                            convergerstate->damp  = DAMP_TOLERANCE ;
                            convergerstate->QDIIS = False ;
                        }
                    }
                }

                /* . DIIS. */
                if ( convergerstate->QDIIS )
                {
                    DIIS_Iterate ( convergerstate, convergerstate->densityp, convergerstate->densityq ) ;
                }
                /* . RCA. */
                else if ( self->QUSERCA && convergerstate->QRCA )
                {
                    RCA_Iterate ( convergerstate, convergerstate->densityp, convergerstate->densityq, &etemp ) ;
                }
                /* . Damping. */
                else if ( ( ! self->QUSERCA ) && convergerstate->QDAMP )
                {
                    Davidson_Damping_Iterate ( convergerstate->iteration, convergerstate->damp, convergerstate->densityp->fock, convergerstate->dfockp1, convergerstate->dfockp2 ) ;
	            if ( convergerstate->densityq != NULL ) Davidson_Damping_Iterate ( convergerstate->iteration, convergerstate->damp, convergerstate->densityq->fock, convergerstate->dfockq1, convergerstate->dfockq2 ) ;
                }
            }

            /* . Save data for RCA or damping. */
            if ( self->QUSERCA )
            {
                RCA_Save ( convergerstate->densityp, convergerstate->rdensityp, convergerstate->rfockp ) ;
                if ( convergerstate->densityq != NULL ) RCA_Save ( convergerstate->densityq, convergerstate->rdensityq, convergerstate->rfockq ) ;
            }
            else
            {
                Davidson_Damping_Fock_Save ( convergerstate->densityp->fock, convergerstate->dfockp1, convergerstate->dfockp2 ) ;
                if ( convergerstate->densityq != NULL ) Davidson_Damping_Fock_Save ( convergerstate->densityq->fock, convergerstate->dfockq1, convergerstate->dfockq2 ) ;
            }

            /* . Update the iteration count only for valid densities. */
            convergerstate->iteration++ ;
        }

        /* . Save energy data. */
        convergerstate->deltaeold = deltae        ;
        convergerstate->diiserror = maximum_error ;
        convergerstate->eold      = etemp         ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make densities.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DIISSCFConverger_MakeDensities ( const DIISSCFConverger *self, DIISSCFConvergerState *convergerstate )
{
    if ( ( self != NULL ) && ( convergerstate != NULL ) )
    {
        auto Real rmsdifferencea = 0.0e+00, rmsdifferenceb = 0.0e+00 ;

        /* . Always build new densities. */
        /* . Create the densities from the Fock matrices - no level-shifting which is often not good with DIIS. */
        rmsdifferencea = QCOnePDM_MakeFromFock ( convergerstate->densityp, convergerstate->orthogonalizer, NULL ) ;
        rmsdifferenceb = QCOnePDM_MakeFromFock ( convergerstate->densityq, convergerstate->orthogonalizer, NULL ) ;

        /* . Save the rms difference. */
        convergerstate->rmsdifference = Maximum ( rmsdifferencea, rmsdifferenceb ) ;
    }
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Save the output Fock matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Davidson_Damping_Fock_Save ( SymmetricMatrix *d0, SymmetricMatrix *d1, SymmetricMatrix *d2 )
{
    SymmetricMatrix_CopyTo ( d1, d2 ) ;
    SymmetricMatrix_CopyTo ( d0, d1 ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a value for the damping factor.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define MAXIMUM_DAMP 256.0e+00
# define SCALE_FACTOR  16.0e+00
static void Davidson_Damping_Factor ( const Real de, const Real dep, Real deavg, Real *damp )
{
   Real biggest ;

   /* . Initialization. */
   biggest = Maximum ( (*damp), deavg ) ;

   /* . Branch on the five cases to be distinguished. */
   if ( de > 0.0e+00 )
   {
      (*damp) = 4.0e+00 * biggest ;
      /* . de > 0.0 , dep > 0.0. */
      if ( dep > 0.0e+00 )
      {
         if      ( de >= ( 4.0e+00  * dep ) ) (*damp)  = SCALE_FACTOR * biggest ;
         else if ( de <= ( 0.25e+00 * dep ) ) (*damp) /= SCALE_FACTOR ;
         else                                 (*damp)  = pow ( de / dep, 2 ) * biggest ;
      }
      /* . de > 0.0 , dep < 0.0. */
      else
      {
         if (   de         >   0.5e+00 * deavg )   (*damp) *= SCALE_FACTOR ;
         if ( ( de - dep ) < ( 0.2e+00 * deavg ) ) (*damp) /= SCALE_FACTOR ;
      }
   }
   else
   {
      /* . de < 0.0 , dep > 0.0. */
      if ( dep > 0.0e+00 )
      {
         (*damp) = 4.0e+00 * biggest ;
         if (   -de         > deavg ) (*damp) *= SCALE_FACTOR ;
         if ( ( -de + dep ) < deavg ) (*damp) /= SCALE_FACTOR ;
      }
      /* . de < 0.0 , dep < 0.0. */
      else
      {
         /* . de > dep. */
         if ( de > dep )
         {
            if ( de <= ( 0.25e+00 * dep ) ) (*damp)  = pow ( de / dep, 2 ) * biggest ;
            else                            (*damp) /= SCALE_FACTOR ;
         }
         /* . de < dep. */
         else
         {
            if      ( abs ( de ) >= ( 2.0e+00 * deavg ) ) (*damp)  = SCALE_FACTOR * biggest ;
            else if ( abs ( de ) <= ( 0.5e+00 * deavg ) ) (*damp) /= SCALE_FACTOR ;
         }
      }
   }

   /* . Limit how big (*damp) can get. */
   if ( (*damp) > MAXIMUM_DAMP ) (*damp) = MAXIMUM_DAMP ;
}
# undef MAXIMUM_DAMP
# undef SCALE_FACTOR

/*----------------------------------------------------------------------------------------------------------------------------------
! . Perform Davidson damping.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Davidson_Damping_Iterate ( const Integer iteration, const Real damp, SymmetricMatrix *d0, SymmetricMatrix *d1, SymmetricMatrix *d2 )
{
   Integer i, n ;
   Real  dum11, dum21, dum22, eps ;

   /* . Always do damping. */
   for ( i = 0 ; i < d0->size ; i++ ) d0->data[i] = ( d0->data[i] + damp * d1->data[i] ) / ( 1.0e+00 + damp ) ;

   /* . Check for extrapolation. */
   if ( ( EXTRAPOLATION_FREQUENCY > 0 ) && ( iteration > 1 ) )
   {
      /* . Check whether any extrapolation needs to be performed. */
      n = iteration - ( iteration / EXTRAPOLATION_FREQUENCY ) * EXTRAPOLATION_FREQUENCY ;
      /* . Perform extrapolation. */
      if ( n == 0 )
      {
         /* . Perform damping and extrapolation for iteration > 1. */
         for ( i = 0, dum11 = 0.0e+00, dum21 = 0.0e+00, dum22 = 0.0e+00 ; i < d0->size ; i++ )
         {
            dum11 += ( d0->data[i] - d1->data[i] ) * ( d0->data[i] - d1->data[i] ) ;
            dum21 += ( d1->data[i] - d2->data[i] ) * ( d0->data[i] - d1->data[i] ) ;
            dum22 += ( d1->data[i] - d2->data[i] ) * ( d1->data[i] - d2->data[i] ) ;
         }
         eps = ( dum11 - dum21 ) / ( dum21 - dum22 ) ;
         if ( dum21 * dum21 < 0.5e+00 * dum11 * dum22 ) eps = 0.0e+00 ;
         if ( eps > 0.5e+00 ) eps = 0.5e+00 ;
         for ( i = 0 ; i < d0->size ; i++ ) d0->data[i] = ( d0->data[i] - eps * d1->data[i] ) / ( 1.0e+00 - eps ) ;
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Save the Fock matrices and calculate the error vectors.
! . Return the maximum absolute value of an element of the error vectors.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real DIIS_Fock_Save ( DIISSCFConvergerState *self, QCOnePDM *densityp, QCOnePDM *densityq, SymmetricMatrix *overlap, Real2DArray *orthogonalizer )
{
    Real maximumError = 0.0e+00 ;
    if ( ( self != NULL ) && ( densityp != NULL ) )
    {
        if ( self->errorp != NULL )
        {
            auto Integer i, m, n, newest ;
            auto Real    f ;

            /* . Update the indexing variables. */
            /* . The newest index is always last in matind and the oldest first. */
            if ( self->matnum < self->ndiis )
            {
               self->matnum++ ;
               newest = self->matind[self->matnum-1] ;
            }
            else
            {
               self->matnum = self->ndiis ;
               newest = self->matind[0] ;
               for ( i = 0 ; i < self->ndiis-1 ; i++ ) self->matind[i] = self->matind[i+1] ;
               self->matind[self->ndiis-1] = newest ;
            }

            /* . First orbital set. */
            /* . Save all data. */
            SymmetricMatrix_CopyTo ( densityp->density, self->xdensityp[newest] ) ;
            SymmetricMatrix_CopyTo ( densityp->fock   , self->xfockp   [newest] ) ;

            /* . Calculate ( f*p*s - s*p*f ) and transform the error vector to the orthogonal basis. */
            if ( ( overlap == NULL ) || ( orthogonalizer == NULL ) )
            {
                AntisymmetricMatrix_CommutatorSS_Fast ( self->errorp[newest] ,
                                                        densityp->fock       ,
                                                        densityp->density    ,
                                                        self->work[0]        ,
                                                        self->work[1]        ,
                                                        self->work[2]        ,
                                                        NULL                 ) ;
            }
            else
            {
                AntisymmetricMatrix_CommutatorTSSST   ( self->errorp[newest] ,
                                                        densityp->fock       ,
                                                        densityp->density    ,
                                                        overlap              ,
                                                        orthogonalizer       ,
                                                        False                , 
                                                        self->work[0]        ,
                                                        self->work[1]        ,
                                                        self->work[2]        ,
                                                        NULL                 ) ;
/*
                AntisymmetricMatrix_CommutatorSSS ( self->asmwork, densityp->fock, densityp->density, overlap,  NULL ) ;
                AntisymmetricMatrix_Transform     ( self->asmwork, orthogonalizer, False, self->errorp[newest], NULL ) ;
*/
            }

            /* . Find the biggest absolute element of the error vector. */
            maximumError = AntisymmetricMatrix_AbsoluteMaximum ( self->errorp[newest] ) ;

# ifdef DIISDEBUG
printf ( "\nMaximum Error = %25.15f\n", maximumError ) ;
# endif

            /* . Evaluate the scalar products between the new and old error vectors. */
            for ( m = 0 ; m < self->matnum ; m++ )
            {
                n = self->matind[m] ;
                f = - AntisymmetricMatrix_Trace2 ( self->errorp[n], self->errorp[newest], NULL ) ;
# ifdef DIISDEBUG
printf ( "Increment%d  %d %d %25.15f %25.15f\n", m, n, newest, f, Real2DArray_Item ( self->bcoeff, n, newest ) ) ;
# endif
               Real2DArray_Item ( self->bcoeff, n, newest ) = f ;
               if ( n != newest ) Real2DArray_Item ( self->bcoeff, newest, n ) = f ;
            }

            /* . Second orbital set. */
            if ( densityq != NULL )
            {
                /* . Save all data. */
                SymmetricMatrix_CopyTo ( densityq->density, self->xdensityq[newest] ) ;
                SymmetricMatrix_CopyTo ( densityq->fock   , self->xfockq   [newest] ) ;

                /* . Calculate ( f*p*s - s*p*f ) and transform the error vector to the orthogonal basis. */
                if ( ( overlap == NULL ) || ( orthogonalizer == NULL ) )
                {
            	    AntisymmetricMatrix_CommutatorSS_Fast ( self->errorq[newest] ,
                                                            densityq->fock       ,
                                                            densityq->density    ,
                                                            self->work[0]        ,
                                                            self->work[1]        ,
                                                            self->work[2]        ,
                                                            NULL                 ) ;
                }
                else
                {
                    AntisymmetricMatrix_CommutatorTSSST   ( self->errorq[newest] ,
                                                            densityq->fock       ,
                                                            densityq->density    ,
                                                            overlap              ,
                                                            orthogonalizer       ,
                                                            False                , 
                                                            self->work[0]        ,
                                                            self->work[1]        ,
                                                            self->work[2]        ,
                                                            NULL                 ) ;
/*
                    AntisymmetricMatrix_CommutatorSSS ( self->asmwork, densityq->fock, densityq->density, overlap,  NULL ) ;
                    AntisymmetricMatrix_Transform     ( self->asmwork, orthogonalizer, False, self->errorq[newest], NULL ) ;
*/
                }

                /* . Find the biggest absolute element of the error vector. */
                maximumError = Maximum ( maximumError, AntisymmetricMatrix_AbsoluteMaximum ( self->errorq[newest] ) ) ;

# ifdef DIISDEBUG
printf ( "\nMaximum Error = %25.15f\n", maximumError ) ;
# endif

                /* . Evaluate the scalar products between the new and old error vectors. */
                for ( m = 0 ; m < self->matnum ; m++ )
                {
                    n = self->matind[m] ;
                    f = - AntisymmetricMatrix_Trace2 ( self->errorq[n], self->errorq[newest], NULL ) ;
# ifdef DIISDEBUG
printf ( "Increment%d  %d %d %25.15f %25.15f\n", m, n, newest, f, Real2DArray_Item ( self->bcoeff, n, newest ) ) ;
# endif
                    Real2DArray_Item ( self->bcoeff, n, newest ) += f ;
                    if ( n != newest ) Real2DArray_Item ( self->bcoeff, newest, n ) += f ;
                }

                /* . Scale the error by 2 to get rough equivalence with the restricted case result. */
                maximumError *= 2.0e+00 ;
            }
        }
    }
    return maximumError ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply the DIIS procedure.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DEVIATION 1.0e-6
static void DIIS_Iterate ( DIISSCFConvergerState *self, QCOnePDM *densityp, QCOnePDM *densityq )
{
   if ( ( self != NULL ) && ( densityp != NULL ) )
   {
      if ( self->errorp != NULL )
      {
         auto Integer              i, j, m, n   ;
         auto Real           deviation, f ;
         auto Real1DArray     *b            ;
         auto Status           localstatus  ;
         auto SymmetricMatrix *a            ;

         /* . Stop if there are too few matrices. */
         if ( self->matnum <= 1 ) return ;

         /* . Set up and solve the linear equations. */
         /* . Loop until there is success. */
         while ( True )
         {
# ifdef DIISDEBUG
printf ( "\nHistory Indices (%d %d):\n", self->matnum, self->ndiis ) ;
for ( i = 0 ; i < self->ndiis ; i++ ) printf ( "%5d", self->matind[i] ) ;
printf ( "\n" ) ;
printf ( "\nDIIS Traces:\n" ) ;
Real2DArray_Print ( self->bcoeff ) ;
# endif
            /* . Allocate the matrix a and the rhs b. */
            Status_Set ( &localstatus, Status_Continue ) ;
            a = SymmetricMatrix_AllocateN ( self->matnum + 1, &localstatus ) ; SymmetricMatrix_Set ( a, 0.0e+00 ) ;
            b = Real1DArray_Allocate      ( self->matnum + 1, &localstatus ) ; Real1DArray_Set     ( b, 0.0e+00 ) ;

            /* . Fill the matrix a. */
            for ( i = 0 ; i < self->matnum ; i++ )
            {
               m = self->matind[i] ;
               for ( j = 0 ; j <= i ; j++ )
               {
                  n = self->matind[j] ;
                  SymmetricMatrix_Set_Component ( a, i, j, Real2DArray_Item ( self->bcoeff, m, n ) ) ;
               }
            }

            /* . Fill the remaining elements of a and b. */
            for ( i = 0 ; i < self->matnum ; i++ )
            {
               SymmetricMatrix_Set_Component ( a, self->matnum, i, -1.0e+00 ) ;
            }
            Real1DArray_Item ( b, self->matnum ) = -1.0e+00 ;

# ifdef DIISDEBUG
printf ( "\nRight Hand Side (%d):\n", self->matnum ) ;
Real1DArray_Print ( b ) ;
printf ( "\nLeft Hand Side:\n" ) ;
SymmetricMatrix_Print ( a ) ;
# endif
            /* . Solve the matrix equation. */
            SymmetricMatrix_SolveLinearEquations ( a, b, &localstatus ) ;

            /* . Determine the sum. */
            Real1DArray_Item ( b, self->matnum ) = -1.0e+00 ;
            deviation = Real1DArray_Sum ( b ) ;

# ifdef DIISDEBUG
printf ( "\nSolution:\n" ) ;
Real1DArray_Print ( b ) ;
# endif
            /* . Success. */
            if ( Status_OK ( &localstatus ) && ( fabs ( deviation ) < DEVIATION ) ) break ;
            /* . An ill-conditioned solution. */
            else
            {
               /* . Free space. */
               Real1DArray_Deallocate     ( &b ) ;
               SymmetricMatrix_Deallocate ( &a ) ;

               /* . Remove the oldest matrix and reset the indexing array. */
               self->matnum-- ;
               for ( i = 0 ; i < self->matnum ; i++ ) self->matind[i] = self->matind[i+1] ;
               for ( i = self->matnum ; i < self->ndiis ; i++ )
               {
                  self->matind[i] = self->matind[i-1] + 1 ;
                  if ( self->matind[i] >= self->ndiis ) self->matind[i] -= self->ndiis ;
               }
               /* . Leave if there are not enough matrices. */
               if ( self->matnum <= 1 ) return ;
            }
         }

         /* . Calculate the new density and Fock matrices. */
         SymmetricMatrix_Set ( densityp->density, 0.0e+00 ) ;
         SymmetricMatrix_Set ( densityp->fock   , 0.0e+00 ) ;
# ifdef DIISDEBUG
printf ( "\nSummation:\n" ) ;
# endif
         for ( m = 0 ; m < self->matnum ; m++ )
         {
            f = Real1DArray_Item ( b, m ) ;
            n = self->matind[m] ;
# ifdef DIISDEBUG
printf ( "m, n, f = %d %d %25.15f\n", m, n, f ) ;
# endif
# ifndef NOUPDATEDENSITY
            SymmetricMatrix_AddScaledMatrix ( densityp->density, f, self->xdensityp[n] ) ;
# endif
            SymmetricMatrix_AddScaledMatrix ( densityp->fock   , f, self->xfockp   [n] ) ;
         }
         if ( densityq != NULL )
         {
            SymmetricMatrix_Set ( densityq->density, 0.0e+00 ) ;
            SymmetricMatrix_Set ( densityq->fock   , 0.0e+00 ) ;
            for ( m = 0 ; m < self->matnum ; m++ )
            {
               f = Real1DArray_Item ( b, m ) ;
               n = self->matind[m] ;
# ifndef NOUPDATEDENSITY
                SymmetricMatrix_AddScaledMatrix ( densityq->density, f, self->xdensityq[n] ) ;
# endif
                SymmetricMatrix_AddScaledMatrix ( densityq->fock   , f, self->xfockq   [n] ) ;
            }
         }

         /* . Free space. */
         Real1DArray_Deallocate     ( &b ) ;
         SymmetricMatrix_Deallocate ( &a ) ;
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Perform the RCA procedure.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void RCA_Iterate ( DIISSCFConvergerState *self, QCOnePDM *densityp, QCOnePDM *densityq, Real *energy )
{
    auto Real a, b, c, d, fac, fmin, mu, mubmin, mut, mut1, mut2, nu ;

    /* . Check factors for spin-unrestricted version! */

    /* . Calculate some traces. */
    a = SymmetricMatrix_Multiply2_Trace ( densityp->fock, densityp->density ) - SymmetricMatrix_Multiply2_Trace ( densityp->fock, self->rdensityp ) ;
    c = SymmetricMatrix_Multiply2_Trace ( self->rfockp,   densityp->density ) - SymmetricMatrix_Multiply2_Trace ( self->rfockp,   self->rdensityp ) ;
    if ( densityq != NULL )
    {
        a += SymmetricMatrix_Multiply2_Trace ( densityq->fock, densityq->density ) - SymmetricMatrix_Multiply2_Trace ( densityq->fock, self->rdensityq ) ;
        c += SymmetricMatrix_Multiply2_Trace ( self->rfockq,   densityq->density ) - SymmetricMatrix_Multiply2_Trace ( self->rfockq,   self->rdensityq ) ;
    }

    /* . Calculate the remaining polynomial coefficients. */
    d  = self->eold ;
    a += c + 2.0e+00 * ( d - (*energy) ) ;
    b  = (*energy) - a - c - d ;

    /* . Find mu by minimizing the cubic polynomial in the range [0,1]. */
    /* . Find the mu at the boundary with the smallest energy (either 0 or 1). */
    fmin   = (*energy) ;
    mubmin = 1.0e+00   ;
    if ( self->eold < (*energy) ) { fmin = self->eold ; mubmin = 0.0e+00 ; }

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
    mu = Maximum ( mu, MINIMUM_MU ) ;

    /* . Form the final matrices. */
    nu = 1.0e+00 - mu ;
# ifndef NOUPDATEDENSITY
    SymmetricMatrix_Scale           ( densityp->density, mu ) ;
    SymmetricMatrix_AddScaledMatrix ( densityp->density, nu, self->rdensityp ) ;
# endif
    SymmetricMatrix_Scale           ( densityp->fock   , mu ) ;
    SymmetricMatrix_AddScaledMatrix ( densityp->fock   , nu, self->rfockp    ) ;
    if ( densityq != NULL )
    {
# ifndef NOUPDATEDENSITY
        SymmetricMatrix_Scale           ( densityq->density, mu ) ;
        SymmetricMatrix_AddScaledMatrix ( densityq->density, nu, self->rdensityq ) ;
# endif
        SymmetricMatrix_Scale           ( densityq->fock   , mu ) ;
        SymmetricMatrix_AddScaledMatrix ( densityq->fock   , nu, self->rfockq    ) ;
    }

    /* . Save the mu parameter. */
    self->rcamu = mu ;

# ifdef ODADEBUG
printf ( "\nODA> a, b, c, d, mu, eold, enew: %25.15f %25.15f %25.15f %25.15f %25.15f %25.15f %25.15f\n", a, b, c, d, mu, self->eold, (*energy) ) ;
# endif

    /* . Reset energy (assuming the interpolation is a reasonable approximation). */
    (*energy) = ( ( a * mu + b ) * mu + c ) * mu + d ;

# ifdef ODADEBUG
printf ( "\nODA> einterp: %25.15f\n", (*energy) ) ;
# endif
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Save data for the RCA procedure.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void RCA_Save ( QCOnePDM *df, SymmetricMatrix *density, SymmetricMatrix *fock )
{
    SymmetricMatrix_CopyTo ( df->density, density ) ;
    SymmetricMatrix_CopyTo ( df->fock   , fock    ) ;
}
