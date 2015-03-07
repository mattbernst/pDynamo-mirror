/*------------------------------------------------------------------------------
! . File      : ADIISSCFConverger.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . An SCF converger combining ADIIS, ODA and DIIS.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "ADIISSCFConverger.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Defaults for some variables. */
# define DEFAULT_DENSITYTOLERANCE   1.0e-8
# define DEFAULT_DIISOFF            0.8e+00
# define DEFAULT_DIISON             0.2e+00
# define DEFAULT_ENERGYTOLERANCE    2.0e-04
# define DEFAULT_MAXIMUMHISTORY     10
# define DEFAULT_MAXIMUMSCFCYCLES   100
# define DEFAULT_MINIMUMCOEFFICIENT 1.0e-03
# define DEFAULT_USEEDIIS           False
# define DEFAULT_USEODA             False

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ADIISSCFConverger *ADIISSCFConverger_Allocate ( void )
{
    ADIISSCFConverger *self = NULL ;
    self = ( ADIISSCFConverger * ) Memory_Allocate ( sizeof ( ADIISSCFConverger ) ) ;
    if ( self != NULL )
    {
        self->densityTolerance   = DEFAULT_DENSITYTOLERANCE   ;
        self->diisOff            = DEFAULT_DIISOFF            ;
        self->diisOn             = DEFAULT_DIISON             ;
        self->energyTolerance    = DEFAULT_ENERGYTOLERANCE    ;
        self->maximumHistory     = DEFAULT_MAXIMUMHISTORY     ;
        self->maximumSCFCycles   = DEFAULT_MAXIMUMSCFCYCLES   ;
        self->minimumCoefficient = DEFAULT_MINIMUMCOEFFICIENT ;
        self->useEDIIS           = DEFAULT_USEEDIIS           ;
        self->useODA             = DEFAULT_USEODA             ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
ADIISSCFConverger *ADIISSCFConverger_Clone ( const ADIISSCFConverger *self )
{
    ADIISSCFConverger *new = NULL ;
    if ( self != NULL )
    {
        new = ADIISSCFConverger_Allocate ( ) ;
        new->densityTolerance   = self->densityTolerance   ;
        new->diisOff            = self->diisOff            ;
        new->diisOn             = self->diisOn             ;
        new->energyTolerance    = self->energyTolerance    ;
        new->maximumHistory     = self->maximumHistory     ;
        new->maximumSCFCycles   = self->maximumSCFCycles   ;
        new->minimumCoefficient = self->minimumCoefficient ;
        new->useEDIIS           = self->useEDIIS           ;
        new->useODA             = self->useODA             ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for continuation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ADIISSCFConverger_Continue ( const ADIISSCFConverger *self, const ADIISSCFConvergerState *state )
{
    Boolean doContinue = False ;
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        doContinue = ( state->iteration < self->maximumSCFCycles ) && ( ! state->isConverged ) ;
    }
    return doContinue ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConverger_Deallocate ( ADIISSCFConverger **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for convergence.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ADIISSCFConverger_IsConverged ( const ADIISSCFConverger *self, ADIISSCFConvergerState *state )
{
    Boolean isConverged = False ;
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        isConverged = ( state->iteration > 0 ) && ( state->rmsDifference <= self->densityTolerance ) && ( state->energyChange <= self->energyTolerance ) ;
        state->isConverged = isConverged ;
    }
    return isConverged ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Start an iteration.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConverger_IterateStart ( const ADIISSCFConverger *self, ADIISSCFConvergerState *state, const Real energy )
{
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        auto Real maximumError = 0.0e+00 ;

        /* . Energies. */
        state->energyChange = energy - state->oldEnergy ;
        state->oldEnergy    = energy ;

        /* . Check for a valid density. */
        if ( state->densityP->isValid )
        {
            /* . Save the input data for the ADIIS and DIIS procedures. */
            maximumError     = ADIISSCFConvergerState_DIISSave ( state, energy ) ;
            state->diisError = maximumError ;

            /* . Start the convergence procedure. */
            if ( state->iteration > 0 )
            {
                /* . Check to see which procedure to do. */
                /* . ODA. */
                if ( self->useODA )
                {
                    if ( state->iterationType == ADIISSCFConvergerIterationType_DIIS )
                    {
                        if ( maximumError >= self->diisOff ) state->iterationType = ADIISSCFConvergerIterationType_ODA ;
                    }
                    else if ( state->iterationType == ADIISSCFConvergerIterationType_ODA )
                    {
                        if ( ( state->iteration > 1 ) && ( maximumError <= self->diisOn ) ) state->iterationType = ADIISSCFConvergerIterationType_DIIS ;
                    }
                    else state->iterationType = ADIISSCFConvergerIterationType_ODA ;
                }
                /* . ADIIS. */
                else
                {
                    if ( state->iterationType == ADIISSCFConvergerIterationType_DIIS )
                    {
                        if ( maximumError >= self->diisOff ) state->iterationType = ADIISSCFConvergerIterationType_ADIIS ;
                    }
                    else if ( state->iterationType == ADIISSCFConvergerIterationType_ADIIS )
                    {
                        if ( ( state->iteration > 1 ) && ( maximumError <= self->diisOn ) ) state->iterationType = ADIISSCFConvergerIterationType_DIIS ;
                    }
                    else state->iterationType = ADIISSCFConvergerIterationType_ADIIS ;
                }

                /* . ADIIS. */
                if ( state->iterationType == ADIISSCFConvergerIterationType_ADIIS )
                {
                    ADIISSCFConvergerState_ADIISSetUp ( state, self->useEDIIS ) ;
                }
                /* . DIIS. */
                else if ( state->iterationType == ADIISSCFConvergerIterationType_DIIS )
                {
                    ADIISSCFConvergerState_DIISIterate ( state ) ;
                }
                /* . ODA. */
                else if ( state->iterationType == ADIISSCFConvergerIterationType_ODA )
                {
                    ADIISSCFConvergerState_ODAIterate ( state, energy, self->minimumCoefficient ) ;
                }
            }

            /* . Save old energy for ODA if an ODA step wasn't done. */
            if ( self->useODA && ( state->iterationType != ADIISSCFConvergerIterationType_ODA ) ) state->odaOldEnergy = energy ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Stop an iteration.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ADIISSCFConverger_IterateStop ( const ADIISSCFConverger *self, ADIISSCFConvergerState *state )
{
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        auto Real rmsDifferenceP = 0.0e+00, rmsDifferenceQ = 0.0e+00 ;

        /* . Check for a valid density. */
        if ( state->densityP->isValid )
        {
            /* . Compute the new density and Fock matrices for ADIIS. */
            if ( ( state->iteration > 0 ) && ( state->iterationType == ADIISSCFConvergerIterationType_ADIIS ) )
            {
                ADIISSCFConvergerState_ADIISUpdate ( state ) ;
            }

            /* . Save data for ODA. */
            if ( self->useODA ) ADIISSCFConvergerState_ODASave ( state ) ;

            /* . Update the iteration count only for valid densities. */
            state->iteration++ ;
        }

        /* . Create the densities from the Fock matrices without level shifting. */
        rmsDifferenceP = QCOnePDM_MakeFromFock ( state->densityP, state->orthogonalizer, NULL ) ;
        rmsDifferenceQ = QCOnePDM_MakeFromFock ( state->densityQ, state->orthogonalizer, NULL ) ;

        /* . Save the rms difference. */
        state->rmsDifference = Maximum ( rmsDifferenceP, rmsDifferenceQ ) ;
    }
}
