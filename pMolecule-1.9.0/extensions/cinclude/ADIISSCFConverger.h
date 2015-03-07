/*------------------------------------------------------------------------------
! . File      : ADIISSCFConverger.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _ADIISSCFCONVERGER
# define _ADIISSCFCONVERGER

# include "ADIISSCFConvergerState.h"
# include "Boolean.h"
# include "Definitions.h"
# include "Integer.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The ADIIS SCF converger type. */
typedef struct {
    Boolean useEDIIS           ;
    Boolean useODA             ;
    Integer maximumHistory     ;
    Integer maximumSCFCycles   ;
    Real    densityTolerance   ;
    Real    diisOff            ;
    Real    diisOn             ;
    Real    energyTolerance    ;
    Real    minimumCoefficient ;
} ADIISSCFConverger ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern ADIISSCFConverger *ADIISSCFConverger_Allocate     ( void ) ;
extern ADIISSCFConverger *ADIISSCFConverger_Clone        ( const ADIISSCFConverger  *self ) ;
extern Boolean            ADIISSCFConverger_Continue     ( const ADIISSCFConverger  *self, const ADIISSCFConvergerState *state ) ;
extern void               ADIISSCFConverger_Deallocate   (       ADIISSCFConverger **self ) ;
extern Boolean            ADIISSCFConverger_IsConverged  ( const ADIISSCFConverger  *self,       ADIISSCFConvergerState *state ) ;
extern void               ADIISSCFConverger_IterateStart ( const ADIISSCFConverger  *self, ADIISSCFConvergerState *state, const Real energy ) ;
extern void               ADIISSCFConverger_IterateStop  ( const ADIISSCFConverger  *self, ADIISSCFConvergerState *state ) ;

# endif
