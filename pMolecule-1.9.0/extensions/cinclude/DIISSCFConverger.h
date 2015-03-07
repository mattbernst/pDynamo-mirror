/*------------------------------------------------------------------------------
! . File      : DIISSCFConverger.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _DIISSCFCONVERGER
# define _DIISSCFCONVERGER

# include "Definitions.h"
# include "DIISSCFConvergerState.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The DIIS SCF converger type. */
typedef struct {
    Boolean QUSERCA          ;
    Integer maximumSCFCycles ;
    Integer ndiis            ;
    Real    densityTolerance ;
    Real    diisonset        ;
    Real    energytolerance  ;
} DIISSCFConverger ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern DIISSCFConverger *DIISSCFConverger_Allocate                ( void ) ;
extern DIISSCFConverger *DIISSCFConverger_Clone                   ( const DIISSCFConverger  *self ) ;
extern Boolean           DIISSCFConverger_Continue                ( const DIISSCFConverger  *self, const DIISSCFConvergerState *convergerstate ) ;
extern Boolean           DIISSCFConverger_Converged               ( const DIISSCFConverger  *self,       DIISSCFConvergerState *convergerstate ) ;
extern void              DIISSCFConverger_Deallocate              (       DIISSCFConverger **self ) ;
extern void              DIISSCFConverger_IterateWithoutDensities ( const DIISSCFConverger  *self, DIISSCFConvergerState *convergerstate, const Real energy ) ;
extern void              DIISSCFConverger_MakeDensities           ( const DIISSCFConverger  *self, DIISSCFConvergerState *convergerstate ) ;

# endif
