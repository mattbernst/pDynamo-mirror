/*------------------------------------------------------------------------------
! . File      : NBModelABFS.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _NBMODELABFS
# define _NBMODELABFS

# include "Boolean.h"
# include "Integer.h"
# include "NBModelABFSState.h"
# include "PairListGenerator.h"
# include "PairwiseInteraction.h"
# include "QCMMLinkAtomCouplingOptions.h"
# include "Real.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The NB model type. */
typedef struct {
    Boolean checkForInverses          ; /* . For testing. */
    Boolean useCentering              ;
    Integer imageExpandFactor         ; /* . For testing. */
    Real    dampingCutoff             ;
    Real    dielectric                ;
    Real    electrostaticScale14      ;
    Real    innerCutoff               ;
    Real    listCutoff                ;
    Real    outerCutoff               ;
    QCMMLinkAtomCoupling qcmmCoupling ;
} NBModelABFS ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern NBModelABFS *NBModelABFS_Allocate       ( void ) ;
extern NBModelABFS *NBModelABFS_Clone          ( const NBModelABFS  *self ) ;
extern void         NBModelABFS_Deallocate     (       NBModelABFS **self ) ;
extern void         NBModelABFS_MMMMEnergy     ( const NBModelABFS  *self, const PairwiseInteractionABFS *mmmmPairwiseInteraction, NBModelABFSState *nbState ) ;
extern void         NBModelABFS_QCMMEnergyLJ   ( const NBModelABFS  *self, const PairwiseInteractionABFS *mmmmPairwiseInteraction, NBModelABFSState *nbState ) ;
extern void         NBModelABFS_QCMMGradients  ( const NBModelABFS  *self, const PairwiseInteractionABFS *qcmmPairwiseInteraction, const PairwiseInteractionABFS *qcqcPairwiseInteraction, NBModelABFSState *nbState ) ;
extern void         NBModelABFS_QCMMPotentials ( const NBModelABFS  *self, const PairwiseInteractionABFS *qcmmPairwiseInteraction, const PairwiseInteractionABFS *qcqcPairwiseInteraction, NBModelABFSState *nbState ) ;
extern Boolean      NBModelABFS_Update         ( const NBModelABFS  *self, const PairListGenerator *generator, NBModelABFSState *nbState, Status *status ) ;

# endif
