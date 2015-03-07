/*------------------------------------------------------------------------------
! . File      : NBModelORCA.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _NBMODELORCA
# define _NBMODELORCA

# include "NBModelORCAState.h"
# include "QCMMLinkAtomCouplingOptions.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The NB model type. */
typedef struct {
    Real                 electrostaticscale14 ;
    QCMMLinkAtomCoupling qcmmcoupling         ;
} NBModelOrca ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern NBModelOrca *NBModelORCA_Allocate       ( void ) ;
extern NBModelOrca *NBModelORCA_Clone          ( const NBModelOrca  *self ) ;
extern void         NBModelORCA_Deallocate     (       NBModelOrca **self ) ;
extern void         NBModelORCA_MMMMEnergy     ( const NBModelOrca  *self, NBModelOrcaState *nbState ) ;
extern void         NBModelORCA_QCMMEnergyLJ   ( const NBModelOrca  *self, NBModelOrcaState *nbState ) ;

# endif
