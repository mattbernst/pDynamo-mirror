/*------------------------------------------------------------------------------
! . File      : NBModelFull.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _NBMODELFULL
# define _NBMODELFULL

# include "NBModelFullState.h"
# include "QCMMLinkAtomCouplingOptions.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The NB model type. */
typedef struct {
    Real                 dielectric           ;
    Real                 electrostaticscale14 ;
    QCMMLinkAtomCoupling qcmmcoupling         ;
} NBModelFull ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern NBModelFull *NBModelFull_Allocate       ( void ) ;
extern NBModelFull *NBModelFull_Clone          ( const NBModelFull  *self ) ;
extern void         NBModelFull_Deallocate     (       NBModelFull **self ) ;
extern void         NBModelFull_MMMMEnergy     ( const NBModelFull  *self, NBModelFullState *nbState ) ;
extern void         NBModelFull_QCMMEnergyLJ   ( const NBModelFull  *self, NBModelFullState *nbState ) ;
extern void         NBModelFull_QCMMGradients  ( const NBModelFull  *self, NBModelFullState *nbState ) ;
extern void         NBModelFull_QCMMPotentials ( const NBModelFull  *self, NBModelFullState *nbState ) ;

# endif
