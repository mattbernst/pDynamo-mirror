/*------------------------------------------------------------------------------
! . File      : NBModelMonteCarlo.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _NBMODELMONTECARLO
# define _NBMODELMONTECARLO

# include "Definitions.h"
# include "NBModelMonteCarloState.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The NB model type. */
typedef struct {
    Real buffer      ;
    Real cutoff      ;
    Real dielectric  ;
    Real underflowel ;
    Real underflowlj ;
} NBModelMonteCarlo ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern NBModelMonteCarlo *NBModelMonteCarlo_Allocate         ( void ) ;
extern NBModelMonteCarlo *NBModelMonteCarlo_Clone            ( const NBModelMonteCarlo  *self ) ;
extern void               NBModelMonteCarlo_Deallocate       (       NBModelMonteCarlo **self ) ;
extern Real               NBModelMonteCarlo_MMMMEnergyFull   ( const NBModelMonteCarlo  *self, NBModelMonteCarloState *nbState ) ;
extern Real               NBModelMonteCarlo_MMMMEnergySingle ( const NBModelMonteCarlo  *self, const Integer isolate, NBModelMonteCarloState *nbState ) ;

# endif
