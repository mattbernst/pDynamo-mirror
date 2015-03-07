/*------------------------------------------------------------------------------
! . File      : NBModelFullState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _NBMODELFULLSTATE
# define _NBMODELFULLSTATE

# include "Coordinates3.h"
# include "Definitions.h"
# include "LJParameterContainer.h"
# include "MMAtomContainer.h"
# include "QCAtomContainer.h"
# include "PairList.h"
# include "Real1DArray.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The NB model state type. */
typedef struct {
    /* . Energies. */
    Real                  emmel          ;
    Real                  emmel14        ;
    Real                  emmlj          ;
    Real                  emmlj14        ;
    Real                  eqcmmlj        ;
    Real                  eqcmmlj14      ;
    /* . Allocated arrays. */
    Boolean              *QE14           ;
    Boolean              *QFREE          ;
    Boolean              *QINCL          ;
    Boolean              *QMM            ;
    Integer              *baindex        ;
    /* . Aliases. */
    Coordinates3         *coordinates3   ;
    Coordinates3         *gradients3     ;
    LJParameterContainer *ljParameters   ;
    LJParameterContainer *ljParameters14 ;
    MMAtomContainer      *mmAtoms        ;
    PairList             *exclusions     ;
    PairList             *interactions14 ;
    QCAtomContainer      *qcAtoms        ;
    Real1DArray          *qcCharges      ;
    Real1DArray          *qcmmPotentials ;
} NBModelFullState ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern NBModelFullState *NBModelFullState_Allocate   ( const Integer n ) ;
extern void              NBModelFullState_Deallocate ( NBModelFullState **self ) ;
extern void              NBModelFullState_Initialize ( NBModelFullState  *self, Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern NBModelFullState *NBModelFullState_Setup      ( MMAtomContainer *mmAtoms, QCAtomContainer *qcAtoms, Selection *fixedAtoms, PairList *exclusions, PairList *interactions14,
                                                                                                        LJParameterContainer *ljParameters, LJParameterContainer *ljParameters14,
                                                                                                                          Real1DArray *qcCharges, Real1DArray *qcmmPotentials ) ;

# endif
