/*------------------------------------------------------------------------------
! . File      : NBModelORCAState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _NBMODELORCASTATE
# define _NBMODELORCASTATE

# include "Coordinates3.h"
# include "Definitions.h"
# include "LJParameterContainer.h"
# include "MMAtomContainer.h"
# include "QCAtomContainer.h"
# include "QCMMLinkAtomCouplingOptions.h"
# include "PairList.h"
# include "Real1DArray.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The NB model state type. */
typedef struct {
    /* . Options. */
    QCMMLinkAtomCoupling     qcmmcoupling     ;
    /* . Energies. */
    Real                     emmel            ;
    Real                     emmel14          ;
    Real                     emmlj            ;
    Real                     emmlj14          ;
    Real                     eqcmmlj          ;
    Real                     eqcmmlj14        ;
    /* . Allocated arrays. */
    Boolean                 *QE14             ;
    Boolean                 *QFREE            ;
    Boolean                 *QINCL            ;
    Boolean                 *QMM              ;
    /* . Aliases. */
    Coordinates3            *coordinates3     ;
    Coordinates3            *gradients3       ;
    LJParameterContainer    *ljParameters     ;
    LJParameterContainer    *ljParameters14   ;
    MMAtomContainer         *mmAtoms          ;
    PairList                *exclusions       ;
    PairList                *interactions14   ;
    QCAtomContainer         *qcAtoms          ;
    /* . MM arrays. */
    Coordinates3            *mmcoordinates3   ;
    Coordinates3            *mmgradients3     ;
    PairList                *mmexclusions     ;
    PairList                *mminteractions14 ;
    Real1DArray             *mmCharges        ;
} NBModelOrcaState ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern NBModelOrcaState *NBModelORCAState_Allocate   ( const Integer n ) ;
extern void              NBModelORCAState_Deallocate ( NBModelOrcaState **self ) ;
extern void              NBModelORCAState_Finalize   ( NBModelOrcaState  *self ) ;
extern void              NBModelORCAState_Initialize ( NBModelOrcaState  *self, Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern NBModelOrcaState *NBModelORCAState_Setup      ( MMAtomContainer *mmAtoms, QCAtomContainer *qcAtoms, Selection *fixedatoms, PairList *exclusions, PairList *interactions14,
                                                                                                        LJParameterContainer *ljParameters, LJParameterContainer *ljParameters14,
                                                                                                                                      const QCMMLinkAtomCoupling qcmmcoupling ) ;

# endif
