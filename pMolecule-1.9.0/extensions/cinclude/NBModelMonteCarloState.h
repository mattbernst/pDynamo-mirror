/*------------------------------------------------------------------------------
! . File      : NBModelMonteCarloState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _NBMODELMONTECARLOSTATE
# define _NBMODELMONTECARLOSTATE

# include "Coordinates3.h"
# include "Definitions.h"
# include "LJParameterContainer.h"
# include "MMAtomContainer.h"
# include "SelectionContainer.h"
# include "SymmetryParameters.h"
# include "Vector3.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The NB model state type. */
typedef struct {
    /* . Dimension. */
    Integer               nisolates          ;
    /* . Scaling. */
    Integer               isolatescale       ;
    Real                  chargeScale        ;
    Real                  epsilonScale       ;
    Real                  sigmaScale         ;
    /* . Energies. */
    Real                  efmmel             ;
    Real                  e1mmel             ;
    Real                  efmmlj             ;
    Real                  e1mmlj             ;
    /* . Allocated arrays. */
    Boolean              *QFREE              ;
    Vector3             **centers            ;
    /* . Aliases. */
    Coordinates3         *coordinates3       ;
    LJParameterContainer *ljParameters       ;
    MMAtomContainer      *mmAtoms            ;
    SelectionContainer   *isolates           ;
    SymmetryParameters   *symmetryParameters ;
} NBModelMonteCarloState ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern NBModelMonteCarloState *NBModelMonteCarloState_Allocate                          ( const Integer n ) ;
extern void                    NBModelMonteCarloState_Deallocate                        ( NBModelMonteCarloState **self ) ;
extern void                    NBModelMonteCarloState_Initialize                        ( NBModelMonteCarloState  *self, Coordinates3 *coordinates3, SymmetryParameters *symmetryParameters ) ;
extern void                    NBModelMonteCarloState_ScaleIsolateInteractionParameters ( NBModelMonteCarloState  *self, const Integer isolate, const Real chargeScale, const Real epsilonScale, const Real sigmaScale ) ;
extern NBModelMonteCarloState *NBModelMonteCarloState_Setup                             ( SelectionContainer *isolates, MMAtomContainer *mmAtoms, LJParameterContainer *ljParameters, Selection *fixedatoms ) ;

# endif
