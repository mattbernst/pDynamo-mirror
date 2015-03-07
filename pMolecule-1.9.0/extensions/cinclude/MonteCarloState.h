/*------------------------------------------------------------------------------
! . File      : MonteCarloState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _MONTECARLOSTATE
# define _MONTECARLOSTATE

# include "Coordinates3.h"
# include "Definitions.h"
# include "Matrix33.h"
# include "NBModelMonteCarlo.h"
# include "NBModelMonteCarloState.h"
# include "SelectionContainer.h"
# include "Status.h"
# include "SymmetryParameters.h"
# include "Vector3.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The Monte Carlo state type. */
typedef struct {
    /* . Counters. */
    Integer                 blocks                ;
    Integer                 moves                 ;
    Integer                 nreject               ;
    Integer                 nrejectm              ;
    Integer                 nrejectt              ;
    Integer                 nrejectv              ;
    Integer                 ntrym                 ;
    Integer                 ntryv                 ;
    /* . Current values. */
    Real                    beta                  ;
    Real                    ecurrent              ;
    Real                    pressure              ;
    Real                    tfact                 ;
    Real                    volume                ;
    /* . Move data. */
    Real                    acceptanceratio       ;
    Real                    rmax                  ;
    Real                    tmax                  ;
    Real                    vmax                  ;
    /* . Statistics. */
    Real                    eav                   ;
    Real                    eav2                  ;
    Real                    etot                  ;
    Real                    etot2                 ;
    Real                    etotb                 ;
    Real                    etotb2                ;
    Real                    hav                   ;
    Real                    hav2                  ;
    Real                    htot                  ;
    Real                    htot2                 ;
    Real                    htotb                 ;
    Real                    htotb2                ;
    Real                    vav                   ;
    Real                    vav2                  ;
    Real                    vtot                  ;
    Real                    vtot2                 ;
    Real                    vtotb                 ;
    Real                    vtotb2                ;
    /* . Arrays to allocate. */
    Real                   *random                ;
    Coordinates3           *oldcoordinates3       ;
    Matrix33               *rotation              ;
    SymmetryParameters     *oldsymmetryParameters ;
    Vector3                *translation           ;
    /* . Aliases. */
    Coordinates3           *coordinates3          ;
    NBModelMonteCarlo      *nbModel               ;
    NBModelMonteCarloState *nbState               ;
    SelectionContainer     *isolates              ;
    SymmetryParameters     *symmetryParameters    ;
} MonteCarloState ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern MonteCarloState *MonteCarloState_Allocate                  ( const Integer nparticles, const Integer nrandom ) ;
extern void             MonteCarloState_Deallocate                ( MonteCarloState **self ) ;
extern void             MonteCarloState_AdjustMoveSizes           ( MonteCarloState  *self ) ;
extern Status           MonteCarloState_MoveIsolate               ( MonteCarloState  *self ) ;
extern Status           MonteCarloState_MoveVolume                ( MonteCarloState  *self ) ;
extern void             MonteCarloState_StatisticsBlockAccumulate ( MonteCarloState  *self ) ;
extern void             MonteCarloState_StatisticsBlockStart      ( MonteCarloState  *self ) ;
extern void             MonteCarloState_StatisticsBlockStop       ( MonteCarloState  *self ) ;
extern void             MonteCarloState_StatisticsStop            ( MonteCarloState  *self ) ;
extern void             MonteCarloState_StatisticsStart           ( MonteCarloState  *self ) ;

# endif
