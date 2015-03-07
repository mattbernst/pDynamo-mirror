/*------------------------------------------------------------------------------
! . File      : HarmonicAngleContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _HARMONICANGLECONTAINER
# define _HARMONICANGLECONTAINER

# include "Definitions.h"
# include "Coordinates3.h"
# include "Selection.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
typedef struct {
    Boolean QACTIVE ;
    Integer atom1   ;
    Integer atom2   ;
    Integer atom3   ;
    Integer type    ;
} HarmonicAngle ;

typedef struct {
    Real eq ;
    Real fc ;
} HarmonicAngleParameter ;

typedef struct {
    Boolean                 QSORTED     ;
    Integer                 nparameters ;
    Integer                 nterms      ;
    HarmonicAngle          *terms       ;
    HarmonicAngleParameter *parameters  ;
} HarmonicAngleContainer ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern void                    HarmonicAngleContainer_ActivateTerms            (       HarmonicAngleContainer  *self ) ;
extern HarmonicAngleContainer *HarmonicAngleContainer_Allocate                 ( const Integer nterms, const Integer nparameters ) ;
extern HarmonicAngleContainer *HarmonicAngleContainer_Clone                    ( const HarmonicAngleContainer  *self ) ;
extern void                    HarmonicAngleContainer_DeactivateFixedAtomTerms (       HarmonicAngleContainer  *self, Selection *fixedatoms ) ;
extern void                    HarmonicAngleContainer_DeactivateQCAtomTerms    (       HarmonicAngleContainer  *self, Selection *qcAtoms, Selection *boundaryatoms ) ;
extern void                    HarmonicAngleContainer_Deallocate               (	    HarmonicAngleContainer **self ) ;
extern Real                    HarmonicAngleContainer_Energy                   ( const HarmonicAngleContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern HarmonicAngleContainer *HarmonicAngleContainer_Merge                    ( const HarmonicAngleContainer  *self, const HarmonicAngleContainer *other, const Integer atomincrement ) ;
extern Integer                 HarmonicAngleContainer_NumberOfInactiveTerms    ( const HarmonicAngleContainer  *self ) ;
extern HarmonicAngleContainer *HarmonicAngleContainer_Prune                    (       HarmonicAngleContainer  *self, Selection *selection ) ;
extern void                    HarmonicAngleContainer_Sort                     (       HarmonicAngleContainer  *self ) ;
extern Integer                 HarmonicAngleContainer_UpperBound               (       HarmonicAngleContainer  *self ) ;

# endif
