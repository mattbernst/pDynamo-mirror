/*------------------------------------------------------------------------------
! . File      : HarmonicImproperContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _HARMONICIMPROPERCONTAINER
# define _HARMONICIMPROPERCONTAINER

# include "Definitions.h"
# include "Coordinates3.h"
# include "Selection.h"

/*------------------------------------------------------------------------------
! . A harmonic improper energy term is defined as follows:
!
!   E = fc * ( phi - phi0 )^2
!
!-----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
typedef struct {
    Boolean QACTIVE ;
    Integer atom1   ;
    Integer atom2   ;
    Integer atom3   ;
    Integer atom4   ;
    Integer type    ;
} HarmonicImproper ;

typedef struct {
    Real eq    ;
    Real fc    ;
    Real coseq ;
    Real sineq ;
} HarmonicImproperParameter ;

typedef struct {
    Boolean                    QSORTED     ;
    Integer                    nparameters ;
    Integer                    nterms      ;
    HarmonicImproper          *terms       ;
    HarmonicImproperParameter *parameters  ;
} HarmonicImproperContainer ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern void                       HarmonicImproperContainer_ActivateTerms            (       HarmonicImproperContainer  *self ) ;
extern HarmonicImproperContainer *HarmonicImproperContainer_Allocate                 ( const Integer nterms, const Integer nparameters ) ;
extern HarmonicImproperContainer *HarmonicImproperContainer_Clone                    ( const HarmonicImproperContainer  *self ) ;
extern void                       HarmonicImproperContainer_DeactivateFixedAtomTerms (       HarmonicImproperContainer  *self, Selection *fixedatoms ) ;
extern void                       HarmonicImproperContainer_DeactivateQCAtomTerms    (       HarmonicImproperContainer  *self, Selection *qcAtoms, Selection *boundaryatoms ) ;
extern void                       HarmonicImproperContainer_Deallocate               (       HarmonicImproperContainer **self ) ;
extern Real                       HarmonicImproperContainer_Energy                   ( const HarmonicImproperContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern void                       HarmonicImproperContainer_FillCosSinValues         (       HarmonicImproperContainer  *self ) ;
extern HarmonicImproperContainer *HarmonicImproperContainer_Merge                    ( const HarmonicImproperContainer  *self, const HarmonicImproperContainer *other, const Integer atomincrement ) ;
extern Integer                    HarmonicImproperContainer_NumberOfInactiveTerms    ( const HarmonicImproperContainer  *self ) ;
extern HarmonicImproperContainer *HarmonicImproperContainer_Prune                    (       HarmonicImproperContainer  *self, Selection *selection ) ;
extern void                       HarmonicImproperContainer_Sort                     (       HarmonicImproperContainer  *self ) ;
extern Integer                    HarmonicImproperContainer_UpperBound               (       HarmonicImproperContainer  *self ) ;

# endif
