/*------------------------------------------------------------------------------
! . File      : HarmonicBondContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _HARMONICBONDCONTAINER
# define _HARMONICBONDCONTAINER

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
    Integer type    ;
} HarmonicBond ;

typedef struct {
    Real eq ;
    Real fc ;
} HarmonicBondParameter ;

typedef struct {
    Boolean                QSORTED     ;
    Integer                nparameters ;
    Integer                nterms      ;
    HarmonicBond          *terms       ;
    HarmonicBondParameter *parameters  ;
} HarmonicBondContainer ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern void                   HarmonicBondContainer_ActivateTerms            (       HarmonicBondContainer  *self ) ;
extern HarmonicBondContainer *HarmonicBondContainer_Allocate                 ( const Integer nterms, const Integer nparameters ) ;
extern HarmonicBondContainer *HarmonicBondContainer_Clone                    ( const HarmonicBondContainer  *self ) ;
extern void                   HarmonicBondContainer_DeactivateFixedAtomTerms (       HarmonicBondContainer  *self, Selection *fixedatoms ) ;
extern void                   HarmonicBondContainer_DeactivateQCAtomTerms    (       HarmonicBondContainer  *self, Selection *qcAtoms, Selection *boundaryatoms ) ;
extern void                   HarmonicBondContainer_Deallocate               (       HarmonicBondContainer **self ) ;
extern Real                   HarmonicBondContainer_Energy                   ( const HarmonicBondContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern Integer                HarmonicBondContainer_IdentifyBoundaryAtoms    (       HarmonicBondContainer  *self, Selection *qcAtoms, Integer **mmboundary, Integer **qcpartners ) ;
extern HarmonicBondContainer *HarmonicBondContainer_Merge                    ( const HarmonicBondContainer  *self, const HarmonicBondContainer *other, const Integer atomincrement ) ;
extern Integer                HarmonicBondContainer_NumberOfInactiveTerms    ( const HarmonicBondContainer  *self ) ;
extern HarmonicBondContainer *HarmonicBondContainer_Prune                    (       HarmonicBondContainer  *self, Selection *selection ) ;
extern void                   HarmonicBondContainer_Sort                     (       HarmonicBondContainer  *self ) ;
extern Integer                HarmonicBondContainer_UpperBound               (       HarmonicBondContainer  *self ) ;

# endif
