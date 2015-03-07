/*------------------------------------------------------------------------------
! . File      : FourierDihedralContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _FOURIERDIHEDRALCONTAINER
# define _FOURIERDIHEDRALCONTAINER

# include "Definitions.h"
# include "Coordinates3.h"
# include "Selection.h"

/*------------------------------------------------------------------------------
! . A fourier dihedral energy term is defined as follows:
!
!   E = fc * ( 1 + cos ( period * phi - delta ) )
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
} FourierDihedral ;

typedef struct {
    Integer period   ;
    Real    fc       ;
    Real    phase    ;
    Real    cosphase ;
    Real    sinphase ;
} FourierDihedralParameter ;

typedef struct {
    Boolean                   QSORTED     ;
    Integer                   nparameters ;
    Integer                   nterms      ;
    FourierDihedral          *terms       ;
    FourierDihedralParameter *parameters  ;
} FourierDihedralContainer ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern void                      FourierDihedralContainer_ActivateTerms            (       FourierDihedralContainer  *self ) ;
extern FourierDihedralContainer *FourierDihedralContainer_Allocate                 ( const Integer nterms, const Integer nparameters ) ;
extern FourierDihedralContainer *FourierDihedralContainer_Clone                    ( const FourierDihedralContainer  *self ) ;
extern void                      FourierDihedralContainer_DeactivateFixedAtomTerms (       FourierDihedralContainer  *self, Selection *fixedatoms ) ;
extern void                      FourierDihedralContainer_DeactivateQCAtomTerms    (       FourierDihedralContainer  *self, Selection *qcAtoms, Selection *boundaryatoms ) ;
extern void                      FourierDihedralContainer_Deallocate               (       FourierDihedralContainer **self ) ;
extern Real                      FourierDihedralContainer_Energy                   ( const FourierDihedralContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern void                      FourierDihedralContainer_FillCosSinPhases         (       FourierDihedralContainer  *self ) ;
extern FourierDihedralContainer *FourierDihedralContainer_Merge                    ( const FourierDihedralContainer  *self, const FourierDihedralContainer *other, const Integer atomincrement ) ;
extern Integer                   FourierDihedralContainer_NumberOfInactiveTerms    ( const FourierDihedralContainer  *self ) ;
extern FourierDihedralContainer *FourierDihedralContainer_Prune                    (       FourierDihedralContainer  *self, Selection *selection ) ;
extern void                      FourierDihedralContainer_Sort                     (       FourierDihedralContainer  *self ) ;
extern Integer                   FourierDihedralContainer_UpperBound               (       FourierDihedralContainer  *self ) ;

# endif
