/*------------------------------------------------------------------------------
! . File      : CMAPDihedralContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _CMAPDIHEDRALCONTAINER
# define _CMAPDIHEDRALCONTAINER

# include "Definitions.h"
# include "BicubicSpline.h"
# include "Boolean.h"
# include "Coordinates3.h"
# include "Integer.h"
# include "Real.h"
# include "Selection.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . A CMAP dihedral term is for CHARMM force fields. It is a table-based energy
! . term that is a function of two dihedral angles.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    Boolean isActive ;
    Integer atom1    ;
    Integer atom2    ;
    Integer atom3    ;
    Integer atom4    ;
    Integer atom5    ;
    Integer atom6    ;
    Integer atom7    ;
    Integer atom8    ;
    Integer type     ;
} CMAPDihedral ;

typedef struct {
    Boolean         isSorted    ;
    Integer         nparameters ;
    Integer         nterms      ;
    CMAPDihedral   *terms       ;
    BicubicSpline **parameters  ;
} CMAPDihedralContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void                   CMAPDihedralContainer_ActivateTerms            (       CMAPDihedralContainer  *self ) ;
extern CMAPDihedralContainer *CMAPDihedralContainer_Allocate                 ( const Integer nterms, const Integer nparameters ) ;
extern CMAPDihedralContainer *CMAPDihedralContainer_Clone                    ( const CMAPDihedralContainer  *self ) ;
extern void                   CMAPDihedralContainer_DeactivateFixedAtomTerms (       CMAPDihedralContainer  *self, Selection *fixedatoms ) ;
extern void                   CMAPDihedralContainer_DeactivateQCAtomTerms    (       CMAPDihedralContainer  *self, Selection *qcAtoms, Selection *boundaryatoms ) ;
extern void                   CMAPDihedralContainer_Deallocate               (       CMAPDihedralContainer **self ) ;
extern Real                   CMAPDihedralContainer_Energy                   ( const CMAPDihedralContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern CMAPDihedralContainer *CMAPDihedralContainer_Merge                    ( const CMAPDihedralContainer  *self, const CMAPDihedralContainer *other, const int atomincrement ) ;
extern Integer                CMAPDihedralContainer_NumberOfInactiveTerms    ( const CMAPDihedralContainer  *self ) ;
extern CMAPDihedralContainer *CMAPDihedralContainer_Prune                    (       CMAPDihedralContainer  *self, Selection *selection ) ;
extern void                   CMAPDihedralContainer_Sort                     (       CMAPDihedralContainer  *self ) ;
extern Integer                CMAPDihedralContainer_UpperBound               (       CMAPDihedralContainer  *self ) ;

# endif
