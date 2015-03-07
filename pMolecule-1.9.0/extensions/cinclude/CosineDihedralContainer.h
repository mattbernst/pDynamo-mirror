/*------------------------------------------------------------------------------
! . File      : CosineDihedralContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _COSINEDIHEDRALCONTAINER
# define _COSINEDIHEDRALCONTAINER

# include "Definitions.h"
# include "Coordinates3.h"
# include "Selection.h"

/*----------------------------------------------------------------------------------------------------------------------------------
!
! . A cosine angle energy term is defined as follows:
!
!   (1) E = Sum_n c_n cos ( n alpha )      ( n = 0, 1, 2 ... )
!
! . or:
!
!   (2) E = Sum_m c_m cos ( m alpha ) +    ( m = 0, 2, 4 ... )
!           Sum_n c_n sin ( n alpha )      ( n = 1, 3, 5 ... )
!
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    Boolean QACTIVE ;
    Integer atom1   ;
    Integer atom2   ;
    Integer atom3   ;
    Integer atom4   ;
    Integer type    ;
} CosineDihedral ;

typedef struct {
    Integer  npowers           ;
    Integer  nterms            ;
    Integer *periods           ;
    Real    *coefficients      ;
    Real    *powercoefficients ;
} CosineDihedralParameter ;

typedef struct {
    Boolean                  QSORTED     ;
    Integer                  nperiods    ; /* . The maximum value. */
    Integer                  nparameters ;
    Integer                  nterms      ;
    CosineDihedral          *terms       ;
    CosineDihedralParameter *parameters  ;
} CosineDihedralContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void                     CosineDihedralContainer_ActivateTerms            (       CosineDihedralContainer *self ) ;
extern CosineDihedralContainer *CosineDihedralContainer_Allocate                 ( const Integer nterms, const Integer nparameters ) ;
extern CosineDihedralContainer *CosineDihedralContainer_Clone                    ( const CosineDihedralContainer  *self ) ;
extern void                     CosineDihedralContainer_DeactivateFixedAtomTerms (       CosineDihedralContainer  *self, Selection *fixedatoms ) ;
extern void                     CosineDihedralContainer_DeactivateQCAtomTerms    (       CosineDihedralContainer  *self, Selection *qcAtoms, Selection *boundaryatoms ) ;
extern void                     CosineDihedralContainer_Deallocate               (       CosineDihedralContainer **self ) ;
extern Real                     CosineDihedralContainer_Energy                   ( const CosineDihedralContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern void                     CosineDihedralContainer_FindMaximumPeriod        (       CosineDihedralContainer  *self ) ;
extern void                     CosineDihedralContainer_MakePowers               (       CosineDihedralContainer  *self ) ;
extern CosineDihedralContainer *CosineDihedralContainer_Merge                    ( const CosineDihedralContainer  *self, const CosineDihedralContainer *other, const Integer atomincrement ) ;
extern Integer                  CosineDihedralContainer_NumberOfInactiveTerms    ( const CosineDihedralContainer  *self ) ;
extern CosineDihedralContainer *CosineDihedralContainer_Prune                    (       CosineDihedralContainer  *self, Selection *selection ) ;
extern void                     CosineDihedralContainer_Sort                     (       CosineDihedralContainer  *self ) ;
extern Integer                  CosineDihedralContainer_UpperBound               (       CosineDihedralContainer  *self ) ;

extern void CosineDihedralParameter_Allocate ( CosineDihedralParameter *self, const Integer nterms ) ;

# endif
