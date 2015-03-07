/*------------------------------------------------------------------------------
! . File      : CosineAngleContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _COSINEANGLECONTAINER
# define _COSINEANGLECONTAINER

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
    Integer type    ;
} CosineAngle ;

typedef struct {
    Integer  npowers           ;
    Integer  nterms            ;
    Integer *periods           ;
    Real    *coefficients      ;
    Real    *powercoefficients ;
} CosineAngleParameter ;

typedef struct {
    Boolean               QSORTED     ;
    Integer               nperiods    ; /* . The maximum value. */
    Integer               nparameters ;
    Integer               nterms      ;
    CosineAngle          *terms       ;
    CosineAngleParameter *parameters  ;
} CosineAngleContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void                  CosineAngleContainer_ActivateTerms            (       CosineAngleContainer *self ) ;
extern CosineAngleContainer *CosineAngleContainer_Allocate                 ( const Integer nterms, const Integer nparameters ) ;
extern CosineAngleContainer *CosineAngleContainer_Clone                    ( const CosineAngleContainer  *self ) ;
extern void                  CosineAngleContainer_DeactivateFixedAtomTerms (       CosineAngleContainer  *self, Selection *fixedatoms ) ;
extern void                  CosineAngleContainer_DeactivateQCAtomTerms    (       CosineAngleContainer  *self, Selection *qcAtoms, Selection *boundaryatoms ) ;
extern void                  CosineAngleContainer_Deallocate               (       CosineAngleContainer **self ) ;
extern Real                  CosineAngleContainer_Energy                   ( const CosineAngleContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern void                  CosineAngleContainer_FindMaximumPeriod        (       CosineAngleContainer  *self ) ;
extern void                  CosineAngleContainer_MakePowers               (       CosineAngleContainer  *self ) ;
extern CosineAngleContainer *CosineAngleContainer_Merge                    ( const CosineAngleContainer  *self, const CosineAngleContainer *other, const Integer atomincrement ) ;
extern Integer               CosineAngleContainer_NumberOfInactiveTerms    ( const CosineAngleContainer  *self ) ;
extern CosineAngleContainer *CosineAngleContainer_Prune                    (       CosineAngleContainer  *self, Selection *selection ) ;
extern void                  CosineAngleContainer_Sort                     (       CosineAngleContainer  *self ) ;
extern Integer               CosineAngleContainer_UpperBound               (       CosineAngleContainer  *self ) ;

extern void CosineAngleParameter_Allocate ( CosineAngleParameter *self, const Integer nterms ) ;

# endif
