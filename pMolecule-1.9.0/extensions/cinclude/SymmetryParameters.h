/*------------------------------------------------------------------------------
! . File      : SymmetryParameters.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SYMMETRYPARAMETERS
# define _SYMMETRYPARAMETERS

# include "Coordinates3.h"
# include "Definitions.h"
# include "Matrix33.h"
# include "Selection.h"
# include "SelectionContainer.h"
# include "Status.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The symmetryParameters type. */
typedef struct {
    Boolean   QM       ;
    Real      a        ;
    Real      b        ;
    Real      c        ;
    Real      alpha    ;
    Real      beta     ;
    Real      gamma    ;
    Matrix33 *inverseM ;
    Matrix33 *M        ;
} SymmetryParameters ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern SymmetryParameters *SymmetryParameters_Allocate                          ( void ) ;
extern Status              SymmetryParameters_CenterCoordinates3ByIndex         ( const SymmetryParameters  *self, const Selection *selection, Coordinates3 *coordinates3 ) ;
extern Status              SymmetryParameters_CenterCoordinates3ByIsolate       ( const SymmetryParameters  *self, const SelectionContainer *isolates, Selection *selection, Coordinates3 *coordinates3 ) ;
extern void                SymmetryParameters_ClearM                            (       SymmetryParameters  *self ) ;
extern void                SymmetryParameters_CopyTo                            ( const SymmetryParameters  *self, SymmetryParameters *other ) ;
extern void                SymmetryParameters_Deallocate                        (       SymmetryParameters **self ) ;
extern void                SymmetryParameters_Displacement                      ( const SymmetryParameters  *self, const Integer a, const Integer b, const Integer c, Vector3 *displacement ) ;
extern void                SymmetryParameters_FindBoxSearchLimits               ( const SymmetryParameters  *self   ,
                                                                                  const Vector3             *lower  ,
                                                                                  const Vector3             *upper  ,
                                                                                  const Vector3             *ilower ,
                                                                                  const Vector3             *iupper ,
                                                                                        Integer             *alow   ,
                                                                                        Integer             *ahigh  ,
                                                                                        Integer             *blow   ,
                                                                                        Integer             *bhigh  ,
                                                                                        Integer             *clow   ,
                                                                                        Integer             *chigh  ) ;
extern Status              SymmetryParameters_FindCenteringTranslation          ( const SymmetryParameters  *self, const Vector3 *point, Vector3 *translation ) ;
extern Boolean             SymmetryParameters_IsMinimumImageConventionSatisfied ( const SymmetryParameters  *self, const Real length ) ;
extern Boolean             SymmetryParameters_IsOrthorhombic                    ( const SymmetryParameters  *self ) ;
extern void                SymmetryParameters_IsotropicScale                    (       SymmetryParameters  *self, const Real scale ) ;
extern void                SymmetryParameters_MakeM                             (       SymmetryParameters  *self ) ;
extern void                SymmetryParameters_MakeMinimumImageVector3           ( const SymmetryParameters  *self, Vector3 *r, Vector3 *dr ) ;
extern void                SymmetryParameters_MakeMinimumImageXYZ               ( const SymmetryParameters  *self, Real *x, Real *y, Real *z, Real *dx, Real *dy, Real *dz ) ;
extern void                SymmetryParameters_MakeOrthorhombicWidths            ( const SymmetryParameters  *self, Vector3 *widths ) ;
extern void                SymmetryParameters_SetCrystalParameters              (       SymmetryParameters  *self  ,
                                                                                  const Real                 a     ,
                                                                                  const Real                 b     ,
                                                                                  const Real                 c     ,
                                                                                  const Real                 alpha ,
                                                                                  const Real                 beta  ,
                                                                                  const Real                 gamma ) ;
extern Real                SymmetryParameters_Volume                            ( const SymmetryParameters  *self ) ;

# endif
