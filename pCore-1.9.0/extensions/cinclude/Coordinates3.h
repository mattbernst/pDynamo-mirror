/*------------------------------------------------------------------------------
! . File      : Coordinates3.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _COORDINATES3
# define _COORDINATES3

# include "Definitions.h"
# include "Matrix33.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "RegularGrid.h"
# include "RegularGridOccupancy.h"
# include "Selection.h"
# include "Status.h"
# include "SymmetricMatrix.h"
# include "Transformation3.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Type. */
# define Coordinates3                     Real2DArray

/* . Procedures. */
# define Coordinates3_1DSlice             Real2DArray_1DSlice
# define Coordinates3_2DSlice             Real2DArray_Slice
# define Coordinates3_AbsoluteMaximum     Real2DArray_AbsoluteMaximum
# define Coordinates3_AddScaledArray      Real2DArray_AddScaledArray
# define Coordinates3_Clone               Real2DArray_Clone
# define Coordinates3_CopyFromArray       Real2DArray_CopyFromArrayByRow
# define Coordinates3_CopyTo              Real2DArray_CopyTo
# define Coordinates3_CopyToArray         Real2DArray_CopyToArrayByRow
# define Coordinates3_Deallocate          Real2DArray_Deallocate
# define Coordinates3_GetItem             Real2DArray_GetItem
# define Coordinates3_Length              Real2DArray_Length
# define Coordinates3_Prune               Real2DArray_PruneByRow
# define Coordinates3_Resize              Real2DArray_Resize
# define Coordinates3_RootMeanSquare      Real2DArray_RootMeanSquare
# define Coordinates3_Scale               Real2DArray_Scale
# define Coordinates3_Set                 Real2DArray_Set
# define Coordinates3_SetByRow            Real2DArray_SetByRow
# define Coordinates3_SetItem             Real2DArray_SetItem
# define Coordinates3_SliceCoordinates3   Real2DArray_Slice
# define Coordinates3_SliceVector3        Real2DArray_1DSlice
# define Coordinates3_ViewOfRaw           Real2DArray_ViewOfRaw

/* . Macros. */
# define Coordinates3_Columns             Real2DArray_Columns
# define Coordinates3_Data                Real2DArray_Data
# define Coordinates3_Item                Real2DArray_Item
# define Coordinates3_ItemIndex           Real2DArray_ItemIndex
# define Coordinates3_ItemPointer         Real2DArray_ItemPointer
# define Coordinates3_RowPointer          Real2DArray_RowPointer
# define Coordinates3_Rows                Real2DArray_Rows

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Coordinates3_DecrementRow( self, i, xij, yij, zij ) \
        { \
	    Coordinates3_Item ( self, i, 0 ) -= xij ; \
	    Coordinates3_Item ( self, i, 1 ) -= yij ; \
	    Coordinates3_Item ( self, i, 2 ) -= zij ; \
        }

# define Coordinates3_DifferenceRow( self, i, j, xij, yij, zij ) \
        { \
	    xij = Coordinates3_Item ( self, i, 0 ) - Coordinates3_Item ( self, j, 0 ) ; \
	    yij = Coordinates3_Item ( self, i, 1 ) - Coordinates3_Item ( self, j, 1 ) ; \
	    zij = Coordinates3_Item ( self, i, 2 ) - Coordinates3_Item ( self, j, 2 ) ; \
        }

# define Coordinates3_GetRow( self, i, x, y, z ) \
        { \
	   x = Coordinates3_Item ( self, i, 0 ) ; \
	   y = Coordinates3_Item ( self, i, 1 ) ; \
	   z = Coordinates3_Item ( self, i, 2 ) ; \
        }

# define Coordinates3_IncrementRow( self, i, xij, yij, zij ) \
        { \
	    Coordinates3_Item ( self, i, 0 ) += xij ; \
	    Coordinates3_Item ( self, i, 1 ) += yij ; \
	    Coordinates3_Item ( self, i, 2 ) += zij ; \
        }

# define Coordinates3_ScaleRow( self, i, value ) \
        { \
	    Coordinates3_Item ( self, i, 0 ) *= value ; \
	    Coordinates3_Item ( self, i, 1 ) *= value ; \
	    Coordinates3_Item ( self, i, 2 ) *= value ; \
        }

# define Coordinates3_SetRow( self, i, x, y, z ) \
        { \
	    Coordinates3_Item ( self, i, 0 ) = x ; \
	    Coordinates3_Item ( self, i, 1 ) = y ; \
	    Coordinates3_Item ( self, i, 2 ) = z ; \
        }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Coordinates3    *Coordinates3_Allocate                                ( const Integer extent ) ;
extern Real             Coordinates3_Angle                                   ( const Coordinates3  *self, const Integer i, const Integer j, const Integer k ) ;
extern Status           Coordinates3_BuildPointFromDistance                  (       Coordinates3  *self, const Integer i, const Integer j, const Real r, const Vector3 *direction ) ;
extern Status           Coordinates3_BuildPointFromDistanceAngle             (       Coordinates3  *self, const Integer i, const Integer j, const Integer k, const Real r, const Real theta, const Vector3 *direction ) ;
extern Status           Coordinates3_BuildPointFromDistanceAngleDihedral     (       Coordinates3  *self, const Integer i, const Integer j, const Integer k, const Integer l, const Real r, const Real theta, const Real phi ) ;
extern Status           Coordinates3_BuildPointFromDistancePlaneAngle        (       Coordinates3  *self, const Integer i, const Integer j, const Integer k, const Integer l, const Real r, const Real planeangle ) ;
extern Status           Coordinates3_BuildPointFromDistanceTetrahedralTripod (       Coordinates3  *self, const Integer i, const Integer j, const Integer k, const Integer l, const Integer m, const Real r ) ;
extern Status           Coordinates3_Center                                  ( const Coordinates3  *self, const Selection *selection, const Real1DArray *weights, Vector3 **center ) ;
extern Real             Coordinates3_Dihedral                                ( const Coordinates3  *self, const Integer i, const Integer j, const Integer k, const Integer l ) ;
extern Real             Coordinates3_Distance                                ( const Coordinates3  *self, const Integer i, const Integer j ) ;
extern void             Coordinates3_EnclosingOrthorhombicBox                ( const Coordinates3  *self, const Selection *selection, const Real1DArray *radii, Vector3 *origin, Vector3 *extents ) ;
extern Status           Coordinates3_FromRegularGrid                         (       Coordinates3 **self, const RegularGrid *grid, Selection *selection ) ;
extern void             Coordinates3_Gather                                  (       Coordinates3  *self, const Coordinates3 *other, const Selection *selection ) ;
extern void             Coordinates3_GatherAddScaledMatrix                   (       Coordinates3  *self, const Real alpha, const Coordinates3 *other, const Selection *selection ) ;
extern Status           Coordinates3_IdentifyOccupiedGridPoints              (       Coordinates3  *self, const RegularGrid *grid             ,
                                                                                                          const Real1DArray *radii            ,
                                                                                                          const Boolean      QMIDPOINTOVERLAP ,
                                                                                                                Selection  **occupied         ) ;
extern SymmetricMatrix *Coordinates3_InertiaMatrix                           ( const Coordinates3  *self, const Selection *selection, const Real1DArray *weights ) ;
extern void             Coordinates3_MakeConformingGrid                      ( const Coordinates3  *self, const Selection    *andSelection   ,
                                                                                                          const RegularGrid  *grid           ,       
                                                                                                          RegularGrid       **conformingGrid ,       
                                                                                                          Integer1DArray    **offSet         ,       
                                                                                                          Status             *status         ) ;
extern void             Coordinates3_MakeConformingGridAndOccupancy          ( const Coordinates3  *self, const Selection       *andSelection   ,                   
                                                                                                          const RegularGrid     *grid           ,                   
                                                                                                          RegularGrid          **conformingGrid ,                   
                                                                                                          RegularGridOccupancy **occupancy      ,                   
                                                                                                          Integer1DArray       **offSet         ,                   
                                                                                                          Status                *status         ) ;
extern RegularGrid     *Coordinates3_MakeGrid                                ( const Coordinates3  *self, const Selection *andSelection, const Real gridSize, Status *status ) ;
extern void             Coordinates3_MakeGridAndOccupancy                    ( const Coordinates3  *self, const Selection       *andSelection ,
                                                                                                          const Real             gridSize     ,
                                                                                                          RegularGrid          **grid         ,
                                                                                                          RegularGridOccupancy **occupancy    ,
                                                                                                          Status                *status       ) ;
extern void             Coordinates3_MakePeriodicGridAndOccupancy            ( const Coordinates3    *self      ,
                                                                               const Vector3         *boxSize   ,
                                                                               const Real             gridSize  ,
                                                                               RegularGrid          **grid      ,
                                                                               RegularGridOccupancy **occupancy ,
                                                                               Status                *status    ) ;
extern void             Coordinates3_MomentsOfInertia                        ( const Coordinates3  *self, const Selection *selection, const Real1DArray *weights, Vector3 *moments, Matrix33 *axes ) ;
extern Real             Coordinates3_RadiusOfGyration                        ( const Coordinates3  *self, const Selection *selection, const Real1DArray *weights ) ;
extern Real             Coordinates3_RMSDeviation                            ( const Coordinates3  *self, const Coordinates3 *other, const Selection *selection, const Real1DArray *weights ) ;
extern void             Coordinates3_Rotate                                  (       Coordinates3  *self, const Matrix33 *rotation, const Selection *selection ) ;
extern Real2DArray     *Coordinates3_RotationTranslationVectors              ( const Coordinates3  *self, Real1DArray *weights, const Boolean QRx, const Boolean QRy, const Boolean QRz,
                                                                                                      const Boolean QTx, const Boolean QTy, const Boolean QTz, const Integer dimension ) ;
extern void             Coordinates3_ScaleRows                               (       Coordinates3  *self, const Real1DArray *rowScalingFactors, Status *status ) ;
extern void             Coordinates3_Scatter                                 ( const Coordinates3  *self,       Coordinates3 *other, const Selection *selection ) ;
extern void             Coordinates3_ScatterAddScaledMatrix                  ( const Coordinates3  *self, const Real alpha, Coordinates3 *other, const Selection *selection ) ;
extern void             Coordinates3_Superimpose                             (       Coordinates3  *self, const Coordinates3 *other, const Selection *selection, const Real1DArray *weights,
                                                                                                                                                  Matrix33 *rotation, Vector3 *translation ) ;
extern void             Coordinates3_ToPrincipalAxes                         (       Coordinates3  *self, const Selection *selection, const Real1DArray *weights ) ;
extern void             Coordinates3_Transform                               (       Coordinates3  *self, const Transformation3 *transformation, const Selection *selection ) ;
extern void             Coordinates3_Translate                               (       Coordinates3  *self, const Vector3 *translation, const Selection *selection ) ;
extern void             Coordinates3_TranslateToCenter                       (       Coordinates3  *self, const Selection *selection, const Real1DArray *weights ) ;

# endif
