/*------------------------------------------------------------------------------
! . File      : PolygonalSurface.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _POLYGONALSURFACE
# define _POLYGONALSURFACE

# include "Boolean.h"
# include "Integer.h"
# include "Integer2DArray.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Status.h"

/* . Simple surface type. */
/* . At the moment only the smallest polygons appropriate for the dimension are supported (i.e. lines for 2D, triangles for 3D). */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The surface type. */
typedef struct {
    Integer2DArray *polygons   ;
    Real2DArray    *normals    ;
    Real2DArray    *vertices   ;
} PolygonalSurface ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern PolygonalSurface *PolygonalSurface_Allocate         ( const Integer ndimensions, const Integer nvertices, const Integer npolygons, Status *status ) ;
extern PolygonalSurface *PolygonalSurface_Clone            ( const PolygonalSurface  *self, Status *status ) ;
extern void              PolygonalSurface_ConstructNormals (       PolygonalSurface  *self ) ;
extern void              PolygonalSurface_CopyTo           ( const PolygonalSurface  *self, PolygonalSurface  *other, Status *status ) ;
extern void              PolygonalSurface_Deallocate       (       PolygonalSurface **self ) ;
extern void              PolygonalSurface_InitializeArrays (       PolygonalSurface  *self ) ;
extern void              PolygonalSurface_NormalizeNormals (       PolygonalSurface  *self ) ;
extern Integer           PolygonalSurface_NumberOfPolygons ( const PolygonalSurface  *self ) ;
extern Integer           PolygonalSurface_NumberOfVertices ( const PolygonalSurface  *self ) ;
extern void              PolygonalSurface_OriginAndExtents ( const PolygonalSurface  *self, Real1DArray *origin, Real1DArray *extents, Status *status ) ;
extern Integer           PolygonalSurface_Rank             ( const PolygonalSurface  *self ) ;
extern void              PolygonalSurface_Resize           (       PolygonalSurface  *self, const Integer nvertices, const Integer npolygons, const Boolean initialize, Status *status ) ;

# endif
