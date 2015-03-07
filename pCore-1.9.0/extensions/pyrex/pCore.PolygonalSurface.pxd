#-------------------------------------------------------------------------------
# . File      : pCore.PolygonalSurface.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions   cimport Boolean, Integer, Real
from pCore.Integer2DArray cimport CInteger2DArray, Integer2DArray, Integer2DArray_Deallocate, Integer2DArray_GetItem, Integer2DArray_SetItem
from pCore.Real1DArray    cimport CReal1DArray, Real1DArray
from pCore.Real2DArray    cimport CReal2DArray, Real2DArray, Real2DArray_Deallocate, Real2DArray_GetItem, Real2DArray_SetItem
from pCore.RealNDArray    cimport CRealNDArray, RealNDArray
from pCore.RegularGrid    cimport CRegularGrid, RegularGrid
from pCore.Status         cimport Status, Status_Continue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PolygonalSurface.h":

    # . The real array type.
    ctypedef struct CPolygonalSurface "PolygonalSurface":
        CInteger2DArray *polygons
        CReal2DArray    *normals
        CReal2DArray    *vertices

    # . Functions.
    cdef CPolygonalSurface *PolygonalSurface_Allocate         ( Integer ndimensions, Integer nvertices, Integer npolygons, Status *status )
    cdef CPolygonalSurface *PolygonalSurface_Clone            ( CPolygonalSurface  *self, Status *status )
    cdef void               PolygonalSurface_ConstructNormals ( CPolygonalSurface  *self )
    cdef void               PolygonalSurface_CopyTo           ( CPolygonalSurface  *self, CPolygonalSurface  *other, Status *status )
    cdef void               PolygonalSurface_Deallocate       ( CPolygonalSurface **self )
    cdef void               PolygonalSurface_InitializeArrays ( CPolygonalSurface  *self )
    cdef void               PolygonalSurface_NormalizeNormals ( CPolygonalSurface  *self )
    cdef Integer            PolygonalSurface_NumberOfPolygons ( CPolygonalSurface  *self )
    cdef Integer            PolygonalSurface_NumberOfVertices ( CPolygonalSurface  *self )
    cdef void               PolygonalSurface_OriginAndExtents ( CPolygonalSurface  *self, CReal1DArray *origin, CReal1DArray *extents, Status *status )
    cdef Integer            PolygonalSurface_Rank             ( CPolygonalSurface  *self )
    cdef void               PolygonalSurface_Resize           ( CPolygonalSurface  *self, Integer nvertices, Integer npolygons, Boolean initialize, Status *status )

cdef extern from "MarchingCubes.h":

    # . Functions.
    cdef CPolygonalSurface *MarchingCubes_GenerateIsosurface ( CRegularGrid *grid, CRealNDArray *data, Real isovalue, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PolygonalSurface:

    cdef CPolygonalSurface *cObject
    cdef public object      isOwner
    cdef public object      label

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PolygonalSurfacePolygonIndexIterator:

    cdef CPolygonalSurface *cObject
    cdef Integer            dimensions
    cdef Integer            length
    cdef Integer            position

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PolygonalSurfacePolygonVertexIterator:

    cdef CPolygonalSurface *cObject
    cdef Integer            dimensions
    cdef Integer            length
    cdef Integer            positiond
    cdef Integer            positionp

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PolygonalSurfaceVertexIterator:

    cdef CPolygonalSurface *cObject
    cdef Integer            dimensions
    cdef Integer            length
    cdef Integer            position
