#-------------------------------------------------------------------------------
# . File      : pCore.PolygonalSurface.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle polygonal surfaces."""

from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PolygonalSurface:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef PolygonalSurface new
        new = self.__class__.WithSizes ( self.rank, self.NumberOfVertices ( ), self.NumberOfPolygons ( ) )
        PolygonalSurface_CopyTo ( self.cObject, new.cObject, NULL )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            PolygonalSurface_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.PolygonalSurface"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer2DArray polygons
        cdef Real2DArray    normals, vertices
        state = { "rank" : self.rank }
        if self.label is not None: state["label"] = self.label
        if self.cObject.polygons != NULL:
            polygons = Integer2DArray.Raw ( )
            polygons.cObject = self.cObject.polygons
            polygons.isOwner = False
            state["polygons"] = polygons
        if self.cObject.normals != NULL:
            normals = Real2DArray.Raw ( )
            normals.cObject = self.cObject.normals
            normals.isOwner = False
            state["normals"] = normals
        if self.cObject.vertices != NULL:
            vertices = Real2DArray.Raw ( )
            vertices.cObject = self.cObject.vertices
            vertices.isOwner = False
            state["vertices"] = vertices
        return state

    def __init__ ( self, Integer rank, Integer polygons, Integer vertices, initialize = False ):
        """Constructor with sizes."""
        self._Initialize ( )
        self._Allocate ( rank, vertices, polygons )
        if initialize: PolygonalSurface_InitializeArrays ( self.cObject )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer2DArray polygons
        cdef Real2DArray    normals, vertices
        # . Basic data.
        ndimensions = state["rank"]
        label       = state.get ( "label", None )
        self.label  = label
        # . Allocate the object.
        self._Allocate ( ndimensions, -1, -1 )
        # . Fill the object.
        if "polygons" in state:
            polygons = state["polygons"]
            Integer2DArray_Deallocate ( &(self.cObject.polygons) )
            self.cObject.polygons = polygons.cObject
            polygons.isOwner      = False
        if "normals" in state:
            normals = state["normals"]
            Real2DArray_Deallocate ( &(self.cObject.normals) )
            self.cObject.normals = normals.cObject
            normals.isOwner      = False
        if "vertices" in state:
            vertices = state["vertices"]
            Real2DArray_Deallocate ( &(self.cObject.vertices) )
            self.cObject.vertices = vertices.cObject
            vertices.isOwner      = False

    def _Allocate ( self, Integer rank, Integer polygons, Integer vertices ):
        """Constructor."""
        cdef Status  status
        status       = Status_Continue
        self.cObject = PolygonalSurface_Allocate ( rank, vertices, polygons, &status )
        if status   != Status_Continue: raise ValueError ( "Object allocation error." )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.label   = None

    def NumberOfPolygons ( self ):
        """Return the number of polygons."""
        return PolygonalSurface_NumberOfPolygons ( self.cObject )

    def NumberOfVertices ( self ):
        """Return the number of vertices."""
        return PolygonalSurface_NumberOfVertices ( self.cObject )

    def OriginAndExtents ( self ):
        """Get the origin and extents of the surface."""
        cdef Real1DArray origin, extents
        origin  = Real1DArray.WithExtent ( self.rank )
        extents = Real1DArray.WithExtent ( self.rank )
        PolygonalSurface_OriginAndExtents ( self.cObject, origin.cObject, extents.cObject, NULL )
        return ( origin, extents )

    def PolygonIndexIterator ( self ):
        """Return an iterator over polygons giving the integer indices of the vertices."""
        return PolygonalSurfacePolygonIndexIterator ( self )

    def PolygonVertexIterator ( self ):
        """Return an iterator over vertices in polygon order."""
        return PolygonalSurfacePolygonVertexIterator ( self )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetNormal ( self, Integer i, data ):
        """Set a normal."""
        cdef Integer d
        cdef Real    value
        if len ( data ) == self.rank:
            if ( i >= 0 ) and ( i < self.NumberOfVertices ( ) ):
                for d from 0 <= d < self.rank:
                    value = data[d]
                    Real2DArray_SetItem ( self.cObject.normals, i, d, value, NULL )
            else: raise IndexError
        else: raise ValueError ( "Invalid data dimension." )

    def SetPolygon ( self, Integer i, data ):
        """Set a vertex."""
        cdef Integer d, value
        if len ( data ) == self.rank:
            if ( i >= 0 ) and ( i < self.NumberOfPolygons ( ) ):
                for d from 0 <= d < self.rank:
                    value = data[d]
                    Integer2DArray_SetItem ( self.cObject.polygons, i, d, value, NULL )
            else: raise IndexError
        else: raise ValueError ( "Invalid data dimension." )

    def SetVertex ( self, Integer i, data ):
        """Set a vertex."""
        cdef Integer d
        cdef Real    value
        if len ( data ) == self.rank:
            if ( i >= 0 ) and ( i < self.NumberOfVertices ( ) ):
                for d from 0 <= d < self.rank:
                    value = data[d]
                    Real2DArray_SetItem ( self.cObject.vertices, i, d, value, NULL )
            else: raise IndexError
        else: raise ValueError ( "Invalid data dimension." )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "Polygonal Surface Summary" )
            summary.Entry ( "Polygons", "{:d}".format ( self.NumberOfPolygons ( ) ) )
            summary.Entry ( "Rank"    , "{:d}".format ( self.rank                 ) )
            summary.Entry ( "Vertices", "{:d}".format ( self.NumberOfVertices ( ) ) )
            summary.Stop ( )

    def VertexIterator ( self ):
        """Return an iterator over vertices in vertex order."""
        return PolygonalSurfaceVertexIterator ( self )

    @classmethod
    def WithSizes ( selfClass, Integer rank, Integer polygons, Integer vertices, initialize = False ):
        """Constructor with sizes."""
        return selfClass ( rank, vertices, polygons, initialize = initialize )

    # . Properties.
    property rank:
        def __get__ ( self ): return PolygonalSurface_Rank ( self.cObject )

#===================================================================================================================================
# . Isosurface constructors.
#===================================================================================================================================
def Isosurface_FromRegularGridData ( RegularGrid grid, RealNDArray data, Real isovalue ):
    """Generate an isosurface from regular grid data."""
    cdef PolygonalSurface self
    cdef Status           status
    self = None
    if ( grid is not None ) and ( data is not None ):
        self         = PolygonalSurface.Raw ( )
        status       = Status_Continue
        self.cObject = MarchingCubes_GenerateIsosurface ( grid.cObject, data.cObject, isovalue, &status )
        self.isOwner = True
        if status   != Status_Continue: raise ValueError ( "Isosurface generation error." )
    return self

#===================================================================================================================================
# . Iterators.
#===================================================================================================================================
cdef class PolygonalSurfacePolygonIndexIterator:
    """An iterator over polygons defined by their indices."""

    def __init__ ( self, PolygonalSurface target ):
        """Constructor."""
        self.cObject    = NULL
        self.dimensions = 0
        self.length     = 0
        self.position   = 0
        if target is not None:
            self.dimensions = target.rank
            self.length     = target.NumberOfPolygons ( )
            if target.cObject != NULL: self.cObject = target.cObject

    def __iter__ ( self ): return self

    def __next__ ( self ):
        """Get the next item."""
        cdef Integer d, p
        if self.position >= self.length: raise StopIteration
        else:
            p     = self.position
            items = []
            for d from 0 <= d < self.dimensions: items.append ( Integer2DArray_GetItem ( self.cObject.polygons, p, d, NULL ) )
            self.position += 1
            return tuple ( items )

cdef class PolygonalSurfacePolygonVertexIterator:
    """An iterator over vertices in the order that they occur in the polygons."""

    def __init__ ( self, PolygonalSurface target ):
        """Constructor."""
        self.cObject    = NULL
        self.dimensions = 0
        self.length     = 0
        self.positiond  = 0
        self.positionp  = 0
        if target is not None:
            self.dimensions = target.rank
            self.length     = target.NumberOfPolygons ( ) * self.dimensions
            if target.cObject != NULL: self.cObject = target.cObject

    def __iter__ ( self ): return self

    def __next__ ( self ):
        """Get the next item."""
        cdef Integer d, p, v
        if ( self.positionp * self.dimensions + self.positiond ) >= self.length: raise StopIteration
        else:
            d     = self.positiond
            p     = self.positionp
            v     = Integer2DArray_GetItem ( self.cObject.polygons, p, d, NULL )
            items = []
            for d from 0 <= d < self.dimensions: items.append ( Real2DArray_GetItem ( self.cObject.vertices, v, d, NULL ) )
            for d from 0 <= d < self.dimensions: items.append ( Real2DArray_GetItem ( self.cObject.normals , v, d, NULL ) )
            self.positiond += 1
            if self.positiond >= self.dimensions:
                self.positiond  = 0
                self.positionp += 1
            return tuple ( items )

cdef class PolygonalSurfaceVertexIterator:
    """An iterator over vertices and their normals."""

    def __init__ ( self, PolygonalSurface target ):
        """Constructor."""
        self.cObject    = NULL
        self.dimensions = 0
        self.length     = 0
        self.position   = 0
        if target is not None:
            self.dimensions = target.rank
            self.length     = target.NumberOfVertices ( )
            if target.cObject != NULL: self.cObject = target.cObject

    def __iter__ ( self ): return self

    def __next__ ( self ):
        """Get the next item."""
        cdef Integer d, p
        if self.position >= self.length: raise StopIteration
        else:
            p     = self.position
            items = []
            for d from 0 <= d < self.dimensions: items.append ( Real2DArray_GetItem ( self.cObject.vertices, p, d, NULL ) )
            for d from 0 <= d < self.dimensions: items.append ( Real2DArray_GetItem ( self.cObject.normals , p, d, NULL ) )
            self.position += 1
            return tuple ( items )
