#-------------------------------------------------------------------------------
# . File      : pCore.Coordinates3.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle 3-D coordinate arrays."""

from CoreObjects   import CLibraryError
from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Coordinates3:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        clone = None
        if self.isOwner:
            clone = self.__deepcopy__ ( None )
        else:
            clone = self.__class__ ( 0 )
            self.owner.MakeCoordinates3View ( clone, self.cObject.offset, self.cObject.length0, self.cObject.stride0, self.cObject.length1, self.cObject.stride1 )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Coordinates3_Deallocate ( &self.cObject )
            self.isOwner = False
        self.owner = None

    def __deepcopy__ ( self, memo ):
        """Copying."""
        cdef Coordinates3 new
        new = self.__class__.WithExtent ( self.rows )
        self.CopyTo ( new )
        return new

    def __getitem__ ( self, indices ):
        """Get an item."""
        cdef Integer i, j
        cdef Real    value
        cdef Status  status
        view = self._CheckSliceArguments ( indices )
        if view is None:
            i          = indices[0]
            j          = indices[1]
            status     = Status_Continue
            value      = Coordinates3_GetItem ( self.cObject, i, j, &status )
            if status != Status_Continue: raise IndexError ( "Indices ({:d},{:d}) out of range.".format ( i, j ) )
            return value
        else:
            return view

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.Coordinates3"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i, j
        cdef Real    value
        extent0 = self.rows
        extent1 = self.columns
        state   = { "shape" : [ extent0, extent1 ], "storage" : "RowMajor" }
        if self.isOwner:
            items = []
            for i from 0 <= i < extent0:
                for j from 0 <= j < extent1:
                    value = Coordinates3_GetItem ( self.cObject, i, j, NULL )
                    items.append ( value )
            state["items"] = items
        else:
            state.update ( { "offset"  : self.cObject.offset ,
                             "owner"   : self.owner          ,
                             "strides" : [ self.cObject.stride0, self.cObject.stride1 ] } )
        if self.numberUndefined > 0:
            indices = list ( self.undefined )
            indices.sort ( )
            state["undefined"] = indices
        return state

    def __init__ ( self, extent ):
        """Constructor with extent."""
        self._Initialize ( )
        self._Allocate ( extent )

    def __len__ ( self ):
        """Return the size of the array."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, indices, Real value ):
        """Set an item."""
        cdef Integer i, j
        cdef Status  status
        view = self._CheckSliceArguments ( indices )
        if view is None:
            i          = indices[0]
            j          = indices[1]
            status     = Status_Continue
            Coordinates3_SetItem ( self.cObject, i, j, value, &status )
            if status != Status_Continue: raise IndexError ( "Indices ({:d},{:d}) out of range.".format ( i, j ) )
        else:
            view.Set ( value )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer i, j, extent0, extent1, stride0, stride1
        cdef Real    value
        cdef Status  status
        ( extent0, extent1 ) = state["shape"]
        items                = state.get ( "items", None )
        if items is not None:
            self._Allocate ( extent0 )
            n = 0
            for i from 0 <= i < extent0:
                for j from 0 <= j < extent1:
                    value = items[n]
                    Coordinates3_SetItem ( self.cObject, i, j, value, NULL )
                    n += 1
        else:
            offset               = state["offset" ]
            owner                = state["owner"  ]
            ( stride0, stride1 ) = state["strides"]
            self._Allocate ( 0 )
            owner.MakeCoordinates3View ( self, offset, extent0, stride0, extent1, stride1 )
        if "undefined" in state:
            self.undefined = set ( state["undefined"] )
        else:
            self.undefined = None

    def _Allocate ( self, Integer extent ):
        """Allocation."""
        if extent >= 0:
            self.cObject = Coordinates3_Allocate ( extent )
            self.isOwner = True
            self.owner   = None
        if ( extent < 0 ) or ( self.cObject == NULL ): raise CLibraryError ( "Memory allocation failure." )

    def _CheckSliceArguments ( self, indices ):
        """Check slice arguments."""
        d    = 0
        isOK = True
        view = None
        if   isinstance ( indices, int   ): view = self._SliceVector3      ( indices )
        elif isinstance ( indices, slice ): view = self._SliceCoordinates3 ( indices )
        elif isinstance ( indices, tuple ) and ( len ( indices ) == 2 ):
            for index in indices:
                if   isinstance ( index, int   ): pass
                elif isinstance ( index, slice ): d+= 1
                else:
                    isOK = False
                    break
        else: isOK = False
        if not isOK: raise TypeError ( "Expecting an integer, a slice or a two-tuple of integers and slices as indices." )
        if   d == 1: view = self._Slice1D ( indices[0], indices[1] )
        elif d == 2: view = self._Slice2D ( indices[0], indices[1] )
        return view

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.owner   = None

    def _Slice1D ( self, i, j ):
        """Create a 1-D slice."""
        cdef Integer start0, start1, stop0, stop1, stride0, stride1
        cdef Real1DArray view
        cdef Status      status
        # . Indices.
        if isinstance ( i, int ):
            start0 = i ; stop0 = i+1 ; stride0 = 1
        else:
            ( start0, stop0, stride0 ) = i.indices ( self.rows )
        if isinstance ( j, int ):
            start1 = j ; stop1 = j+1 ; stride1 = 1
        else:
            ( start1, stop1, stride1 ) = j.indices ( self.columns )
        # . View.
        view = Real1DArray ( 0 )
        view.isOwner = False
        view.owner   = self
        status       = Status_Continue
        Coordinates3_1DSlice ( self.cObject, start0, stop0, stride0, start1, stop1, stride1, view.cObject, &status )
        if status   != Status_Continue: raise CLibraryError ( "Invalid 1-D slice." )
        return view

    def _Slice2D ( self, i, j ):
        """Create a 2-D slice."""
        cdef Integer start0, start1, stop0, stop1, stride0, stride1
        cdef Real2DArray view
        cdef Status      status
        # . Indices.
        ( start0, stop0, stride0 ) = i.indices ( self.rows    )
        ( start1, stop1, stride1 ) = j.indices ( self.columns )
        # . View.
        view = Real2DArray ( 0, 0 )
        view.isOwner = False
        view.owner   = self
        status       = Status_Continue
        Coordinates3_2DSlice ( self.cObject, start0, stop0, stride0, start1, stop1, stride1, view.cObject, &status )
        if status   != Status_Continue: raise CLibraryError ( "Invalid 2-D slice." )
        return view

    def _SliceCoordinates3 ( self, i ):
        """Create a Coordinates3 slice."""
        cdef Coordinates3 view
        cdef Integer      start, stop, stride
        cdef Status       status
        ( start, stop, stride ) = i.indices ( self.rows )
        view         = self.__class__ ( 0 )
        view.isOwner = False
        view.owner   = self
        status       = Status_Continue
        Coordinates3_SliceCoordinates3 ( self.cObject, start, stop, stride, 0, 3, 1, view.cObject, &status )
        if status   != Status_Continue: raise CLibraryError ( "Invalid Coordinates3 slice." )
        return view

    def _SliceVector3 ( self, Integer i ):
        """Create a Vector3 slice."""
        cdef Vector3 view
        cdef Status  status
        view = Vector3 ( )
        view.isOwner = False
        view.owner   = self
        status       = Status_Continue
        Coordinates3_SliceVector3 ( self.cObject, i, i+1, 1, 0, 3, 1, view.cObject, &status )
        if status   != Status_Continue: raise CLibraryError ( "Invalid Vector3 slice." )
        return view

    def AbsoluteMaximum ( self ):
        """Return the maximum absolute value in the array."""
        return Coordinates3_AbsoluteMaximum ( self.cObject )

    def AddScaledMatrix ( self, Real alpha, Coordinates3 other ):
        """Add a scaled matrix."""
        Coordinates3_AddScaledArray ( self.cObject, alpha, other.cObject, NULL )

    def AddScaledVector3 ( self, Integer i, Real alpha, Vector3 v ):
        """Increment a row by a scaled vector."""
        self[i,0] = self[i,0] + alpha * v[0]
        self[i,1] = self[i,1] + alpha * v[1]
        self[i,2] = self[i,2] + alpha * v[2]

    def Angle ( self, Integer i, Integer j, Integer k ):
        """Calculate the angle between three points."""
        return Coordinates3_Angle ( self.cObject, i, j, k )

    def BuildPointFromDistance ( self, Integer i, Integer j, Real r, Vector3 direction ):
        """Build point |i| from point |j|, a distance |r| and a direction |direction|."""
        cdef Status status
        status = Coordinates3_BuildPointFromDistance ( self.cObject, i, j, r, direction.cObject )
        if status == Status_Success:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def BuildPointFromDistanceAngle ( self, Integer i, Integer j, Integer k, Real r, Real theta, Vector3 direction ):
        """Build point |i| from points |j| and |k|, a distance |r|, an angle |theta| and a direction |direction|."""
        cdef Status status
        status = Coordinates3_BuildPointFromDistanceAngle ( self.cObject, i, j, k, r, theta, direction.cObject )
        if status == Status_Success:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def BuildPointFromDistanceAngleDihedral ( self, Integer i, Integer j, Integer k, Integer l, Real r, Real theta, Real phi ):
        """Build point |i| from points |j|, |k| and |l|, a distance |r|, an angle |theta| and a dihedral |phi|."""
        cdef Status status
        status = Coordinates3_BuildPointFromDistanceAngleDihedral ( self.cObject, i, j, k, l, r, theta, phi )
        if status == Status_Success:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def BuildPointFromDistancePlaneAngle ( self, Integer i, Integer j, Integer k, Integer l, Real r, Real planeangle ):
        """Build point |i| from points |j|, |k| and |l|, a distance |r| and a plane angle |planeangle|."""
        cdef Status status
        status = Coordinates3_BuildPointFromDistancePlaneAngle ( self.cObject, i, j, k, l, r, planeangle )
        if status == Status_Success:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def BuildPointFromDistanceTetrahedralTripod ( self, Integer i, Integer j, Integer k, Integer l, Integer m, Real r ):
        """Build point |i| from points |j|, |k|, |l| and |m| and a distance |r| using a tetrahedral tripod."""
        cdef Status status
        status = Coordinates3_BuildPointFromDistanceTetrahedralTripod ( self.cObject, i, j, k, l, m, r )
        if status == Status_Success:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def Center ( self, Selection selection = None, Real1DArray weights = None ):
        """Determine the center of the matrix."""
        cdef Vector3 center
        center         = Vector3.Raw ( )
        center.isOwner = True
        if selection is None:
            if weights is None: Coordinates3_Center ( self.cObject, NULL,              NULL,            &(center.cObject) )
            else:               Coordinates3_Center ( self.cObject, NULL,              weights.cObject, &(center.cObject) )
        else:
            if weights is None: Coordinates3_Center ( self.cObject, selection.cObject, NULL,            &(center.cObject) )
            else:               Coordinates3_Center ( self.cObject, selection.cObject, weights.cObject, &(center.cObject) )
        return center

    def CopyFromArray ( self, Real1DArray vector, Selection selection = None ):
        """Fill the matrix from a vector."""
        if selection is None: Coordinates3_CopyFromArray ( self.cObject, vector.cObject, NULL              )
        else:                 Coordinates3_CopyFromArray ( self.cObject, vector.cObject, selection.cObject )

    def CopyTo ( self, Coordinates3 other ):
        """Copying."""
        Coordinates3_CopyTo ( self.cObject, other.cObject, NULL )

    def CopyToArray ( self, Real1DArray vector, Selection selection = None ):
        """Empty the matrix to a vector."""
        if selection is None: Coordinates3_CopyToArray ( self.cObject, vector.cObject, NULL              )
        else:                 Coordinates3_CopyToArray ( self.cObject, vector.cObject, selection.cObject )

    def Decrement ( self, Integer i, Vector3 v ):
        """Decrement a row of the matrix."""
        self[i,0] = self[i,0] - v[0]
        self[i,1] = self[i,1] - v[1]
        self[i,2] = self[i,2] - v[2]

    def Dihedral ( self, Integer i, Integer j, Integer k, Integer l ):
        """Calculate the dihedral between four points."""
        return Coordinates3_Dihedral ( self.cObject, i, j, k, l )

    def Displacement ( self, Integer i, Integer j, dr = None ):
        """Displacement between two points i and j."""
        cdef Vector3 vij
        if dr is None: vij = Vector3.Uninitialized ( )
        else:          vij = dr
        vij[0] = self[i,0] - self[j,0]
        vij[1] = self[i,1] - self[j,1]
        vij[2] = self[i,2] - self[j,2]
        return vij

    def Distance ( self, Integer i, Integer j ):
        """Calculate the distance between two points."""
        return Coordinates3_Distance ( self.cObject, i, j )

    def EnclosingOrthorhombicBox ( self, Selection selection = None, Real1DArray radii = None ):
        """Find the enclosing box around the matrix."""
        cdef Vector3 extents, origin
        extents = Vector3.Uninitialized ( )
        origin  = Vector3.Uninitialized ( )
        if selection is None:
            if radii is None: Coordinates3_EnclosingOrthorhombicBox ( self.cObject, NULL, NULL,          origin.cObject, extents.cObject )
            else:             Coordinates3_EnclosingOrthorhombicBox ( self.cObject, NULL, radii.cObject, origin.cObject, extents.cObject )
        else:
            if radii is None: Coordinates3_EnclosingOrthorhombicBox ( self.cObject, selection.cObject, NULL,          origin.cObject, extents.cObject )
            else:             Coordinates3_EnclosingOrthorhombicBox ( self.cObject, selection.cObject, radii.cObject, origin.cObject, extents.cObject )
        return ( origin, extents )

    def FlagCoordinateAsDefined ( self, i ):
        """Flag a coordinate as being defined."""
        if self.undefined is not None: self.undefined.discard ( i )

    def FlagCoordinateAsUndefined ( self, i ):
        """Flag a coordinate as being undefined."""
        if self.undefined is None: self.undefined = set ( )
        self.undefined.add ( i )

    def Gather ( self, Coordinates3 other, Selection selection = None ):
        """Gather items from a sparse other, indexed by selection, to a compact self."""
        if selection is None: Coordinates3_CopyTo ( other.cObject, self.cObject , NULL              )
        else:                 Coordinates3_Gather ( self.cObject , other.cObject, selection.cObject )

    def GatherAddScaledMatrix ( self, Real alpha, Coordinates3 other, Selection selection = None ):
        """Gather items, with adding and scaling, from a sparse other, indexed by selection, to a compact self."""
        if selection is None: Coordinates3_AddScaledArray        ( self.cObject, alpha, other.cObject, NULL              )
        else:                 Coordinates3_GatherAddScaledMatrix ( self.cObject, alpha, other.cObject, selection.cObject )

    def GetRow ( self, row ):
        """Get a row as a vector."""
        cdef Real x, y, z
        cdef Integer    i
        try:    i = int ( row )
        except: raise TypeError ( "Invalid row index type." )
        if ( i < 0 ) or ( i >= self.rows ): raise IndexError ( "Row index - {:d} - out of range.".format ( i ) )
        else:
            x = Coordinates3_GetItem ( self.cObject, i, 0, NULL )
            y = Coordinates3_GetItem ( self.cObject, i, 1, NULL )
            z = Coordinates3_GetItem ( self.cObject, i, 2, NULL )
            return Vector3.WithValues ( x, y, z )

    def IdentifyOccupiedGridPoints ( self, RegularGrid grid, Real1DArray radii ):
        """Identify occupied grid points."""
        cdef Selection occupied
        occupied         = Selection.Raw ( )
        occupied.isOwner = True
        Coordinates3_IdentifyOccupiedGridPoints ( self.cObject, grid.cObject, radii.cObject, CFalse, &occupied.cObject )
        return occupied

    def IdentifyUnoccupiedGridPoints ( self, RegularGrid grid, Real1DArray radii ):
        """Identify unoccupied grid points."""
        occupied = self.IdentifyOccupiedGridPoints ( grid, radii )
        return occupied.Complement ( upperBound = grid.NumberOfGridPoints ( ) )

    def Increment ( self, Integer i, Vector3 v ):
        """Increment a row of the matrix."""
        self[i,0] = self[i,0] + v[0]
        self[i,1] = self[i,1] + v[1]
        self[i,2] = self[i,2] + v[2]

    def InertiaMatrix ( self, Selection selection = None, Real1DArray weights = None  ):
        """Determine the inertia matrix."""
        cdef SymmetricMatrix inertia
        inertia         = SymmetricMatrix.Raw ( )
        inertia.isOwner = True
        if selection is None:
            if weights is None: inertia.cObject = Coordinates3_InertiaMatrix ( self.cObject, NULL, NULL            )
            else:               inertia.cObject = Coordinates3_InertiaMatrix ( self.cObject, NULL, weights.cObject )
        else:
            if weights is None: inertia.cObject = Coordinates3_InertiaMatrix ( self.cObject, selection.cObject, NULL            )
            else:               inertia.cObject = Coordinates3_InertiaMatrix ( self.cObject, selection.cObject, weights.cObject )
        return inertia

    def MakeCoordinates3View ( self, Coordinates3 view, Integer offset, Integer extent0, Integer stride0, Integer extent1, Integer stride1 ):
        """Make a Coordinates3 view."""
        cdef Status status
        status = Status_Continue
        Coordinates3_ViewOfRaw ( view.cObject, offset, extent0, stride0, extent1, stride1, self.cObject.data, self.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        view.isOwner = False
        view.owner   = self

    def MakeVector3View ( self, Vector3 view, Integer offset, Integer extent, Integer stride ):
        """Make a Vector3 view."""
        cdef Status status
        status = Status_Continue
        Vector3_ViewOfRaw ( view.cObject, offset, extent, stride, self.cObject.data, self.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        view.isOwner = False
        view.owner   = self

    def Make1DView ( self, Real1DArray view, Integer offset, Integer extent, Integer stride ):
        """Make a 1D view."""
        cdef Status status
        status = Status_Continue
        Real1DArray_ViewOfRaw ( view.cObject, offset, extent, stride, self.cObject.data, self.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        view.isOwner = False
        view.owner   = self

    def Make2DView ( self, Real2DArray view, Integer offset, Integer extent0, Integer stride0, Integer extent1, Integer stride1 ):
        """Make a 2D view."""
        cdef Status status
        status = Status_Continue
        Real2DArray_ViewOfRaw ( view.cObject, offset, extent0, stride0, extent1, stride1, self.cObject.data, self.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        view.isOwner = False
        view.owner   = self

    def Merge ( self, items, information = {} ):
        """Merging."""
        cdef Real   value
        cdef Integer      i, j, n, n0
        cdef Coordinates3 item, new
        # . Find the size of the merged item.
        n0 = 0
        for item in [ self ] + items:
            n0 = n0 + item.rows
        # . Allocate the object.
        new = self.__class__.WithExtent ( n0 )
        # . Fill the object and assign references.
        n0 = 0
        for item in [ self ] + items:
            n = item.rows
            for i from 0 <= i < n:
                for j from 0 <= j < self.columns:
                    value = Coordinates3_GetItem ( item.cObject, i, j, NULL )
                    Coordinates3_SetItem ( new.cObject, i + n0, j, value, NULL )
            n0 = n0 + n
        return new

    def MomentsOfInertia ( self, Selection selection = None, Real1DArray weights = None  ):
        """Calculate the moments of inertia."""
        cdef Matrix33 axes
        cdef Vector3  moments
        axes    = Matrix33.Null ( )
        moments = Vector3.Null  ( )
        if selection is None:
            if weights is None: Coordinates3_MomentsOfInertia ( self.cObject, NULL             , NULL           , moments.cObject, axes.cObject )
            else:               Coordinates3_MomentsOfInertia ( self.cObject, NULL             , weights.cObject, moments.cObject, axes.cObject )
        else:
            if weights is None: Coordinates3_MomentsOfInertia ( self.cObject, selection.cObject, NULL           , moments.cObject, axes.cObject )
            else:               Coordinates3_MomentsOfInertia ( self.cObject, selection.cObject, weights.cObject, moments.cObject, axes.cObject )
        return ( moments, axes )

    def Print ( self, indexWidth =  6, itemFormat = "{:18.8f}", itemWidth  = 19, log = logFile, title = None ):
        """Printing."""
        cdef Real value
        cdef Integer    i, irow
        if LogFileActive ( log ):
            table = log.GetTable ( columns = [ indexWidth, itemWidth, itemWidth, itemWidth ] )
            table.Start ( )
            if title is not None: table.Title ( title )
            table.Heading ( None )
            table.Heading ( "0" )
            table.Heading ( "1" )
            table.Heading ( "2" )
            for irow in range ( self.rows ):
                table.Entry ( "{:d}".format ( irow ) )
                for i in range ( self.columns ):
                    value = Coordinates3_GetItem ( self.cObject, irow, i, NULL )
                    table.Entry ( itemFormat.format ( value ) )
            table.Stop ( )

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef Coordinates3 new
        new         = self.__class__.Raw ( )
        new.cObject = Coordinates3_Prune ( self.cObject, selection.cObject )
        new.isOwner = True
        return new

    def RadiusOfGyration ( self, Selection selection = None, Real1DArray weights = None ):
        """Determine the radius of gyration."""
        cdef Real value
        if selection is None:
            if weights is None: value = Coordinates3_RadiusOfGyration ( self.cObject, NULL, NULL            )
            else:               value = Coordinates3_RadiusOfGyration ( self.cObject, NULL, weights.cObject )
        else:
            if weights is None: value = Coordinates3_RadiusOfGyration ( self.cObject, selection.cObject, NULL            )
            else:               value = Coordinates3_RadiusOfGyration ( self.cObject, selection.cObject, weights.cObject )
        return value

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def RMSDeviation ( self, Coordinates3 other, Selection selection = None, Real1DArray weights = None ):
        """Determine the RMS deviation between two matrices."""
        if selection is None:
            if weights is None: rms = Coordinates3_RMSDeviation ( self.cObject, other.cObject, NULL, NULL            )
            else:               rms = Coordinates3_RMSDeviation ( self.cObject, other.cObject, NULL, weights.cObject )
        else:
            if weights is None: rms = Coordinates3_RMSDeviation ( self.cObject, other.cObject, selection.cObject, NULL            )
            else:               rms = Coordinates3_RMSDeviation ( self.cObject, other.cObject, selection.cObject, weights.cObject )
        return rms

    def RMSValue ( self ):
        """Determine the RMS value of the matrix elements."""
        return Coordinates3_RootMeanSquare ( self.cObject )

    def Rotate ( self, Matrix33 rotation, Selection selection = None ):
        """Rotation."""
        if selection is None: Coordinates3_Rotate ( self.cObject, rotation.cObject, NULL              )
        else:                 Coordinates3_Rotate ( self.cObject, rotation.cObject, selection.cObject )

    def RotationTranslationVectors ( self, qrotation, qtranslation, Integer dimension = 0, Real1DArray weights = None ):
        """Construct a set of rotation and translation vectors."""
        cdef Boolean     qr, qt
        cdef Real1DArray reference, scalars
        cdef Real2DArray vectors
        qr = CFalse
        qt = CFalse
        if qrotation:    qr = CTrue
        if qtranslation: qt = CTrue
        vectors         = Real2DArray.Raw ( )
        vectors.isOwner = True
        if weights is None: vectors.cObject = Coordinates3_RotationTranslationVectors ( self.cObject, NULL,            qr, qr, qr, qt, qt, qt, dimension )
        else:               vectors.cObject = Coordinates3_RotationTranslationVectors ( self.cObject, weights.cObject, qr, qr, qr, qt, qt, qt, dimension )
        reference = Real1DArray.WithExtent ( vectors.rows    )
        scalars   = Real1DArray.WithExtent ( vectors.columns )
        reference.Set ( 0.0 )
        self.CopyToArray ( reference )
        Real2DArray_VectorMultiply ( CTrue, 1.0, vectors.cObject, reference.cObject, 0.0, scalars.cObject, NULL )
        return ( vectors, scalars )

    def Scale ( self, Real value ):
        """Scaling."""
        Coordinates3_Scale ( self.cObject, value )

    def Scatter ( self, Coordinates3 other, Selection selection = None ):
        """Scatter items from a compact self to a sparse other, indexed by selection."""
        if selection is None: Coordinates3_CopyTo  ( self.cObject, other.cObject, NULL              )
        else:                 Coordinates3_Scatter ( self.cObject, other.cObject, selection.cObject )

    def ScatterAddScaledMatrix ( self, Real alpha, Coordinates3 other, Selection selection = None ):
        """Scatter items, with adding and scaling, from a compact self to a sparse other, indexed by selection."""
        if selection is None: Coordinates3_AddScaledArray         ( other.cObject, alpha, self.cObject , NULL              )
        else:                 Coordinates3_ScatterAddScaledMatrix ( self.cObject , alpha, other.cObject, selection.cObject )

    def Set ( self, Real value ):
        """Set all the elements of a matrix."""
        Coordinates3_Set ( self.cObject, value )

    def SetRowSelection ( self, Selection selection, Real value ):
        """Set the values of selected rows."""
        Coordinates3_SetByRow ( self.cObject, selection.cObject, value )

    def Superimpose ( self, Coordinates3 other, Selection selection = None, Real1DArray weights = None, returnTransformation = False ):
        """Superimpose the matrix onto a other matrix."""
        cdef CMatrix33 *cRotation
        cdef CVector3  *cTranslation
        cdef Matrix33   rotation
        cdef Vector3    translation
        if returnTransformation:
            rotation     = Matrix33.Null ( )
            translation  = Vector3.Null  ( )
            cRotation    = rotation.cObject
            cTranslation = translation.cObject
        else:
            cRotation    = NULL
            cTranslation = NULL
        if selection is None:
            if weights is None: Coordinates3_Superimpose ( self.cObject, other.cObject, NULL, NULL           , cRotation, cTranslation )
            else:               Coordinates3_Superimpose ( self.cObject, other.cObject, NULL, weights.cObject, cRotation, cTranslation )
        else:
            if weights is None: Coordinates3_Superimpose ( self.cObject, other.cObject, selection.cObject, NULL           , cRotation, cTranslation )
            else:               Coordinates3_Superimpose ( self.cObject, other.cObject, selection.cObject, weights.cObject, cRotation, cTranslation )
        if returnTransformation: return ( rotation, translation )

    def ToPrincipalAxes ( self, Selection selection = None, Real1DArray weights = None ):
        """Do a principal axis transformation."""
        if selection is None:
            if weights is None: Coordinates3_ToPrincipalAxes ( self.cObject, NULL, NULL            )
            else:               Coordinates3_ToPrincipalAxes ( self.cObject, NULL, weights.cObject )
        else:
            if weights is None: Coordinates3_ToPrincipalAxes ( self.cObject, selection.cObject, NULL            )
            else:               Coordinates3_ToPrincipalAxes ( self.cObject, selection.cObject, weights.cObject )

    def Transform ( self, Transformation3 transformation3, Selection selection = None ):
        """Transformation."""
        if selection is None: Coordinates3_Transform ( self.cObject, transformation3.cObject, NULL              )
        else:                 Coordinates3_Transform ( self.cObject, transformation3.cObject, selection.cObject )

    def Translate ( self, Vector3 translation, Selection selection = None ):
        """Translation."""
        if selection is None: Coordinates3_Translate ( self.cObject, translation.cObject, NULL              )
        else:                 Coordinates3_Translate ( self.cObject, translation.cObject, selection.cObject )

    def TranslateToCenter ( self, Selection selection = None, Real1DArray weights = None ):
        """Translate the matrix to its center."""
        if selection is None:
            if weights is None: Coordinates3_TranslateToCenter ( self.cObject, NULL, NULL            )
            else:               Coordinates3_TranslateToCenter ( self.cObject, NULL, weights.cObject )
        else:
            if weights is None: Coordinates3_TranslateToCenter ( self.cObject, selection.cObject, NULL            )
            else:               Coordinates3_TranslateToCenter ( self.cObject, selection.cObject, weights.cObject )

    @classmethod
    def WithExtent ( selfClass, extent ):
        """Constructor with extent."""
        return selfClass ( extent )

    # . Properties.
    property columns:
        def __get__ ( self ): return Coordinates3_Length ( self.cObject, 1 )

    property numberUndefined:
        def __get__ ( self ):
            if self.undefined is None: return 0
            else:                      return len ( self.undefined )

    property rank:
        def __get__ ( self ): return 2

    property rows:
        def __get__ ( self ): return Coordinates3_Length ( self.cObject, 0 )

    property shape:
        def __get__ ( self ):
            return [ self.rows, self.columns ]

    property size:
        def __get__ ( self ): return Coordinates3_Length ( self.cObject, -1 )

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def Coordinates3_FromGrid ( RegularGrid grid, Selection selection = None ):
    """Constructor from grid."""
    cdef Coordinates3 self
    self         = Coordinates3.Raw ( )
    self.isOwner = True
    if selection is None: CCoordinates3_FromRegularGrid ( &self.cObject, grid.cObject, NULL              )
    else:                 CCoordinates3_FromRegularGrid ( &self.cObject, grid.cObject, selection.cObject )
    return self
