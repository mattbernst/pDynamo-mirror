#-------------------------------------------------------------------------------
# . File      : pCore.RealNDArray.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle real N-D arrays."""

from CoreObjects   import CLibraryError
from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealNDArray:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        cdef CArrayView *cView
        cdef RealNDArray clone
        clone = None
        if self.isOwner:
            clone = self.__deepcopy__ ( None )
        else:
            clone = self.__class__.Raw ( )
            clone._Allocate ( self.shape, doCObject = False )
            clone.isOwner = False
            clone.owner   = self.owner
            cView         = ArrayView_Clone ( self.cObject.view, NULL )
            if   self.owner.rank == 1: clone._ViewOf1DArray ( clone.owner, cView )
            elif self.owner.rank == 2: clone._ViewOf2DArray ( clone.owner, cView )
            else:                      clone._ViewOfNDArray ( clone.owner, cView )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        Memory_Deallocate_Integer ( &self.indices )
        if self.isOwner:
            RealNDArray_Deallocate ( &self.cObject )
            self.isOwner = False
        self.owner = None

    def __deepcopy__ ( self, memo ):
        """Copying."""
        cdef RealNDArray new
        new = self.__class__.WithExtents ( *self.shape )
        self.CopyTo ( new )
        return new

    def __getitem__ ( self, indices ):
        """Set an item."""
        cdef Integer d
        cdef Real    value
        cdef Status  status
        view = self._CheckSliceArguments ( indices )
        if view is None:
            for d from 0 <= d < self.rank:
                self.indices[d] = indices[d]
            status     = Status_Continue
            value      = RealNDArray_GetItem ( self.cObject, self.indices, &status )
            if status != Status_Continue: raise IndexError ( "Indices out of range." )
            return value
        else:
            return view

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.RealNDArray"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer d, n
        cdef Real    value
        cdef CArrayViewItemIterator *iterator
        if self.isOwner:
            iterator = ArrayViewItemIterator_Allocate ( self.cObject.view, NULL )
            ArrayViewItemIterator_Initialize ( iterator )
            items = []
            while True:
                n = ArrayViewItemIterator_Next ( iterator )
                if n < 0: break
                value = self.cObject.data[n]
                items.append ( value )
            ArrayViewItemIterator_Deallocate ( &iterator )
            return { "items" : items, "shape" : self.shape, "storage" : "FirstDimensionMajor" }
        else:
            shape   = []
            strides = []
            for d from 0 <= d < self.rank:
                shape.append   ( self.cObject.view.extents[d] )
                strides.append ( self.cObject.view.strides[d] )
            return { "offset"  : self.cObject.view.offset ,
                     "owner"   : self.owner               ,
                     "shape"   : shape                    ,
                     "storage" : "FirstDimensionMajor"    ,
                     "strides" : strides                  }

    def __init__ ( self, *extents ):
        """Constructor with extents."""
        self._Initialize ( )
        self._Allocate ( extents )

    def __len__ ( RealNDArray self ):
        """Return the size of the array."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, indices, Real value ):
        """Set an item."""
        cdef Integer d
        cdef Status  status
        view = self._CheckSliceArguments ( indices )
        if view is None:
            for d from 0 <= d < self.rank:
                self.indices[d] = indices[d]
            status     = Status_Continue
            RealNDArray_SetItem ( self.cObject, self.indices, value, &status )
            if status != Status_Continue: raise IndexError ( "Indices out of range." )
        else:
#            print "View> ", type ( view ), view.rank, view.shape
            view.Set ( value )

    def __setstate__ ( self, state ):
        """Return the state."""
        cdef Integer d, n
        cdef Real    value
        cdef CArrayViewItemIterator *iterator
        items = state.get ( "items", None )
        shape = state["shape"]
        if items is not None:
            self._Allocate ( shape )
            iterator = ArrayViewItemIterator_Allocate ( self.cObject.view, NULL )
            ArrayViewItemIterator_Initialize ( iterator )
            while True:
                n = ArrayViewItemIterator_Next ( iterator )
                if n < 0: break
                value = items[n]
                self.cObject.data[n] = value
            ArrayViewItemIterator_Deallocate ( &iterator )
        else:
            offset  = state["offset" ]
            owner   = state["owner"  ]
            strides = state["strides"]
            self._Allocate ( shape, doCObject = False )
            self.isOwner = False
            self.owner   = owner
            cView        = ArrayView_Allocate ( self.rank, NULL )
            n = 1
            for d from 0 <= d < self.rank:
                cView.extents[d] = shape  [d]
                cView.strides[d] = strides[d]
                n *= shape[d]
            cView.offset = offset
            cView.size   = n
            if   owner.rank == 1: self._ViewOf1DArray ( owner, cView )
            elif owner.rank == 2: self._ViewOf2DArray ( owner, cView )
            else:                 self._ViewOfNDArray ( owner, cView )

    def _Allocate ( self, shape, doCObject = True ):
        """Constructor."""
        cdef Integer d
        cdef Status  status
        status    = Status_Continue
        self.rank = len ( shape )
        if self.rank > 0:
            self.indices = Memory_Allocate_Array_Integer ( self.rank )
            if self.indices == NULL: status = Status_MemoryAllocationFailure
        if doCObject and ( status == Status_Continue ):
            for d from 0 <= d < self.rank:
                self.indices[d] = shape[d]
            self.cObject = RealNDArray_Allocate ( self.rank, self.indices, &status )
            self.isOwner = True
            self.owner   = None
        if status != Status_Continue: raise ValueError ( "Object allocation error." )

    def _CheckSliceArguments ( self, indices ):
        """Check slice arguments."""
        cdef CArrayView     *cView       = NULL
        cdef CMultiSliceX   *cMultiSlice = NULL
        cdef CSliceX        *cSlice      = NULL
        cdef CStatusHandler  handler
        # . Basic checks.
        d    = 0
        isOK = True
        view = None
        if isinstance ( indices, tuple ) and ( len ( indices ) == self.rank ):
            for ( i, index ) in enumerate ( indices ):
                if   isinstance ( index, int   ): pass
                elif isinstance ( index, slice ): d += 1
                else:
                    isOK = False
                    break
        else: isOK = False
        if not isOK: raise TypeError ( "Expecting an {:d}-tuple of integers and slices as indices.".format ( self.rank ) )
        # . Process a generalized slice expression.
        if d > 0:
            StatusHandler_Initialize ( &handler )
            StatusHandler_Record ( &handler, MultiSliceX_Allocate ( &cMultiSlice, self.rank ) )
            if StatusHandler_Continue ( &handler ) == CTrue:
                for ( i, ( index, extent ) ) in enumerate ( zip ( indices, self.shape ) ):
                    cSlice = &(cMultiSlice.items[i])
                    if isinstance ( index, int ):
                        StatusHandler_Record ( &handler, SliceX_SetFromScalar ( cSlice, index, extent ) )
#                        print "\nInteger> ", index, extent, cSlice.isScalar, cSlice.extent, cSlice.start, cSlice.stop, cSlice.stride
                    elif isinstance ( index, slice ):
                        ( start, stop, stride ) = index.indices ( extent )
                        StatusHandler_Record ( &handler, SliceX_SetFromSlice ( cSlice, start, stop, stride, extent ) )
#                        print "\nSlice> ", start, stop, stride, extent, cSlice.isScalar, cSlice.extent, cSlice.start, cSlice.stop, cSlice.stride
                MultiSliceX_SetRank ( cMultiSlice )
                StatusHandler_Record ( &handler, ArrayView_Slice ( self.cObject.view, cMultiSlice, &cView ) )
#                print "\nView> ", cView.offset, cView.rank, cView.size
#                for i in range ( cView.rank ):
#                    print i, cView.extents[i], cView.strides[i]
                MultiSliceX_Deallocate ( &cMultiSlice )
                if StatusHandler_Continue ( &handler ) == CTrue:
                    if   d == 1: view = self._Slice1D ( cView )
                    elif d == 2: view = self._Slice2D ( cView )
                    else:        view = self._SliceND ( cView )
                else: raise CLibraryError ( "Invalid array slice." )
        return view

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.indices = NULL
        self.isOwner = False
        self.owner   = None
        self.rank    = -1

    cdef _Slice1D ( self, CArrayView *cView ):
        """Create a 1-D slice."""
        cdef Real1DArray view
        view = Real1DArray ( 0 )
        view.isOwner = False
        view.owner   = self
        RealNDArray_Slice1D ( self.cObject, cView, view.cObject, NULL )
        return view

    cdef _Slice2D ( self, CArrayView *cView ):
        """Create a 2-D slice."""
        cdef Real2DArray view
        view = Real2DArray ( 0, 0 )
        view.isOwner = False
        view.owner   = self
        RealNDArray_Slice2D ( self.cObject, cView, view.cObject, NULL )
        return view

    cdef _SliceND ( self, CArrayView *cView ):
        """Create a N-D slice."""
        cdef RealNDArray view
        view = RealNDArray.Raw ( )
        view._Allocate ( cView.rank * [ 1 ], doCObject = False )
        view.isOwner = False
        view.owner   = self
        view.cObject = RealNDArray_SliceND ( self.cObject, cView, NULL )
#        print "ND Slice> ", self, view
#        if self.cObject == NULL : print "Self is Null"
#        if view.cObject == NULL : print "View is Null"
#        if cView == NULL : print "cView is Null"
        return view

    cdef _ViewOf1DArray ( self, Real1DArray other, CArrayView *cView ):
        """Make a ND view of a 1D array."""
        cdef Status status
        status       = Status_Continue
        self.cObject = RealNDArray_ViewOfRaw ( cView, other.cObject.data, other.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        self.isOwner = False
        self.owner   = other

    cdef _ViewOf2DArray ( self, Real2DArray other, CArrayView *cView ):
        """Make a ND view of a 2D array."""
        cdef Status status
        status       = Status_Continue
        self.cObject = RealNDArray_ViewOfRaw ( cView, other.cObject.data, other.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        self.isOwner = False
        self.owner   = other

    cdef _ViewOfNDArray ( self, RealNDArray other, CArrayView *cView ):
        """Make a ND view of a ND array."""
        cdef Status status
        status       = Status_Continue
        self.cObject = RealNDArray_ViewOfRaw ( cView, other.cObject.data, other.cObject.capacity, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        self.isOwner = False
        self.owner   = other

    def CopyTo ( self, RealNDArray other ):
        """Copying."""
        cdef Status status
        status = Status_Continue
        RealNDArray_CopyTo ( self.cObject, other.cObject, &status )
        if status != Status_Continue: raise CLibraryError ( "Non-conformable arrays." )

    def Extent ( RealNDArray self, Integer dimension ):
        """Return the extent of a dimension."""
        return RealNDArray_Extent ( self.cObject, dimension, NULL )

    def Print ( self, indexWidth =  4 , itemFormat = "{:18.8f}" , itemsPerRow = 2    ,
                      itemWidth  = 19 , log        = logFile    , title       = None ):
        """Printing."""
        cdef Integer d, n
        cdef CArrayViewItemIterator *iterator
        if LogFileActive ( log ) and ( self.size > 0 ):
            # . Get the columns.
            columns = itemsPerRow * ( self.rank * [ indexWidth ] + [ itemWidth ] )
            # . Initialization.
            iterator = ArrayViewItemIterator_Allocate ( self.cObject.view, NULL )
            ArrayViewItemIterator_Initialize ( iterator )
            # . Output the table.
            table = log.GetTable ( columns = columns )
            table.Start ( )
            if title is not None: table.Title ( title )
            while True:
                n = ArrayViewItemIterator_NextWithIndices ( iterator, self.indices )
                if n < 0: break
                for d from 0 <= d < self.rank:
                    table.Entry ( "{:d}".format ( self.indices[d] ) )
                table.Entry ( itemFormat.format ( self.cObject.data[n] ) )
            table.Stop ( )
            # . Finish up.
            ArrayViewItemIterator_Deallocate ( &iterator )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Set ( RealNDArray self, Real value ):
        """Set all the elements of a real N-D array."""
        RealNDArray_Set ( self.cObject, value )

    @classmethod
    def WithExtents ( selfClass, *extents ):
        """Constructor with extents."""
        return selfClass ( *extents )

    # . Properties.
    # . Rank already an attribute.
    property shape:
        def __get__ ( self ):
            """Return the shape of the array."""
            cdef Integer d
            extents = []
            for d from 0 <= d < RealNDArray_Rank ( self.cObject ):
                extents.append ( RealNDArray_Extent ( self.cObject, d, NULL ) )
            return extents

    property size:
        def __get__ ( self ): return RealNDArray_Size ( self.cObject )
