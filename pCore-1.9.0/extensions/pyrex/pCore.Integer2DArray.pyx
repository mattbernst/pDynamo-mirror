#-------------------------------------------------------------------------------
# . File      : pCore.Integer2DArray.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle integer 2-D arrays."""

from CoreObjects   import CLibraryError
from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Integer2DArray:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        clone = None
        if self.isOwner:
            clone = self.__deepcopy__ ( None )
        else:
            clone = self.__class__ ( 0, 0 )
            if self.owner.rank >= self.rank: self.owner.Make2DView ( clone     , self.cObject.offset, self.cObject.length0, self.cObject.stride0, self.cObject.length1, self.cObject.stride1 )
            else:                            clone._ViewOf1DArray  ( self.owner, self.cObject.offset, self.cObject.length0, self.cObject.stride0, self.cObject.length1, self.cObject.stride1 )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Integer2DArray_Deallocate ( &self.cObject )
            self.isOwner = False
        self.owner = None

    def __deepcopy__ ( self, memo ):
        """Copying."""
        cdef Integer2DArray new
        new = self.__class__.WithExtents ( self.rows, self.columns )
        self.CopyTo ( new )
        return new

    def __getitem__ ( self, indices ):
        """Get an item."""
        cdef Integer i, j, value
        cdef Status  status
        view = self._CheckSliceArguments ( indices )
        if view is None:
            i          = indices[0]
            j          = indices[1]
            status     = Status_Continue
            value      = Integer2DArray_GetItem ( self.cObject, i, j, &status )
            if status != Status_Continue: raise IndexError ( "Indices ({:d},{:d}) out of range.".format ( i, j ) )
            return value
        else:
            return view

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.Integer2DArray"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i, j
        extent0 = self.rows
        extent1 = self.columns
        if self.isOwner:
            items = []
            for i from 0 <= i < extent0:
                for j from 0 <= j < extent1:
                    items.append ( Integer2DArray_GetItem ( self.cObject, i, j, NULL ) )
            return { "items" : items, "shape" : [ extent0, extent1 ], "storage" : "RowMajor" }
        else:
            return { "offset"  : self.cObject.offset ,
                     "owner"   : self.owner          ,
                     "shape"   : [ extent0, extent1 ],
                     "storage" : "RowMajor"          ,
                     "strides" : [ self.cObject.stride0, self.cObject.stride1 ] }

    def __init__ ( self, extent0, extent1 ):
        """Constructor with extents."""
        self._Initialize ( )
        self._Allocate ( extent0, extent1 )

    def __len__ ( self ):
        """Return the size of the array."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, indices, Integer value ):
        """Set an item."""
        cdef Integer i, j
        cdef Status  status
        view = self._CheckSliceArguments ( indices )
        if view is None:
            i          = int ( indices[0] )
            j          = int ( indices[1] )
            status     = Status_Continue
            Integer2DArray_SetItem ( self.cObject, i, j, value, &status )
            if status != Status_Continue: raise IndexError ( "Indices ({:d},{:d}) out of range.".format ( i, j ) )
        else:
            view.Set ( value )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer i, j, extent0, extent1, stride0, stride1, value
        cdef Status  status
        ( extent0, extent1 ) = state["shape"]
        items                = state.get ( "items", None )
        if items is not None:
            self._Allocate ( extent0, extent1 )
            n = 0
            for i from 0 <= i < extent0:
                for j from 0 <= j < extent1:
                    value = items[n]
                    Integer2DArray_SetItem ( self.cObject, i, j, value, NULL )
                    n += 1
        else:
            offset               = state["offset" ]
            owner                = state["owner"  ]
            ( stride0, stride1 ) = state["strides"]
            self._Allocate ( 0, 0 )
            if owner.rank >= self.rank: owner.Make2DView    ( self , offset, extent0, stride0, extent1, stride1 )
            else:                       self._ViewOf1DArray ( owner, offset, extent0, stride0, extent1, stride1 )

    def _Allocate ( self, extent0, extent1 ):
        """Allocation."""
        cdef Status  status
        status       = Status_Continue
        self.cObject = Integer2DArray_Allocate ( extent0, extent1, &status )
        if status   != Status_Continue: raise ValueError ( "Object allocation error." )
        self.isOwner = True
        self.owner   = None

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.owner   = None

    def _CheckSliceArguments ( self, indices ):
        """Check slice arguments."""
        d    = 0
        isOK = True
        view = None
        if isinstance ( indices, tuple ) and ( len ( indices ) == 2 ):
            for index in indices:
                if   isinstance ( index, int   ): pass
                elif isinstance ( index, slice ): d+= 1
                else:
                    isOK = False
                    break
        else: isOK = False
        if not isOK: raise TypeError ( "Expecting a two-tuple of integers and slices as indices." )
        if   d == 1: view = self._Slice1D ( indices[0], indices[1] )
        elif d == 2: view = self._Slice2D ( indices[0], indices[1] )
        return view

    def _Slice1D ( self, i, j ):
        """Create a 1-D slice."""
        cdef Integer        start0, start1, stop0, stop1, stride0, stride1
        cdef Integer1DArray view
        cdef Status         status
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
        view = Integer1DArray ( 0 )
        view.isOwner = False
        view.owner   = self
        status       = Status_Continue
        Integer2DArray_1DSlice ( self.cObject, start0, stop0, stride0, start1, stop1, stride1, view.cObject, &status )
        if status   != Status_Continue: raise CLibraryError ( "Invalid 1-D slice." )
        return view

    def _Slice2D ( self, i, j ):
        """Create a 2-D slice."""
        cdef Integer        start0, start1, stop0, stop1, stride0, stride1
        cdef Integer2DArray view
        cdef Status         status
        # . Indices.
        ( start0, stop0, stride0 ) = i.indices ( self.rows    )
        ( start1, stop1, stride1 ) = j.indices ( self.columns )
        # . View.
        view = self.__class__ ( 0, 0 )
        view.isOwner = False
        view.owner   = self
        status       = Status_Continue
        Integer2DArray_Slice ( self.cObject, start0, stop0, stride0, start1, stop1, stride1, view.cObject, &status )
        if status   != Status_Continue: raise CLibraryError ( "Invalid 2-D slice." )
        return view

    def _ViewOf1DArray ( self, Integer1DArray other, Integer offset, Integer extent0, Integer stride0, Integer extent1, Integer stride1 ):
        """Make a 2D view of a 1D array."""
        cdef Status status
        status = Status_Continue
        Integer2DArray_ViewOfRaw ( self.cObject, offset, extent0, stride0, extent1, stride1, other.cObject.data, other.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        self.isOwner = False
        self.owner   = other

    def CopyTo ( self, Integer2DArray other ):
        """Copying."""
        Integer2DArray_CopyTo ( self.cObject, other.cObject, NULL )

    def Extent ( self, Integer dimension ):
        """Return the extent of a dimension."""
        return Integer2DArray_Length ( self.cObject, dimension )

    def Make1DView ( self, Integer1DArray view, Integer offset, Integer extent, Integer stride ):
        """Make a 1D view."""
        cdef Status status
        status = Status_Continue
        Integer1DArray_ViewOfRaw ( view.cObject, offset, extent, stride, self.cObject.data, self.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        view.isOwner = False
        view.owner   = self

    def Make2DView ( self, Integer2DArray view, Integer offset, Integer extent0, Integer stride0, Integer extent1, Integer stride1 ):
        """Make a 2D view."""
        cdef Status status
        status = Status_Continue
        Integer2DArray_ViewOfRaw ( view.cObject, offset, extent0, stride0, extent1, stride1, self.cObject.data, self.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        view.isOwner = False
        view.owner   = self

    def Print ( self, columnLabels = None , indexWidth = 6    , itemFormat = "{:10d}" , itemsPerRow =  10   ,
                      itemWidth    = 11   , rowLabels  = None , log        = logFile  , title       = None ):
        """Printing."""
        if LogFileActive ( log ):
            ( numberRows, numberColumns ) = self.shape
            widthIndex = 0
            # . Column labels - column numbers.
            labels = []
            for i in range ( numberColumns ): labels.append ( "{:d}".format ( i ) )
            columnLabelSets = [ ( None, labels, "center" ) ]
            if columnLabels is not None:
                for ( tag, items ) in columnLabels:
                    alignment  = "right"
                    widthIndex = max ( widthIndex, len ( tag ) )
                    labels     = []
                    for i in range ( len ( items ) ):
                        item = items[i]
                        if   isinstance ( item, float ): labels.append ( itemFormat.format ( item ) )
                        elif isinstance ( item, int   ): labels.append (     "{:d}".format ( item ) )
                        else:
                            alignment = "center"
                            labels.append ( item )
                    columnLabelSets.append ( ( tag, labels, alignment ) )
            # . Row labels.
            if ( rowLabels is None ) or ( len ( rowLabels ) != numberRows ):
                rowLabels = []
                for i in range ( numberRows ): rowLabels.append ( "{:d}".format ( i ) )
                widthIndex = max ( widthIndex, len ( "{:d}".format ( numberRows ) ) )
            else:
                for label in rowLabels: widthIndex = max ( len ( label ), widthIndex )
            widthIndex = max ( indexWidth, widthIndex + 2 )
            # . Tables.
            ( ntables, nlast ) = divmod ( numberColumns, itemsPerRow )
            if nlast > 0:
                ntables = ntables + 1
            for itable in range ( ntables ):
                if ( itable == ntables - 1 ) and ( nlast > 0 ): nelements = nlast
                else:                                           nelements = itemsPerRow
                columns = [ widthIndex ]
                for i in range ( nelements ): columns.append ( itemWidth )
                table = log.GetTable ( columns = columns )
                table.Start ( )
                if title is not None: table.Title ( title )
                for ( tag, labels, alignment ) in columnLabelSets:
                    table.Entry ( tag, alignment = "left" )
                    for i in range ( nelements ): table.Entry ( labels[i + itable * itemsPerRow], alignment = alignment )
                for irow in range ( numberRows ):
                    table.Entry ( rowLabels[irow], alignment = "left" )
                    for i in range ( nelements ): table.Entry ( itemFormat.format ( Integer2DArray_GetItem ( self.cObject, irow, i + itable * itemsPerRow, NULL ) ) )
                table.Stop ( )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Set ( self, Integer value ):
        """Set all the items."""
        Integer2DArray_Set ( self.cObject, value )

    @classmethod
    def WithExtents ( selfClass, extent0, extent1 ):
        """Constructor."""
        cdef Integer2DArray self
        return selfClass ( extent0, extent1 )

    # . Properties.
    property columns:
        def __get__ ( self ): return Integer2DArray_Length ( self.cObject, 1 )

    property isSquare:
        def __get__ ( self ): return ( self.rows == self.columns )

    property rank:
        def __get__ ( self ): return 2

    property rows:
        def __get__ ( self ): return Integer2DArray_Length ( self.cObject, 0 )

    property shape:
        def __get__ ( self ): return [ self.rows, self.columns ]

    property size:
        def __get__ ( self ): return Integer2DArray_Length ( self.cObject, -1 )
