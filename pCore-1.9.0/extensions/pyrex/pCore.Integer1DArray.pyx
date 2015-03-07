#-------------------------------------------------------------------------------
# . File      : pCore.Integer1DArray.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle integer 1-D arrays."""

from CoreObjects   import CLibraryError
from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Integer1DArray:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        clone = None
        if self.isOwner:
            clone = self.__deepcopy__ ( None )
        else:
            clone = self.__class__ ( 0 )
            self.owner.Make1DView ( clone, self.cObject.offset, self.cObject.length, self.cObject.stride )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Integer1DArray_Deallocate ( &self.cObject )
            self.isOwner = False
        self.owner = None

    def __deepcopy__ ( self, memo ):
        """Copying."""
        cdef Integer1DArray new
        new = self.__class__.WithExtent ( self.size )
        self.CopyTo ( new )
        return new

    def __getitem__ ( self, i ):
        """Get an item."""
        if isinstance ( i, int ):
            if i < 0: i += len ( self )
            if ( i < 0 ) or ( i >= len ( self ) ): raise IndexError ( "Index {:d} out of range.".format ( i ) )
            return Integer1DArray_GetItem ( self.cObject, i, NULL )
        elif isinstance ( i, slice ):
            return self._Slice ( *i.indices ( len ( self ) ) )
        else: raise TypeError ( "Expecting integer or slice not {!r}.".format ( type ( i ) ) )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.Integer1DArray"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i
        if self.isOwner:
            items = []
            for i from 0 <= i < self.size:
                items.append ( Integer1DArray_GetItem ( self.cObject, i, NULL ) )
            return { "items" : items }
        else:
            return { "extent" : self.cObject.length ,
                     "offset" : self.cObject.offset ,
                     "owner"  : self.owner          ,
                     "stride" : self.cObject.stride }

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

    def __setitem__ ( self, i, Integer value ):
        """Set an item."""
        if isinstance ( i, int ):
            if i < 0: i += len ( self )
            if ( i < 0 ) or ( i >= len ( self ) ): raise IndexError ( "Index {:d} out of range.".format ( i ) )
            Integer1DArray_SetItem ( self.cObject, i, value, NULL )
        elif isinstance ( i, slice ):
            view = self._Slice ( *i.indices ( len ( self ) ) )
            view.Set ( value )
        else: raise TypeError ( "Expecting integer or slice not {!r}.".format ( type ( i ) ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer i, extent, offset, size, stride, value
        items = state.get ( "items", None )
        if items is not None:
            size = len ( items )
            self._Allocate ( size )
            for i from 0 <= i < size:
                value = items[i]
                Integer1DArray_SetItem ( self.cObject, i, value, NULL )
        else:
            extent = state["extent"]
            offset = state["offset"]
            owner  = state["owner" ]
            stride = state["stride"]
            self._Allocate ( 0 )
            owner.Make1DView ( self, offset, extent, stride )

    def _Allocate ( self, size ):
        """Allocation."""
        cdef Status  status
        status       = Status_Continue
        self.cObject = Integer1DArray_Allocate ( size, &status )
        if status   != Status_Continue: raise ValueError ( "Object allocation error." )
        self.isOwner = True
        self.owner   = None

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.owner   = None

    def _Slice ( self, Integer start, Integer stop, Integer step ):
        """Create a slice."""
        cdef Integer1DArray view
        cdef Status      status
        view = self.__class__ ( 0 )
        view.isOwner = False
        view.owner   = self
        status       = Status_Continue
        Integer1DArray_Slice ( self.cObject, start, stop, step, view.cObject, &status )
        if status   != Status_Continue: raise CLibraryError ( "Invalid slice." )
        return view

    def CopyTo ( self, Integer1DArray other ):
        """Copying."""
        Integer1DArray_CopyTo ( self.cObject, other.cObject, NULL )

    def Make1DView ( self, Integer1DArray view, Integer offset, Integer extent, Integer stride ):
        """Make a 1D view."""
        cdef Status status
        status = Status_Continue
        Integer1DArray_ViewOfRaw ( view.cObject, offset, extent, stride, self.cObject.data, self.cObject.size, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating view." )
        view.isOwner = False
        view.owner   = self

    def Print ( self, itemFormat = "{:10d}", itemsPerRow = 10, itemWidth = 11, log = logFile, title = None ):
        """Printing."""
        dimension = len ( self )
        if LogFileActive ( log ) and ( dimension > 0 ):
            # . Find the number of rows.
            ( numberOfRows, numberRemaining ) = divmod ( dimension, itemsPerRow )
            if numberRemaining > 0: numberOfRows += 1
            # . Initialize the table.
            width = len ( "{:d}".format ( dimension - 1 ) ) + 1
            if itemsPerRow > 1: columns = [ width, 2, width ]
            else:               columns = [ width ]
            columns += min ( dimension, itemsPerRow ) * [ itemWidth ]
            table    = log.GetTable ( columns = columns )
            table.Start ( )
            if title is not None: table.Title ( title )
            # . Row output.
            start = 0
            for iRow in range ( numberOfRows ):
                if ( iRow == numberOfRows - 1 ) and ( numberRemaining > 0 ): numberOfItems = numberRemaining
                else:                                                        numberOfItems = itemsPerRow
                table.Entry ( "{:d}".format ( start ) )
                if itemsPerRow > 1:
                    table.Entry ( " -" )
                    table.Entry ( "{:d}".format ( start + numberOfItems - 1 ) )
                for i in range ( numberOfItems ):
                    value = Integer1DArray_GetItem ( self.cObject, i + start, NULL )
                    table.Entry ( itemFormat.format ( value ) )
                start += numberOfItems
            # . Finish up.
            table.Stop ( )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Set ( self, Integer value ):
        """Set all the items."""
        Integer1DArray_Set ( self.cObject, value )

    @classmethod
    def WithExtent ( selfClass, extent ):
        """Constructor."""
        return selfClass ( extent )

    # . Properties.
    property rank:
        def __get__ ( self ): return 1

    property shape:
        def __get__ ( self ): return [ self.size ]

    property size:
        def __get__ ( self ): return Integer1DArray_Length ( self.cObject )
