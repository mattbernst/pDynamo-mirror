#-------------------------------------------------------------------------------
# . File      : pCore.SelectionContainer.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle selection containers."""

from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultItemName = "Selection"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SelectionContainer:

    def __copy__ ( self ):
        """Copying."""
        cdef SelectionContainer new
        new         = self.__class__.Raw ( )
        new.cObject = SelectionContainer_Clone ( self.cObject )
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            SelectionContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.SelectionContainer"

    def __getitem__ ( self, Integer index ):
        """Get an item."""
        cdef Selection new
        if ( index < 0 ) or ( index >= len ( self ) ):
            raise IndexError
        else:
            new         = Selection.Raw ( )
            new.cObject = self.cObject.items[index]
            new.isOwner = False
            return new

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i, s
        cdef CSelection *selection
        items = []
        n    = len ( self )
        if n > 0:
            for i from 0 <= i < n:
                selection = self.cObject.items[i]
                indices   = []
                for s from 0 <= s < selection.nindices:
                    indices.append ( selection.indices[s] )
                items.append ( indices )
        return { "items" : items }

    def __init__ ( self, size, itemName = _DefaultItemName ):
        """Constructor with size."""
        self._Initialize ( )
        self._Allocate ( size )
        self.itemName = itemName

    def __len__ ( self ):
        """Length."""
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nitems

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer i, s
        cdef CSelection *selection
        items = state["items"]
        size  = len ( items )
        self._Allocate ( size )
        if size > 0:
            for ( i, indices ) in enumerate ( items ):
                selection = Selection_Allocate ( len ( indices ) )
                for ( s, index ) in enumerate ( indices ):
                    selection.indices[s] = index
                Selection_Sort ( selection )
                self.cObject.items[i] = selection

    def _Allocate ( self, size ):
        """Allocation."""
        self.cObject = SelectionContainer_Allocate ( size )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject  = NULL
        self.isOwner  = False
        self.itemName = _DefaultItemName

    def GetAllIndices ( self, Selection indices ):
        """Get the union of all the container selections which contain any of the indices in |indices|."""
        cdef Selection new
        new         = Selection.Raw ( )
        new.cObject = SelectionContainer_GetAllIndices ( self.cObject, indices.cObject )
        new.isOwner = True
        return new

    def MergeIsolates ( self, Selection tomerge ):
        """Merge isolates within the container."""
        SelectionContainer_MergeIsolates ( self.cObject, tomerge.cObject )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SummaryEntry ( self, summary ):
        """Summary entry."""
        if summary is not None: summary.Entry ( "Number of " + self.itemName + "s", "{:d}".format ( len ( self ) ) )
