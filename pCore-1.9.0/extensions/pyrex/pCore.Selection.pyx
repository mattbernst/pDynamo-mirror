#-------------------------------------------------------------------------------
# . File      : pCore.Selection.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle selections."""

from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Selection:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef Selection new
        new         = self.__class__.Raw ( )
        new.cObject = Selection_Clone ( self.cObject )
        new.isOwner = True
        return new

    def __contains__ ( self, Integer value ):
       """Membership."""
       if self.cObject == NULL:
           return False
       else:
           Selection_MakeFlags ( self.cObject, -1 )
           if ( value >= 0 ) and ( value < self.cObject.nflags ):
               return ( self.cObject.flags[value] == CTrue )
           else:
               return False

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Selection_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getitem__ ( self, Integer index ):
        """Get an item."""
        if ( index < 0 ) or ( index >= len ( self ) ): raise IndexError
        else:                                          return self.cObject.indices[index]

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.Selection"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i
        items = []
        for i from 0 <= i < self.size:
            items.append ( self.cObject.indices[i] )
        return { "items" : items }

    def __init__ ( self, iterable ):
        """Constructor given an iterable."""
        cdef Integer i, v
        size = len ( iterable )
        self._Initialize ( )
        self._Allocate ( size )
        for ( i, v ) in enumerate ( iterable ):
            if v < 0: raise ValueError ( "Invalid item in selection - {:d}.".format ( v ) )
            self.cObject.indices[i] = v
        Selection_Sort ( self.cObject )

    def __len__ ( self ):
        """Length."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

#    def __setitem__ ( self, Integer index, Integer value ):
#        """Set an item - this should not be allowed and will be removed."""
#        if ( index < 0 ) or ( index >= len ( self ) ): raise IndexError
#        else:
#            self.cObject.indices[index] = value
#            self.cObject.QSORTED        = CFalse

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer i
        items = state["items"]
        self._Allocate ( len ( items ) )
        for i from 0 <= i < self.size:
            self.cObject.indices[i] = items[i]
        self.Sort ( )

    def _Allocate ( self, size ):
        """Allocation."""
        self.cObject = Selection_Allocate ( size )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def Complement ( self, Integer upperBound = -1 ):
        """Return a complementary selection."""
        cdef Selection new
        new         = Selection.Raw ( )
        new.cObject = Selection_Complement ( self.cObject, upperBound )
        new.isOwner = True
        return new

    def Exclude ( self, object indices ):
        """Exclude indices from the selection - this is done in place."""
        cdef Integer i, n, upperBound
        if len ( indices ) > 0:
            upperBound = Selection_UpperBound ( self.cObject )
            Selection_MakeFlags ( self.cObject, upperBound )
            for i in indices:
                if ( i >= 0 ) and ( i < upperBound ): self.cObject.flags[i] = CFalse
            n = 0
            for i from 0 <= i < upperBound:
                if ( self.cObject.flags[i] == CTrue ):
                    self.cObject.indices[n] = i
                    n = n + 1
            self.cObject.nindices = n
            Selection_ClearPositions ( self.cObject )
            self.cObject.QSORTED = CTrue

    @classmethod
    def FromIterable ( selfClass, iterable ):
        """Constructor given an iterable."""
        return selfClass ( iterable )

    def Increment ( self, Integer increment ):
        """Increment the indices in the selection."""
        cdef Integer i, lower
        if len ( self ) > 0:
            lower = self.cObject.indices[0]
            if ( lower + increment ) < 0: raise ValueError ( "Invalid selection increment - {:d}.".format ( increment ) )
            Selection_ClearRepresentations ( self.cObject )
            for i from 0 <= i < self.cObject.nindices:
                self.cObject.indices[i] = self.cObject.indices[i] + increment

    def Merge ( self, items, information = {} ):
        """Merging - slow version."""
        # . Initialization.
        new = None
        # . Get the selection increments.
        increments = information.get ( "indexIncrements", None )
        # . Do the merge.
        if increments is not None:
            # . Get the old indices.
            items = []
            for ( increment, item ) in zip ( increments, [ self ] + items ):
                state = item.__getstate__ ( )
                old   = state["items"]
                for i in old: items.append ( i + increment )
            # . Allocate the object.
            new = self.__class__.FromIterable ( items )
        return new

    def Position ( self, Integer value ):
       """Return the position of value in the selection."""
       if self.cObject == NULL:
           return -1
       else:
           Selection_MakePositions ( self.cObject, -1 )
           if ( value >= 0 ) and ( value < self.cObject.npositions ): return self.cObject.positions[value]
           else:                                                      return -1

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef Selection new
        new         = Selection.Raw ( )
        new.cObject = Selection_Prune ( self.cObject, selection.cObject )
        new.isOwner = True
        return new

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Sort ( self ):
        """Sorting."""
        Selection_Sort ( self.cObject )

    def UpperBound ( self ):
        """Return the upperbound for the selection."""
        return Selection_UpperBound ( self.cObject )

    property size:
        def __get__ ( self ): return Selection_Size ( self.cObject )

#===================================================================================================================================
# . Static methods.
#===================================================================================================================================
def Selection_SortComparison ( Selection self, Selection other ):
    """A comparison of two selections."""
    return Selection_Compare ( self.cObject, other.cObject )
