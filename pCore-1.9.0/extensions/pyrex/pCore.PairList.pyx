#-------------------------------------------------------------------------------
# . File      : pCore.PairList.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle cross- and self-pairlists."""

from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CrossPairList:

    # . Public methods.
    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            PairList_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.PairList"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer *cIndices, i, n
        state = None
        n     = len ( self )
        data  = []
        if ( n > 0 ):
            cIndices = PairList_ToIntegerPairArray ( self.cObject )
            if ( cIndices != NULL ):
                for i from 0 <= i < 2 * n: data.append ( cIndices[i] )
                Memory_Deallocate_Integer ( &cIndices )
                state = { "items" : data, "shape" : [ n, 2 ] }
                if self.label is not None: state["label"] = self.label
        return state

    def __init__ ( self, label = "Cross Pair List" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )
        self.label = label

    def __iter__ ( self ):
        """Return an iterator."""
        return CrossPairListIterator ( self )

    def __len__ ( self ):
        """The number of pairs."""
        return PairList_Length ( self.cObject )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer *cIndices, i, j, n
        if state is not None:
            indices  = state["items"]
            n        = len ( indices )
            cIndices = Memory_Allocate_Array_Integer ( n )
            for i from 0 <= i < n: cIndices[i] = indices[i]
            self.cObject = PairList_FromIntegerPairArray ( CFalse, n // 2, cIndices )
            self.isOwner = True
            Memory_Deallocate_Integer ( &cIndices )
            if "label" in state: self.label = state["label"]

    def _Allocate ( self ):
        """Constructor."""
        self.cObject = PairList_Allocate ( CFalse )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.label   = None

    def Merge ( self, items, information = {} ):
        """Merging - slow version."""
        # . Initialization.
        new = None
        # . Get the selection increments.
        increments = information.get ( "indexIncrements", None )
        # . Do the merge.
        if increments is not None:
            # . Get the old indices.
            data = []
            for ( increment, item ) in zip ( increments, [ self ] + items ):
                state = item.__getstate__ ( )
                if state is not None:
                    old   = state["items"]
                    for i in old: data.append ( i + increment )
            # . Allocate the object.
            n0  = len ( data ) // 2
            new = self.__class__.Raw ( )
            new.__setstate__ ( { "items" : data, "label" : self.label, "shape" : [ n0, 2 ] } )
        return new

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SummaryEntry ( self, summary ):
        """Summary entry."""
        if summary is not None: summary.Entry ( self.label, "{:d}".format ( len ( self ) ) )

#===================================================================================================================================
# . Class methods.
#===================================================================================================================================
def CrossPairList_FromIntegerPairs ( indices ):
    """Create a cross-pairlist from a list of integer indices."""
    cdef Integer *cIndices, i, n
    cdef CrossPairList self
    n = len ( indices ) // 2
    if ( n > 0 ):
        cIndices = Memory_Allocate_Array_Integer ( 2 * n )
        for i from 0 <= i < 2 * n: cIndices[i] = indices[i]
        self         = CrossPairList.Raw ( )
        self.cObject = PairList_FromIntegerPairArray ( CFalse, n, cIndices )
        self.isOwner = True
        Memory_Deallocate_Integer ( &cIndices )
    else:
        self = CrossPairList.Raw ( )
    return self

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CrossPairListIterator:

    # . Public methods.
    def __init__ ( self, CrossPairList target ):
        """Constructor."""
        self.cObject = NULL
        if target is not None:
            if target.cObject != NULL:
                self.cObject = target.cObject
                List_Iterate_Initialize ( self.cObject.pairs )
                self.indexedSelection = PairList_Iterate ( self.cObject )
                self.position         = 0

    def __next__ ( self ):
        """Get the next pair."""
        cdef Integer i, j
        # . Get the next pair if necessary.
        if self.indexedSelection != NULL:
            if self.position         >= self.indexedSelection.nindices:
                self.indexedSelection = PairList_Iterate ( self.cObject )
                self.position         = 0
        # . Now return something.
        if self.indexedSelection == NULL: raise StopIteration
        else:
            i = self.indexedSelection.index
            j = self.indexedSelection.indices[self.position]
            self.position = self.position + 1
            return ( i, j )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SelfPairList:

    # . Public methods.
    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            PairList_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.PairList"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer *cIndices, i, n
        state = None
        n     = len ( self )
        data  = []
        if ( n > 0 ):
            cIndices = PairList_ToIntegerPairArray ( self.cObject )
            if ( cIndices != NULL ):
                for i from 0 <= i < 2 * n: data.append ( cIndices[i] )
                Memory_Deallocate_Integer ( &cIndices )
                state = { "items" : data, "shape" : [ n, 2 ] }
                if self.label is not None: state["label"] = self.label
        return state

    def __init__ ( self, label = "Self Pair List" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )
        self.label = label

    def __iter__ ( self ):
        """Return an iterator."""
        return SelfPairListIterator ( self )

    def __len__ ( self ):
        """The number of pairs."""
        return PairList_Length ( self.cObject )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer *cIndices, i, j, n
        if state is not None:
            indices  = state["items"]
            n        = len ( indices )
            cIndices = Memory_Allocate_Array_Integer ( n )
            for i from 0 <= i < n: cIndices[i] = indices[i]
            self.cObject = PairList_FromIntegerPairArray ( CTrue, n // 2, cIndices )
            self.isOwner = True
            Memory_Deallocate_Integer ( &cIndices )
            if "label" in state: self.label = state["label"]

    def _Allocate ( self ):
        """Constructor."""
        self.cObject = PairList_Allocate ( CTrue )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.label   = None

    def GetIsolates ( self, Integer upperBound = -1 ):
        """Get a selection container with isolates from the pairlist."""
        cdef SelectionContainer new
        new         = SelectionContainer.Raw ( )
        new.cObject = SelfPairList_ToIsolateSelectionContainer ( self.cObject, upperBound )
        new.isOwner = True
        return new

    def Merge ( self, items, information = {} ):
        """Merging - slow version."""
        # . Initialization.
        new = None
        # . Get the selection increments.
        increments = information.get ( "indexIncrements", None )
        # . Do the merge.
        if increments is not None:
            # . Get the old indices.
            data = []
            for ( increment, item ) in zip ( increments, [ self ] + items ):
                state = item.__getstate__ ( )
                if state is not None:
                    old   = state["items"]
                    for i in old: data.append ( i + increment )
            # . Allocate the object.
            n0  = len ( data ) // 2
            new = self.__class__.Raw ( )
            new.__setstate__ ( { "items" : data, "label" : self.label, "shape" : [ n0, 2 ] } )
        return new

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef SelfPairList new
        new         = SelfPairList.Raw ( )
        new.cObject = SelfPairList_FromSelfPairList ( self.cObject, selection.cObject, NULL, CTrue )
        new.isOwner = True
        new.label   = self.label 
        return new

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SummaryEntry ( self, summary ):
        """Summary entry."""
        if summary is not None: summary.Entry ( self.label, "{:d}".format ( len ( self ) ) )

    def UpperBound ( self ):
        """Return the upperbound for the selection."""
        return SelfPairList_UpperBound ( self.cObject )

#===================================================================================================================================
# . Class methods.
#===================================================================================================================================
def SelfPairList_FromIntegerPairs ( indices ):
    """Create a self-pairlist from a list of integer indices."""
    cdef Integer *cIndices, i, n
    cdef SelfPairList self
    n = len ( indices ) // 2
    if ( n > 0 ):
        cIndices = Memory_Allocate_Array_Integer ( 2 * n )
        for i from 0 <= i < 2 * n: cIndices[i] = indices[i]
        self = SelfPairList.Raw ( )
        self.cObject = PairList_FromIntegerPairArray ( CTrue, n, cIndices )
        Memory_Deallocate_Integer ( &cIndices )
    else:
        self = SelfPairList.Raw ( )
    return self

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SelfPairListIterator:

    # . Public methods.
    def __init__ ( self, SelfPairList target ):
        """Constructor."""
        self.cObject = NULL
        if target is not None:
            if target.cObject != NULL:
                self.cObject = target.cObject
                List_Iterate_Initialize ( self.cObject.pairs )
                self.indexedSelection = PairList_Iterate ( self.cObject )
                self.position         = 0

    def __next__ ( self ):
        """Get the next pair."""
        cdef Integer i, j
        # . Get the next pair if necessary.
        if self.indexedSelection != NULL:
            if self.position         >= self.indexedSelection.nindices:
                self.indexedSelection = PairList_Iterate ( self.cObject )
                self.position         = 0
        # . Now return something.
        if self.indexedSelection == NULL: raise StopIteration
        else:
            i = self.indexedSelection.index
            j = self.indexedSelection.indices[self.position]
            self.position = self.position + 1
            return ( i, j )
