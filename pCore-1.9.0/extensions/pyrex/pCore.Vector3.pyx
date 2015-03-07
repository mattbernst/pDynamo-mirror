#-------------------------------------------------------------------------------
# . File      : pCore.Vector3.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle vectors of length 3."""

from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

# . For the moment this does not support slicing. Worth doing?
# . Currently, these objects can only be views of coordinates3 objects.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Vector3:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        clone = None
        if self.isOwner:
            clone = self.__deepcopy__ ( None )
        else:
            clone = self.__class__ ( )
            self.owner.MakeVector3View ( clone, self.cObject.offset, self.cObject.length, self.cObject.stride )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Vector3_Deallocate ( &self.cObject )
            self.isOwner = False
        self.owner = None

    def __deepcopy__ ( self, memo ):
        """Copying."""
        cdef Vector3 new
        new = self.__class__.Uninitialized ( )
        self.CopyTo ( new )
        return new

    def __getitem__ ( self, i ):
        """Get an item."""
        if isinstance ( i, int ):
            if i < 0: i += len ( self )
            if ( i < 0 ) or ( i >= len ( self ) ): raise IndexError ( "Index {:d} out of range.".format ( i ) )
            return Vector3_GetItem ( self.cObject, i, NULL )
        else: raise TypeError ( "Expecting integer not {!r}.".format ( type ( i ) ) )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.Vector3"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i
        if self.isOwner:
            items = []
            for i from 0 <= i < self.size:
                items.append ( Vector3_GetItem ( self.cObject, i, NULL ) )
            return { "items" : items }
        else:
            return { "extent" : self.cObject.length ,
                     "offset" : self.cObject.offset ,
                     "owner"  : self.owner          ,
                     "stride" : self.cObject.stride }

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def __len__ ( Vector3 self ):
        """Return the size of the vector."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, i, Real value ):
        """Set an item."""
        if isinstance ( i, int ):
            if i < 0: i += len ( self )
            if ( i < 0 ) or ( i >= len ( self ) ): raise IndexError ( "Index {:d} out of range.".format ( i ) )
            Vector3_SetItem ( self.cObject, i, value, NULL )
        else: raise TypeError ( "Expecting integer not {!r}.".format ( type ( i ) ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer i, extent, offset, size, stride
        cdef Real    value
        items = state.get ( "items", None )
        if items is not None:
            size = len ( items )
            self._Allocate ( )
            for i from 0 <= i < size:
                value = items[i]
                Vector3_SetItem ( self.cObject, i, value, NULL )
        else:
            extent = state["extent"]
            offset = state["offset"]
            owner  = state["owner" ]
            stride = state["stride"]
            self._Allocate ( )
            owner.MakeVector3View ( self, offset, extent, stride )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = Vector3_Allocate ( )
        self.isOwner = True
        self.owner   = None

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.owner   = None

    def AbsoluteMaximum ( Vector3 self ):
        """Return the maximum absolute value in the vector."""
        return Vector3_AbsoluteMaximum ( self.cObject )

    def AddScalar ( Vector3 self, Real value ):
        """Add a scalar."""
        Vector3_AddScalar ( self.cObject, value )

    def AddScaledVector3 ( Vector3 self, Real value, Vector3 other ):
        """Add a scaled vector."""
        Vector3_AddScaledVector ( self.cObject, value, other.cObject, NULL )

# . Better term?
#    def BasisVector ( self, Integer i ):
#        """Make a basis vector."""
#        Vector3_Set ( self.cObject, 0.0 )
#        self[i] = 1.0

    def CopyTo ( Vector3 self, Vector3 other ):
        """In-place copy of one vector to another."""
        Vector3_CopyTo ( self.cObject, other.cObject, NULL )

    def Cross ( Vector3 self, Vector3 other ):
        """In-place cross-product."""
        Vector3_CrossProduct ( self.cObject, other.cObject )

    def Divide ( self, Vector3 other ):
        """Divide two arrays element wise."""
        Vector3_Divide ( self.cObject, other.cObject, NULL )

    def Dot ( Vector3 self, Vector3 other ):
        """Return the dot product of two vectors."""
        cdef Real result
        result = Vector3_Dot ( self.cObject, other.cObject, NULL )
        return result

    def Multiply ( self, Vector3 other ):
        """Multiply two arrays element wise."""
        Vector3_Multiply ( self.cObject, other.cObject, NULL )

    def Norm2 ( Vector3 self ):
        """Return the norm of the vector."""
        return Vector3_Norm2 ( self.cObject )

    def Normalize ( Vector3 self, tolerance = None ):
        """Normalization."""
        cdef Real tol
        if tolerance is None:
            Vector3_Normalize ( self.cObject, NULL, NULL )
        else:
            tol = float ( tolerance )
            Vector3_Normalize ( self.cObject, &tol, NULL )

    @classmethod
    def Null ( selfClass ):
        """Constructor."""
        cdef Vector3 self
        self = selfClass.Uninitialized ( )
        self.Set ( 0.0 )
        return self

    def Print ( self, itemFormat = "{:18.8f}", itemWidth = 19, log = logFile, title = None ):
        """Printing."""
        cdef Real value
        cdef Integer    i
        if LogFileActive ( log ):
            table = log.GetTable ( columns = 3 * [ itemWidth ] )
            table.Start ( )
            if title is not None: table.Title   ( title )
            for i in range ( 3 ): table.Heading ( "{:d}".format ( i ) )
            for i in range ( 3 ):
                value = Vector3_GetItem ( self.cObject, i, NULL )
                table.Entry   ( itemFormat.format ( value ) )
            table.Stop ( )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def RMSValue ( Vector3 self ):
        """Determine the RMS value of the elements."""
        return Vector3_RootMeanSquare ( self.cObject )

    def Scale ( Vector3 self, Real value ):
        """Scale all the elements of a vector."""
        Vector3_Scale ( self.cObject, value )

    def Set ( Vector3 self, Real value ):
        """Set all the elements of a vector."""
        Vector3_Set ( self.cObject, value )

    def Translate ( self, Vector3 translation ):
        """Translation."""
        Vector3_AddScaledVector ( self.cObject, 1.0, translation.cObject, NULL )

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( )

    @classmethod
    def WithValues ( selfClass, a, b, c ):
        """Constructor."""
        cdef Vector3 self
        self = selfClass.Uninitialized ( )
        self[0] = a ; self[1] = b ; self[2] = c
        return self

    # . Properties.
    property rank:
        def __get__ ( self ): return 1

    property shape:
        def __get__ ( self ): return [ self.size ]

    property size:
        def __get__ ( self ): return Vector3_Length ( self.cObject )
