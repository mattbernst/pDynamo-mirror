#-------------------------------------------------------------------------------
# . File      : pCore.Matrix33.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle matrices of size 3 x 3."""

import math

from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

# . For the moment this does not support slicing. Worth doing?
# . Currently, these objects are never views.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Matrix33:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef Matrix33 new
        new = self.__class__.Uninitialized ( )
        self.CopyTo ( new )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Matrix33_Deallocate ( &self.cObject )
            self.isOwner = False
        self.owner = None

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getitem__ ( self, indices ):
        """Get an item."""
        cdef Integer i, j
        cdef Real    value
        cdef Status  status
        try:    ( i, j ) = indices
        except: raise TypeError ( "Expecting two-tuple of integers as indices." )
        status = Status_Continue
        value  = Matrix33_GetItem ( self.cObject, i, j, &status )
        if status != Status_Continue: raise IndexError ( "Indices ({:d},{:d}) out of range.".format ( i, j ) )
        return value

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.Matrix33"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i, j
        cdef Real    value
        items = []
        for i from 0 <= i < 3:
            for j from 0 <= j < 3:
                value = Matrix33_GetItem ( self.cObject, i, j, NULL )
                items.append ( value )
        return { "items" : items, "shape" : [ 3, 3 ], "storage" : "RowMajor" }

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def __len__ ( Matrix33 self ):
        """Return the size of the matrix."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, indices, Real value ):
        """Set an item."""
        cdef Integer i, j
        cdef Status  status
        try:    ( i, j ) = indices
        except: raise TypeError ( "Expecting two-tuple of integers as indices." )
        status = Status_Continue
        Matrix33_SetItem ( self.cObject, i, j, value, &status )
        if status != Status_Continue: raise IndexError ( "Indices ({:d},{:d}) out of range.".format ( i, j ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        items = state["items"]
        self._Allocate ( )
        n = 0
        for i from 0 <= i < 3:
            for j from 0 <= j < 3:
                Matrix33_SetItem ( self.cObject, i, j, items[n], NULL )
                n = n + 1

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = Matrix33_Allocate ( )
        self.isOwner = True
        self.owner   = None

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.owner   = None

    def ApplyTo ( self, Vector3 vector3 ):
        """Apply the matrix to another object - usually a vector."""
        Matrix33_ApplyToVector3 ( self.cObject, vector3.cObject )

    def CopyTo ( self, Matrix33 other ):
        """Copying."""
        Matrix33_CopyTo ( self.cObject, other.cObject, NULL )

    @classmethod
    def Identity ( selfClass ):
        """Constructor."""
        cdef Matrix33 self
        self = selfClass.Null ( )
        self[0,0] = 1.0 ; self[1,1] = 1.0 ; self[2,2] = 1.0
        return self

    def Invert ( self, Matrix33 other ):
        """Invert the matrix and put the result in other."""
        Matrix33_Invert ( other.cObject, self.cObject )

# . Make properties?

    def IsIdentity ( self ):
        """Test for the identity."""
        return ( Matrix33_IsIdentity ( self.cObject ) == CTrue )

    def IsImproperRotation ( self ):
        """Test for an improper rotation."""
        return ( Matrix33_IsImproperRotation ( self.cObject ) == CTrue )

    def IsOrthogonal ( self ):
        """Test for orthogonality."""
        return ( Matrix33_IsOrthogonal ( self.cObject ) == CTrue )

    def IsProperRotation ( self ):
        """Test for a proper rotation."""
        return ( Matrix33_IsProperRotation ( self.cObject ) == CTrue )

    @classmethod
    def Null ( selfClass ):
        """Constructor."""
        cdef Matrix33 self
        self = selfClass.Uninitialized ( )
        self.Set ( 0.0 )
        return self

    def PostMultiplyBy ( self, Matrix33 other ):
        """Postmultiply by another matrix."""
        Matrix33_PostMultiplyBy ( self.cObject, other.cObject )

    def PreMultiplyBy ( self, Matrix33 other ):
        """Premultiply by another matrix."""
        Matrix33_PreMultiplyBy ( self.cObject, other.cObject )

    def Print ( self, indexWidth = 6, itemFormat = "{:18.8f}", itemWidth = 19, log = logFile, title = None ):
        """Printing."""
        cdef Real value
        if LogFileActive ( log ):
            table = log.GetTable ( columns = [ indexWidth, itemWidth, itemWidth, itemWidth ] )
            table.Start ( )
            if title is not None: table.Title ( title )
            table.Heading ( None )
            for i in range ( 3 ): table.Heading ( "{:d}".format ( i ) )
            for irow in range ( 3 ):
                table.Entry ( "{:d}".format ( irow ) )
                for i in range ( 3 ):
                    value = Matrix33_GetItem ( self.cObject, irow, i, NULL )
                    table.Entry (  itemFormat.format ( value ) )
            table.Stop ( )

    def RandomRotation ( self, randomNumberGenerator, tolerance = 1.0e-20 ):
        """Fill the matrix with a random rotation."""
        ralpha = 2.0 * math.pi * ( randomNumberGenerator.NextReal ( ) - 0.5 )
        raxis  = Vector3.Null ( )
        norm2  = -1.0
        while norm2 < tolerance:
            for i in range ( 3 ):
                raxis[i] = 2.0 * ( randomNumberGenerator.NextReal ( ) - 0.5 )
            norm2 = raxis.Norm2 ( )
        raxis.Normalize ( )
        self.RotationAboutAxis ( ralpha, raxis )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Reflection ( self, Vector3 normal ):
        """Fill the matrix with a reflection through the origin."""
        CMatrix33_Reflection ( &self.cObject, normal.cObject )

    def RotationAboutAxis ( self, Real angle, Vector3 axis ):
        """Fill the matrix with a rotation about a normalized axis."""
        cdef Real x, y, z
        x = Vector3_GetItem ( axis.cObject, 0, NULL )
        y = Vector3_GetItem ( axis.cObject, 1, NULL )
        z = Vector3_GetItem ( axis.cObject, 2, NULL )
        CMatrix33_RotationAboutAxis ( &self.cObject, angle, x, y, z )

    def Scale ( self, Real value ):
        """Scaling."""
        Matrix33_Scale ( self.cObject, value )

    def Set ( Matrix33 self, Real value ):
        """Set all the elements of a matrix."""
        Matrix33_Set ( self.cObject, value )

    def Transpose ( Matrix33 self ):
        """Transpose the matrix."""
        Matrix33_Transpose ( self.cObject, NULL )

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( )

    # . Properties.
    property columns:
        def __get__ ( self ): return Matrix33_Length ( self.cObject, 1 )

    property isSquare:
        def __get__ ( self ): return ( self.rows == self.columns )

    property rank:
        def __get__ ( self ): return 2

    property rows:
        def __get__ ( self ): return Matrix33_Length ( self.cObject, 0 )

    property shape:
        def __get__ ( self ): return [ self.rows, self.columns ]

    property size:
        def __get__ ( self ): return Matrix33_Length ( self.cObject, -1 )

#===================================================================================================================================
# . Class methods - can't have same name as instance methods.
#===================================================================================================================================
def Matrix33_RandomRotation ( randomNumberGenerator, tolerance = 1.0e-20 ):
    """Construct a random rotation."""
    cdef Matrix33 self
    self = Matrix33.Uninitialized ( )
    self.RandomRotation ( randomNumberGenerator, tolerance = tolerance )
    return self

def Matrix33_RotationAboutAxis ( Real angle, Vector3 axis ):
    """Construct a rotation corresponding to a rotation about a normalized axis."""
    cdef Matrix33 self
    self = Matrix33.Uninitialized ( )
    self.RotationAboutAxis ( angle, axis )
    return self
