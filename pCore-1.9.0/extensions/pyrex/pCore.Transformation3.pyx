#-------------------------------------------------------------------------------
# . File      : pCore.Transformation3.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle transformations of dimension 3."""

from Clone         import Clone
from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Transformation3:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef Transformation3 new
        new = self.__class__.Uninitialized ( )
        Transformation3_Copy ( new.cObject, self.cObject )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Transformation3_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.Transformation3"

    def __getstate__ ( self ):
        """Return the state."""
        return { "rotation" : Clone ( self.rotation ), "translation" : Clone ( self.translation ) }

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        rotation    = state["rotation"   ]
        translation = state["translation"]
        rotation.CopyTo    ( self.rotation    )
        translation.CopyTo ( self.translation )

    def _Allocate ( self ):
        """Allocation."""
        cdef Matrix33 r
        cdef Vector3  t
        # . C object.
        self.cObject = Transformation3_Allocate ( )
        self.isOwner = True
        # . Python views.
        if self.cObject != NULL:
            r                = Matrix33.Raw ( )
            r.cObject        = self.cObject.rotation
            r.isOwner        = False
            r.owner          = self
            t                = Vector3.Raw  ( )
            t.cObject        = self.cObject.translation
            t.isOwner        = False
            t.owner          = self
            self.rotation    = r
            self.translation = t

    def _Initialize ( self ):
        """Initialization."""
        self.cObject     = NULL
        self.isOwner     = False
        self.rotation    = None
        self.translation = None

    def IsIdentity ( self ):
        """Is this transformation the identity?"""
        return ( Transformation3_IsIdentity ( self.cObject ) == CTrue )

    @classmethod
    def Null ( selfClass ):
        """Constructor."""
        cdef Transformation3 self
        self = selfClass.Uninitialized ( )
        self.rotation.Set    ( 0.0 )
        self.translation.Set ( 0.0 )
        return self

    def Orthogonalize ( self, Matrix33 A, Matrix33 B ):
        """Orthogonalization."""
        Transformation3_Orthogonalize ( self.cObject, A.cObject, B.cObject )

    def Print ( self, log = logFile, title = None ):
        """Printing."""
        if LogFileActive ( log ):
            if title is not None: log.Heading ( title, QBLANKLINE = True )
            self.rotation.Print    ( log = log, title = "Rotation"    )
            self.translation.Print ( log = log, title = "Translation" )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( )

#===================================================================================================================================
# . Class methods.
#===================================================================================================================================
def Transformation3_FromSymmetryOperationString ( ostring ):
    """Construct a transformation from a string encoding a crystallographic symmetry operation."""
    cdef Transformation3 self
    self = None
    if isinstance ( ostring, basestring ):
        self = Transformation3.Null ( )
        # . Remove spaces and enclosing parentheses.
        t = ostring.replace ( " ", "" )
        l = list ( t )
        if l[ 0] == "(" : l.pop ( 0 )
        if l[-1] == ")" : l.pop (   )
        t   = "".join ( l )
        # . Loop over the three different specifications.
        for ( i, s ) in enumerate ( t.split ( "," ) ):
            # . Ensure all numbers end in a decimal point.
            ns = [ s[0:1] ]
            for ( j, c ) in enumerate ( s[1:] ):
                if s[j:j+1].isdigit ( ) and not ( c.isdigit ( ) or c == "." ): ns.append ( "." )
                ns.append ( c.lower ( ) )
            if ns[-1].isdigit ( ): ns.append ( "." )
            s = "".join ( ns )
            # . Determine the transformation.
            x  = 0.0
            y  = 0.0
            z  = 0.0
            t  = eval ( s, {}, { "x" : x, "y" : y, "z" : z } )
            x  = 1.0
            rx = eval ( s, {}, { "x" : x, "y" : y, "z" : z } ) - t
            y  = 1.0
            ry = eval ( s, {}, { "x" : x, "y" : y, "z" : z } ) - t - rx
            z  = 1.0
            rz = eval ( s, {}, { "x" : x, "y" : y, "z" : z } ) - t - rx - ry
            self.translation[i] = t
            self.rotation [i,0] = rx
            self.rotation [i,1] = ry
            self.rotation [i,2] = rz
    return self
