#-------------------------------------------------------------------------------
# . File      : pMolecule.SymmetryParameters.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Interface to the symmetry parameters C type."""

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetryParameters:
    """Define a set of symmetry parameters."""

    def __copy__ ( self ):
        """Copying."""
        state = self.__getstate__ ( )
        new   = self.__class__.Raw ( )
        new.__setstate__ ( state )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            SymmetryParameters_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.SymmetryParameters"

    def __getstate__ ( self ):
        """Return the state."""
        return { "a" : self.cObject.a, "b" : self.cObject.b, "c" : self.cObject.c, "alpha" : self.cObject.alpha, "beta" : self.cObject.beta, "gamma" : self.cObject.gamma }

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
        a     = state["a"]
        b     = state["b"]
        c     = state["c"]
        alpha = state["alpha"]
        beta  = state["beta" ]
        gamma = state["gamma"]
        self.SetCrystalParameters ( a, b, c, alpha, beta, gamma )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = SymmetryParameters_Allocate ( )
        self.isOwner = True
        self.AllocatePythonObjects ( )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject  = NULL
        self.isOwner  = True
        self.M        = None
        self.inverseM = None

    def AllocatePythonObjects ( self ):
        """Allocate the Python objects."""
        cdef Matrix33 inverseM, M
        M                = Matrix33.Raw ( )
        M.isOwner        = False
        M.cObject        = NULL
        self.M           = M
        inverseM         = Matrix33.Raw ( )
        inverseM.isOwner = False
        inverseM.cObject = NULL
        self.inverseM    = inverseM

    def CenterCoordinatesByAtom ( self, Coordinates3 coordinates3, Selection selection = None ):
        """Center coordinates by atom."""
        cdef CSelection *cSelection
        if selection is None: cSelection = NULL
        else:                 cSelection = selection.cObject
        SymmetryParameters_CenterCoordinates3ByIndex ( self.cObject, cSelection, coordinates3.cObject )

    def CenterCoordinatesByIsolate ( self, Coordinates3 coordinates3, SelectionContainer isolates, Selection selection = None ):
        """Center coordinates by isolate."""
        cdef CSelection *cSelection
        if selection is None: cSelection = NULL
        else:                 cSelection = selection.cObject
        SymmetryParameters_CenterCoordinates3ByIsolate ( self.cObject, isolates.cObject, cSelection, coordinates3.cObject )

    def MakeMinimumImageVector3 ( self, Vector3 r, Vector3 dr = None ):
        """Apply the minimum image convention to a vector."""
        cdef CVector3 *cdr
        if dr is None: cdr = NULL
        else:          cdr = dr.cObject
        SymmetryParameters_MakeMinimumImageVector3 ( self.cObject, r.cObject, cdr )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Set the crystal parameters."""
        cdef Matrix33 inverseM, M
        if self.cObject == NULL: self._Allocate ( )
        SymmetryParameters_SetCrystalParameters ( self.cObject, a, b, c, alpha, beta, gamma )
        # . Python representations of component objects.
        M                = self.M
        M.cObject        = self.cObject.M
        inverseM         = self.inverseM
        inverseM.cObject = self.cObject.inverseM

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( )

    # . Properties.
    property a:
        def __get__ ( self ): return self.cObject.a
    property b:
        def __get__ ( self ): return self.cObject.b
    property c:
        def __get__ ( self ): return self.cObject.c
    property alpha:
        def __get__ ( self ): return self.cObject.alpha
    property beta:
        def __get__ ( self ): return self.cObject.beta
    property gamma:
        def __get__ ( self ): return self.cObject.gamma
    property volume:
        def __get__ ( self ): return SymmetryParameters_Volume ( self.cObject )
