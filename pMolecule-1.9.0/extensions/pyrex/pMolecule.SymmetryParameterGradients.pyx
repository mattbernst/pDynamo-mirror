#-------------------------------------------------------------------------------
# . File      : pMolecule.SymmetryParameterGradients.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Interface to the symmetry parameter gradients C type."""

from pCore import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetryParameterGradients:
    """Define the gradients corresponding to a set of symmetry parameters."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            SymmetryParameterGradients_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.SymmetryParameterGradients"

    def __getstate__ ( self ):
        """Return the state."""
        return { "dEdA"     : self.cObject.dEda     ,
                 "dEdB"     : self.cObject.dEdb     ,
                 "dEdC"     : self.cObject.dEdc     ,
                 "dEdAlpha" : self.cObject.dEdalpha ,
                 "dEdBeta"  : self.cObject.dEdbeta  ,
                 "dEdGamma" : self.cObject.dEdgamma ,
                 "dEdM"     : self.dEdM             }

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
        self.cObject.dEda     = state["dEdA"    ]
        self.cObject.dEdb     = state["dEdB"    ]
        self.cObject.dEdc     = state["dEdC"    ]
        self.cObject.dEdalpha = state["dEdAlpha"]
        self.cObject.dEdbeta  = state["dEdBeta" ]
        self.cObject.dEdgamma = state["dEdGamma"]
        state["dEdM"].CopyTo ( self.dEdM )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = SymmetryParameterGradients_Allocate ( )
        self.isOwner = True
        self.AllocatePythonObjects ( )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = True
        self.dEdM    = None

    def AllocatePythonObjects ( self ):
        """Allocate the Python objects."""
        cdef Matrix33 dEdM
        dEdM         = Matrix33.Raw ( )
        dEdM.isOwner = False
        self.dEdM    = dEdM
        if self.cObject != NULL: dEdM.cObject = self.cObject.dEdM
        else:                    dEdM.cObject = NULL

    def MakeCrystalDerivatives ( self, SymmetryParameters symmetryParameters ):
        """Make the crystal derivatives from the lattice derivatives."""
        SymmetryParameterGradients_CrystalDerivatives ( self.cObject, symmetryParameters.cObject )

    def MakeFractionalDerivatives ( self, SymmetryParameters symmetryParameters, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Convert [r,M] derivatives to [f,M] derivatives."""
#        print "M:"
#        symmetryParameters.M.Print ( )
#        print "Inverse M:"
#        symmetryParameters.inverseM.Print ( )
        SymmetryParameterGradients_FractionalDerivatives ( self.cObject, symmetryParameters.cObject, coordinates3.cObject, gradients3.cObject )

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

    # . Properties.
    property dEdA:
        def __get__ ( self ): return self.cObject.dEda
    property dEdB:
        def __get__ ( self ): return self.cObject.dEdb
    property dEdC:
        def __get__ ( self ): return self.cObject.dEdc
    property dEdAlpha:
        def __get__ ( self ): return self.cObject.dEdalpha
    property dEdBeta:
        def __get__ ( self ): return self.cObject.dEdbeta
    property dEdGamma:
        def __get__ ( self ): return self.cObject.dEdgamma
