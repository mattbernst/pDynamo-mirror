#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModel.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the non-bonding model for MM/MM and QC/MM interactions.

The model is defined to be independent of any system data.
"""

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Default names.
_STATENAMES = [ "nbState", "qcmmstate" ]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModel:
    """The base class for non-bonding models.

    This class should not be used directly.
    """

    def __copy__ ( self ):
        """Copying."""
        options = self.__getstate__ ( )
        new     = self.__class__ ( **options )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        pass

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModel"

    def __getstate__ ( self ):
        """Return the state."""
        return { }

    def __init__ ( self, **options ):
        """Constructor with options."""
        self._Initialize ( )
        self._Allocate   ( )
        self.SetOptions ( **options )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        self.SetOptions ( **state )

    def _Allocate ( self ):
        """Allocation."""
        pass

    def _Initialize ( self ):
        """Initialization."""
        pass

    def Clear ( self, configuration ):
        """Clear up temporary data."""
        if configuration is not None:
            for name in _STATENAMES:
                if hasattr ( configuration, name ): delattr ( configuration, name )

    def Energy ( self, configuration ):
        """Energy and gradients."""
        pass

    @classmethod
    def FromOptions ( selfClass, **options ):
        """Constructor from options."""
        return selfClass ( **options )

    def QCMMGradients ( self, configuration ):
        """Calculate the QC/MM electrostatic gradients."""
        pass

    def QCMMPotentials ( self, configuration ):
        """Calculate the QC/MM electrostatic potentials."""
        pass

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        pass

    def SetUp ( self, MMAtomContainer mmAtoms, QCAtomContainer qcAtoms, LJParameterContainer ljParameters, LJParameterContainer ljParameters14, Selection fixedatoms, SelfPairList interactions14, SelfPairList exclusions, symmetry, isolates, configuration, log = logFile ):
        """Set up the energy calculation."""
        pass

    def Summary ( self, log = logFile ):
        """Summary."""
        pass
