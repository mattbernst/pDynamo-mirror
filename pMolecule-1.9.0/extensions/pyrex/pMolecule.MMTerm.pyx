#-------------------------------------------------------------------------------
# . File      : pMolecule.MMTerm.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Base class for MM terms."""

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MMTerm:
    """The base class for MM terms.

    This class should not be used directly.
    """

    # . Public methods.
    def __dealloc__ ( self ):
        """Finalization."""
        pass

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.MMTerm"

    def __init__ ( self, numberOfParameters, numberOfTerms, label = "MM Term" ):
        """Constructor."""
        pass

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def _Allocate ( self, numberOfParameters, numberOfTerms ):
        """Allocation."""
        pass

    def _Initialize ( self ):
        """Initialization."""
        pass

    def MergeKeys ( self, parameterKeys, parameters ):
        """Merge parameter keys."""
        newKeys = set ( parameterKeys )
        if len ( newKeys ) < len ( parameterKeys ):
            # . Keys.
            newKeys = list ( newKeys )
            newKeys.sort ( )
            oldToNew = []
            for oldKey in parameterKeys:
                oldToNew.append ( newKeys.index ( oldKey ) )
            # . Parameters.
            newParameters = []
            for newKey in newKeys:
                newParameters.append ( parameters[parameterKeys.index ( newKey )] )
            # . Finish up.
            return ( True , oldToNew, newKeys      , newParameters )
        else:
            return ( False, None    , parameterKeys, parameters    )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self
