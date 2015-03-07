#-------------------------------------------------------------------------------
# . File      : MMParameterContainer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Base class for MM parameter containers."""

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . The character for parameter keys.
_ParameterKeyFieldSeparator = ":"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMParameterContainer ( object ):
    """Base class for MM parameter containers."""

    defaultAttributes = { "items"               : None      ,
                          "keySeparator"        : _ParameterKeyFieldSeparator ,
                          "label"               : None      ,
                          "properties"          : None      ,
                          "rawItems"            : None      ,
                          "termLabel"           : None      ,
                          "useStrictAssignment" : True      }

    defaultStrictAssignment = True
    defaultTermLabel        = "MM Term"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        return {}

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        self.termLabel           = self.__class__.defaultTermLabel
        self.useStrictAssignment = self.__class__.defaultStrictAssignment

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        pass

    @staticmethod
    def MakeKey ( *arguments ):
        """Make a key."""
        return tuple ( arguments )

    def MakeMMTerms ( self, atomTypeLabels, termIndices ):
        """Make the appropriate MM terms given a list of term indices."""
        missingParameters = set ( )
        mmTerms           = []
        return ( mmTerms, missingParameters )

    def MakeMMTermsFromConnectivity ( self, atomTypeLabels, connectivity ):
        """Make the appropriate MM terms given a connectivity."""
        termIndices = connectivity.GetTermIndices ( ) # . Doesn't exist.
        return self.MakeMMTerms ( atomTypeLabels, termIndices )

    def ProcessRawItems ( self ):
        """Process the raw items."""
        pass

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
