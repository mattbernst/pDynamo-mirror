#-------------------------------------------------------------------------------
# . File      : Symmetry.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Basic symmetry class."""

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Symmetry ( object ):
    """Hold symmetry information."""

    # . Default attributes.
    defaultAttributes = { "crystalClass"    : None ,
                          "transformations" : None }

    def __getstate__ ( self ):
        state = {}
        for key in self.__class__.defaultAttributes:
            value = self.__dict__.get ( key, None )
            if value is not None: state[key] = value
        return state

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ):
            setattr ( self, key, value )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        for key in self.__class__.defaultAttributes:
            value = state.get ( key, None )
            if value is not None: self.__dict__[key] = value

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            if self.crystalClass is not None:
                summary = log.GetSummary ( )
                summary.Start ( "Symmetry Summary" )
                summary.Entry ( "Crystal Class", self.crystalClass.Label ( ) )
                summary.Stop ( )
            if self.transformations is not None:
                self.transformations.Summary ( log = log )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
