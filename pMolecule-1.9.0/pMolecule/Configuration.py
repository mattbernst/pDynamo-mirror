#-------------------------------------------------------------------------------
# . File      : Configuration.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Class for holding the configuration of a system.

A configuration holds mutable or temporary data such as coordinates, a
wavefunction or velocities.

Attributes are either permanent (e.g. coordinates3 and symmetryParameters),
persistent (e.g. NB and QC states) or temporary (e.g. things that are
re-evaluated at each energy calculation).
"""

from pCore import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Configuration ( object ):
    """Define a configuration."""

    permanentAttributes = [ "coordinates3", "symmetryParameters" ]

    def __getstate__ ( self ):
        """Return a state consisting of permanent attributes only."""
        mapping = {}
        for attribute in self.__class__.permanentAttributes:
            item = getattr ( self, attribute, None )
            if item is not None: mapping[attribute] = item
        return mapping

    def __init__ ( self ): pass

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ): self.__dict__.update ( state )

    def Clear ( self ):
        """Clear all non-permanent attributes."""
        for key in self.__dict__.keys ( ):
            if ( key not in self.__class__.permanentAttributes ): delattr ( self, key )
        self.temporaryAttributes.clear ( )

    def ClearTemporaryAttributes ( self ):
        """Clear temporary attributes."""
        for key in self.temporaryAttributes: delattr ( self, key )
        self.temporaryAttributes.clear ( )

    def DeleteTemporaryAttribute ( self, label ):
        """Delete a temporary attribute."""
        if label in self.temporaryAttributes:
            delatt ( self, label )
            self.temporaryAttributes.remove ( label )
        else: raise AttributeError ( "Temporary attribute \"" + label + "\" not found in configuration." )

    def FlagAttributeAsTemporary ( self, label ):
        """Flag an attribute as temporary."""
        if ( label not in self.__class__.permanentAttributes ) and ( label in self.__dict__ ): self.temporaryAttributes.add ( label )
        else: raise AttributeError ( "Attribute \"" + label + "\" not found in configuration." )

    def Merge ( self, others, information = {} ):
        """Merging."""
        # . Gather coordinates3.
        toMergeItems = []
        for item in [ self ] + list ( others ):
            coordinates3 = getattr ( item, "coordinates3", None )
            if coordinates3 is not None: toMergeItems.append ( coordinates3 )
        # . Merging.
        if len ( toMergeItems ) == ( len ( others ) + 1 ):
            merged = self.__class__ ( )
            merged.coordinates3 = toMergeItems[0].Merge ( toMergeItems[1:], information )
        else:
            merged = None
        # . Finish up.
        return merged

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        coordinates3 = getattr ( self, "coordinates3", None )
        if coordinates3 is None:
            pruned = None
        else:
            pruned = self.__class__ ( )
            pruned.coordinates3 = coordinates3.Prune ( selection )
        return pruned

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass ( )
        return self

    def SetTemporaryAttribute ( self, label, value ):
        """Set a temporary attribute."""
        if ( label not in self.__class__.permanentAttributes ):
            setattr ( self, label, value )
            self.temporaryAttributes.add ( label )
        else:
            raise AttributeError ( "Cannot make the attribute \"" + label + "\" temporary." )

    @property
    def temporaryAttributes ( self ):
        if getattr ( self, "_temporaryAttributes", None ) is None:
            self._temporaryAttributes = set ( )
        return self._temporaryAttributes

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
