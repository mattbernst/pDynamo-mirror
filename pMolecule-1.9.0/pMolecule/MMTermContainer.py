#-------------------------------------------------------------------------------
# . File      : MMTermContainer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A class for a container holding MM terms."""

from pCore  import RawObjectConstructor, SingleObjectContainer
from MMTerm import MMTerm

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMTermContainer ( SingleObjectContainer ):
    """A container class for MM terms."""

    def __getstate__ ( self ):
        state = {}
        for item in self.items: state[item.label] = item
        return state

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, items ):
        self.__dict__["items"] = items.values ( )

    def ItemClass ( self ): return MMTerm

    def ItemName ( self ): return "MM Term"

    def Merge ( self, others, information = {} ):
        """Merging."""
        # . Initialization.
        merged       = None
        toMergeItems = [ self ] + list ( others )
        # . Need increments.
        increments = information.get ( "atomIncrements", None )
        if ( increments is not None ) and ( len ( increments ) == len ( toMergeItems ) ):
            # . Get unique keys and states.
            keys   = set ( )
            states = []
            for toMerge in toMergeItems:
                state = toMerge.__getstate__ ( )
                keys.update ( state.keys ( ) )
                states.append ( state )
            keys = list ( keys )
            keys.sort ( )
            # . Do merging.
            if len ( keys ) > 0:
                mergedItems = []
                for key in keys:
                    items         = []
                    newIncrements = []
                    for ( state, increment ) in zip ( states, increments ):
                        item = state.get ( key, None )
                        if item is not None:
                            items.append         ( item      )
                            newIncrements.append ( increment )
                    if len ( items ) > 0:
                        mergedItems.append ( items[0].Merge ( items[1:], information = { "atomIncrements" : newIncrements } ) )
                # . Create the object.
                if len ( mergedItems ) > 0: merged = self.__class__ ( mergedItems )
        # . Finish up.
        return merged

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        # . Initialization.
        pruned = None
        # . Prune the items.
        prunedItems = []
        for item in self.items:
            if hasattr ( item, "Prune" ):
                prunedItem = item.Prune ( selection, information = information )
                if prunedItem is not None: prunedItems.append ( prunedItem )
        # . Construct the object.
        if len ( prunedItems ) > 0: pruned = self.__class__ ( prunedItems )
        # . Finish up.
        return pruned

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass ( )
        return self

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
