#-------------------------------------------------------------------------------
# . File      : SoftConstraintContainer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Class for holding soft constraints."""

from pCore          import logFile, LogFileActive, RawObjectConstructor
from SoftConstraint import SoftConstraint

#===============================================================================
# . Class.
#===============================================================================
class SoftConstraintContainer ( dict ):
    """A class defining a container for soft constraints."""

    def __getstate__ ( self ): return dict ( self )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ): self.update ( state )

    def __setitem__ ( self, k, value ):
        """Set an item."""
        if isinstance ( k, basestring ) and isinstance ( value, SoftConstraint ):
            super ( SoftConstraintContainer, self ).__setitem__ ( k, value )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy and, optionally, the gradients."""
        scstate = {}
        for sc in self: scstate[sc] = self[sc].Energy ( coordinates3, gradients3 )
        esc     = 0.0
        for sc in scstate: esc += scstate[sc][0]
        return ( ( "Soft Constraint", esc ), scstate )

    def Merge ( self, others, information = {} ):
        """Merging."""
        # . No checking for duplicate keys.
        # . Initialization.
        merged        = None
        toMergeItems  = [ self ] + list ( others )
        numberToMerge = len ( toMergeItems )
        # . Need increments.
        increments = information.get ( "atomIncrements", None )
        if ( increments is not None ) and ( len ( increments ) == numberToMerge ):
            mergedItems = {}
            for ( toMerge, increment ) in zip ( toMergeItems, increments ):
                for ( key, constraint ) in toMerge.iteritems ( ):
                    mergedItems[key] = constraint.Merge ( atomIncrement )
            # . Construct the object.
            if len ( mergedItems ) > 0: merged = self.__class__ ( mergedItems )
        # . Finish up.
        return merged

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        # . Initialization.
        pruned = None
        # . Prune the items.
        prunedItems = {}
        for ( key, item ) in self.iteritems ( ):
            if hasattr ( item, "Prune" ):
                prunedItem = item.Prune ( selection )
                if prunedItem is not None: prunedItems[key] = prunedItem
        # . Construct the object.
        if len ( prunedItems ) > 0: pruned = self.__class__ ( prunedItems )
        # . Finish up.
        return pruned

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            frequencies = {}
            for key in self:
                sclabel = self[key].Label ( )
                if sclabel in frequencies: frequencies[sclabel] += 1
                else:                      frequencies[sclabel]  = 1
            keys = frequencies.keys ( )
            keys.sort ( )
            summary = log.GetSummary ( )
            summary.Start ( "SoftConstraint Container Summary" )
            for key in keys: summary.Entry ( key, "{:d}".format ( frequencies[key] ) )
            summary.Stop ( )

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :

    pass
