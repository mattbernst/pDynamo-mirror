#-------------------------------------------------------------------------------
# . File      : Prune.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Functions for pruning chemical objects.

Easier to use the Prune method of the appropriate object directly.
"""

#===================================================================================================================================
# . Class.
#===================================================================================================================================
def PruneByAtom ( item, selection ):
    """Prune |item| by |selection|."""
    if hasattr ( item, "Prune" ): return item.Prune ( selection )
    else: raise ValueError ( "Cannot prune object of type {!r}.".format ( type ( object ) ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
