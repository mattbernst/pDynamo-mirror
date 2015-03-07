#-------------------------------------------------------------------------------
# . File      : Merge.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Functions for merging chemical objects."""

#===================================================================================================================================
# . Class.
#===================================================================================================================================
def MergeByAtom ( *arguments ):
    """Merge the arguments."""
    merged = None
    if len ( arguments ) > 0:
        if isinstance ( arguments[0], ( list, tuple ) ): items = arguments[0]
        else:                                            items = arguments
        if hasattr ( items[0], "Merge" ): return items[0].Merge ( items[1:] )
        else: raise ValueError ( "Cannot merge objects of type {!r}.".format ( type ( object ) ) )
    return merged

#===================================================================================================================================
# . Class.
#===================================================================================================================================
def MergeRepeatByAtom ( item, repeat, **keywordArguments ):
    """Return a merged object consisting of |repeat| copies of |item|."""
    merged = None
    if repeat <= 0:
        raise ValueError ( "Invalid repeat argument." )
    elif repeat == 1:
        merged = item
    elif hasattr ( item, "MergeRepeat" ):
        merged = item.MergeRepeat ( repeat, information = keywordArguments )
    elif hasattr ( item, "Merge" ):
        merged = item.Merge ( ( repeat - 1 ) * [ item ], information = keywordArguments )
    else:
        raise ValueError ( "Cannot merge objects of type {!r}.".format ( type ( object ) ) )
    return merged

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
