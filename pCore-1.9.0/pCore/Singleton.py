#-------------------------------------------------------------------------------
# . File      : Singleton.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A singleton class from the Python website and a null object class modified from the Python cookbook (originally by Dinu C. Gherman)."""

#===============================================================================
# . Class.
#===============================================================================
class Singleton ( object ):
    """A singleton class.

    Here singletons are used as definition objects whose internal state does not change.
    """

    def __init__( self ):
        pass

    def __new__ ( selfClass, *arguments, **keywordArguments ):
        """Allocate the singleton."""
        it = selfClass.__dict__.get ( "__it__" )
        if it is None:
            selfClass.__it__ = it = object.__new__ ( selfClass )
            it.Initialize ( *arguments, **keywordArguments )
        return it

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( self.__class__, ( ) )

    def Initialize ( self, *arguments, **keywordArguments ):
        """Initialization - can use to set variables."""
        pass

#===============================================================================
# . Class.
#===============================================================================
# . This is not used for the moment.
class Null ( Singleton ):
    """A class implementing the null object pattern."""

    # . More methods need to be added.
    def __bool__ ( self ):
        return False

    def __call__ ( self, *arguments, **keywordArguments ):
        return self

    def __delattr__ ( self, name ):
        return self

    def __eq__ ( self, other ):
        if isinstance ( other, Null ): return True
        else:                          return False

    def __ge__ ( self, other ): return True

    def __getattr__ ( self, name ):
        return self

    def __getnewarguments__ ( self, ): return ( )

    # . This needed for pickling but it should be unnecessary.
    def __getmodule__ ( self ):
        return "pCore.Singleton"

    def __gt__ ( self, other ):
        if isinstance ( other, Null ): return False
        else:                          return True

    # . This needs looking at.
    def __hash__ ( self ): return -1

    def __init__ ( self, *arguments, **keywordArguments ):
        super ( Null, self ).__init__ ( *arguments, **keywordArguments )
#        return None

    def __iter__ ( self ):
        return self

    def __le__ ( self, other ):
        if isinstance ( other, Null ): return True
        else:                          return False

    def __len__ ( self ):
        return 0

    def __lt__ ( self, other ): return False

    def __ne__ ( self, other ): return ( self == other )

    def __next__ ( self ): return self.next ( )

    __nonzero__ = __bool__

    def __repr__ ( self ):
        return "<Null>"

    def __setattr__ ( self, name, value ):
        return self

    def __str__ ( self ):
        return "Null"

    def next ( self ): raise StopIteration

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
