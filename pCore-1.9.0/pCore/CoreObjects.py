#-------------------------------------------------------------------------------
# . File      : CoreObjects.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Definitions of the core object types used in pCore and dependent packages."""

#===================================================================================================================================
# . Classes for lists containing restricted types of object.
#===================================================================================================================================
class CLibraryError ( Exception ): pass

#===================================================================================================================================
# . Classes for properties.
#===================================================================================================================================
class pObjectProperty ( property ):
    """Properties for objects of arbitrary type."""

    def __init__( self, name, objtype, default = None, doc = "" ):
        """Constructor.

        |objtype| can be a single class or a tuple of classes.
        """
        self.default = default
        self.name    = name
        self.objtype = objtype
        self.__doc__ = doc

    def __delete__( self, obj ):
        """Deleting."""
        obj.__dict__.pop ( self.name, None )

    def __get__ ( self, obj, objtype = None ):
        """Getting."""
        try:    return obj.__dict__[self.name]
        except: return self.default
#        print obj, self.name, objtype
#        item = obj.__dict__.get ( self.name )
#        if item is None: return self.default
#        else:            return item

    def __set__( self, obj, value ):
        """Setting."""
        if isinstance ( value, self.objtype ): obj.__dict__[self.name] = value
        else: raise AttributeError ( "Wrong object type. Expecting {!r} but got {!r}.".format ( self.objtype, type ( value ) ) )

#===================================================================================================================================
# . Classes for standard container types.
#===================================================================================================================================
class SingleObjectContainer ( object ):
    """A basic container class whose contents can be set only once and which only contains one class of object."""

    def __delitem__ ( self, i ):
        """Delete an item."""
        raise AttributeError ( "Cannot delete a container item." )

    def __getitem__ ( self, i ):
        """Get an item."""
        if isinstance ( i, int ):
            try   : return self.items[i]
            except: raise IndexError ( "Container index out of range." )
        else: raise TypeError ( "Container index must be integer." )

    def __init__ ( self, *arguments, **keywordArguments ):
        """Constructor."""
        # . Items.
        self.__dict__["items"] = None
        # . Initialization from sequence.
        if len ( arguments ) > 0: self.items = arguments[0]

    def __len__ ( self ):
        """Length."""
        return len ( self.items )

    def __setitem__ ( self, i, value ):
        """Set an item."""
        raise AttributeError ( "Cannot reset a container item." )

    def ItemClass ( self ): return object

    def ItemName ( self ): return "Object"

    def SummaryEntry ( self, summary ):
        """Summary entry."""
        if summary is not None: summary.Entry ( "Number of " + self.ItemName ( ) + "s", "{:d}".format ( len ( self ) ) )

    # . Properties.
    # . Items.
    def __DelItems ( self ):
        """Deletion not possible."""
        raise AttributeError ( "Deletion of container items is not possible." )
    def __GetItems ( self ):
        """Return a list."""
        items = self.__dict__["items"]
        if items is None: items = []
        return items
    def __SetItems ( self, items ):
        """The item list can only be set once."""
        if self.__dict__["items"] is None:
            try:
                for item in items:
                    if not isinstance ( item, self.ItemClass ( ) ):
                        raise TypeError ( "Wrong container item type. Expecting {!r} but got {!r}.".format ( self.ItemClass ( ), type ( item ) ) )
                self.__dict__["items"] = list ( items )
            except:
                raise TypeError ( "A non-iterable object is being used to set container items." )
        else:
            raise AttributeError ( "Container items can only be set once." )
    items = property ( __GetItems, __SetItems, __DelItems, "Container Items." )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
