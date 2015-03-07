#-------------------------------------------------------------------------------
# . File      : MMAtomTypeContainer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MM atom types."""

from MMModelError import MMModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMAtomType ( object ):
    """An MM atom type."""

    defaultAttributes = { "atomicNumber"  :   -1,
                          "charge"        :  0.0,
                          "hydrogenType"  : None,
                          "label"         : None,
                          "properties"    : None }

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    @classmethod
    def FromKeysValues ( selfClass, keys, values ):
        """Constructor from keys and values."""
        self = selfClass ( )
        self.properties = {}
        for ( key, value ) in zip ( keys, values ): self.properties[key] = value
        self.ProcessProperties ( )
        return self

    def ProcessProperties ( self ):
        """Process properties for needed attributes."""
        if self.properties is not None:
            self.atomicNumber = self.properties.get ( "Atomic Number" ,   -1 )
            self.charge       = self.properties.get ( "Charge"        ,  0.0 )
            self.label        = self.properties.get ( "Label"         , None )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMAtomTypeContainer ( object ):
    """A container for MM atom types."""

    defaultAttributes = { "items"           :            None,
                          "parameterFields" :            None,
                          "properties"      :            None,
                          "label"           :            None,
                          "rawItems"        :            None,
                          "termLabel"       : "MM Atom Type" }

    #yaml_tag = "!MMAtomTypeContainer"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Columns.
        columnKeys                  = self.parameterFields
        mapping["Parameter Fields"] = columnKeys
        # . Rows.
        rowKeys = self.rawItems.keys ( )
        rowKeys.sort ( )
        rows = []
        for rowKey in rowKeys:
            row = []
            for columnKey in columnKeys:
                row.append ( self.rawItems[rowKey].properties[columnKey] )
            rows.append ( row )
        mapping["Parameter Values"] = rows
        # . Other data.
        if self.label is not None: mapping["Label"] = self.label
        return mapping

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Basic construction.
            self.__init__ ( )
            columns              = mapping.pop ( "Parameter Fields" )
            rows                 = mapping.pop ( "Parameter Values" )
            self.label           = mapping.pop ( "Label", None )
            self.parameterFields = columns
            self.properties      = dict ( mapping )
            # . Get token indices.
            l = columns.index ( "Label" )
            # . Create the raw items.
            self.rawItems = {}
            for row in rows:
                key = row[l]
                self.rawItems[key] = MMAtomType.FromKeysValues ( columns, row )
            if len ( self.rawItems ) != len ( rows ): raise
            # . Create the items that are to be used.
            self.ProcessRawItems ( )
#        except Exception as e:
#            print e[0]
        except:
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def GetItem ( self, label ):
        """Get an item given its label."""
        return self.items.get ( label, None )

    def ProcessRawItems ( self ):
        """Process the raw items."""
        self.items = self.rawItems
        if self.items is not None:
            if "Hydrogen Type" in self.parameterFields:
                for item in self.items.values ( ):
                    hydrogenLabel     = item.properties["Hydrogen Type"]
                    item.hydrogenType = self.GetItem ( hydrogenLabel )

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
