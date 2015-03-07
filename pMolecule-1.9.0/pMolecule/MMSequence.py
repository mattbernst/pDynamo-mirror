#-------------------------------------------------------------------------------
# . File      : MMSequence.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MM sequence classes."""

from MMModelError import MMModelError

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Key separators.
_KeyTokenSeparator = "_"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMSequenceAtom ( object ):
    """An MM atom type."""

    attributeLabelMappings = { "Charge" : "charge"    ,
                               "Label"  : "label"     ,
                               "Type"   : "typeLabel" }

    defaultAttributes = { "charge"    : None,
                          "label"     : None,
                          "typeLabel" : None }

    def __init__ ( self, **keywordArguments ):
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

    @classmethod
    def DefaultAttributeLabels ( selfClass ):
        return [ "Label", "Type", "Charge" ]

    @classmethod
    def FromList ( selfClass, attributeList, attributes ):
        """Constructor from list."""
        keywordArguments = {}
        for ( key, value ) in zip ( attributes, attributeList ):
            if value is not None: keywordArguments[selfClass.attributeLabelMappings[key]] = value
        self = selfClass ( **keywordArguments )
        return self

    def ToList ( self, attributes ):
        """Get a list representation of the object."""
        selfList = []
        for attribute in attributes:
            selfList.append ( getattr ( self, self.__class__.attributeLabelMappings[attribute], None ) )
        return selfList

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMSequenceComponent ( object ):
    """An MM sequence component."""

    attributeLabelMappings = { "Label" : "label" }

    defaultAttributes = { "atomDictionary" : None ,
                          "atomFields"     : None ,
                          "atoms"          : None ,
                          "label"          : None }

    #yaml_tag = "!MMSequenceComponent"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Atom and bonds.
        for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atoms", "atomFields", "Atoms", "Atom Fields", MMSequenceAtom ), ):
            items = getattr ( self, tag, None )
            if items is not None:
                listItems = []
                keys      = getattr ( self, keyTag, None )
                if keys is None: keys = itemObject.DefaultAttributeLabels ( )
                for item in items:
                    listItems.append ( item.ToList ( keys ) )
                mapping[label   ] = listItems
                mapping[keyLabel] = keys
        # . Other attributes.
        for ( fullKey, key ) in self.__class__.attributeLabelMappings.iteritems ( ):
            attribute = getattr ( self, key, None )
            if attribute is not None: mapping[fullKey] = attribute
        # . Finish up.
        return mapping

    def __init__ ( self, **keywordArguments ):
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            self.__init__ ( )
            # . Atom and bonds.
            for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atoms", "atomFields", "Atoms", "Atom Fields", MMSequenceAtom ), ):
                listItems = mapping.get ( label, None )
                if listItems is not None:
                    items = []
                    keys  = mapping[keyLabel]
                    for listItem in listItems:
                        items.append ( itemObject.FromList ( listItem, keys ) )
                    setattr ( self, tag   , items )
                    setattr ( self, keyTag, keys  )
            # . Other attributes.
            for ( fullKey, key ) in self.__class__.attributeLabelMappings.iteritems ( ):
                attribute = mapping.pop ( fullKey, None )
                if attribute is not None: setattr ( self, key, attribute )
#        except Exception as e:
#            print e[0]
        except:
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def MakeAtomDictionary ( self ):
        """Make an atom dictionary."""
        if self.atomDictionary is None:
            index = {}
            for atom in self.atoms:
                index[atom.label] = atom
            self.atomDictionary = index

    @staticmethod
    def MakeKey ( label ):
        """Make a key."""
        tokens = label.lower ( ).split ( )
        return _KeyTokenSeparator.join ( tokens )

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

    def TypeSequenceComponentAtoms ( self, component, atomCharges = None, atomTypes = None ):
        """Type the atoms of a sequence component."""
        if atomCharges is None: atomCharges = {}
        if atomTypes   is None: atomTypes   = {}
        self.MakeAtomDictionary ( )
        for atom in component.children:
            mmAtom = self.atomDictionary.get ( atom.label, None )
            if mmAtom is not None:
                index = atom.index
                atomCharges[index] = mmAtom.charge
                atomTypes  [index] = mmAtom.typeLabel
        return ( atomTypes, atomCharges )

    # . Key property.
    def __GetKey ( self ):
        """Get a key."""
        if self.__dict__.get ( "key", None ) is None:
            self.__dict__["key"] = self.MakeKey ( self.label )
        return self.__dict__["key"]
    key = property ( __GetKey, None, None, "Key." )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMSequenceLink ( object ):
    """An MM sequence component link."""

    attributeLabelMappings = { "Label"         : "label"        ,
                               "Left Variant"  : "leftVariant"  ,
                               "Right Variant" : "rightVariant" }

    defaultAttributes = { "label"        : None ,
                          "leftVariant"  : None , 
                          "rightVariant" : None }

    #yaml_tag = "!MMSequenceLink"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping         = {}
        previousState   = None
        previousVariant = None
        for ( fullKey, key ) in self.__class__.attributeLabelMappings.iteritems ( ):
            attribute = getattr ( self, key, None )
            if attribute is not None:
                if key.endswith ( "Variant" ):
                    if attribute is previousVariant:
                        state           = previousState
                    else:
                        state           = attribute.__getstate__ ( )
                        previousState   = state
                        previousVariant = attribute
                    attribute = state
                mapping[fullKey] = attribute
        return mapping

    def __init__ ( self, **keywordArguments ):
        self.fullKeyAttributes = {}
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Build the object.
            self.__init__ ( )
            previousState   = None
            previousVariant = None
            for ( fullKey, key ) in self.__class__.attributeLabelMappings.iteritems ( ):
                attribute = mapping.pop ( fullKey, None )
                if attribute is not None:
                    if key.endswith ( "Variant" ):
                        if attribute is previousState:
                            variant = previousVariant
                        else:
                            variant = MMSequenceVariant ( )
                            variant.__setstate__ ( attribute )
                            previousState   = attribute
                            previousVariant = variant
                        attribute = variant
                    setattr ( self, key, attribute )
#        except Exception as e:
#            print e[0]
        except:
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    @staticmethod
    def MakeKeys ( linkLabel, leftComponentLabel, rightComponentLabel ):
        """Make keys."""
        # . Get the labels.
        label1 = _KeyTokenSeparator.join ( leftComponentLabel.lower  ( ).split ( ) )
        label2 = _KeyTokenSeparator.join ( rightComponentLabel.lower ( ).split ( ) )
        # . Get the keys.
        key00 = _KeyTokenSeparator.join ( linkLabel.lower ( ).split ( ) )
        keyX0 = _KeyTokenSeparator.join ( [ key00, label1         ] )
        key0Y = _KeyTokenSeparator.join ( [ key00, ""    , label2 ] )
        key   = _KeyTokenSeparator.join ( [ keyX0,         label2 ] )
        # . Most specific to most general.
        return ( key, keyX0, key0Y, key00 )

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

    # . Key property.
    def __GetKey ( self ):
        """Get a key."""
        if self.__dict__.get ( "key", None ) is None:
            tokens = self.label.lower ( ).split ( )
            label1 = getattr ( self.leftVariant , "componentLabel", None )
            label2 = getattr ( self.rightVariant, "componentLabel", None )
            if label1 is not None: tokens1 = label1.lower ( ).split ( )
            if label2 is not None: tokens2 = label2.lower ( ).split ( )
            if label1 is None:
                if label2 is not None: tokens.extend ( [ "" ] + tokens2 )
            else:
                tokens.extend ( tokens1 )
                if label2 is not None: tokens.extend ( tokens2 )
            self.__dict__["key"] = _KeyTokenSeparator.join ( tokens )
        return self.__dict__["key"]
    key = property ( __GetKey, None, None, "Key." )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMSequenceVariant ( MMSequenceComponent ):
    """An MM sequence component variant."""

    attributeLabelMappings = dict ( MMSequenceComponent.attributeLabelMappings )
    attributeLabelMappings.update ( { "Component" : "componentLabel" } )

    defaultAttributes = dict ( MMSequenceComponent.defaultAttributes )
    defaultAttributes.update ( { "componentLabel" : None } )

    #yaml_tag = "!MMSequenceVariant"

    @staticmethod
    def MakeKeys ( variantLabel, componentLabel ):
        """Make keys."""
        tokens0 = variantLabel.lower   ( ).split ( )
        tokens  = tokens0 + componentLabel.lower ( ).split ( )
        key0    = _KeyTokenSeparator.join ( tokens0 )
        key     = _KeyTokenSeparator.join ( tokens  )
        # . Most specific to most general.
        return ( key, key0 )

    # . Key property.
    def __GetKey ( self ):
        """Get a key."""
        if self.__dict__.get ( "key", None ) is None:
            tokens = self.label.lower ( ).split ( )
            label  = getattr ( self, "componentLabel", None )
            if label is not None: tokens.extend ( label.lower ( ).split ( ) )
            self.__dict__["key"] = _KeyTokenSeparator.join ( tokens )
        return self.__dict__["key"]
    key = property ( __GetKey, None, None, "Key." )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
