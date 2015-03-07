#-------------------------------------------------------------------------------
# . File      : PDBComponent.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""PDB Components."""

import os, os.path

from pCore     import logFile, LogFileActive
from pMolecule import AromaticSingleBond, DoubleBond, SingleBond, TripleBond, UndefinedBond,  \
                      MMSequenceAtom, MMSequenceComponent, MMSequenceLink, MMSequenceVariant, \
                      Sequence, System

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Bond label/object mappings.
_BondLabelObjectMapping = { "Aromatic Single" : AromaticSingleBond ( ) ,
                            "Double"          : DoubleBond         ( ) ,
                            "Single"          : SingleBond         ( ) ,
                            "Triple"          : TripleBond         ( ) ,
                            "Undefined"       : UndefinedBond      ( ) }
_BondObjectLabelMapping = { AromaticSingleBond ( ) : "Aromatic Single" ,
                            DoubleBond         ( ) : "Double"          ,
                            SingleBond         ( ) : "Single"          ,
                            TripleBond         ( ) : "Triple"          ,
                            UndefinedBond      ( ) : "Undefined"       }

# . Key separators.
_KeyTokenSeparator = "_"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponentError ( Exception ):
    pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponent ( object ):
    """A PDB component."""

    attributeLabelMappings = { "Component Class"      : "componentClass"     ,
                               "Formal Charge"        : "formalCharge"       ,
                               "Label"                : "label"              ,
                               "Left Link"            : "leftLink"           ,
                               "Left Termination"     : "leftTermination"    ,
                               "Name"                 : "name"               ,
                               "PDB Class"            : "pdbClass"           ,
                               "Is In Chain"          : "isInChain"          ,
                               "Is Chain Terminating" : "isChainTerminating" ,
                               "Is Heteroatom"        : "isHeteroatom"       ,
                               "Right Link"           : "rightLink"          ,
                               "Right Termination"    : "rightTermination"   ,
                               "Variants"             : "variants"           }

    defaultAttributes = { "atomFields"         : None ,
                          "atoms"              : None ,
                          "bondFields"         : None ,
                          "bonds"              : None ,
                          "componentClass"     : None ,   
                          "formalCharge"       : None ,   
                          "label"              : None ,  
                          "leftLink"           : None ,   
                          "leftTermination"    : None ,   
                          "name"               : None ,   
                          "pdbClass"           : None ,   
                          "isInChain"          : None ,  
                          "isChainTerminating" : None ,  
                          "isHeteroatom"       : None ,   
                          "rightLink"          : None ,   
                          "rightTermination"   : None ,   
                          "variants"           : None }

    #yaml_tag = "!PDBComponent"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Atom and bonds.
        for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atoms", "atomFields", "Atoms", "Atom Fields", PDBComponentAtom ),
                                                              ( "bonds", "bondFields", "Bonds", "Bond Fields", PDBComponentBond ) ):
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
            for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atoms", "atomFields", "Atoms", "Atom Fields", PDBComponentAtom ),
                                                                  ( "bonds", "bondFields", "Bonds", "Bond Fields", PDBComponentBond ) ):
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
            raise PDBComponentError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    @classmethod
    def FromSystem ( selfClass, system, label = None ):
        """Generate a component from a system."""
        # . Atoms.
        atoms        = []
        formalCharge = 0
        labelIndices = {}
        for ( i, atom ) in enumerate ( system.atoms ):
            formalCharge += atom.formalCharge
            atomLabel     = atom.label
            if len ( atomLabel ) > 3: pdbAlign = 0
            else:                     pdbAlign = 1
            labelIndices[i] = atomLabel
            atoms.append ( PDBComponentAtom ( atomicNumber = atom.atomicNumber, formalCharge = atom.formalCharge, label = atomLabel, pdbAlign = pdbAlign ) )
        # . Bonds.
        bonds = []
        if system.connectivity is not None:
            for bond in system.connectivity.bonds:
                bonds.append ( PDBComponentBond ( atomLabel1 = labelIndices[bond.i], atomLabel2 = labelIndices[bond.j], bondType = bond.type ) )
        # . Create the component.
        self = selfClass ( atoms = atoms, bonds = bonds, formalCharge = formalCharge, label = label )
        return self

    @staticmethod
    def MakeKey ( label ):
        """Make a key."""
        tokens = label.lower ( ).split ( )
        return _KeyTokenSeparator.join ( tokens )

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            if self.label is None: summary.Start ( "PDB Component Summary" )
            else:                  summary.Start ( "Summary for PDB Component " + self.label.upper ( ) )
            summary.Entry ( "Number of Atoms", "{:d}".format ( len ( self.atoms ) ) )
            summary.Entry ( "Number of Bonds", "{:d}".format ( len ( self.bonds ) ) )
            summary.Entry ( "Formal Charge",   "{:d}".format ( self.formalCharge  ) )
            summary.Stop ( )

    def ToMMSequenceObject ( self ):
        """Make an MM sequence object."""
        mmAtoms = []
        for atom in self.atoms:
            mmAtoms.append ( MMSequenceAtom ( charge = 0.0, label = atom.label, typeLabel = atom.label ) ) 
        mmItem = MMSequenceComponent ( atoms = mmAtoms, label = self.label )
        return mmItem

    def ToSystem ( self ):
        """Generate a system from the component."""
        # . Bonds.
        labelIndices = {}
        for ( i, atom ) in enumerate ( self.atoms ): labelIndices[atom.label] = i
        bonds = []
        for bond in self.bonds: bonds.append ( ( labelIndices[bond.atomLabel1], labelIndices[bond.atomLabel2], bond.bondType ) )
        # . System from atoms and bonds.
        system = System.FromAtoms ( self.atoms, bonds = bonds, withSequence = True )
        # . Set additional atom data.
        for ( pAtom, sAtom ) in zip ( self.atoms, system.atoms ):
            sAtom.formalCharge = pAtom.formalCharge
            sAtom.label        = pAtom.label
        # . Label.
        system.label = "PDB Component"
        if self.label is not None: system.label += " " + self.label
        # . Define component labels.
        for entity in system.sequence.children:
            for component in entity.children:
                component.label = self.label
        # . Finish up.
        return system

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
class PDBComponentAtom ( object ):
    """Atom in a PDB component."""

    # . To follow is position in component. Only useful for links and variants which can involve changes in order.
    # . Also cannot be used to put atom at beginning of a residue unless use special notation - e.g. toFollow = "start"?

    attributeLabelMappings = { "Atomic Number" : "atomicNumber" ,
                               "Formal Charge" : "formalCharge" ,
                               "Label"         : "label"        ,
                               "PDB Align"     : "pdbAlign"     ,
                               "To Follow"     : "toFollow"     }

    defaultAttributes = { "atomicNumber" :   -1 ,
                          "formalCharge" :    0 ,
                          "label"        : None ,
                          "pdbAlign"     :    0 ,
                          "toFollow"     : None }

    def __init__ ( self, **keywordArguments ):
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

    @classmethod
    def DefaultAttributeLabels ( selfClass ):
        return [ "Atomic Number", "Formal Charge", "Label" ]

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
class PDBComponentBond ( object ):
    """Bond in a PDB component."""

    attributeLabelMappings = { "Atom 1" : "atomLabel1"    ,
                               "Atom 2" : "atomLabel2"    ,
                               "Type"   : "bondTypeLabel" }

    defaultAttributes = { "atomLabel1"    : None ,
                          "atomLabel2"    : None ,
                          "bondType"      : None ,
                          "bondTypeLabel" : None }

    def __init__ ( self, **keywordArguments ):
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )
        self.CheckBond ( )

    def CheckBond ( self ):
        """Check the bond attributes."""
        if   self.bondType      is not None: self.bondTypeLabel = _BondObjectLabelMapping[self.bondType     ]
        elif self.bondTypeLabel is not None: self.bondType      = _BondLabelObjectMapping[self.bondTypeLabel]

    @classmethod
    def DefaultAttributeLabels ( selfClass ):
        return [ "Atom 1", "Atom 2", "Type" ]

    @classmethod
    def FromList ( selfClass, attributeList, attributes ):
        """Constructor from list."""
        keywordArguments = {}
        for ( key, value ) in zip ( attributes, attributeList ):
            if value is not None: keywordArguments[selfClass.attributeLabelMappings[key]] = value
        self = selfClass ( **keywordArguments )
#        if self.bondType is None:
#            print self.atomLabel1, self.atomLabel2, self.bondTypeLabel, keywordArguments
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
class PDBComponentLink ( object ):
    """Non-directional link between two PDB components."""

    attributeLabelMappings = { "Bond Type"     : "bondTypeLabel"  ,
                               "Label"         : "label"          ,
                               "Left Atom"     : "leftAtomLabel"  ,
                               "Left Variant"  : "leftVariant"    ,
                               "Right Atom"    : "rightAtomLabel" ,
                               "Right Variant" : "rightVariant"   }

    defaultAttributes = { "bondType"       : None ,
                          "bondTypeLabel"  : None ,
                          "label"          : None ,
                          "leftAtomLabel"  : None ,
                          "leftVariant"    : None ,
                          "rightAtomLabel" : None ,
                          "rightVariant"   : None }

    #yaml_tag = "!PDBComponentLink"

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
        self.CheckBond  ( )

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
                            variant = PDBComponentVariant ( )
                            variant.__setstate__ ( attribute )
                            previousState   = attribute
                            previousVariant = variant
                        attribute = variant
                    setattr ( self, key, attribute )
            # . Checks.
            self.CheckBond ( )
#        except Exception as e:
#            print e[0]
        except:
            raise PDBComponentError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def CheckBond ( self ):
        """Check the bond attributes."""
        if   self.bondType      is not None: self.bondTypeLabel = _BondObjectLabelMapping[self.bondType     ]
        elif self.bondTypeLabel is not None: self.bondType      = _BondLabelObjectMapping[self.bondTypeLabel]

    @staticmethod
    def MakeKeys ( linkLabel, leftComponentLabel, rightComponentLabel ):
        """Make keys."""
        # . Get the labels.
        label1 = _KeyTokenSeparator.join ( leftComponentLabel.lower  ( ).split ( ) )
        label2 = _KeyTokenSeparator.join ( rightComponentLabel.lower ( ).split ( ) )
        # . Get the keys.
        key00 = _KeyTokenSeparator.join ( linkLabel.lower ( ).split ( ) )
        keyX0 = _KeyTokenSeparator.join ( [ key00, label1        ] )
        key0Y = _KeyTokenSeparator.join ( [ key00, ""   , label2 ] )
        keyXY = _KeyTokenSeparator.join ( [ keyX0,        label2 ] )
        # . Most specific to most general.
        return ( keyXY, keyX0, key0Y, key00 )

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

    def ToMMSequenceObject ( self ):
        """Make an MM sequence object."""
        leftVariant = self.leftVariant.ToMMSequenceObject ( )
        if self.leftVariant is self.rightVariant: rightVariant = leftVariant
        else: rightVariant = self.rightVariant.ToMMSequenceObject ( )
        mmItem = MMSequenceLink ( label = self.label, leftVariant = leftVariant, rightVariant = rightVariant )
        return mmItem

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
class PDBComponentVariant ( object ):
    """Variant of a PDB component."""

    attributeLabelMappings = { "Atoms To Delete"      : "atomsToDelete"      ,
                               "Bonds To Delete"      : "bondsToDelete"      ,
                               "Component"            : "componentLabel"     ,
                               "Formal Charges"       : "formalCharges"      ,
                               "Is Chain Terminating" : "isChainTerminating" ,
                               "Label"                : "label"              }

    defaultAttributes = { "atomsToAddFields"   : None ,
                          "atomsToAdd"         : None , # . List of PDBComponentAtom objects.
                          "atomsToDelete"      : None , # . List of atom labels.
                          "bondsToAddFields"   : None ,
                          "bondsToAdd"         : None , # . List of PDBComponentBond objects.
                          "bondsToDelete"      : None , # . List of lists of atom label pairs.
                          "bondTypesFields"    : None ,
                          "bondTypes"          : None , # . List of PDBComponentBond objects.
                          "componentLabel"     : None ,
                          "formalCharges"      : None , # . Mapping of formal charges to atom labels.
                          "isChainTerminating" : None ,
                          "label"              : None }

    #yaml_tag = "!PDBComponentVariant"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Atom and bonds.
        for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atomsToAdd", "atomsToAddFields", "Atoms To Add", "Atom To Add Fields", PDBComponentAtom ),
                                                              ( "bondsToAdd", "bondsToAddFields", "Bonds To Add", "Bond To Add Fields", PDBComponentBond ),
                                                              ( "bondTypes" , "bondTypesFields" , "Bond Types"  , "Bond Type Fields"  , PDBComponentBond ) ):
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
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            self.__init__ ( )
            # . Atom and bonds.
            for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atomsToAdd", "atomsToAddFields", "Atoms To Add", "Atom To Add Fields", PDBComponentAtom ) ,
                                                                  ( "bondsToAdd", "bondsToAddFields", "Bonds To Add", "Bond To Add Fields", PDBComponentBond ) ,
                                                                  ( "bondTypes" , "bondTypesFields" , "Bond Types"  , "Bond Type Fields"  , PDBComponentBond ) ):
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
            raise PDBComponentError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    @staticmethod
    def MakeKeys ( variantLabel, componentLabel ):
        """Make keys."""
        tokens0 = variantLabel.lower   ( ).split ( )
        tokens  = tokens0 + componentLabel.lower ( ).split ( )
        key0    = _KeyTokenSeparator.join ( tokens0 )
        key     = _KeyTokenSeparator.join ( tokens  )
        # . Most specific to most general.
        return ( key, key0 )

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

    def ToMMSequenceObject ( self ):
        """Make an MM sequence object."""
        mmAtoms = []
        if self.atomsToAdd is not None:
            for atom in self.atomsToAdd:
                mmAtoms.append ( MMSequenceAtom ( charge = 0.0, label = atom.label, typeLabel = atom.label ) ) 
        mmItem = MMSequenceVariant ( atoms = mmAtoms, componentLabel = self.componentLabel, label = self.label )
        return mmItem

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
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
