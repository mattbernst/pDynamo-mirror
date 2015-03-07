#-------------------------------------------------------------------------------
# . File      : PDBModel.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""PDB model classes."""

import operator, os.path

from pCore               import Coordinates3, logFile, LogFileActive, Real1DArray, Selection, TreeBranchNode, TreeLeafNode, TreeRootNode, YAMLPickle, YAMLUnpickle
from pMolecule           import Atom, Bond, Sequence, SequenceComponent, SequenceEntity, SequenceLinearPolymer, SequenceLink, SequenceVariant, System
from PDBComponentLibrary import PDBComponentLibrary

#
# . Notes:
#
#   Proper pickle format for model (for pDesktop)?
#
#   Symmetry operations could be added as a fourth level, e.g. s:e:c:a.
#   It appears that SegID is no longer present in the official PDB format.
#
#   Use table format for linear polymer, link and variant data in model file?
#

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . AltLoc data.
_DefaultAltLoc = ""

# . Link data.
_DefaultLinkLabel = "GenericSingle"

# . Value for undefined coordinate entries.
_PDBAtomDataUndefined = 9999.0

#===================================================================================================================================
# . PDB Atoms.
#===================================================================================================================================
class PDBModelAtom ( TreeBranchNode ):
    """A class to represent a PDB atom."""

    defaultAttributes = { "atomicNumber" : -1 ,
                          "formalCharge" :  0 }
    defaultAttributes.update ( TreeBranchNode.defaultAttributes )

    @classmethod
    def FromPDBComponentAtom ( selfClass, atom, **keywordArguments ):
        """Constructor from a PDBComponentAtom."""
        self = selfClass ( **keywordArguments )
        for attribute in ( "atomicNumber", "formalCharge", "label" ): setattr ( self, attribute, getattr ( atom, attribute ) )
        return self

    def GetData ( self, altLoc = None ):
        """Return data for an atom."""
        defaultData = self.childIndex.get ( _DefaultAltLoc, None )
        if altLoc is None:
            # . Get data of the highest occupancy.
            if defaultData is None:
                maximumDatum     = None
                maximumOccupancy = - _PDBAtomDataUndefined
                for datum in self.children:
                    if datum.occupancy > maximumOccupancy:
                        maximumDatum     = datum
                        maximumOccupancy = datum.occupancy
                return maximumDatum
            # . Return the default data.
            else:
                return defaultData
        else:
            # . Return the data with the given altLoc or, if absent, the default data.
            return self.childIndex.get ( altLoc, defaultData )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBModelAtomData ( TreeLeafNode ):
    """Store coordinate data about an atom."""

    defaultAttributes = { "charge"            : _PDBAtomDataUndefined ,
                          "occupancy"         : _PDBAtomDataUndefined ,
                          "radius"            : _PDBAtomDataUndefined ,
                          "temperatureFactor" : _PDBAtomDataUndefined ,
                          "u00"               : _PDBAtomDataUndefined ,
                          "u11"               : _PDBAtomDataUndefined ,
                          "u22"               : _PDBAtomDataUndefined ,
                          "u01"               : _PDBAtomDataUndefined ,
                          "u02"               : _PDBAtomDataUndefined ,
                          "u12"               : _PDBAtomDataUndefined ,
                          "x"                 : _PDBAtomDataUndefined ,
                          "y"                 : _PDBAtomDataUndefined ,
                          "z"                 : _PDBAtomDataUndefined }
    defaultAttributes.update ( TreeLeafNode.defaultAttributes )

    def CopyTo ( self, other ):
        """Copy data to another datum."""
        for attribute in self.__class__.defaultAttributes:
            setattr ( other, attribute, getattr ( self, attribute, _PDBAtomDataUndefined ) )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBModelBond ( object ):
    """A class to represent a PDB bond."""

    defaultAttributes = { "atom1"    : None ,
                          "atom2"    : None ,
                          "bondType" : None ,
                          "_key"     : None }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

    @classmethod
    def FromPDBComponentBond ( selfClass, bond, component1, component2 ):
        """Constructor from a PDBComponentBond."""
        atom1 = component1.childIndex.get ( bond.atomLabel1, None )
        atom2 = component2.childIndex.get ( bond.atomLabel2, None )
        if ( atom1 is None ) or ( atom2 is None ):
            raise IndexError ( "Unknown atom in bond: " + bond.atomLabel1 + "-" + bond.atomLabel2 + "." )
        self  = selfClass ( atom1 = atom1, atom2 = atom2, bondType = bond.bondType )
        return self

    @staticmethod
    def MakeKey ( label1, label2 ):
        """Get a bond key from labels."""
        return ( max ( label1, label2 ), min ( label1, label2 ) )

    @property
    def key ( self ):
        """The bond key."""
        if self._key is None:
            self._key = PDBModelBond.MakeKey ( self.atom1.path, self.atom2.path )
        return self._key

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBModelComponent ( TreeBranchNode ):
    """A class to represent a PDB model component."""

    defaultAttributes = { "bondIndex" : None ,
                          "bonds"     : None }
    defaultAttributes.update ( TreeBranchNode.defaultAttributes )

    def AddBond ( self, bond ):
        """Add a bond to the component."""
        if bond is not None:
            if isinstance ( bond, PDBModelBond ): self.bonds.append ( bond )
            else:                                 raise TypeError  ( "Invalid bond being appended to PDB component: {!r}.".format ( type ( bond ) ) )
            if bond.key in self.bondIndex:        raise ValueError ( "Duplicate bond in component: " + bond.atom1.label + "-" + bond.atom2.label + "." )
            else:                                 self.bondIndex[bond.key] = bond

    def ApplyLibraryVariant ( self, libraryVariant ):
        """Apply a library variant."""
        # . Bonds to delete.
        items = getattr ( libraryVariant, "bondsToDelete", None )
        if items is not None:
            for ( atomLabel1, atomLabel2 ) in items: self.DeleteBond ( atomLabel1, atomLabel2 )
        # . Atoms to delete.
        items = getattr ( libraryVariant, "atomsToDelete", None )
        if items is not None:
            for label in items: self.DeleteAtom ( label )
        # . Atoms to add.
        items = getattr ( libraryVariant, "atomsToAdd", None )
        if items is not None:
            for atom in items: self.AddChild ( PDBModelAtom.FromPDBComponentAtom ( atom ), toFollow = atom.toFollow )
        # . Bonds to add.
        items = getattr ( libraryVariant, "bondsToAdd", None )
        if items is not None:
            for bond in items: self.AddBond ( PDBModelBond.FromPDBComponentBond ( bond, self, self ) )
        # . Formal charges to change.
        items = getattr ( libraryVariant, "formalCharges", None )
        if items is not None:
            for ( atomLabel, formalCharge ) in items.iteritems ( ):
                atom = self.childIndex.get ( atomLabel, None )
                if atom is not None: atom.formalCharge = formalCharge
        # . Bond orders to change.
        items = getattr ( libraryVariant, "bondTypes", None )
        if items is not None:
            for item in items:
                bond = self.GetBondFromAtomLabels ( item.atomLabel1, item.atomLabel2 )
                if bond is not None: bond.bondType = item.bondType

    def ClearAtoms ( self ):
        """Clear atoms and bonds from the component."""
        self.childIndex = {}
        self.children   = []
        self.bondIndex  = {}
        self.bonds      = []

    def DeleteAtom ( self, atomLabel ):
        """Delete an atom and its associated bonds given a atomLabel."""
        atom = self.childIndex.get ( atomLabel, None )
        if atom is not None:
            self.children.remove ( atom )
            del self.childIndex[atomLabel]
            self.DeleteBondsToAtom ( atom )

    def DeleteBond ( self, atomLabel1, atomLabel2 ):
        """Delete a bond given two atom labels."""
        bond = self.GetBondFromAtomLabels ( atomLabel1, atomLabel2 )
        if bond is not None:
            self.bonds.remove ( bond )
            del self.bondIndex[bond.key]

    def DeleteBondsToAtom ( self, atom ):
        """Delete bonds involving an atom."""
        toRemove = []
        for bond in self.bonds:
            if ( bond.atom1 is atom ) or ( bond.atom2 is atom ):
                toRemove.append ( bond )
                del self.bondIndex[bond.key]
        for bond in toRemove: self.bonds.remove ( bond )

    def GetBondFromAtomLabels ( self, atomLabel1, atomLabel2 ):
        """Get a bond given two atom labels."""
        bond = None
        if len ( self.bonds ) > 0:
            atom1 = self.childIndex.get ( atomLabel1, None )
            atom2 = self.childIndex.get ( atomLabel2, None )
            if ( atom1 is not None ) and ( atom2 is not None ):
                bond = self.bondIndex.get ( PDBModelBond.MakeKey ( atom1.path, atom2.path ), None )
        return bond

    def Initialize ( self ):
        """Initialization."""
        super ( PDBModelComponent, self ).Initialize ( )
        self.bondIndex = {}
        self.bonds     = []

    def OrderHydrogens ( self, embedded = False ):
        """Order hydrogens either after the atoms to which they are bound or at end of the component. Unbound hydrogens are always placed at the end."""
        # . Find an index of hydrogens - needed for hydrogens with zero or two or more bonds to heavy atoms.
        hydrogenConnections  = {}
        hydrogenIndex        = {}
        for atom in self.children:
            if atom.atomicNumber == 1: hydrogenIndex      [atom.label] = atom
            else:                      hydrogenConnections[atom.label] = []
        # . Find hydrogens bound to each non-hydrogen atom. Bonds between hydrogens are ignored.
        for bond in self.bonds:
            atomicNumber1 = bond.atom1.atomicNumber
            atomicNumber2 = bond.atom2.atomicNumber
            if   ( atomicNumber1 == 1 ) and ( atomicNumber2 != 1 ): hydrogenConnections[bond.atom2.label].append ( bond.atom1.label )
            elif ( atomicNumber1 != 1 ) and ( atomicNumber2 == 1 ): hydrogenConnections[bond.atom1.label].append ( bond.atom2.label )
        # . Create a new atom list.
        # . Heavy atoms first.
        atoms = []
        if embedded: allAtoms = atoms
        else:        allAtoms = []
        for atom in self.children:
            if atom.atomicNumber != 1:
                atoms.append ( atom )
                hydrogenLabels = hydrogenConnections.pop ( atom.label )
                hydrogenLabels.sort ( )
                for hydrogenLabel in hydrogenLabels:
                    if hydrogenLabel in hydrogenIndex:
                        allAtoms.append ( hydrogenIndex.pop ( hydrogenLabel ) )
        # . Remaining hydrogens.
        hydrogenLabels = hydrogenIndex.keys ( )
        hydrogenLabels.sort ( )
        for hydrogenLabel in hydrogenLabels: allAtoms.append ( hydrogenIndex[hydrogenLabel] )
        # . Append the two lists for the non-embedded option.
        if not embedded: atoms.extend ( allAtoms )
        # . Check everything is OK.
        if len ( atoms ) == len ( self.children ): self.children = atoms
        else: raise ValueError ( "Logic error in hydrogen ordering." )

#===================================================================================================================================
# . Class for entities.
#===================================================================================================================================
class PDBModelEntity ( TreeBranchNode ):
    """A set of PDB model components."""

    defaultAttributes = { "linkIndex"          : None ,
                          "nonDefaultLinks"    : None ,
                          "nonDefaultVariants" : None ,
                          "variantIndex"       : None }
    defaultAttributes.update ( TreeBranchNode.defaultAttributes )

    def AddLink ( self, link ):
        """Add a link to the entity."""
        if isinstance ( link, PDBModelLink ) and ( not link.isExtraEntity ):
            if link.key not in self.linkIndex:
                self.linkIndex[link.key] = link
                self.nonDefaultLinks.append ( link )

    def AddVariant ( self, variant ):
        """Add a variant to the entity."""
        if isinstance ( variant, PDBModelVariant ):
            if variant.key not in self.variantIndex:
                self.variantIndex[variant.key] = variant
                self.nonDefaultVariants.append ( variant )

    @classmethod
    def FromModelFileMapping ( selfClass, mapping, root ):
        """Constructor from a model file mapping."""
        label = mapping.get ( "Label", root.defaultLabel )
        self  = selfClass ( label = label )
        root.AddChild ( self ) # . Entity needs to be rooted before calculate paths.
        # . Components.
        for componentLabel in mapping.get ( "Components", [] ):
            genericLabel = root.ParseLabel ( componentLabel )[0]
            component    = PDBModelComponent ( genericLabel = genericLabel, label = componentLabel )
            self.AddChild ( component )
        # . Non-default links.
        for subMapping in mapping.get ( "Links", [] ):
            self.AddLink ( PDBModelLink.FromModelFileMapping ( subMapping, self, root, usePath = False ) )
        # . Non-default variants.
        for subMapping in mapping.get ( "Variants", [] ):
            self.AddVariant ( PDBModelVariant.FromModelFileMapping ( subMapping, self ) )
        # . Finish up.
        return self

    def GatherLinks ( self, library, missingItems ):
        """Gather the links for the entity."""
        return list ( self.nonDefaultLinks )

    def GatherVariants ( self, library, missingItems ):
        """Gather the variants for the entity."""
        # . Get the non-default variants by component.
        variantDictionary = {}
        for variant in self.nonDefaultVariants:
            toAdd = variantDictionary.pop ( variant.component, [] )
            toAdd.append ( variant )
            variantDictionary[variant.component] = toAdd
        # . Loop over components.
        variants = []
        for component in self.children:
            toAdd = variantDictionary.pop ( component, None )
            if toAdd is None:
                libraryComponent = library.GetComponent ( component.genericLabel, missingItems = missingItems )
                if ( libraryComponent is not None ) and ( libraryComponent.variants is not None ):
                    for libraryVariant in libraryComponent.variants:
                        variants.append ( PDBModelVariant ( component, label = libraryVariant ) )
            else: variants.extend ( toAdd )
        return variants

    def Initialize ( self ):
        """Initialization."""
        super ( PDBModelEntity, self ).Initialize ( )
        self.linkIndex          = {}
        self.nonDefaultLinks    = []
        self.nonDefaultVariants = []
        self.variantIndex       = {}

    def MakeAtomicModelFromComponentLibrary ( self, library, usedComponents, missingItems ):
        """Make an atomic model for the entity using PDB components."""
        # . Loop over components.
        for component in self.children:
            libraryComponent = library.GetComponent ( component.genericLabel, missingItems = missingItems )
            if libraryComponent is not None:
                usedComponents.add ( libraryComponent.key )
                for atom in libraryComponent.atoms: component.AddChild ( PDBModelAtom.FromPDBComponentAtom ( atom ) )
                bonds = libraryComponent.bonds
                if bonds is not None:
                    for bond in bonds: component.AddBond ( PDBModelBond.FromPDBComponentBond ( bond, component, component ) )

    def PrintComponentFrequency ( self, log = logFile, title = None ):
        """Print the component sequence of the entity."""
        if LogFileActive ( log ):
            # . Frequencies.
            frequencies = {}
            for component in self.children:
                label = component.genericLabel
                if label in frequencies: frequencies[library] += 1
                else:                    frequencies[library]  = 1
            keys = frequencies.keys ( )
            keys.sort ( )
            # . Printing.
            length = min ( 10, len ( frequencies ) )
            table  = log.GetTable ( columns = length * [ 6, 6 ] )
            table.Start ( )
            if title is None: table.Title ( self.Label ( ) + " Component Frequency" )
            else:             table.Title ( title                                 )
            for key in keys:
                table.Entry ( key                                , alignment = "r" )
                table.Entry ( "{:d}".format ( frequencies[key] ) , alignment = "r" )
            table.Stop ( )

    def PrintComponentSequence ( self, log = logFile, title = None ):
        """Print the component sequence of the entity."""
        if LogFileActive ( log ):
            length = min ( 8, len ( self.children ) )
            table  = log.GetTable ( columns = length * [ 15 ] )
            table.Start ( )
            if title is None: table.Title ( self.Label ( ) + " Component Sequence" )
            else:             table.Title ( title )
            table.Title ( title )
            for component in self.children:
                table.Entry ( component.label, alignment = "r" )
            table.Stop ( )

    def ToModelFileMapping ( self ):
        """Get a mapping suitable for a model file."""
        mapping = {}
        # . Label.
        if ( self.label is not None ) and ( len ( self.label ) > 0 ): mapping["Label"] = self.label
        # . Components.
        items = []
        for item in self.children: items.append ( item.label )
        if len ( items ) > 0: mapping["Components"] = items
        # . Non-default links.
        items = []
        for item in self.nonDefaultLinks: items.append ( item.ToModelFileMapping ( usePath = False ) )
        if len ( items ) > 0: mapping["Links"] = items
        # . Non-default variants.
        items = []
        for item in self.nonDefaultVariants: items.append ( item.ToModelFileMapping ( ) )
        if len ( items ) > 0: mapping["Variants"] = items
        return mapping

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBModelLinearPolymer ( object ):
    """A class to represent a linear polymer."""

    defaultAttributes = { "isCyclic"               : False ,
                          "leftTerminalComponent"  : None  ,
                          "leftTermination"        : None  ,
                          "linkIndex"              : None  ,
                          "nonDefaultLinks"        : None  ,
                          "rightTerminalComponent" : None  ,
                          "rightTermination"       : None  }

    def __init__ ( self, **keywordArguments ):
        """Constuctor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )
        self.CheckTerminations ( )
        self.linkIndex       = {}
        self.nonDefaultLinks = []

    def AddLink ( self, link ):
        """Add a non-default link to the polymer."""
        if isinstance ( link, PDBModelLink ):
            if link.componentKey not in self.linkIndex:
                self.linkIndex[link.componentKey] = link
                self.nonDefaultLinks.append ( link )

    def CheckTerminations ( self ):
        """Check the terminations."""
        # . Checks.
        if ( self.leftTerminalComponent is None ) or ( self.rightTerminalComponent is None ):
            raise ValueError ( "Linear polymer has undefined terminal components." )
        leftEntity  = self.leftTerminalComponent.parent
        rightEntity = self.rightTerminalComponent.parent
        if ( leftEntity is None ) or ( rightEntity is None ):
            raise ValueError ( "Linear polymer terminal components with undefined entities." )
        if ( leftEntity is not rightEntity ):
            raise ValueError ( "Linear polymer cannot be defined straddling different entities." )
        leftIndex  = leftEntity.children.index  ( self.leftTerminalComponent  )
        rightIndex = rightEntity.children.index ( self.rightTerminalComponent )
        if leftIndex > rightIndex:
            raise ValueError ( "Linear polymer with its left terminal component to the right of its right terminal component." )
        # . Save data.
        self._entity     = leftEntity
        self._leftIndex  = leftIndex
        self._rightIndex = rightIndex

    @classmethod
    def FromModelFileMapping ( selfClass, mapping, root ):
        """Constructor from a model file mapping."""
        leftTerminalComponent  = root.GetDescendantFromPath ( mapping["Left Terminal Component" ] )
        rightTerminalComponent = root.GetDescendantFromPath ( mapping["Right Terminal Component"] )
        # . Terminations.
        isCyclic         = mapping.get ( "Cyclic Termination", False )
        leftTermination  = mapping.get ( "Left Termination"  , None  )
        rightTermination = mapping.get ( "Right Termination" , None  )
        # . Object.
        self = selfClass ( isCyclic = isCyclic, leftTerminalComponent  = leftTerminalComponent , leftTermination  = leftTermination  ,
                                                rightTerminalComponent = rightTerminalComponent, rightTermination = rightTermination )
        # . Non-default links.
        for subMapping in mapping.get ( "Polymer Links", [] ):
            self.AddLink ( PDBModelLink.FromModelFileMapping ( subMapping, None, root, usePath = True ) )
        # . Finish up.
        return self

    def GetLink ( self, leftComponent, rightComponent, library, missingItems ):
        """Get a polymer link between neighboring components."""
        # . Initialization.
        link = None
        # . Get components.
        leftLibraryComponent  = library.GetComponent ( leftComponent.genericLabel.lower  ( ), missingItems = missingItems )
        rightLibraryComponent = library.GetComponent ( rightComponent.genericLabel.lower ( ), missingItems = missingItems )
        if ( leftLibraryComponent is not None ) and ( rightLibraryComponent is not None ):
            # . Try for a non-default link.
            link = self.linkIndex.get ( PDBModelLink.MakeKey ( leftComponent.label, rightComponent.label ), None )
            if link is None:
                rightLink = leftLibraryComponent.rightLink
                leftLink  =  rightLibraryComponent.leftLink
            # . Incompatible or unknown links.
            if ( leftLink is None ) or ( rightLink is None ) or ( leftLink != rightLink ):
                if leftLink  is None: leftLink  = "Unknown"
                if rightLink is None: rightLink = "Unknown"
                raise ValueError ( "Incompatible or unknown polymer link types: " + leftLibraryComponent.path  + ":" + rightLinkLabel + \
                                                                         " with " + rightLibraryComponent.path + ":" + leftLinkLabel  + "." )
            # . Links the same.
            else:
                libraryLink = library.GetLink ( leftLink, leftComponent.genericLabel, rightComponent.genericLabel, missingItems = missingItems )
                if libraryLink is not None: link = PDBModelLink ( leftComponent, rightComponent, label = leftLink )
        return link

    def GatherLinks ( self, library, missingItems ):
        """Gather the links for the polymer."""
        links  = []
        # . Polymer links.
        entity = self._entity
        for i in range ( self._leftIndex, self._rightIndex ):
            link = self.GetLink ( entity.children[i], entity.children[i+1], library, missingItems )
            if link is not None: links.append ( link )
        # . Cyclic termination.
        if self.isCyclic:
            link = self.GetLink ( self.rightTerminalComponent, self.leftTerminalComponent, library, missingItems )
            if link is not None: links.append ( link )
        return links

    def GatherVariants ( self, library, missingItems ):
        """Gather the variants for the polymer."""
        # . Process left and right terminations.
        variants = []
        if not self.isCyclic:
            leftLibraryComponent  = library.GetComponent ( self.leftTerminalComponent.genericLabel , missingItems = missingItems )
            if self.leftTermination  is None: self.leftTermination  = leftLibraryComponent.leftTermination
            rightLibraryComponent = library.GetComponent ( self.rightTerminalComponent.genericLabel, missingItems = missingItems )
            if self.rightTermination is None: self.rightTermination = rightLibraryComponent.rightTermination
            if self.leftTermination  is not None: variants.append ( PDBModelVariant ( self.leftTerminalComponent , label = self.leftTermination  ) )
            if self.rightTermination is not None: variants.append ( PDBModelVariant ( self.rightTerminalComponent, label = self.rightTermination ) )
        return variants

    def ToModelFileMapping ( self ):
        """Get a mapping suitable for a model file."""
        mapping = { "Left Terminal Component"  : self.leftTerminalComponent.path  ,
                    "Right Terminal Component" : self.rightTerminalComponent.path }
        # . Terminations.
        if self.isCyclic:
            mapping["Cyclic Termination"] = True
        else:
            if self.leftTermination  is not None: mapping["Left Termination" ] = self.leftTermination
            if self.rightTermination is not None: mapping["Right Termination"] = self.rightTermination
        # . Non-default links.
        links = []
        for link in self.nonDefaultLinks: links.append ( link.ToModelFileMapping ( ) )
        if len ( links ) > 0: mapping["Polymer Links"] = links
        # . Finish up.
        return mapping

    def ToSequenceObject ( self, mapping ):
        """Make a sequence object."""
        return SequenceLinearPolymer ( isCyclic = self.isCyclic, leftTerminalComponent  = mapping[self.leftTerminalComponent ] ,
                                                                 rightTerminalComponent = mapping[self.rightTerminalComponent] )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBModelLink ( object ):
    """A class to represent a LINK or SSBOND from a PDB file."""

    defaultAttributes = { "label"          : _DefaultLinkLabel ,
                          "leftComponent"  : None ,
                          "rightComponent" : None ,
                          "_componentKey"  : None ,
                          "_key"           : None }

    def __init__ ( self, leftComponent, rightComponent, **keywordArguments ):
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )
        self.leftComponent  = leftComponent
        self.rightComponent = rightComponent

    @classmethod
    def FromModelFileMapping ( selfClass, mapping, parent, root, usePath = True ):
        """Constructor from a model file mapping."""
        label = mapping.get ( "Label", root.defaultLabel )
        if usePath:
            leftComponent  = root.GetDescendantFromPath ( mapping["Left Component" ] )
            rightComponent = root.GetDescendantFromPath ( mapping["Right Component"] )
        else:
            leftComponent  = parent.childIndex.get ( mapping["Left Component" ], None )
            rightComponent = parent.childIndex.get ( mapping["Right Component"], None )
        return selfClass ( leftComponent, rightComponent, label = label )

    def GetAndApply ( self, library, usedLinks, missingItems ):
        """Get the link from the PDB component library and then apply it."""
        libraryLink = library.GetLink ( self.label, self.leftComponent.genericLabel, self.rightComponent.genericLabel, missingItems = missingItems )
        if libraryLink is not None:
            usedLinks.add ( libraryLink.key )
            self.leftComponent.ApplyLibraryVariant  ( libraryLink.leftVariant  )
            self.rightComponent.ApplyLibraryVariant ( libraryLink.rightVariant )
            return PDBModelBond ( atom1 = self.leftComponent.childIndex.get  ( libraryLink.leftAtomLabel , None ),
                                  atom2 = self.rightComponent.childIndex.get ( libraryLink.rightAtomLabel, None ), bondType = libraryLink.bondType )
        else:
            return None

    @staticmethod
    def MakeKey ( *labels ):
        """Make a link key."""
        return tuple ( labels )

    def ToModelFileMapping ( self, usePath = True ):
        """Get a mapping suitable for a model file."""
        # . Links here are always global.
        if usePath: mapping = { "Left Component" : self.leftComponent.path , "Right Component" : self.rightComponent.path  }
        else:       mapping = { "Left Component" : self.leftComponent.label, "Right Component" : self.rightComponent.label }
        if self.label is not None: mapping["Label"] = self.label
        return mapping

    def ToSequenceObject ( self, mapping ):
        """Make a sequence object."""
        return SequenceLink ( label = self.label, leftComponent = mapping[self.leftComponent], rightComponent = mapping[self.rightComponent] )

    @property
    def componentKey ( self ):
        """Link key with components only."""
        if self._componentKey is None: self._componentKey = PDBModelLink.MakeKey ( self.leftComponent.path, self.rightComponent.path )
        return self._componentKey

    @property
    def isExtraEntity ( self ):
        """Is the link between entities?"""
        return ( self.leftComponent.parent is not self.rightComponent.parent )

    @property
    def isIntraComponent ( self ):
        """Is the link within a component."""
        return ( self.leftComponent is self.rightComponent )

    @property
    def key ( self ):
        """Link key."""
        if self._key is None: self._key = PDBModelLink.MakeKey ( self.label, self.leftComponent.path, self.rightComponent.path )
        return self._key

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBModelVariant ( object ):
    """A class to represent a PDB model variant."""

    defaultAttributes = { "label"     : _DefaultLinkLabel ,
                          "component" :              None ,
                          "_key"      :              None }

    def __init__ ( self, component, **keywordArguments ):
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )
        self.component = component

    @classmethod
    def FromModelFileMapping ( selfClass, mapping, parent ):
        """Constructor from a model file mapping."""
        label          = mapping["Label"]
        componentLabel = mapping["Component"]
        component      = parent.childIndex[componentLabel]
        return selfClass ( component = component, label = label )

    def GetAndApply ( self, library, usedVariants, missingItems ):
        """Get the variant from the PDB component library and then apply it."""
        libraryVariant = library.GetVariant ( self.label, self.component.genericLabel, missingItems = missingItems )
        if libraryVariant is not None:
            usedVariants.add ( libraryVariant.key )
            self.component.ApplyLibraryVariant ( libraryVariant )

    def ToModelFileMapping ( self ):
        """Get a mapping suitable for a model file."""
        # . Variants here are always local to an entity.
        return { "Label" : self.label, "Component" : self.component.label }

    def ToSequenceObject ( self, mapping ):
        """Make a sequence object."""
        return SequenceVariant ( label = self.label, component = mapping[self.component] )

    @staticmethod
    def MakeKey ( linkLabel, componentLabel ):
        """Make a link key."""
        return ( linkLabel, componentLabel )

    @property
    def key ( self ):
        """Link key."""
        if self._key is None: self._key = PDBModelVariant.MakeKey ( self.label, self.component.path )
        return self._key

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBModel ( TreeRootNode ):
    """A PDB model."""

    defaultAttributes = { "altLocs"             : None ,
                          "assemblies"          : None ,
                          "bondIndex"           : None ,   
                          "bonds"               : None ,
                          "linearPolymers"      : None ,
                          "links"               : None ,
                          "nonDefaultLinkIndex" : None,
                          "nonDefaultLinks"     : None ,
                          "numberOfUnknowns"    : 0    ,
                          "variants"            : None }
    defaultAttributes.update ( TreeRootNode.defaultAttributes )

    def AddBond ( self, bond ):
        """Add a bond to the model."""
        if bond is not None:
            if isinstance ( bond, PDBModelBond ): self.bonds.append ( bond )
            else:                                 raise TypeError  ( "Invalid bond being appended to PDB model: {!r}.".format ( type ( bond ) ) )
            if bond.key in self.bondIndex:        raise ValueError ( "Duplicate bond in model: " + bond.atom1.label + "-" + bond.atom2.label + "." )
            else:                                 self.bondIndex[bond.key] = bond

    def AddChild ( self, child, toFollow = None ):
        """Add a child (with sorting)."""
        isOK = super ( PDBModel, self ).AddChild ( child, toFollow = toFollow )
        if isOK: self.children.sort ( key = operator.attrgetter ( "label" ) )
        return isOK

    # . Check for overlapping polymers?
    # . For moment restrict to one per entity.
    def AddLinearPolymer ( self, linearPolymer ):
        """Add a linear polymer."""
        isOK = False
        if isinstance ( linearPolymer, PDBModelLinearPolymer ):
            self.linearPolymers.append ( linearPolymer )
            isOK = True
        return isOK

    def AddLink ( self, link ):
        """Add a link to a component, an entity or to the model itself."""
        if isinstance ( link, PDBModelLink ):
            # . Put the link in the correct place.
            # . Model.
            if link.isExtraEntity:
                if link.key not in self.nonDefaultLinkIndex:
                    # . Append the link.
                    self.nonDefaultLinkIndex[link.key] = link
                    self.nonDefaultLinks.append ( link )
            # . Entity.
            else:
                link.leftComponent.parent.AddLink ( link )
        else: raise TypeError ( "Invalid PDB link type being added to model: {!r}.".format ( type ( link ) ) )

    # . It would be nice to make this general enough to put in Tree.
    def AtomAliasesAssign ( self, aliases ):
        """Assign atom aliases."""
        # . Default labeling.
        for entity in self.children:
            for component in entity.children:
                for atom in component.children:
                    atom.aliases = [ ( entity.label, component.label, atom.label ) ]
        # . Check for aliases.
        if ( aliases is not None ) and ( len ( aliases ) > 0 ):
            # . Loop over aliases in reverse order of priority.
            for ( newPattern, oldPattern ) in reversed ( aliases ):
                # . Parsing.
                newTags = self.ParsePathToFields ( newPattern, fields = [ 1, 3, 1 ] )
                oldTags = self.ParsePathToFields ( oldPattern, fields = [ 1, 3, 1 ] )
                # . Check wild-card patterns.
                for ( newTag, oldTag ) in zip ( newTags, oldTags ):
                    if ( ( ( newTag == self.wildCard ) and ( oldTag != self.wildCard ) ) or \
                         ( ( newTag != self.wildCard ) and ( oldTag == self.wildCard ) ) ):
                        raise TypeError ( "Aliases do not have the same wildcard pattern: " + newPattern + " and " + oldPattern + "." )
                # . Unpack newTags.
                ( entityLabel, resName, resSeq, iCode, atomLabel ) = newTags
                # . Identify atoms conforming to newpattern.
                # . Entity match.
                for entity in self.children:
                    label = entity.label
                    if ( entityLabel == self.wildCard ) or ( entityLabel == label ):
                        # . Entity part of the alias.
                        if entityLabel == self.wildCard: aliasE = [ label      ]
                        else:                            aliasE = [ oldTags[0] ]
                        # . Component match.
                        for component in entity.children:
                            ( componentResName, componentResSeq, componentICode ) = self.ParseLabel ( component.label, fields = 3 )
                            if ( ( ( resName == self.wildCard ) or ( resName == componentResName ) ) and \
                                 ( ( resSeq  == self.wildCard ) or ( resSeq  == componentResSeq  ) ) and \
                                 ( ( iCode   == self.wildCard ) or ( iCode   == componentICode   ) ) ):
                                # . Component part of the alias.
                                if resName == self.wildCard: cName = componentResName
                                else:                        cName = oldTags[1]
                                if resSeq  == self.wildCard: cSeq  = componentResSeq
                                else:                        cSeq  = oldTags[2]
                                if iCode   == self.wildCard: cCode = componentICode
                                else:                        cCode = oldTags[3]
                                aliasEC = aliasE + [ self.MakeLabel ( cName, cSeq, cCode ) ]
                                # . Atom match.
                                for atom in component.children:
                                    if ( atomLabel == self.wildCard ) or ( atomLabel == atom.label ):
                                        # . Set the alias.
                                        if atomLabel == self.wildCard: aliasECA = aliasEC + [ atom.label ]
                                        else:                          aliasECA = aliasEC + [ oldTags[4] ] 
                                        atom.aliases.insert ( 0, tuple ( aliasECA ) )

    def AtomAliasesRemove ( self ):
        """Remove atom aliases."""
        for entity in self.children:
            for component in entity.children:
                for atom in component.children: delattr ( atom, "aliases" )

    def ClearAtoms ( self ):
        """Clear the atoms (and bonds) from a model."""
        # . Component data.
        for entity in self.children:
            for component in entity.children:
                component.ClearAtoms ( )
        # . Model data.
        self.bondIndex = {}
        self.bonds     = []
        self.links     = []
        self.variants  = []

    def ExtractAtoms ( self, other, log = logFile ):
        """Extract atoms from another model."""
        # . Aliases are not put here for the moment as extraction is done by component.
        # . To use aliases:
        # . * Create an intermediate model using ExtractAtoms.
        # . * Create the final model by using aliases on the intermediate model.
        if isinstance ( other, PDBModel ):
            if self.NumberOfAtoms ( ) > 0: raise ValueError ( "The model extracting atoms already has atoms." )
            # . Loop over components.
            numberExtracted = 0
            for entity in self.children:
                for component in entity.children:
                    otherComponent = other.GetDescendantFromLabels ( entity.label, component.label )
                    if otherComponent is not None:
                        for otherAtom in otherComponent.children:
                            component.AddChild ( otherAtom )
                        numberExtracted += len ( otherComponent.children )
            # . Finish up.
            if numberExtracted > 0: self.Verify ( )
            # . Do some printing.
            if LogFileActive ( log ): log.Paragraph ( "Number of atoms extracted from model = {:d}.".format ( numberExtracted ) )
        else: raise TypeError ( "Invalid PDB model: {!r}.".format ( type ( other ) ) )

    def ExtractAtomData ( self, other, aliases = [], log = logFile ):
        """Extract atom data from the atoms of another model."""
        if isinstance ( other, PDBModel ):
            # . Assign atom aliases.
            self.AtomAliasesAssign ( aliases )
            # . Loop over atoms.
            numberExtracted = 0
            for entity in self.children:
                for component in entity.children:
                    for atom in component.children:
                        for ( entityLabel, componentLabel, atomLabel ) in atom.aliases:
                            otherAtom = other.GetDescendantFromLabels ( entityLabel, componentLabel, atomLabel )
                            if otherAtom is not None:
                                for datum in otherAtom.children:
                                    newDatum = PDBModelAtomData ( )
                                    datum.CopyTo ( newDatum )
                                    atom.AddChild ( newDatum )
                                numberExtracted += 1
                                break
            # . Finish up.
            if numberExtracted > 0: self.IdentifyAltLocs ( )
            # . Remove atom aliases.
            self.AtomAliasesRemove ( )
            # . Do some printing.
            if LogFileActive ( log ): log.Paragraph ( "Number of atom data extracted from model = {:d}.".format ( numberExtracted ) )
        else: raise TypeError ( "Invalid PDB model: {!r}.".format ( type ( other ) ) )

    @classmethod
    def FromModelFileMapping ( selfClass, mapping ):
        """Constructor from a model file mapping."""
        label = mapping.get ( "Label", None )
        self  = selfClass ( label = label )
        # . Entities.
        for subMapping in mapping.get ( "Entities", [] ):
            self.AddChild ( PDBModelEntity.FromModelFileMapping ( subMapping, self ) )
        # . Linear polymers.
        for subMapping in mapping.get ( "Linear Polymers", [] ):
            self.AddLinearPolymer ( PDBModelLinearPolymer.FromModelFileMapping ( subMapping, self ) )
        # . Links.
        for subMapping in mapping.get ( "Links", [] ):
            self.AddLink ( PDBModelLink.FromModelFileMapping ( subMapping, self, self, usePath = True ) )
        # . Finish up.
        return self

    def GatherLinks ( self, library, missingItems ):
        """Gather the links for the model."""
        self.links = list ( self.nonDefaultLinks )
        for entity in self.children:
            self.links.extend ( entity.GatherLinks  ( library, missingItems ) )
        for polymer in self.linearPolymers:
            self.links.extend ( polymer.GatherLinks ( library, missingItems ) )

    def GatherVariants ( self, library, missingItems ):
        """Gather the variants for the model."""
        self.variants = []
        for entity in self.children:
            self.variants.extend ( entity.GatherVariants  ( library, missingItems ) )
        for polymer in self.linearPolymers:
            self.variants.extend ( polymer.GatherVariants ( library, missingItems ) )

    def IdentifyAltLocs ( self ):
        """Get a mapping of non-default AltLoc labels to atom paths."""
        self.altLocs = {}
        for entity in self.children:
            for component in entity.children:
                for atom in component.children:
                    for altLoc in atom.childIndex:
                        if altLoc != _DefaultAltLoc:
                            items = self.altLocs.get ( altLoc, [] )
                            items.append ( atom.path )
                            self.altLocs[altLoc] = items

    def IdentifyUnknownAtoms ( self ):
        """Identify unknown atoms in the model."""
        n = 0
        for entity in self.children:
            for component in entity.children:
                for atom in component.children:
                    if atom.atomicNumber < 0: n += 1
        self.numberOfUnknowns = n
        return n

    def Initialize ( self ):
        """Initialization."""
        super ( PDBModel, self ).Initialize ( )
        self.altLocs             = {}
        self.bondIndex           = {}
        self.bonds               = []
        self.linearPolymers      = []
        self.links               = []
        self.nonDefaultLinkIndex = {}
        self.nonDefaultLinks     = []
        self.variants            = []

    def MakeAtomicModelFromComponentLibrary ( self, libraryPaths = None, log = logFile ):
        """Make an atomic model from the model and the PDB component library."""
        # . Clear existing atoms.
        self.ClearAtoms ( )
        # . Get the library.
        library = PDBComponentLibrary ( paths = libraryPaths )
        # . Initialization.
        usedComponents = set ( )
        usedLinks      = set ( )
        usedVariants   = set ( )
        missingItems   = set ( )
        # . Build atoms and bonds for the components.
        for entity in self.children:
            entity.MakeAtomicModelFromComponentLibrary ( library, usedComponents, missingItems )
        # . Gather and apply variants.
        self.GatherVariants ( library, missingItems )
        for variant in self.variants:
            variant.GetAndApply ( library, usedVariants, missingItems )
# . Gather atoms and bonds?
# . Problem is coordinating atom and bond lists with AddChild, etc.
        # . Gather and apply links.
        self.GatherLinks ( library, missingItems )
        for link in self.links:
            self.AddBond ( link.GetAndApply ( library, usedLinks, missingItems ) )
        # . Print items that were employed.
        if LogFileActive ( log ):
            table = log.GetTable ( columns = [ 12, 24 ] )
            table.Start ( )
            table.Title ( "Employed PDB Components" )
            table.Heading ( "Type"  )
            table.Heading ( "Label" )
            for ( items, label ) in ( ( usedComponents, "Component" ), ( usedLinks, "Link" ), ( usedVariants, "Variant" ) ):
                keys = list ( items )
                keys.sort ( )
                for key in keys:
                    table.Entry ( label, alignment = "l" )
                    table.Entry ( key,   alignment = "l" )
            table.Stop ( )
        # . Check for missing items.
        if len ( missingItems ) > 0:
            if LogFileActive ( log ):
                missingItems = list ( missingItems )
                missingItems.sort ( )
                table = log.GetTable ( columns = [ 12, 24 ] )
                table.Start ( )
                table.Title ( "Missing PDB Components" )
                table.Heading ( "Type"  )
                table.Heading ( "Label" )
                for ( label, key ) in missingItems:
                    table.Entry ( label, alignment = "l" )
                    table.Entry ( key,   alignment = "l" )
                table.Stop ( )
            raise ValueError ( "There are {:d} missing items.".format ( len ( missingItems ) ) )

    def MakeAtomicNumbers ( self ):
        """Return a list of the atomic numbers of the atoms in the model."""
        atomicNumbers = []
        for entity in self.children:
            for component in entity.children:
                for atom in component.children:
                    atomicNumbers.append ( atom.atomicNumber )
        return atomicNumbers

    def MakeBonds ( self, atomMapping ):
        """Make a list of bonds."""
        bonds = []
        for modelBond in self.bonds:
            bonds.append ( Bond ( i = atomMapping[modelBond.atom1], j = atomMapping[modelBond.atom2], type = modelBond.bondType ) )
        for entity in self.children:
            for component in entity.children:
                for modelBond in component.bonds:
                    bonds.append ( Bond ( i = atomMapping[modelBond.atom1], j = atomMapping[modelBond.atom2], type = modelBond.bondType ) )
        return bonds

    def MakeCoordinates3 ( self, altLoc ):
        """Make a coordinates3 object from the model.

        If altLoc is not specified, the coordinates of the most occupied location for the atom are taken.
        """
        coordinates3 = Coordinates3.WithExtent ( self.NumberOfAtoms ( ) )
        coordinates3.Set ( _PDBAtomDataUndefined )
        n = 0
        for entity in self.children:
            for component in entity.children:
                for atom in component.children:
                    data = atom.GetData ( altLoc )
                    if data is None:
                        coordinates3.FlagCoordinateAsUndefined ( n )
                    else:
                        coordinates3[n,0] = data.x
                        coordinates3[n,1] = data.y
                        coordinates3[n,2] = data.z
                    n += 1
        return coordinates3

    def MakeData ( self, attribute, altLoc ):
        """Return an array with the appropriate atom data."""
        values = Real1DArray.WithExtent ( self.NumberOfAtoms ( ) )
        n      = 0
        for entity in self.children:
            for component in entity.children:
                for atom in component.children:
                    data      = atom.GetData ( altLoc )
                    values[n] = getattr ( data, attribute, _PDBAtomDataUndefined )
                    n += 1
        return values

    def MakeSequence ( self ):
        """Make a sequence from the model."""
        # . Basic sequence.
        atomMapping = {}
        index       = 0
        sequence    = Sequence ( )
        for modelEntity in self.children:
            entity = SequenceEntity ( label = modelEntity.label )
            sequence.AddChild ( entity )
            for modelComponent in modelEntity.children:
                component = SequenceComponent ( genericLabel = modelComponent.genericLabel, label = modelComponent.label )
                entity.AddChild ( component )
                for modelAtom in modelComponent.children:
                    atom = Atom ( atomicNumber = modelAtom.atomicNumber, formalCharge = modelAtom.formalCharge, index = index, label = modelAtom.label )
                    component.AddChild ( atom )
                    atomMapping[modelAtom] = index
                    index += 1
        # . Check for links and variants - done here until paths, etc. sorted out.
        # . Make link and variant lists.
        doLinearPolymers = len ( self.linearPolymers ) > 0
        doLinks          = len ( self.links          ) > 0
        doVariants       = len ( self.variants       ) > 0
        # . Make component index.
        if doLinearPolymers or doLinks or doVariants:
            componentMapping = self.MakeSequenceComponentMapping ( sequence )
        # . Linear polymers.
        if doLinearPolymers > 0:
            items = []
            for item in self.linearPolymers:
                items.append ( item.ToSequenceObject ( componentMapping ) )
            if len ( items ) > 0:
                sequence.__dict__["linearPolymers"] = items
        # . Links.
        if doLinks:
            items = []
            for item in self.links:
                items.append ( item.ToSequenceObject ( componentMapping ) )
            if len ( items ) > 0:
                sequence.__dict__["links"] = items
        # . Variants.
        if doVariants:
            items = []
            for item in self.variants:
                items.append ( item.ToSequenceObject ( componentMapping ) )
            if len ( items ) > 0:
                sequence.__dict__["variants"] = items
        return ( sequence, atomMapping )

    def MakeSequenceComponentMapping ( self, sequence ):
        """Make a mapping between model and sequence components."""
        # . Model components.
        modelComponentMapping = {}
        for entity in self.children:
            for modelComponent in entity.children:
                modelComponentMapping[modelComponent.path] = modelComponent
        # . Sequence components.
        sequenceComponentMapping = {}
        for entity in sequence.children:
            for component in entity.children:
                sequenceComponentMapping[modelComponentMapping[component.path]] = component 
        return sequenceComponentMapping

    def MakeSystem ( self, altLoc = None, label = None ):
        """Make a system from the model."""
        ( sequence, atomMapping ) = self.MakeSequence ( )
        bonds                     = self.MakeBonds        ( atomMapping )
        coordinates3              = self.MakeCoordinates3 ( altLoc      )
        if sequence is not None:
            system                = System.FromSequence ( sequence, bonds = bonds )
            system.coordinates3   = coordinates3
            if        label is not None: system.label = label
            elif self.label is not None: system.label = self.label
        return system

    def NumberOfAtoms ( self ):
        """Return the number of atoms."""
        n = 0
        for entity in self.children:
            for component in entity.children: n += len ( component.children )
        return n

    def NumberOfBonds ( self ):
        """Return the number of bonds."""
        n = len ( self.bonds )
        for entity in self.children:
            for component in entity.children: n += len ( component.bonds )
        return n

    def NumberOfComponents ( self ):
        """Return the number of components."""
        n = 0
        for entity in self.children: n += len ( entity.children )
        return n

    def OrderHydrogens ( self, embedded = False ):
        """Order hydrogens within each component.

        Ordering is done by inserting hydrogens after the heavy atoms to which they
        are bound (embedded = True) or at end of each component (embedded = False)."""
        for entity in self.children:
            for component in entity.children: component.OrderHydrogens ( embedded = embedded )

    def PrintComponentFrequency ( self, log = logFile, title = "PDB Model Component Frequency" ):
        """Print the component frequencies."""
        if LogFileActive ( log ):
            log.Heading ( title )
            for entity in self.children: entity.PrintComponentFrequency ( log = log )

    def PrintComponentSequence ( self, log = logFile, title = "PDB Model Component Sequence" ):
        """Print the component sequences."""
        if LogFileActive ( log ):
            log.Heading ( title )
            for entity in self.children: entity.PrintComponentSequence ( log = log )

    def Summary ( self, log = logFile, title = None ):
        """Summary."""
        if LogFileActive ( log ):
            numberOfBonds = self.NumberOfBonds ( )
            summary       = log.GetSummary ( )
            if title is None: summary.Start ( "PDB Model Summary" )
            else:             summary.Start ( title )
            summary.Entry ( "Atoms"          , "{:d}".format ( self.NumberOfAtoms      ( ) ) )
            summary.Entry ( "Components"     , "{:d}".format ( self.NumberOfComponents ( ) ) )
            summary.Entry ( "Entities  "     , "{:d}".format ( len ( self.children       ) ) )
            summary.Entry ( "Linear Polymers", "{:d}".format ( len ( self.linearPolymers ) ) )
            if len ( self.altLocs    ) > 0: summary.Entry ( "Non-blank AltLoc IDs", "{:d}".format ( len ( self.altLocs     ) ) )
            if ( self.assemblies is not None ) and ( len ( self.assemblies ) > 0 ):
                                            summary.Entry ( "Assemblies"          , "{:d}".format ( len ( self.assemblies  ) ) )
            if numberOfBonds           > 0: summary.Entry ( "Bonds"               , "{:d}".format ( numberOfBonds            ) )
            if len ( self.links      ) > 0: summary.Entry ( "Links"               , "{:d}".format ( len ( self.links       ) ) )
            if self.numberOfUnknowns   > 0: summary.Entry ( "Unknown Elements"    , "{:d}".format ( self.numberOfUnknowns    ) )
            if len ( self.variants   ) > 0: summary.Entry ( "Variants"            , "{:d}".format ( len ( self.variants    ) ) )
            summary.Stop ( )
 
    def ToModelFileMapping ( self ):
        """Get a mapping suitable for a model file."""
        mapping = {}
        # . Label.
        if self.label is not None: mapping["Label"] = self.label
        # . Entities.
        items = []
        for item in self.children: items.append ( item.ToModelFileMapping ( ) )
        if len ( items ) > 0: mapping["Entities"] = items
        # . Linear polymers.
        items = []
        for item in self.linearPolymers: items.append ( item.ToModelFileMapping ( ) )
        if len ( items ) > 0: mapping["Linear Polymers"] = items
        # . Links.
        items = []
        for item in self.nonDefaultLinks: items.append ( item.ToModelFileMapping ( usePath = True ) )
        if len ( items ) > 0: mapping["Links"] = items
        # . Finish up.
        return mapping

    def Verify ( self ):
        """Verify the model."""
        self.IdentifyAltLocs      ( )
        self.IdentifyUnknownAtoms ( )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def PDBModel_FromModelFile ( path, log = logFile ):
    """Construct a model from a model file."""
    mapping = YAMLUnpickle ( path )
    return PDBModel.FromModelFileMapping ( mapping )

def PDBModel_ToModelFile ( path, model, log = logFile ):
    """Write a model to a model file."""
    mapping = model.ToModelFileMapping ( )
    YAMLPickle ( path, mapping )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
