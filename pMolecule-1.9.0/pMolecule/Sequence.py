#-------------------------------------------------------------------------------
# . File      : Sequence.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A data structure for handling macromolecular sequence information."""

"""
A Sequence is structured as a tree in a three-level hierarchy that
consists of entities, components and atoms. Each item in the tree has
a path which is unique. The path is constructed from the item's label
and the labels of its antecedents (or parents) in the tree. Thus, the
path of an atom is:

               entityLabel:componentLabel:atomLabel

Labels in a path are separated by ":" by default, although this is
customizable.

The items in a sequence are ordered in the sense that the children of
a given node in the tree are stored in an ordered fashion.

To accommodate the PDB naming system, the labels of items can consist of
a concatenation of fields separated (optionally) by ".". For example, a
component label in a PDB sequence could be 'ASP.1.A".
"""

import copy

from pCore   import logFile, LogFileActive, Selection, TreeBranchNode, TreeLeafNode, TreeRootNode

from Atom    import Atom
from Element import PeriodicTable

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Sequence ( TreeRootNode ):
    """A sequence."""

    defaultAttributes = { "linearPolymers" : None ,
                          "links"          : None ,
                          "variants"       : None }
    defaultAttributes.update ( TreeRootNode.defaultAttributes )

    def AtomIndex ( self, atomPath ):
        """Get the index of the atom with a given path."""
        # . Initialization.
        index = -1
        # . Parse the path.
        ( entityLabel, componentLabel, atomLabel ) = self.ParsePath ( atomPath, labels = 3 )
        # . Entity/component/atom match.
        for entity in self.children:
            if entityLabel == entity.label:
                for component in entity.children:
                    if componentLabel == component.label:
                        for atom in component.children:
                            if atomLabel == atom.label:
                                index = atom.index
                                break
        # . Finish up.
        if index < 0: raise KeyError ( "Atom path - {:s} - not found.".format ( atomPath ) )
        return index

    def FieldsLabelMatch ( self, pFields, label ):
        """Match a set of fields versus a label and check for a match."""
        matched = True
        tokens  = self.ParseLabel ( label )
        if len ( tokens ) <= len ( pFields ):
            for ( field, token ) in zip ( pFields, tokens ):
                if ( field != self.wildCard ) and ( field != token.strip ( ) ):
                    matched = False
                    break
            if matched:
                for field in pFields[len ( tokens ):]:
                    if ( field != "" ) and ( field != self.wildCard ):
                        matched = False
                        break
        elif pFields[-1] != self.wildCard:
            matched = False
        return matched

    @classmethod
    def FromAtomContainer ( selfClass, atoms, componentLabel = "UNK.1", entityLabel = "A" ):
        """Create a default sequence given an atom container."""
        # . Sequence tree.
        if componentLabel == None: componentLabel = selfClass.defaultAttributes["defaultLabel"]
        if    entityLabel == None:    entityLabel = selfClass.defaultAttributes["defaultLabel"]
        self      = selfClass ( )
        entity    = SequenceEntity    ( label = entityLabel    )
        component = SequenceComponent ( label = componentLabel )
        self.AddChild   ( entity    )
        entity.AddChild ( component )
        # . Are the atom labels unique?
        labels = set ( )
        for atom in atoms:
            if atom.label is not None: labels.add ( atom.label )
        reassignLabels = ( len ( atoms ) != len ( labels ) )
        # . Assign atoms.
        for ( i, atom ) in enumerate ( atoms ):
            atom.index = i
            if reassignLabels:
                atom.label = PeriodicTable.Symbol ( atom.atomicNumber, index = i )
            component.AddChild ( atom )
        # . Set path to label if the component and entity labels are not set.
        if ( len ( componentLabel ) <= 0 ) and ( len ( entityLabel ) <= 0 ):
            for atom in atoms: atom._path = atom.label
        # . There is no associated data.
        return self

    def GatherAtoms ( self ):
        """Gather all atoms."""
        atoms = []
        for entity in self.children:
            for component in entity.children:
                for atom in component.children: atoms.append ( atom )
        return atoms

    @classmethod
    def FromAtomPaths ( selfClass, atomPaths, atoms = None ):
        """Create a sequence from a list of atom paths."""
        # . Argument check.
        if ( atoms is not None ) and ( len ( atoms ) != len ( atomPaths ) ):
            raise ValueError ( "The \"atomPaths\" and \"atoms\" arguments are incompatible." )
        createAtoms = ( atoms is None )
        # . Initialization.
        self = selfClass ( )
        # . Process the paths.
        duplicateAtoms = 0
        index          = 0
        for atomPath in atomPaths:
            ( entityLabel, componentLabel, atomLabel ) = self.ParsePath ( atomPath, labels = 3 )
            ( entity, component, atom ) = self.GetDescendantsFromLabels ( entityLabel, componentLabel, atomLabel )
            if entity is None:
                entity = SequenceEntity ( label = entityLabel )
                self.AddChild ( entity )
            if component is None:
                genericLabel = self.ParseLabel ( componentLabel, fields = 1 )[0]
                component    = SequenceComponent ( genericLabel = genericLabel, label = componentLabel )
                entity.AddChild ( component )
            if atom is None:
                if createAtoms:
                    atom = Atom ( index = index, label = atomLabel )
                else:
                    atom       = atoms[index]
                    atom.index = index
                    atom.label = atomLabel
                index += 1
                component.AddChild ( atom )
            else:
                duplicateAtoms += 1
        if duplicateAtoms > 0:
            raise ValueError ( "There were {:d} duplicate atom paths in the sequence.".format ( duplicateAtoms ) )
        # . Finish up.
        return self

    @classmethod
    def FromMapping ( selfClass, mapping ):
        """Constructor from a mapping."""
        self = selfClass ( )
        # . Label.
        label = mapping.get ( "label", None )
        if label is not None: self.label = label
        # . Entities.
        componentMapping = {}
        for subMapping in mapping.get ( "entities", [] ):
            self.AddChild ( SequenceEntity.FromMapping ( subMapping, self, componentMapping ) )
        # . Linear polymers.
        for ( leftPath, rightPath, isCyclic ) in mapping.get ( "linearPolymers", [] ):
            self.linearPolymers.append ( SequenceLinearPolymer ( leftTerminalComponent  = componentMapping[leftPath ] ,
                                                                 rightTerminalComponent = componentMapping[rightPath] ,
                                                                 isCyclic = isCyclic ) )
        # . Links.
        for ( leftPath, rightPath, label ) in mapping.get ( "links", [] ):
            self.links.append ( SequenceLink ( leftComponent  = componentMapping[leftPath ] ,
                                               rightComponent = componentMapping[rightPath] ,
                                               label = label ) )
        # . Variants.
        for ( path, labels ) in mapping.get ( "variants", {} ).iteritems ( ):
            component = componentMapping[path]
            for label in labels:
                self.variants.append ( SequenceVariant ( component = component, label = label ) )
        # . Finish up.
        return self

    def Initialize ( self ):
        """Initialization."""
        super ( Sequence, self ).Initialize ( )
        self.linearPolymers = []
        self.links          = []
        self.variants       = []

    def Merge ( self, others, information = {} ):
        """Merging."""
        # . Initialization.
        entityLabels = {}
        merged       = self.__class__ ( )
        # . Check to see whether the merged atoms already exist.
        atomMapping      = information.get ( "atomMapping", None )
        mergedAtomsExist = ( atomMapping is not None )
        # . Loop over all sequences.
        index = 0
        for sequence in [ self ] + others:
            # . Tree items.
            mapping = {}
            for oldEntity in sequence.children:
                # . Get the entity label.
                oldLabel = oldEntity.label
                if oldLabel in entityLabels:
                    n         = entityLabels[oldLabel]
                    newLabel  = oldLabel + merged.fieldSeparator + repr ( n )
                    n        += 1
                else:
                    newLabel = oldLabel
                    n        = 0
                entityLabels[oldLabel] = n
                # . Entity.
                newEntity = SequenceEntity ( label = newLabel )
                merged.AddChild ( newEntity )
                # . Components.
                for oldComponent in oldEntity.children:
                    newComponent = SequenceComponent ( genericLabel = oldComponent.genericLabel, label = oldComponent.label )
                    newEntity.AddChild ( newComponent )
                    mapping[oldComponent] = newComponent
                    # . Atoms.
                    for oldAtom in oldComponent.children:
                        if mergedAtomsExist: newAtom = atomMapping[oldAtom]
                        else:                newAtom = copy.copy ( oldAtom )
                        newAtom.index = index
                        newComponent.AddChild ( newAtom )
                        index += 1
            # . Non-tree items.
            # . Linear polymers.
            for item in sequence.linearPolymers:
                merged.linearPolymers.append ( SequenceLinearPolymer ( isCyclic = item.isCyclic, leftTerminalComponent  = mapping[item.leftTerminalComponent ] , \
                                                                                                 rightTerminalComponent = mapping[item.rightTerminalComponent] ) )
            # . Links.
            for item in sequence.links:
                merged.links.append ( SequenceLink ( label = item.label, leftComponent = mapping[item.leftComponent], rightComponent = mapping[item.rightComponent] ) )
            # . Variants.
            for item in sequence.variants:
                merged.variants.append ( SequenceVariant ( component = mapping[item.component], label = item.label ) )
        # . Finish up.
        return merged

    def NumberOfAtoms ( self ):
        """Return the number of atoms."""
        n = 0
        for entity in self.children:
            for component in entity.children: n += len ( component.children )
        return n

    def NumberOfComponents ( self ):
        """Return the number of components."""
        n = 0
        for entity in self.children: n += len ( entity.children )
        return n

    # . It would be nice to make this more general too (in Tree).
    def ParseAtomPattern ( self, atomPattern ):
        """Parse an atom pattern."""
        # . Split the pattern into labels.
        labels = self.ParsePath ( atomPattern, labels = 3 )
        # . Split the labels into fields.
        pFields = [ [] for n in range ( 3 ) ]
        for ( i, label ) in enumerate ( labels ):
            # . Stripped tokens.
            tokens = self.ParseLabel ( label )
            fields = []
            for token in tokens: fields.append ( token.strip ( ) )
            # . Remove all trailing empty fields.
            n = len ( fields )
            for field in reversed ( fields ):
                if field == self.defaultLabel: n -= 1
                else:                          break
            fields = fields[0:n]
            # . Condense multiple trailing wildcards to a single one.
            n = 0
            for field in reversed ( fields ):
                if field == self.wildCard: n += 1
                else: break
            if n > 1: fields = fields[0:-(n-1)]
            # . Save.
            pFields[i].extend ( fields )
        # . At a minimum, a field must consist of a default label.
        for pField in pFields:
            if len ( pField ) <= 0:
                pField.append ( self.defaultLabel )
        return ( pFields[0], pFields[1], pFields[2] )

    def PrintComponentFrequency ( self, log = logFile, title = "Sequence Component Frequency" ):
        """Print the component frequencies."""
        if LogFileActive ( log ):
            log.Heading ( title )
            for entity in self.children: entity.PrintComponentFrequency ( log = log )

    def PrintComponentSequence ( self, log = logFile, title = "Sequence Component Sequence" ):
        """Print the component sequences."""
        if LogFileActive ( log ):
            log.Heading ( title )
            for entity in self.children: entity.PrintComponentSequence ( log = log )

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        # . Find atoms to keep.
        atoms  = self.GatherAtoms ( )
        toKeep = []
        for s in selection: toKeep.append ( atoms[s] )
        # . Check to see whether the pruned atoms already exist.
        atomMapping      = information.get ( "atomMapping", None )
        prunedAtomsExist = ( atomMapping is not None )
        # . Construct the pruned tree.
        mapping = {}
        pruned  = Sequence ( )
        for ( index, atom ) in enumerate ( toKeep ):
            component      = atom.parent
            componentLabel = component.label
            entityLabel    = component.parent.label
            newEntity      = pruned.childIndex.get ( entityLabel, None )
            if newEntity is None:
                newEntity = SequenceEntity ( label = entityLabel )
                pruned.AddChild ( newEntity )
            newComponent = newEntity.childIndex.get ( componentLabel, None )
            if newComponent is None:
                newComponent = SequenceComponent ( genericLabel = component.genericLabel, label = componentLabel )
                newEntity.AddChild ( newComponent )
                mapping[component] = newComponent
            if prunedAtomsExist: newAtom = atomMapping[atom]
            else:                newAtom = copy.copy ( atom )
            newAtom.index = index
            newComponent.AddChild ( newAtom )
        # . Non-tree items.
        # . Linear polymers.
        for item in self.linearPolymers:
            if ( item.leftTerminalComponent in mapping ) and ( item.rightTerminalComponent in mapping ):
                pruned.linearPolymers.append ( SequenceLinearPolymer ( isCyclic = item.isCyclic, leftTerminalComponent  = mapping[item.leftTerminalComponent ] , \
                                                                                                 rightTerminalComponent = mapping[item.rightTerminalComponent] ) )
        # . Links.
        for item in self.links:
            if ( item.leftComponent in mapping ) and ( item.rightComponent in mapping ):
                pruned.links.append ( SequenceLink ( label = item.label, leftComponent = mapping[item.leftComponent], rightComponent = mapping[item.rightComponent] ) )
        # . Variants.
        for item in self.variants:
            if item.component in mapping:
                pruned.variants.append ( SequenceVariant ( component = mapping[item.component], label = item.label ) )
        # . Finish up.
        return pruned

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( valueWidth = 10 )
            summary.Start ( "Sequence Summary" )
            summary.Entry ( "Number of Atoms"          , "{:d}".format ( self.NumberOfAtoms      ( ) ) )
            summary.Entry ( "Number of Components"     , "{:d}".format ( self.NumberOfComponents ( ) ) )
            summary.Entry ( "Number of Entities"       , "{:d}".format ( len ( self.children       ) ) )
            summary.Entry ( "Number of Linear Polymers", "{:d}".format ( len ( self.linearPolymers ) ) )
            summary.Entry ( "Number of Links"          , "{:d}".format ( len ( self.links          ) ) )
            summary.Entry ( "Number of Variants"       , "{:d}".format ( len ( self.variants       ) ) )
            summary.Stop ( )

    def ToMapping ( self ):
        """Get a mapping suitable for serialization."""
        mapping = {}
        # . Label.
        if self.label is not None: mapping["label"] = self.label
        # . Entities.
        items = []
        for item in self.children: items.append ( item.ToMapping ( ) )
        if len ( items ) > 0: mapping["entities"] = items
        # . Other data.
        # . Linear polymers.
        items = []
        for item in self.linearPolymers:
            items.append ( [ item.leftTerminalComponent.path, item.rightTerminalComponent.path, item.isCyclic ] )
        if len ( items ) > 0:
            mapping["linearPolymerFields"] = [ "leftTerminalComponent", "rightTerminalComponent", "isCyclic" ] 
            mapping["linearPolymers"     ] = items
        # . Links.
        items = []
        for item in self.links:
            items.append ( [ item.leftComponent.path, item.rightComponent.path, item.label ] )
        if len ( items ) > 0:
            mapping["linkFields"] = [ "leftComponent", "rightComponent", "label" ] 
            mapping["links"     ] = items
        # . Variants.
        items = {}
        for item in self.variants:
            path = item.component.path
            data = items.get ( path, None )
            if data is None:
                data = []
                items[path] = data
            data.append ( item.label )
        if len ( items ) > 0:
            mapping["variantFields"] = [ "component", "labels" ]
            mapping["variants"     ] = items
        # . Finish up.
        return mapping

    @property
    def linearPolymerIndex ( self ):
        """Get an index of the components in linear polymers."""
        if "_linearPolymerIndex" not in self.__dict__:
            index = {}
            for ( n, item ) in enumerate ( self.linearPolymers ):
                for component in item.ComponentIterator ( ):
                    index[component] = n
            self.__dict__["_linearPolymerIndex"] = index
        return self.__dict__["_linearPolymerIndex"]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SequenceComponent ( TreeBranchNode ):
    """A sequence component."""

    defaultAttributes = { "genericLabel" : None }
    defaultAttributes.update ( TreeBranchNode.defaultAttributes )

#===================================================================================================================================
# . Class for entities.
#===================================================================================================================================
class SequenceEntity ( TreeBranchNode ):
    """A sequence entity."""

    defaultAttributes = { }
    defaultAttributes.update ( TreeBranchNode.defaultAttributes )

    @classmethod
    def FromMapping ( selfClass, mapping, root, componentMapping ):
        """Constructor from a model file mapping."""
        label = mapping.get ( "label", root.defaultLabel )
        self  = selfClass ( label = label )
        root.AddChild ( self ) # . Entity needs to be rooted before calculate paths.
        # . Components.
        for ( label, genericLabel, atomLabels ) in mapping.get ( "components", [] ):
            component = SequenceComponent ( genericLabel = genericLabel, label = label )
            self.AddChild ( component )
            for atomLabel in atomLabels:
                atom = Atom ( label = atomLabel )
                component.AddChild ( atom )
            componentMapping[component.path] = component
        # . Finish up.
        return self

    def Label ( self ):
        """Return a label."""
        if self.label == "": return "Unnamed Entity"
        else:                return "Entity " + self.label

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
            else:             table.Title ( title                        )
            for key in keys:
                table.Entry ( key                               , alignment = "r" )
                table.Entry ( "{:d}".format ( frequencies[key] ), alignment = "r" )
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

    def ToMapping ( self ):
        """Return a mapping suitable for serialization."""
        mapping = {}
        if ( self.label is not None ) and ( len ( self.label ) > 0 ): mapping["label"] = self.label 
        components = []
        for component in self.children:
            atoms = []
            for atom in component.children: atoms.append ( atom.label )
            components.append ( [ component.label, component.genericLabel, atoms ] )
        if len ( components ) > 0:
            mapping["componentFields"] = [ "label", "genericLabel", "atoms" ]
            mapping["components"     ] = components
        return mapping

#===================================================================================================================================
# . Class for linear polymers.
#===================================================================================================================================
class SequenceLinearPolymer ( object ):
    """A sequence linear polymer."""

    defaultAttributes = { "isCyclic"               : False ,
                          "leftTerminalComponent"  : None  ,
                          "rightTerminalComponent" : None  }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

    def ComponentIterator ( self ): return SequenceLinearPolymerComponentIterator ( self )

#===================================================================================================================================
# . Class for iteration over the components of a linear polymer.
#===================================================================================================================================
class SequenceLinearPolymerComponentIterator ( object ):
    """Iterator over the components of a linear polymer."""

    def __init__ ( self, polymer ):
        """Constructor."""
        leftTerminus       = polymer.leftTerminalComponent
        self.components    = leftTerminus.parent.children
        self.isFinished    = False
        self.polymer       = polymer
        self.position      = self.components.index ( leftTerminus )
        self.rightTerminus = polymer.rightTerminalComponent

    def __iter__ ( self ): return self

    def __next__ ( self ): return self.next ( )

    def next ( self ):
        """Next item."""
        if self.isFinished: raise StopIteration
        component       = self.components[self.position]
        self.isFinished = ( component is self.rightTerminus )
        self.position  += 1
        return component

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SequenceLink ( object ):
    """A sequence link."""

    defaultAttributes = { "label"          : None ,
                          "leftComponent"  : None , 
                          "rightComponent" : None }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SequenceVariant ( object ):
    """A sequence variant."""

    defaultAttributes = { "component" : None ,
                          "label"     : None }
    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

#===================================================================================================================================
# . Test the module.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
