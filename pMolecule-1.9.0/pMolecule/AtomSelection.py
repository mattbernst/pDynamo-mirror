#-------------------------------------------------------------------------------
# . File      : AtomSelection.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Atom selection."""

from pCore  import CrossPairList_FromSingleCoordinates3, logFile, LogFileActive, Selection
from Atom   import Atom
from System import System

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AtomSelectionError ( Exception ): pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AtomSelection ( object ):
    """An atom selection."""

    defaultAttributes = { "selection" : None ,
                          "system"    : None }

    def __and__ ( self, other ):
        """Intersection."""
        self.CheckSystems ( other )
        return self.__class__ ( selection = ( self.selection & other.selection ), system = self.system )

    def __contains__ ( self, item ):
        """Containing using Atom instance, integer index or string path."""
        if   isinstance ( item, Atom ): return ( ( item.index in self.selection ) and ( item is self.system.atoms[item.index] ) )
        elif isinstance ( item, int  ): return ( item in self.selection )
        elif isinstance ( item, str  ): return ( self.system.sequence.AtomIndex ( item ) in self.selection )

    def __iter__ ( self ): return iter ( self.selection )

    def __iand__ ( self, other ):
        """In-place intersection."""
        self.CheckSystems ( other )
        self.selection &= other.selection
        return self

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )
        self.CheckAttributes ( )

    def __invert__ ( self ):
        """Complement."""
        return self.Complement ( )

    def __ior__ ( self, other ):
        """In-place union."""
        self.CheckSystems ( other )
        self.selection |= other.selection
        return self

    def __isub__ ( self, other ):
        """In-place difference."""
        self.CheckSystems ( other )
        self.selection -= other.selection
        return self

    def __ixor__ ( self, other ):
        """In-place symmetric difference."""
        self.CheckSystems ( other )
        self.selection ^= other.selection
        return self

    def __len__ ( self ): return len ( self.selection )

    def __or__ ( self, other ):
        """Union."""
        self.CheckSystems ( other )
        return self.__class__ ( selection = ( self.selection | other.selection ), system = self.system )

    def __sub__ ( self, other ):
        """Difference."""
        self.CheckSystems ( other )
        return self.__class__ ( selection = ( self.selection - other.selection ), system = self.system )

    def __xor__ ( self, other ):
        """Symmetric difference."""
        self.CheckSystems ( other )
        return self.__class__ ( selection = ( self.selection ^ other.selection ), system = self.system )

    def Add ( self, index ):
        """Add an index to the selection."""
        if ( index < 0 ) or ( index >= len ( self.system.atoms ) ): raise AtomSelectionError ( "Index - {:d} - out of range [0,{:d}].".format ( index, len ( self.system.atoms ) ) )
        self.selection.add ( index )

    def ByBondedNeighbor ( self, iterations = 1 ):
        """Expand the selection by bonded neighbors."""
        if self.system.connectivity.HasFullConnectivity ( ) and ( iterations > 0 ):
            # . Get all connected atoms.
            bonds            = self.system.connectivity.bonds
            currentSelection = self.selection
            for iteration in range ( iterations ):
                updatedSelection = set ( currentSelection )
                for index in currentSelection:
                    updatedSelection.update ( bonds.GetConnectedAtoms ( index ) )
                currentSelection = updatedSelection
            # . Finish up.
            return self.__class__ ( selection = currentSelection, system = self.system )
        else: return self

    def ByComponent ( self ):
        """Expand the selection by component."""
        # . Get all the unique components.
        unique = set ( )
        for index in self.selection:
            unique.add ( self.system.atoms[index].parent )
        # . Get the indices of all atoms.
        selection = set ( )
        for component in unique:
            for atom in component.children: selection.add ( atom.index )
        # . Finish up.
        return self.__class__ ( selection = selection, system = self.system )

    def ByEntity ( self ):
        """Expand the selection by entity."""
        # . Get all the unique entities.
        unique = set ( )
        for index in self.selection:
            unique.add ( self.system.atoms[index].parent.parent )
        # . Get the indices of all atoms.
        selection = set ( )
        for entity in unique:
            for component in entity.children:
                for atom in component.children: selection.add ( atom.index )
        # . Finish up.
        return self.__class__ ( selection = selection, system = self.system )

    def ByIsolate ( self ):
        """Expand the selection by isolate."""
        if self.system.connectivity.HasFullConnectivity ( ):
            # . Get all the unique isolates.
            isolateIndex = self.system.connectivity.isolateIndex
            isolates     = self.system.connectivity.isolates
            unique       = set ( )
            for index in self.selection:
                unique.add ( isolateIndex[index] )
            # . Get the indices of all atoms.
            selection = set ( )
            for index in unique:
                for i in isolates[index]: selection.add ( i )
            # . Finish up.
            return self.__class__ ( selection = selection, system = self.system )
        else: return self

    def ByLinearPolymer ( self ):
        """Expand the selection by linear polymer."""
        # . Get all the unique linear polymers.
        linearPolymerIndex = self.system.sequence.linearPolymerIndex
        unique             = set ( )
        for index in self.selection:
            polymer = linearPolymerIndex.get ( self.system.atoms[index].parent, None )
            if polymer is not None: unique.add ( polymer )
        # . Get the indices of all atoms.
        selection = set ( )
        for index in unique:
            polymer = self.system.sequence.linearPolymers[index]
            for component in polymer.ComponentIterator ( ):
                for atom in component.children: selection.add ( atom.index )
        # . Finish up.
        return self.__class__ ( selection = selection, system = self.system )

    def ByRingSet ( self ):
        """Expand the selection by ring set."""
        if self.system.connectivity.HasFullConnectivity ( ):
            # . Get all the unique ring sets.
            ringSetIndex = self.system.connectivity.ringSetIndex
            ringSets     = self.system.connectivity.ringSets
            unique       = set ( )
            for index in self.selection:
                ringSet = ringSetIndex.get ( index, -1 )
                if ringSet >= 0: unique.add ( ringSet )
            # . Get the indices of all atoms.
            selection = set ( )
            for index in unique:
                for ring in ringSets[index]:
                    for i in ring: selection.add ( i )
            # . Finish up.
            return self.__class__ ( selection = selection, system = self.system )
        else: return self

    def CheckAttributes ( self ):
        """Check attributes."""
        try:
            if not isinstance ( self.selection, set    ): self.selection = set ( self.selection )
            if not isinstance ( self.system   , System ): raise
            if len ( self.selection ) > 0:
                if ( min ( self.selection ) < 0 ) or ( max ( self.selection ) >= len ( self.system.atoms ) ): raise
        except:
            raise AtomSelectionError ( "Invalid selection attributes." )

    def CheckSystems ( self, other ):
        """Check for system compatibility between two selections."""
        if self.system is not other.system: raise AtomSelectionError ( "Incompatible systems in selections." )

    def Complement ( self ):
        """Complement of the selection."""
        return self.__class__ ( selection = ( set ( range ( len ( self.system.atoms ) ) ) - self.selection ), system = self.system )

    @classmethod
    def FromAtomAttributes ( selfClass, system, attribute, value ):
        """Return a selection of atoms for which an attribute has a certain value."""
        data      = system.atoms.GetItemAttributes ( attribute )
        selection = set ( )
        for ( index, datum ) in enumerate ( data ):
            if datum == value: selection.add ( index )
        return selfClass ( selection = selection, system = system )

    @classmethod
    def FromAtomPattern ( selfClass, system, atomPattern ):
        """Get a selection given an atom pattern."""
        # . Initialization.
        selection = set ( )
        sequence  = system.sequence
        # . Parse the pattern to get the label selection fields.
        ( eSFields, cSFields, aSFields ) = sequence.ParseAtomPattern ( atomPattern )
        # . Check for matches.
        for entity in sequence.children:
            if sequence.FieldsLabelMatch ( eSFields, entity.label ):
                for component in entity.children:
                    if sequence.FieldsLabelMatch ( cSFields, component.label ):
                        for atom in component.children:
                            if sequence.FieldsLabelMatch ( aSFields, atom.label ): selection.add ( atom.index )
        # . Finish up.
        return selfClass ( selection = selection, system = system )

    def Summary ( self, log = logFile ):
        """Summarizing."""
        if LogFileActive ( log ):
            log.Paragraph ( "Number of selected atoms = {:d}".format ( len ( self.selection ) ) )

    def Within ( self, distance, excludeSelf = False ):
        """Selection of atoms within a given distance."""
        # . Find all atoms within distance of the current selection.
        selfSelection = self.selection
        system        = self.system
        pairlist      = CrossPairList_FromSingleCoordinates3 ( system.coordinates3, selection1 = Selection.FromIterable ( selfSelection ), safety = distance )
        selection     = set ( selfSelection )
        for ( i, j ) in pairlist: selection.add ( j )
        # . Exclude self?
        if excludeSelf: selection -= self.selection
        # . Finish up.
        return self.__class__ ( selection = selection, system = self.system )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
