#-------------------------------------------------------------------------------
# . File      : Atom.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes for handling atoms."""

import copy

from StringIO import StringIO

from pCore   import logFile, LogFileActive, Real1DArray, Selection, SingleObjectContainer, TreeLeafNode
from Element import Element, PeriodicTable

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Atom ( TreeLeafNode ):
    """An atom."""

    defaultAttributes = { "atomicNumber" :    -1 ,
                          "connections"  :     0 ,
                          "formalCharge" :     0 ,
                          "hydrogens"    :     0 ,
                          "index"        :    -1 ,
                          "isAromatic"   : False ,
                          "valence"      :     0 }
    defaultAttributes.update ( TreeLeafNode.defaultAttributes )

    def __getattribute__ ( self, name ):
        """Check for an atom attribute and then search for an equivalent element attribute if nothing found."""
        try:
            return object.__getattribute__ ( self, name )
        except AttributeError:
            atomicNumber = self.__dict__.get ( "atomicNumber", -1 )
            if hasattr ( PeriodicTable.Element ( atomicNumber ), name ): return PeriodicTable.Element ( atomicNumber ).__dict__[name]
            else: raise AttributeError ( "Unknown atom or elemental attribute: " + name + "." )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AtomContainer ( SingleObjectContainer ):
    """A container class for atoms."""

    def __init__ ( self ):
        """Constructor."""
        super ( AtomContainer, self ).__init__ ( )

    def _SetItemAttributes ( self, mapping ):
        """Set item attributes given a mapping."""
        for ( attribute, data ) in mapping.iteritems ( ):
            if isinstance ( data, dict ):
                for ( index, value ) in data.iteritems ( ):
                    setattr ( self.items[index], attribute, value )
            else:
                for ( index, value ) in enumerate ( data ):
                    setattr ( self.items[index], attribute, value )

    def _SetItemsFromIterable ( self, iterable ):
        """The iterable sequence must consist of Atom objects, atomic numbers or symbols.
        """
        # . Test for an iterable sequence.
        try:    items = iter ( iterable )
        except: TypeError ( "Argument initializer must be an iterable sequence." )
        # . Create the atoms by either using an atom directly or searching for an atomic number.
        atoms     = []
        unlabeled = 0
        for item in items:
            if isinstance ( item, Atom ): atom = item
            else:
                if   isinstance ( item, float ):       atomicNumber = PeriodicTable.AtomicNumberFromMass ( item )
                elif isinstance ( item, int   ):       atomicNumber = item
                elif isinstance ( item, str   ):       atomicNumber = PeriodicTable.AtomicNumber ( item )
                elif hasattr ( item, "atomicNumber" ): atomicNumber = item.atomicNumber
                else: raise TypeError ( "Unrecognized atom initializer: " + str ( item ) + "." )
                atom = Atom ( atomicNumber = atomicNumber )
            atoms.append ( atom )
            if atom.label is None: unlabeled += 1
        self.items = atoms
        # . Make labels and reindex.
        if unlabeled > 0: self.MakeLabels ( )
        self.Reindex ( )

    def ElementDecomposition ( self, selection = None ):
        """Return a dictionary containing the indexes of the atoms of each element."""
        if selection is None: indices = range ( len ( self ) )
        else:                 indices = selection
        decomposition = { }
        for index in indices:
            atom         = self.items[index]
            atomicNumber = atom.atomicNumber
            found        = decomposition.get ( atomicNumber, [] )
            found.append ( index )
            decomposition[atomicNumber] = found
        return decomposition

    def ElementFrequencies ( self, selection = None ):
        """Return a dictionary containing the element frequencies.

        A selection is optional.
        """
        if selection is None: indices = range ( len ( self ) )
        else:                 indices = selection
        frequencies = { }
        for index in indices:
            atom         = self.items[index]
            atomicNumber = atom.atomicNumber
            if atomicNumber in frequencies: frequencies[atomicNumber] = frequencies[atomicNumber] + 1
            else:                           frequencies[atomicNumber] = 1
        return frequencies

    def FormulaString ( self ):
        """Return a formula string."""
        frequencies = self.ElementFrequencies ( )
        keys        = frequencies.keys ( )
        keys.sort ( )
        outstring   = StringIO ( )
        for key in keys: outstring.write ( "{:s}{:d}".format ( PeriodicTable.Symbol ( key ), frequencies[key] ) )
        formula = outstring.getvalue ( )
        outstring.close ( )
        return formula

    @classmethod
    def FromIterable ( selfClass, iterable, attributes = None ):
        """Constructor from iterable."""
        self = selfClass ( )
        self._SetItemsFromIterable ( iterable )
        if attributes is not None: self._SetItemAttributes ( attributes )
        return self

    @classmethod
    def FromMapping ( selfClass, mapping ):
        """Constructor from mapping."""
        atomicNumbers = mapping.pop ( "atomicNumber" )
        return selfClass.FromIterable ( atomicNumbers, attributes = mapping )

    def GetItemAttributes ( self, attributeLabel, asDictionary = False, selection = None ):
        """Return a sequence of item attributes."""
        # . Initialization.
        data = None
        if ( attributeLabel in Atom.defaultAttributes ) or ( attributeLabel in Element.defaultAttributes ):
            if selection is None: indices = range ( len ( self ) )
            else:                 indices = selection
            if indices is not None:
                # . Use a dictionary for non-default values only.
                if asDictionary:
                    data         = {}
                    if attributeLabel in Atom.defaultAttributes: defaultValue =    Atom.defaultAttributes[attributeLabel]
                    else:                                        defaultValue = Element.defaultAttributes[attributeLabel]
                    for i in indices:
                        value = getattr ( self[i], attributeLabel )
                        if value != defaultValue: data[i] = value
                # . Use a sequence for all values.
                else:
                    attribute = getattr ( self[0], attributeLabel )
                    # . Float attribute.
                    if isinstance ( attribute, ( float ) ):
                        data = Real1DArray.WithExtent ( len ( indices ) )
                        for ( i, index ) in enumerate ( indices ): data[i] = getattr ( self[index], attributeLabel )
                    # . Other attribute.
                    else:
                        data = []
                        for i in indices: data.append ( getattr ( self[i], attributeLabel ) )
        return data

    def ItemClass ( self ): return Atom

    def ItemName ( self ): return "Atom"

    def MakeLabels ( self ):
        """Make default atom labels."""
        for ( i, atom ) in enumerate ( self.items ):
            atom.label = PeriodicTable.Symbol ( atom.atomicNumber, index = i )

    def Merge ( self, others, information = {} ):
        """Merging."""
        # . Create the new instance.
	merged =  self.__class__ ( )
        # . Initialization.
        atomIncrements = []
        atomMapping    = {}
        atoms          = []
        numberAtoms    = 0
        # . Loop over the items.
        for item in ( [ self ] + list ( others ) ):
            if isinstance ( item, AtomContainer ):
                # . Create the atom increment.
                atomIncrements.append ( numberAtoms )
                numberAtoms += len ( item )
                # . Copy the atoms.
                for atom in item:
                    newAtom = copy.copy ( atom )
                    atoms.append ( newAtom )
                    atomMapping[atom] = newAtom
            else:
                raise ValueError ( "Invalid atom container in merge." )
        # . Save the atoms and information.
        merged.__dict__["items"] = atoms
        information["atomContainer"  ] = merged
        information["atomIncrements" ] = atomIncrements
        information["atomMapping"    ] = atomMapping
        information["indexIncrements"] = atomIncrements
        # . Finish up.
        return merged

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        atomMapping = {}
        pruned      = self.__class__ ( )
        for ( key, attribute ) in self.__dict__.iteritems ( ):
            if key == "items":
                items = []
                for i in selection:
                    oldAtom = self.items[i]
                    newAtom = copy.copy ( oldAtom )
                    items.append ( newAtom )
                    atomMapping[oldAtom] = newAtom
                pruned.items = items
            else:
                if hasattr ( attribute, "Prune" ):
                    item = attribute.Prune ( selection )
                    if item is not None: pruned.__dict__[key] = item
        information["atomContainer"] = pruned
        information["atomMapping"  ] = atomMapping
        return pruned

    def Reindex ( self ):
        """Reindex the container items."""
        for ( i, atom ) in enumerate ( self.items ):
            atom.index = i

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            frequencies = self.ElementFrequencies ( )
            nhydro = 0
            nother = 0
            for key in frequencies.keys ( ):
                if   key == 1: nhydro = frequencies[key]
                elif key <= 0: nother = nother + frequencies[key]
            natoms = len ( self )
            nheavy = natoms - nhydro - nother
            summary = log.GetSummary ( )
            summary.Start ( "Atom Container Summary" )
            summary.Entry ( "Number of Atoms"       , "{:d}".format ( natoms ) )
            summary.Entry ( "Number of Heavy Atoms" , "{:d}".format ( nheavy ) )
            summary.Entry ( "Number of Hydrogens"   , "{:d}".format ( nhydro ) )
            summary.Entry ( "Number of Unknowns"    , "{:d}".format ( nother ) )
            summary.Stop ( )

    def ToMapping ( self ):
        """Return a mapping for serialization."""
        mapping      = { "atomicNumber" : self.GetItemAttributes ( "atomicNumber" ) }
        formalCharge = self.GetItemAttributes ( "formalCharge", asDictionary = True )
        if len ( formalCharge ) > 0: mapping["formalCharge"] = formalCharge
        return mapping

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
