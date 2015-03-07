#-------------------------------------------------------------------------------
# . File      : Bond.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes for handling bonds in a connectivity."""

import operator

from pCore import Clone, logFile, LogFileActive, SelfPairList_FromCoordinates3, SelfPairList_FromIntegerPairs, SingleObjectContainer, Singleton

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
class BondDefinition ( Singleton ):
    """Base class for bond definitions."""

    defaultattributes = { "bondOrder"  :          -1 ,
                          "isAromatic" :       False ,
                          "label"      : "Undefined" }

    def Initialize ( self ):
        """Initialization."""
        self.__dict__.update ( self.__class__.defaultattributes )

# . Derived bond classes.
class UndefinedBond ( BondDefinition ):

    pass

class NullBond ( BondDefinition ):

    defaultattributes = { "bondOrder"  :      0 ,
                          "isAromatic" :  False ,
                          "label"      : "Null" }

class SingleBond ( BondDefinition ):

    defaultattributes = { "bondOrder"  :        1 ,
                          "isAromatic" :    False ,
                          "label"      : "Single" }

class DoubleBond ( BondDefinition ):

    defaultattributes = { "bondOrder"  :        2 ,
                          "isAromatic" :    False ,
                          "label"      : "Double" }

class TripleBond ( BondDefinition ):

    defaultattributes = { "bondOrder"  :        3 ,
                          "isAromatic" :    False ,
                          "label"      : "Triple" }

class AromaticSingleBond ( BondDefinition ):

    defaultattributes = { "bondOrder"  :                1 ,
                          "isAromatic" :             True ,
                          "label"      : "AromaticSingle" }

class AromaticDoubleBond ( BondDefinition ):

    defaultattributes = { "bondOrder"  :                2 ,
                          "isAromatic" :             True ,
                          "label"      : "AromaticDouble" }

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Bond labels to definition type mapping.
_BondLabelsToType = { "Undefined"      : UndefinedBond      ( ) ,
                      "Null"           : NullBond           ( ) ,
                      "Single"         : SingleBond         ( ) ,
                      "Double"         : DoubleBond         ( ) ,
                      "Triple"         : TripleBond         ( ) ,
                      "AromaticSingle" : AromaticSingleBond ( ) ,
                      "AromaticDouble" : AromaticDoubleBond ( ) }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Bond ( object ):
    """A class to represent a bond."""

    def __init__ ( self, i = -1, j = -1, type = UndefinedBond ( ) ):
        """Constructor."""
        self.i    = i
        self.j    = j
        self.type = type

    def Merge ( self, atomIncrement ):
        """Merging."""
        new = Clone ( self )
        new.i += atomIncrement
        new.j += atomIncrement
        return new

    def Order ( self ):
        """Ordering."""
        if self.i < self.j: self.i, self.j = self.j, self.i

    def Other ( self, i ):
        """Return the index of the other atom."""
        if   self.i == i: return self.j
        elif self.j == i: return self.i
        else:             return -1

    def Prune ( self, selection ):
        """Pruning."""
        new = None
        if ( self.i in selection ) and ( self.j in selection ):
            new   = Clone ( self )
            new.i = selection.Position ( self.i )
            new.j = selection.Position ( self.j )
        return new

    @property
    def sortKey ( self ): return ( self.i, self.j )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class BondContainer ( SingleObjectContainer ):
    """A container class for bonds."""

    def __init__ ( self ):
        """Constructor."""
        super ( BondContainer, self ).__init__ ( )
        self.ClearRepresentations ( )
        self.QSORTED = False

    def _SetItemsFromIterable ( self, iterable ):
        """The iterable sequence must consist of Bond objects or index numbers of atoms and, optionally, a bond type."""
        # . Initialization.
        bonds = []
        # . Dictionary.
        if isinstance ( iterable, dict ):
            for ( key, values ) in iterable.iteritems ( ):
                bondType = _BondLabelsToType[key]
                for ( i, j ) in values:
                    bonds.append ( Bond ( i = i, j = j, type = bondType ) )
        else:
            # . Test for an iterable sequence.
            try:    items = iter ( iterable )
            except: TypeError ( "Argument initializer must be an iterable sequence." )
            # . Create the bonds.
            for item in items:
                # . Use a bond object directly.
                if isinstance ( item, Bond ): bonds.append ( item )
                # . Construct a new object.
                else:
                    new = Bond ( i = item[0], j = item[1] )
                    if len ( item ) > 2:
                        bondType = item[2]
                        if   isinstance ( bondType, BondDefinition ): new.type = bondType
                        elif isinstance ( bondType, basestring     ): new.type = _BondLabelsToType[bondType]
                    bonds.append ( new )
        # . Finish up.
        self.items = bonds
        self.Sort ( )

    def ClearRepresentations ( self ):
        """Clear all representations."""
        self.atomconnections = None
        self.bondconnections = None
        self.upperBound      = -1

    @classmethod
    def FromCoordinates3 ( selfClass, coordinates3, radii = None, safety = 0.45 ):
        """Estimate bonds from coordinates."""
        pairlist = SelfPairList_FromCoordinates3 ( coordinates3, radii = radii, safety = safety )
        self     = BondContainer.FromPairList ( pairlist )
        return self

    @classmethod
    def FromIterable ( selfClass, iterable ):
        """Constructor from iterable."""
        self = selfClass ( )
        self._SetItemsFromIterable ( iterable )
        return self

    @classmethod
    def FromPairList ( selfClass, pairlist ):
        """Constructor from a pairlist."""
        items = []
        for ( i, j ) in pairlist:
            items.append ( Bond ( i = i, j = j ) )
        self       = selfClass ( )
        self.items = items
        self.Sort ( )
        return self

    def GetBond ( self, i, j ):
        """Get a bond given the indices of two atoms."""
        bond = None
        for c in self.GetConnections ( i ):
            connection = self.items[c]
            if connection.Other ( i ) == j:
                bond = connection
                break
        return bond

    def GetConnectedAtoms ( self, i ):
        """Get the atoms connected to an atom."""
        if ( self.atomconnections is not None ):
            try:    iatom = int ( i )
            except: raise TypeError ( "Atom index must be an integer." )
            if ( iatom >= 0 ) and ( iatom < self.upperBound ):
                try:    return self.atomconnections[iatom]
                except: return []
            else: raise IndexError ( "Connection index out of range: {:d}.".format ( iatom ) )
        else: raise AttributeError ( "Connection representation does not exist." )

    def GetConnections ( self, i ):
        """Get the connections for an atom."""
        if ( self.bondconnections is not None ):
            try:    iatom = int ( i )
            except: raise TypeError ( "Atom index must be an integer." )
            if ( iatom >= 0 ) and ( iatom < self.upperBound ):
                try:    return self.bondconnections[iatom]
                except: return []
            else: raise IndexError ( "Connection index out of range." )
        else: raise AttributeError ( "Connection representation does not exist." )

    def IdentifyBoundaryAtoms ( self, selection, results ):
        """Return a dictionary of boundary atoms and their boundary atom, out-of-selection and in-selection partners.

        |results| is a dictionary that must exist on entry.
        """
        # . Initialization.
        if len ( selection ) > 0:
            # . Make sure the connection representation exists with a sufficient upperBound.
            self.MakeConnections ( upperBound = max ( selection.UpperBound ( ), self.UpperBound ( ) ) )
            # . Loop to find the boundary atoms.
            boundaryatoms = set ( )
            for i in selection:
                for j in self.GetConnectedAtoms ( i ):
                    if j not in selection: boundaryatoms.add ( j ) # . Boundary atoms are outside the selection.
            # . Boundary atoms exist.
            if len ( boundaryatoms ) > 0:
                # . Find the partners of the boundary atoms.
                for b in boundaryatoms:
                    data = results.get ( b, [ set ( ), set ( ), set ( ) ] )
                    for i in self.GetConnectedAtoms ( b ):
                        if   ( i in boundaryatoms ): data[0].add ( i ) # . Boundary atom partners.
                        elif ( i not in selection ): data[1].add ( i ) # . Out-of-selection partners.
                        else:                        data[2].add ( i ) # . In-selection partners.
                    results[b] = data

    def ItemClass ( self ): return Bond

    def ItemName ( self ): return "Bond"

    def MakeConnections ( self, upperBound = -1 ):
        """Make the connections representation of the bond container."""
        # . Make the representation.
        if ( self.atomconnections is None ) or ( self.bondconnections is None ):
            self.ClearRepresentations ( )
            self.Sort ( )
            # . Set the upper bound of the representation.
            self.upperBound = max ( self.UpperBound ( ), upperBound )
            # . Initialization.
            self.atomconnections = {}
            self.bondconnections = {}
            for ( b, bond ) in enumerate ( self.items ):
                i = bond.i
                j = bond.j
                if i not in self.atomconnections: self.atomconnections[i] = []
                if j not in self.atomconnections: self.atomconnections[j] = []
                if i not in self.bondconnections: self.bondconnections[i] = []
                if j not in self.bondconnections: self.bondconnections[j] = []
                self.atomconnections[i].append ( j )
                self.atomconnections[j].append ( i )
                self.bondconnections[i].append ( b )
                self.bondconnections[j].append ( b )
            # . Sorting.
            for value in self.atomconnections.values ( ): value.sort ( )
            for value in self.bondconnections.values ( ): value.sort ( )
        # . Check the upperbound.
        else:
            if upperBound > self.upperBound: self.upperBound = upperBound

    def Merge ( self, others, information = {} ):
        """Merging."""
        # . Initialization.
        merged = None
        # . Get the increments.
        increments = information.get ( "atomIncrements", None )
        if increments is not None:
            # . Create the new instance.
            merged = self.__class__ ( )
            # . Loop over the items.
            members = []
            for ( item, increment ) in zip ( ( [ self ] + list ( others ) ), increments ):
                if isinstance ( item, BondContainer ):
                    # . Copy the atoms.
                    for member in item:
                        newMember = member.Merge ( increment )
                        members.append ( newMember )
            # . Save the members.
            merged.__dict__["items"] = members
            merged.ClearRepresentations ( )
            merged.Sort ( )
        # . Finish up.
        return merged

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        new = self.__class__ ( )
        for ( key, attribute ) in self.__dict__.iteritems ( ):
            if key == "items":
                items = []
                for item in self.items:
                    newitem = item.Prune ( selection )
                    if newitem is not None: items.append ( newitem )
                new.items = items
            else:
                if hasattr ( attribute, "Prune" ):
                    item = attribute.Prune ( selection, information = information )
                    if item is not None: new.__dict__[key] = item
        new.QSORTED = self.QSORTED
        return new

    def Sort ( self ):
        """Sorting."""
        QSORTED = self.__dict__.get ( "QSORTED", False )
        if not QSORTED:
            items = self.items
            for item in items: item.Order ( )
            items.sort ( key = operator.attrgetter ( "sortKey" ) )
            for i in range ( len ( self ), 1, -1 ):
                if items[i-1].sortKey == items[i-2].sortKey: del items[i-1]
            self.QSORTED = True

    def ToIndices ( self ):
        """Return a list of bond indices."""
        indices = []
        for bond in self.items:
            indices.append ( bond.i )
            indices.append ( bond.j )
        return indices

    def ToMapping ( self ):
        """Return a mapping for serialization."""
        mapping = {}
        for bond in self.items:
            label = bond.type.label
            items = mapping.get ( label, None )
            if items is None:
                items = []
                mapping[label] = items
            items.append ( [ bond.i, bond.j ] )
        return mapping

    def ToPairList ( self ):
        """Return a pairlist of bond indices."""
        indices = self.ToIndices ( )
        return SelfPairList_FromIntegerPairs ( indices )

    def UpperBound ( self ):
        """Return the upperbound of the container (the highest index + 1)."""
        self.Sort ( )
        if len ( self ) > 0: return self.items[-1].i + 1
        else:                return 0

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
