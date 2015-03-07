#-------------------------------------------------------------------------------
# . File      : Connectivity.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A "connectivity" is a set of atoms and bonds and related information such as angles, dihedrals, molecules and ringsets."""

import itertools

from pCore       import FindRingSets, logFile, LogFileActive, SelfPairList_FromIntegerPairs

from Angle       import AngleContainer
from Aromaticity import DetermineAromaticAtoms # . Can eventually move to RingSetContainer?
from Atom        import AtomContainer
from Bond        import BondContainer
from Dihedral    import DihedralContainer
from Element     import PeriodicTable

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Connectivity ( object ):
    """A class to store connectivity data."""

    # . Attributes that can be derived from atom and bond information.
    derivableAttributes = [ "_angles"   , \
                            "_dihedrals", \
                            "_isolates" , \
                            "_ringSets"   ]

    # . Local attributes.
    localAttributes = [ "_atoms", "_bonds" ] + derivableAttributes

    def __init__ ( self ):
        """Constructor."""
        self.Initialize ( )

    def __len__ ( self ):
        """Return the number of atoms in the connectivity."""
        if self.atoms is None: return 0
        else:                  return len ( self.atoms )

    def BondsFromCoordinates ( self, coordinates, radii, safety ):
        """Estimate bonds from coordinates using a distance search."""
        self.__dict__["_bonds"] = BondContainer.FromCoordinates3 ( coordinates, radii = radii, safety = safety )
        self.CompleteConnectivity ( )

    def ClearConnectivity ( self, clearBonds = False ):
        """Clear the connectivity except atom and (optionally) bond information."""
        if clearBonds: self.__dict__.pop ( "_bonds", None )
        for attribute in Connectivity.derivableAttributes:
            self.__dict__.pop ( attribute, None )

    def CompleteConnectivity ( self ):
        """Build the connectivity from atom and bond information."""
        if ( self.atoms is not None ) and ( self.bonds is not None ):
            self.bonds.MakeConnections ( upperBound = len ( self.atoms ) )
            bondpairlist = self.bonds.ToPairList ( )
#            print "Bond Pair List:"
#            for ( i, j ) in bondpairlist: print i, j
            # . Quantities derived from bonds.
            self.__dict__["_angles"   ] = AngleContainer.FromBondContainer    ( self.bonds )
            self.__dict__["_dihedrals"] = DihedralContainer.FromBondContainer ( self.bonds )
            self.__dict__["_isolates" ] = bondpairlist.GetIsolates ( upperBound = len ( self.atoms ) )
            self.isolates.itemName      = "Isolate"
# . A fix until a better ring-finding algorithm is available.
            if self.__dict__.get ( "_ringSets", None ) is None: self.__dict__["_ringSets"] = self.FindRingSets ( )
            # . Set atom attributes.
            if len ( self.bonds ) > 0:
                # . Connections and valence.
                for ( i, atom ) in enumerate ( self.atoms ):
                    atom.connections = len ( self.bonds.GetConnectedAtoms ( i ) )
                    atom.hydrogens   = len ( self.GetConnectedHydrogens   ( i ) )
                    atom.valence     = 0
                    for b in self.bonds.GetConnections ( i ): atom.valence += self.bonds[b].type.bondOrder
                # . Aromaticity.
                DetermineAromaticAtoms ( self )
            else:
                for atom in self.atoms:
                    atom.connections = 0
                    atom.hydrogens   = 0
                    atom.valence     = 0

    def FindRingSets ( self ):
        """Find ring sets."""
        # . Connections are cloned in FindRingSets so will remain unchanged.
        bonds    = self.bonds
        ringSets = None
        if ( bonds.atomconnections is not None ):
            isolates = self.__dict__.get ( "_isolates", None )
            # . Slow version.
            if isolates is None:
                ringSets = FindRingSets ( bonds.atomconnections )
            # . Less slow version if isolates are present.
            else:
                connections = {}
                ringSets    = []
                for isolate in isolates:
                    ntotal =  0
                    for i in isolate:
                        localBonds = bonds.atomconnections.get ( i, None )
                        if localBonds is not None:
                            connections[i] = localBonds
                            ntotal        += len ( localBonds )
                    # . For an isolate there are no rings if the (number of connections)/2 is less than the number of atoms.
                    if ( ntotal // 2 ) >= len ( isolate ):
                        newSets = FindRingSets ( connections )
                        if len ( newSets ) > 0: ringSets.extend ( newSets )
                    # . Clear up.
                    connections.clear ( )
        return ringSets

    @classmethod
    def FromAtomContainer ( selfClass, atoms, bonds = None, ringSets = None ):
        """Constructor an atom container and optionally, bonds ."""
        self = selfClass ( )
        self.__dict__["_atoms"] = atoms
        if bonds is not None:
            self.__dict__["_bonds"] = BondContainer.FromIterable ( bonds )
        if ringSets is not None: self.__dict__["_ringSets"] = ringSets # . Temporary.
        self.CompleteConnectivity ( )
        return self

    def Get12Interactions ( self ):
        """Get a list of 12 interactions."""
        i12 = None
        if self.HasFullConnectivity ( ):
            # . 12 lists.
            i12 = set ( )
            for item in self.bonds:
                i = item.i
                j = item.j
                i12.add ( ( max ( i, j ), min ( i, j ) ) )
        return i12

    def Get13Interactions ( self ):
        """Get a list of 13 interactions.

        It is necessary to exclude possible 12 interactions (as in 3-membered rings).
        """
        i13 = None
        if self.HasFullConnectivity ( ):
            # . 12 interactions.
            i12 = self.Get12Interactions ( )
            # . Possible 13 interactions.
            i13p = set ( )
            for item in self.angles:
                i = item.i
                k = item.k
                i13p.add ( ( max ( i, k ), min ( i, k ) ) )
            # . Get the final list.
            i13 = i13p.difference ( i12 )
        return i13

    def GetAngleIndices ( self ):
        """Get the angle indices."""
        ijk = []
        if self.angles is not None:
            for angle in self.angles:
                ijk.append ( ( angle.i, angle.j, angle.k ) )
        return ijk

    def GetConnectedHydrogens ( self, i ):
        """Get the hydrogens connected to an atom."""
        h = []
        for j in self.bonds.GetConnectedAtoms ( i ):
            if self.atoms[j].atomicNumber == 1: h.append ( j )
        return h

    def GetDihedralIndices ( self ):
        """Get the dihedral angle indices."""
        ijkl = []
        if self.dihedrals is not None:
            for dihedral in self.dihedrals:
                ijkl.append ( ( dihedral.i, dihedral.j, dihedral.k, dihedral.l ) )
        return ijkl

    def GetPlanarAtomIndices ( self ):
        """Get the indices of the planar atoms and surrounding atoms."""
        ijkl = []
        # . Get the number of atoms with 3 connections.
        atomsc3 = []
        for ( i, atom ) in enumerate ( self.atoms ):
            if atom.connections == 3: atomsc3.append ( i )
        # . Assume that these are planar.
        if len ( atomsc3 ) > 0:
            for i in atomsc3:
                jkl = self.bonds.GetConnectedAtoms ( i )
                ijkl.append ( [ i ] + jkl )
        return ijkl

    def HasFullConnectivity ( self ):
        """Return a flag indicating whether the connectivity is fully defined."""
        QOK = ( self.atoms     is not None ) and ( self.bonds    is not None ) and ( self.angles   is not None ) and \
              ( self.dihedrals is not None ) and ( self.isolates is not None ) and ( self.ringSets is not None )
        return QOK

    def Initialize ( self ):
        """Initialization."""
        for key in self.__class__.localAttributes: setattr ( self, key, None )

    def Make1234And14PairLists ( self ):
        """Make 1234 and 14 pairlists.

        It is necessary to ensure that the 14 list does not include any 123 interactions.
        """
        i1234 = None
        i14   = None
        if self.HasFullConnectivity ( ):
            # . 123 lists.
            i123 = set ( )
            for item in self.bonds:
                i = item.i
                j = item.j
                i123.add ( ( max ( i, j ), min ( i, j ) ) )
            for item in self.angles:
                i = item.i
                k = item.k
                i123.add ( ( max ( i, k ), min ( i, k ) ) )
            # . Possible 14 lists.
            i14p = set ( )
            for item in self.dihedrals:
                i = item.i
                l = item.l
                i14p.add ( ( max ( i, l ), min ( i, l ) ) )
            # . Get the final lists.
            i1234 = i123.union ( i14p )
            i14   = i14p.difference ( i123 )
            # . Generate the lists even if there are no interactions.
            i1234 = list ( itertools.chain ( * list ( i1234 ) ) )
            i1234 = SelfPairList_FromIntegerPairs ( i1234 )
            i14   = list ( itertools.chain ( * list ( i14   ) ) )
            i14   = SelfPairList_FromIntegerPairs ( i14   )
        return ( i1234, i14 )

    def Merge ( self, others, information = {} ):
        """Merging."""
        toMergeItems  = [ self ] + list ( others )
        merged        = self.__class__ ( )
        atomContainer = information.get ( "atomContainer", None )
        if atomContainer is None:
            atoms         = [ getattr ( item, "atoms" ) for item in toMergeItems ]
            atomContainer = atoms[0].Merge ( atoms[1:], information = information )
        merged.__dict__["_atoms"] = atomContainer
        if self.HasFullConnectivity ( ):
            bonds    = [ getattr ( item, "bonds"    ) for item in toMergeItems ]
            ringSets = [ getattr ( item, "ringSets" ) for item in toMergeItems ]
            merged.__dict__["_bonds"   ] = bonds[0].Merge     ( bonds[1:], information = information )
            merged.__dict__["_ringSets"] = self.MergeRingSets ( ringSets, information.get ( "atomIncrements", None ) )
            merged.CompleteConnectivity ( )
        return merged

    def MergeRingSets ( self, ringSets, increments ):
        """Merge ring sets."""
        merged = None
        if increments is not None:
            merged = []
            for ( oldSets, increment ) in zip ( ringSets, increments ):
                for oldSet in oldSets:
                    newSet = []
                    for oldRing in oldSet:
                        newRing = []
                        for i in oldRing:
                            newRing.append ( i + increment )
                        newSet.append ( newRing )
                    merged.append ( newSet )
                    #print oldSet, increment, newSet, merged
        return merged

    def NumberOfAngles ( self ):
        """The number of angles."""
        if self.angles is None: return 0
        else:                   return len ( self.angles )

    def NumberOfAtoms ( self ):
        """The number of atoms."""
        if self.atoms is None: return 0
        else:                  return len ( self.atoms )

    def NumberOfBonds ( self ):
        """The number of bonds."""
        if self.bonds is None: return 0
        else:                  return len ( self.bonds )

    def NumberOfDihedrals ( self ):
        """The number of dihedrals."""
        if self.dihedrals is None: return 0
        else:                      return len ( self.dihedrals )

    def NumberOfIsolates ( self ):
        """The number of isolates."""
        if self.isolates is None: return 0
        else:                     return len ( self.isolates )

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        pruned        = self.__class__ ( )
        atomContainer = information.get ( "atomContainer", None )
        if atomContainer is None: atomContainer = self.atoms.Prune ( selection, information = information )
        pruned.__dict__["_atoms"] = atomContainer
        if self.HasFullConnectivity ( ):
            pruned.__dict__["_bonds"   ] = self.bonds.Prune   ( selection, information = information )
            pruned.__dict__["_ringSets"] = self.PruneRingSets ( selection, information = information )
            pruned.CompleteConnectivity ( )
        return pruned

    def PruneRingSets ( self, selection, information = {} ):
        """Prune ring sets."""
        pruned = []
        for oldSet in self.ringSets:
            newSet = []
            for oldRing in oldSet:
                include = True
                for i in oldRing:
                    if i not in selection:
                        include = False
                        break
                if include:
                    newRing = []
                    for i in oldRing:
                        newRing.append ( selection.Position ( i ) )
                    newSet.append ( newRing )
            if len ( newSet ) > 0: pruned.append ( newSet )
        return pruned

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "Connectivity Summary" )
            if self.atoms     is not None: summary.Entry ( "Atoms"     , "{:d}".format ( len ( self.atoms     ) ) )
            if self.bonds     is not None: summary.Entry ( "Bonds"     , "{:d}".format ( len ( self.bonds     ) ) )
            if self.angles    is not None: summary.Entry ( "Angles"    , "{:d}".format ( len ( self.angles    ) ) )
            if self.dihedrals is not None: summary.Entry ( "Dihedrals" , "{:d}".format ( len ( self.dihedrals ) ) )
            if self.isolates  is not None: summary.Entry ( "Isolates"  , "{:d}".format ( len ( self.isolates  ) ) )
            if self.ringSets  is not None: summary.Entry ( "Ring Sets" , "{:d}".format ( len ( self.ringSets  ) ) )
            summary.Stop ( )

    def ToMapping ( self ):
        """Return a minimal mapping for serialization."""
        mapping = { "atoms" : self.atoms.ToMapping ( ) }
        if self.bonds    is not None: mapping["bonds"   ] = self.bonds.ToMapping ( )
        if self.ringSets is not None: mapping["ringSets"] = self.ringSets
        return mapping

    # . Properties.
    @property
    def angles    ( self ): return self.__dict__["_angles"   ]
    @property
    def atoms     ( self ): return self.__dict__["_atoms"    ]
    @property
    def bonds     ( self ): return self.__dict__["_bonds"    ]
    @property
    def dihedrals ( self ): return self.__dict__["_dihedrals"]
    @property
    def isolates  ( self ): return self.__dict__["_isolates" ]
    @property
    def ringSets  ( self ): return self.__dict__["_ringSets" ]

    @property
    def isolateIndex ( self ):
        if "_isolateIndex" not in self.__dict__:
            index = [ -1 for i in range ( len ( self.atoms ) ) ]
            for ( n, isolate ) in enumerate ( self.isolates ):
                for i in isolate: index[i] = n
            self.__dict__["_isolateIndex"] = index
        return self.__dict__["_isolateIndex"]
    @property
    def ringSetIndex ( self ):
        if "_ringSetIndex" not in self.__dict__:
            index = {}
            for ( n, ringSet ) in enumerate ( self.ringSets ):
                for ring in ringSet:
                    for i in ring: index[i] = n
            self.__dict__["_ringSetIndex"] = index
        return self.__dict__["_ringSetIndex"]

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
