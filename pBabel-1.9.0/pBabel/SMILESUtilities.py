#-------------------------------------------------------------------------------
# . File      : SMILESUtilities.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Miscellaneous SMILES classes and parameters."""

from pCore   import FindRingSets, logFile, LogFileActive
from pMolecule import AromaticSingleBond, AromaticDoubleBond, DoubleBond, NullBond, SingleBond, TripleBond, System

#===============================================================================
# . SMILES parameters.
#===============================================================================
# . Chirality data.
# . The allowed number of connections for each class.
CHIRALITYCLASSCONNECTIONS = { "AL" : 2, "OH" : 6, "SP" : 4, "TB" : 5, "TH" : 4 }

# . The default class for each connectivity.
CHIRALITYDEFAULTCLASSES = { 2 : "AL", 4 : "TH", 5 : "TB", 6 : "OH" }

# . The maximum number permitted for each class.
CHIRALITYMAXIMUMNUMBERS = { "AL" : 2, "OH" : 30, "SP" : 3, "TB" : 20, "TH" : 2 }

# . Element data.
# . Atomic numbers of aromatic elements.
ELEMENTSAROMATIC = ( 6, 7, 8, 15, 16, 33, 34 )

# . Organic elements.
ELEMENTSORGANIC  = ( 5, 6, 7, 8, 9, 15, 16, 17, 35, 53 )

# . The maximum possible valencies for aromatic atoms.
VALENCIESAROMATIC = { 6 : 4, 7 : 5, 8 : 2, 15 : 5, 16 : 6, 33 : 5, 34 : 6 }

# . The possible valencies for atoms in the organic subset when connected hydrogens are specified implicitly.
# . Note that the minimum number of hydrogens are added so as to be consistent with explicitly-specified bonds.
VALENCIESORGANIC = { 5: ( 3, ), 6: ( 4, ), 7: ( 3, 5 ), 8: ( 2, ), 9: ( 1, ), 15: ( 3, 5 ), 16: ( 2, 4, 6 ), 17: ( 1, ), 35: ( 1, ), 53: ( 1, ) }

#===============================================================================
# . SMILES atom class.
#===============================================================================
class SMILESAtom ( object ):
    """A class to represent a SMILES atom."""

    # . Default attributes.
    defaultattributes = { "isAromatic"         : False,
                          "QREDUCED"          : False,
                          "QRING"             : False,
                          "aromaticelectrons" :     0,
                          "aromaticvalence"   :     0,
                          "atomicNumber"      :    -1,
                          "chiralityclass"    :  None,
                          "chiralitynumber"   :     0,
                          "connections"       :     0,
                          "formalCharge"      :     0,
                          "implicithydrogens" :     0,
                          "index"             :    -1,
                          "isotope"           :     0,
                          "valence"           :     0 }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultattributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in keywordArguments.iteritems ( ):                 setattr ( self, key, value )

    def AddImplicitHydrogens ( self ):
        """Add implicit hydrogens to an atom of the organic subset."""
        if ( self.QREDUCED ) and ( self.atomicNumber in ELEMENTSORGANIC ):
            # . Get the existing connections and valence minus implicit hydrogens.
            connections = self.connections - self.implicithydrogens
            valence     = self.valence     - self.implicithydrogens
            # . Determine the minimum number of hydrogens that it is necessary to add to satisfy the atom's valence.
            for v in VALENCIESORGANIC[self.atomicNumber]:
                hcount = v - valence
                if hcount >= 0: break
            # . Set the hydrogen count (and modify other data) only if there is a non-negative hydrogen count.
            if hcount >= 0:
                self.implicithydrogens = hcount
                self.connections       = connections + hcount
                self.valence           = valence     + hcount
        self.QREDUCED = False

    def IsPossibleAromatic ( self ):
        """Check to see if the atom is a possible aromatic.

        This procedure does (and should) not require the presence of implicit
        hydrogens to work.

        All aromatic atoms obey the following rules:

        1. It must be one of the possible aromatic elements.
        2. It cannot have a connectivity of more than 3 (including implicit
           hydrogens - although these are not important for determining
           aromaticity).
        3. Carbons must be sp or sp2 hybridized (i.e. have attached multiple
           bonds or be only 3-connected) but there are no such limitations on
           other atoms.
        4. It must be in a ring and bound to at least two other possible
           aromatic elements that are also in rings.

        Element-specific rules:

        C: The third connection can be a single bond to any other atom (I),
           although double bonds to non-carbon, non-ring atoms are allowed
           too (as long as the total valence is 4) (II).

        N, P, As: These atoms can be 3-coordinate and 3-valence (I) or
                  2-coordinate or 3-coordinate and 4-valence in which case
                  the double bond must be to an extra-cyclic oxygen (II).

        O: Oxygen can only ever have two connections (I).

        S, Se: These atoms can be 2-coordinate or 3-coordinate and
               4-valence with the third connection a double bond to an
               extra-cyclic oxygen (I).

        Donated electrons and extra valencies:

        Element            Case     Donated Electrons      Extra Valence

        C                    I              1                    1
                             II             0                    0
        N, P, As             I              2                    0
                             II             1                    1
        O                    I              2                    0
        S, Se                I              2                    0

        The number of donated electrons is augmented or diminished by the
        charge on the atom.

        User-specified aromatics must obey the same rules. If they do not an
        error is produced. A user-specified aromatic does, however, indicate
        that a carbon is sp/sp2 hybridized even it only has single bonds. It
        also changes slightly the way in which N, P and As are handled.
        """
        # . Initialization.
        aromaticelectrons = 0
        QPOSSIBLE         = False

        # . The atom must be one of the resonant elements.
        if self.atomicNumber in ELEMENTSAROMATIC:

            # . Determine some information about neighbouring atoms and bonds.
            # . Initialization.
            naromatic  = 0
            ndbatoms   = 0
            ndboxygens = 0
            valence    = self.implicithydrogens

            # . Loop over the bonds.
            for bond in self.bonds:
                other    = bond.OtherAtom ( self )
                valence += bond.type.bondOrder
                if ( other.atomicNumber in ELEMENTSAROMATIC ) and ( other.QRING ): naromatic += 1
                if ( other.atomicNumber != 6 ) and ( other.QRING ) and ( bond.type.bondOrder == 2 ):
                    ndbatoms += 1
                    if ( other.atomicNumber == 8 ): ndboxygens += 1

            # . Set the connections for the self.
            connections = len ( self.bonds ) + self.implicithydrogens

            # . In the general case, the number of attached in-ring aromatic atoms must be at least two,
            # . the number of connections (including explicitly-attached hydrogens) must not exceed three
            # . and the valence must not exceed the maximum possible.
            if ( self.QRING ) and ( naromatic >= 2 ) and ( connections <= 3 ) and ( valence <= VALENCIESAROMATIC[self.atomicNumber] ):

                # . Process rules for specific atoms.
                # . C.
                if self.atomicNumber == 6:
                    # . sp/sp2 hybridized (or existing aromatic) and a maximum valence of 4.
                    if ( ( valence > connections ) or ( self.isAromatic ) ) and ( valence <= 4 ):
                        # . The atom can be aromatic.
                        QPOSSIBLE = True
                        # . Extra-cyclic double-bond.
                        if ndbatoms > 0: aromaticelectrons = 0
                        # . Normal case.
                        else:            aromaticelectrons = 1
                # . N, P, As.
                elif self.atomicNumber in ( 7, 15, 33 ):
                    # . Existing aromatic.
                    if self.isAromatic:
                        # . 3-coordinate/3-valence.
                        if ( connections == 3 ) and ( valence == 3 ):
                            aromaticelectrons  = 2
                            QPOSSIBLE = True
                        # . 2-coordinate/2-valence or 3-coordinate/4-valence with an extracyclic doubly-bond oxygen.
                        elif ( ( connections == 2 ) and ( valence == 2 ) ) or ( ( connections == 3 ) and ( valence == 4 ) and ( ndboxygens > 0 ) ):
                            aromaticelectrons = 1
                            QPOSSIBLE         = True
                    # . Not specified as aromatic
                    else:
                        # . 3-coordinate/3-valence (with or without an implicit hydrogen).
                        if ( ( connections == 2 ) and ( valence == 2 ) ) or ( ( connections == 3 ) and ( valence == 3 ) ):
                            aromaticelectrons = 2
                            QPOSSIBLE         = True
                        # . 2-coordinate/2-valence or 3-coordinate/4-valence with an extracyclic doubly-bond oxygen.
                        elif ( ( connections == 2 ) and ( valence == 3 ) ) or ( ( connections == 3 ) and ( valence == 5 ) and ( ndboxygens > 0 ) ):
                            aromaticelectrons  = 1
                            QPOSSIBLE          = True
                # . O.
                elif self.atomicNumber == 8:
                    # . 2-coordinate/2-valence.
                    if ( connections == 2 ) and ( valence == 2 ):
                        aromaticelectrons = 2
                        QPOSSIBLE         = True
                # . S, Se.
                elif self.atomicNumber in ( 16, 34 ):
                    # . 2-coordinate/2-valence or 3-coordinate/4-valence with an extra-cyclic oxygen.
                    if ( ( connections == 2 ) and ( valence == 2 ) ) or ( ( connections == 3 ) and ( valence == 4 ) and ( ndboxygens > 0 ) ):
                        aromaticelectrons = 2
                        QPOSSIBLE         = True

            # . Make some final checks for possible aromatic atoms.
            if QPOSSIBLE:
                # . Adjust the number of aromatic electrons by the atom's formal charge.
                aromaticelectrons -= self.formalCharge
                # . To be aromatic the number of electrons must be 0, 1 or 2.
                if aromaticelectrons in ( 0, 1, 2 ):
                    # . Determine the aromatic valence.
                    if aromaticelectrons == 1: aromaticvalence = 1
                    else:                      aromaticvalence = 0
                    # . Set the variables unless the maximum valence has been exceeded for existing aromatic atoms.
                    if self.isAromatic and ( aromaticvalence + valence > VALENCIESAROMATIC[self.atomicNumber] ):
                        QPOSSIBLE = False
                    else:
                        self.aromaticelectrons = aromaticelectrons
                        self.aromaticvalence   = aromaticvalence
                else:
                    QPOSSIBLE = False

        # . Return the flag.
        return QPOSSIBLE

    def SetValence ( self ):
        """Determine the connections and valence for the atom."""
        self.connections = self.implicithydrogens
        self.valence     = self.implicithydrogens
        if self.isAromatic: self.valence += self.aromaticvalence
        if len ( self.bonds ) > 0:
            for bond in self.bonds:
                self.connections += 1
                self.valence     += bond.type.bondOrder

    def VerifyAromaticity ( self ):
        """Verify an aromatic atom using local rules."""
        if ( self.isAromatic ) and ( not self.IsPossibleAromatic ( ) ): raise SMILESConnectivityError ( "Invalid aromatic atom.", atom = self )

    def VerifyChirality ( self ):
        """Verify the chirality of an atom but do nothing with it."""
        if self.chiralityclass is not None:
            # . Reduced chirality class.
            if self.chiralityclass == "??":
                self.chiralityclass = CHIRALITYDEFAULTCLASSES.get ( self.connections, None )
                if self.chiralityclass is None: raise SMILESConnectivityError ( "An atom has a number of connections for which there is no default chirality class.", atom = self )
           # . Check the number of connections and the number of the class.
            if self.connections != CHIRALITYCLASSCONNECTIONS[self.chiralityclass]:
                raise SMILESConnectivityError ( "An atom has a number of connections which is incompatible with the chirality class.", atom = self )
            if ( self.chiralitynumber <= 0 ) or ( self.chiralitynumber > CHIRALITYMAXIMUMNUMBERS[self.chiralityclass] ):
                raise SMILESConnectivityError ( "The number of a specified chirality class is invalid.", atom = self )

#===============================================================================
# . SMILES bond class.
#===============================================================================
class SMILESBond ( object ):
    """A class to represent a SMILES bond."""

    def __init__ ( self, **keywordArguments ):
        for ( key, value ) in keywordArguments.iteritems ( ):
            setattr ( self, key, value )

    def OtherAtom ( self, atom ):
        """Return the other atom of the bond."""
        if self.atom1 is atom: return self.atom2
        else:                  return self.atom1

#===============================================================================
# . SMILES connectivity class.
#===============================================================================
class SMILESConnectivity ( object ):
    """A class to represent a SMILES connectivity."""

    def __init__ ( self ):
        """Constructor."""
        self.atoms    = []
        self.bonds    = []
        self.ringsets = []

    def __len__ ( self ):
        """The number of atoms."""
        return len ( self.atoms )

    def AddAtom ( self, atom ):
        """Add an atom."""
        if not hasattr ( atom, "bonds" ): atom.bonds = []
        self.atoms.append ( atom )

    def AddBond ( self, bond ):
        """Add a bond."""
        # . Check whether a bond of the same type already exists.
        QADD = True
        for oldbond in bond.atom1.bonds:
            if bond.atom2 is oldbond.OtherAtom ( bond.atom1 ):
                if bond.type is oldbond.type: QADD = False
                else: raise SMILESConnectivityError ( "Adding a second bond of different type between two atoms.", bond = oldbond )
        # . Add the bond.
        if QADD:
            self.bonds.append ( bond )
            bond.atom1.bonds.append ( bond )
            bond.atom2.bonds.append ( bond )

    def Aromaticize ( self ):
        """Aromaticize the connectivity."""
        pass

    def Finalize ( self ):
        """Create a final connectivity.

        To be used once all AddAtom and AddBond calls have been made.
        """
        # . Define the rings and flag ring atoms.
        self.FindRings                ( )
        # . Atom checks.
        for atom in self.atoms:
            atom.VerifyAromaticity    ( )
            atom.SetValence           ( )
            atom.AddImplicitHydrogens ( )
            atom.VerifyChirality      ( )
        # . Bond checks - change non-aromatic to aromatic bonds between atoms.
        for bond in self.bonds:
            if bond.atom1.isAromatic and bond.atom2.isAromatic and ( not bond.type.isAromatic ):
                if   bond.type is SingleBond ( ): bond.type = AromaticSingleBond ( )
                elif bond.type is DoubleBond ( ): bond.type = AromaticDoubleBond ( )
                else: raise SMILESConnectivityError ( "Invalid bond type between aromatic atoms.", bond = bond )
        # . Global checks.
        self.VerifyAromaticity        ( )
        self.RemoveExplicitHydrogens  ( )
        self.IndexAtoms               ( )

    def FindRings ( self ):
        """Find ringsets and flag ring atoms."""
        self.ringsets = []
        if len ( self.atoms ) > 0:
            for atom in self.atoms: atom.QRING = False
            nodetable = self.IntegerRepresentation ( )
            ringsetsi = FindRingSets ( nodetable )
            if len ( ringsetsi ) > 0:
                for ringseti in ringsetsi:
                    ringset = []
                    for ringi in ringseti:
                        ring = []
                        for i in ringi: ring.append ( self.atoms[i] )
                        ringset.append ( ring )
                    self.ringsets.append ( ringset )
                for ringset in self.ringsets:
                    for ring in ringset:
                        for atom in ring: atom.QRING = True

    @classmethod
    def FromSystem ( selfClass, system ):
        """Create a connectivity from a system."""
        self = selfClass ( )
        for atom in system.atoms:
            self.AddAtom ( SMILESAtom ( **atom.__dict__ ) )
        for bond in system.connectivity.bonds:
            self.AddBond ( SMILESBond ( atom1 = self.atoms[bond.i], atom2 = self.atoms[bond.j], type = bond.type ) )
        self.Finalize ( )
        return self

    def IndexAtoms ( self ):
        """Index the atoms."""
        for ( i, atom ) in enumerate ( self.atoms ): atom.index = i

    def IntegerRepresentation ( self ):
        """Convert the connectivity to an integer representation."""
        nodetable = {}
        if len ( self.atoms ) > 0:
            self.IndexAtoms ( )
            for atom in self.atoms:
                table = [ ]
                for bond in atom.bonds:
                    other = bond.OtherAtom ( atom )
                    table.append ( other.index )
                if len ( table ) > 0: nodetable[atom.index] = table
        return nodetable

    def Kekularize ( self ):
        """De-aromaticize the connectivity."""
        pass

    def RemoveExplicitHydrogens ( self ):
        """Remove explicitly-specified hydrogens if possible.

        A hydrogen is removable if there are no attached hydrogens, no isotope
        specification, only a single valence and the attached atom has no chiral
        specification. Nothing is done if the atom's bond lists do not exist.
        Note that in the present version, chirality is not handled (as it is not
        stored).
        """
        hydrogens = []
        for atom in self.atoms:
            if ( atom.atomicNumber == 1 ) and ( atom.implicithydrogens == 0 ) and ( atom.isotope == 0 ) and ( atom.valence == 1 ) and ( len ( atom.bonds ) == 1 ):
                bond = atom.bonds[0]
                host = bond.OtherAtom ( atom )
                if host.atomicNumber != 1:
                    host.formalCharge      += atom.formalCharge
                    host.implicithydrogens += 1
                    host.bonds.remove ( bond )
                    self.bonds.remove ( bond )
                    hydrogens.append  ( atom )
        for hydrogen in hydrogens:
            self.atoms.remove ( hydrogen )

    def Summary ( self, log = logFile ):
        """Write a summary."""
        if LogFileActive ( log ):
            # . Find some atom statistics.
            charge    = 0
            naromatic = 0
            nexplicit = 0
            nimplicit = 0
            nring     = 0
            nunknown  = 0
            for atom in self.atoms:
                charge    += atom.formalCharge
                nimplicit += atom.implicithydrogens
                if atom.isAromatic: naromatic += 1
                if atom.QRING    : nring     += 1
                if   atom.atomicNumber <= 0: nunknown  += 1
                elif atom.atomicNumber == 1: nexplicit += 1
            natoms = len ( self.atoms ) + nimplicit
            nbonds = len ( self.bonds ) + nimplicit
            nheavy = natoms - nexplicit - nimplicit - nunknown
            # . Ring statistics.
            nrings = 0
            for ringset in self.ringsets: nrings += len ( ringset )
            # . Summary.
            summary = log.GetSummary ( )
            summary.Start ( "SMILES Connectivity Summary" )
            summary.Entry ( "Number of Atoms",     "{:d}".format ( natoms                ) )
            summary.Entry ( "Number of Bonds",     "{:d}".format ( nbonds                ) )
            summary.Entry ( "Number of Ring Sets", "{:d}".format ( len ( self.ringsets ) ) )
            summary.Entry ( "Number of Rings",     "{:d}".format ( nrings                ) )
            summary.Entry ( "Charge",              "{:d}".format ( charge                ) )
            summary.Entry ( "Heavy Atoms",         "{:d}".format ( nheavy                ) )
            summary.Entry ( "Explicit Hydrogens",  "{:d}".format ( nexplicit             ) )
            summary.Entry ( "Implicit Hydrogens",  "{:d}".format ( nimplicit             ) )
            summary.Entry ( "Unknown Atoms",       "{:d}".format ( nunknown              ) )
            summary.Entry ( "Aromatic Atoms",      "{:d}".format ( naromatic             ) )
            summary.Entry ( "Ring Atoms",          "{:d}".format ( nring                 ) )
            summary.Stop ( )

    def ToSystem ( self, **keywordArguments ):
        """Return a system.

        Implicit hydrogens are added unless |QOMITIMPLICITHYDROGENS| is True.
        The hydrogens are placed after their host atom.
        """
        # . Options.
        QOMIT = keywordArguments.get ( "QOMITIMPLICITHYDROGENS", False )
        # . Atoms.
        if QOMIT:
            atoms = self.atoms
        else:
            atoms = []
            index = 0
            for atom in self.atoms:
                atom.index = index
                atoms.append ( atom )
                for i in range ( atom.implicithydrogens ):
                    atoms.append ( 1 )
                index += atom.implicithydrogens + 1
        # . Bonds.
        bonds = []
        for bond in self.bonds:
            bonds.append ( ( bond.atom1.index, bond.atom2.index, bond.type ) )
        if not QOMIT:
            for atom in self.atoms:
                for i in range ( atom.implicithydrogens ):
                    bonds.append ( ( atom.index, atom.index + i + 1, SingleBond ( ) ) )
            self.IndexAtoms ( )
        # . Make the system.
        system = System.FromAtoms ( atoms, bonds = bonds )
        # . Set the aromaticity and formalCharge of the new atoms.
        for atom in self.atoms:
            index = atom.index
            system.atoms[index].formalCharge = atom.formalCharge
            system.atoms[index].isAromatic   = atom.isAromatic
        # . Finish up.
        return system

    def VerifyAromaticity ( self ):
        """Verify the existing aromaticity of the connectivity."""
        # . Find aromatic atoms.
        aromatics = set ( )
        for atom in self.atoms:
            if atom.isAromatic: aromatics.add ( atom )
        # . Continue only if there are aromatics.
        if len ( aromatics ) > 0:
            # . Do each ringset separately.
            for ringset in self.ringsets:
                # . Find aromatic rings.
                aromaticrings = []
                for ring in ringset:
                    isAromatic = True
                    for atom in ring:
                        if not atom.isAromatic:
                            isAromatic = False
                            break
                    if isAromatic: aromaticrings.append ( ring )
                # . Check each ring separately - eventually combinations of rings should be done.
                for ring in aromaticrings:
                    # . Determine the number of electrons donated to the ring.
                    nelectrons = 0
                    for atom in ring: nelectrons += atom.aromaticelectrons
                    # . Does this ring obey Huckel's 4N+2 rule?
                    n = ( nelectrons - 2 ) // 4
                    if ( n > 0 ) and ( 4 * n + 2 == nelectrons ):
                        for atom in ring: aromatics.discard ( atom )
            # . There should not be any unused aromatics.
            if len ( aromatics ) > 0: raise SMILESConnectivityError ( "Some aromatic atoms have been found not to be aromatic.", atoms = aromatics )

#===============================================================================
# . Error classes.
#===============================================================================
class SMILESConnectivityError ( Exception ):

    def __init__ ( self, *arguments, **keywordArguments ):
        """Constructor."""
        super ( SMILESConnectivityError, self ).__init__ ( *arguments )
        if len ( arguments ) == 0: self.args = ( "SMILES connectivity error.", )
        self.atom  = keywordArguments.get ( "atom" , None )
        self.bond  = keywordArguments.get ( "bond" , None )
        self.atoms = keywordArguments.get ( "atoms", None )

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
