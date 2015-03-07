#-------------------------------------------------------------------------------
# . File      : MOL2FileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
#===================================================================================================================================
# . Classes and functions to read MOL2 files.
#===================================================================================================================================

from pCore        import Coordinates3, logFile, LogFileActive, Real1DArray, TextFileReader
from ExportImport import _Importer
from pMolecule    import CrystalClass_FromSpaceGroupNumber, PeriodicTable, System, SymmetryParameters, \
                         AromaticSingleBond, DoubleBond, SingleBond, TripleBond, UndefinedBond

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Bond type definitions.
_MOL2BONDTYPES = { "1"  : SingleBond         ( ) ,
                   "2"  : DoubleBond         ( ) ,
                   "3"  : TripleBond         ( ) ,
                   "am" : SingleBond         ( ) ,
                   "ar" : AromaticSingleBond ( ) ,
                   "du" : UndefinedBond      ( ) ,
                   "un" : UndefinedBond      ( ) ,
                   "nc" : None                   }

# . Generic (non-element) atom types.
_GENERICATOMTYPES = ( "Any", "Du", "Hal", "Het", "Hev", "LP" )

# . Section header.
_RTISTRING  = "@<TRIPOS>"

#===================================================================================================================================
# . MOL2 file reader class.
#===================================================================================================================================
class MOL2FileReader ( TextFileReader ):
    """MOL2FileReader is the class for MOL2 files that are to be read."""

    def AtomicNumberFromAtomType ( self, atomtype ):
        """Get an atomic number from atom type."""
        atomicNumber = -1
        token        = atomtype.split ( ".", 1 )[0]
        if token not in _GENERICATOMTYPES: atomicNumber = PeriodicTable.AtomicNumber ( token )
        return atomicNumber

    def GetLine ( self, QWARNING = False ):
        """Get a line."""
        try:
            while True:
                line = next ( self.file ).strip ( )
                self.nlines += 1
                if ( len ( line ) > 0 ) and ( not line.startswith ( "#" ) ): break
            return line
        except:
            if QWARNING: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            # . Parse the data.
            try:
                # . Parse all entries.
                QCONTINUE = True
                while QCONTINUE:
                    line = self.GetLine ( )
#                    print line
                    if   line.startswith ( _RTISTRING + "ATOM"     ): self.ParseAtomSection     ( )
                    elif line.startswith ( _RTISTRING + "BOND"     ): self.ParseBondSection     ( )
                    elif line.startswith ( _RTISTRING + "CRYSIN"   ): self.ParseCrysinSection   ( )
                    elif line.startswith ( _RTISTRING + "MOLECULE" ): self.ParseMoleculeSection ( )
#                    elif line.startswith ( _RTISTRING ):              self.ParseUnknownSection  ( )
            except EOFError:
                pass
	    # . Close the file.
            self.WarningStop ( )
	    self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.QPARSED = True

    def ParseAtomSection ( self ):
        """Parse the ATOM section."""
        if hasattr ( self, "natoms" ):
            atomicCharges = Real1DArray.WithExtent ( self.natoms )
            atomicCharges.Set ( 0.0 )
            atomicNumbers = []
            atomNames     = []
            xyz           = Coordinates3.WithExtent ( self.natoms )
            xyz.Set ( 0.0 )
            for i in range ( self.natoms ):
                items = self.GetTokens ( converters = [ int, None, float, float, float, None, None, None, float ] )
                if len ( items ) < 6:
                    self.Warning ( "Invalid ATOM line.", True )
                else:
                    atomicNumber = PeriodicTable.AtomicNumber ( items[1] )
                    if atomicNumber <= 0: atomicNumber = self.AtomicNumberFromAtomType ( items[5] )
                    atomicNumbers.append ( atomicNumber )
                    atomNames.append ( items[1] )
                    xyz[i,0] = items[2]
                    xyz[i,1] = items[3]
                    xyz[i,2] = items[4]
                    if len ( items ) >= 9: atomicCharges[i] = items[8]
            self.atomicCharges = atomicCharges
            self.atomicNumbers = atomicNumbers
            self.atomNames     = atomNames
            self.xyz           = xyz
        else:
            self.Warning ( "Unknown number of atoms in molecule.", True )

    def ParseBondSection ( self ):
        """Parse the BOND section."""
        if hasattr ( self, "nbonds" ):
            bonds  = []
            natoms = getattr ( self, "natoms", -1 )
            for i in range ( self.nbonds ):
                items = self.GetTokens ( converters = [ int, int, int, None ] )
                if len ( items ) < 4:
                    self.Warning ( "Invalid BOND line.", True )
                else:
                    atom1    = items[1] - 1
                    atom2    = items[2] - 1
                    bondtype = _MOL2BONDTYPES.get ( items[3], UndefinedBond ( ) )
                    if ( atom1 < 0 ) or ( atom1 >= self.natoms ) or ( atom2 < 0 ) or ( atom2 >= self.natoms ): self.Warning ( "Bond atom indices out of range: {:d}, {:d}.".format ( atom1, atom2 ), True )
                    if bondtype is not None: bonds.append ( ( max ( atom1, atom2 ), min ( atom1, atom2 ), bondtype ) )
            bonds.sort ( )
            self.bonds = bonds
        else:
            self.Warning ( "Unknown number of bonds in molecule.", True )

    def ParseCrysinSection ( self ):
        """Parse the CRYSIN section."""
        # . Items are a, b, c, alpha, beta, gamma, space group number and crystal setting.
        items  = self.GetTokens ( converters = 6 * [ float ] + 2 * [ int ] )
        try:
            setting    = items.pop ( -1 )
            spaceGroup = items.pop ( -1 )
            self.crystalClass       = CrystalClass_FromSpaceGroupNumber ( spaceGroup )
            self.symmetryParameters = SymmetryParameters ( )
            self.symmetryParameters.SetCrystalParameters ( *items )
        except:
            self.Warning ( "Invalid CRYSIN line.", False )

    def ParseMoleculeSection ( self ):
        """Parse the MOLECULE section."""
        # . Just parse the first two lines for the moment.
        self.label = self.GetLine ( )                             # . Molecule name.
        items      = self.GetTokens ( converters = [ int, int ] ) # . Number of atoms and bonds.
        if len ( items ) > 0: self.natoms = items[0]
        if len ( items ) > 1: self.nbonds = items[1]

    def ParseUnknownSection ( self ):
        """Parse an unknown section."""
        while True:
            line = self.GetLine ( )
            if line.startswith ( _RTISTRING ): break

    def ToAtomNames ( self ):
        """Return atom names."""
        atomNames = None
        if self.QPARSED: return getattr ( self, "atomNames", None )
        return atomNames

    def ToCharges ( self ):
        """Return charges."""
        charges = None
        if self.QPARSED: return getattr ( self, "atomicCharges", None )
        return charges

    def ToCoordinates3 ( self ):
        """Return a coordinates3 object."""
        xyz = None
        if self.QPARSED: return getattr ( self, "xyz", None )
        return xyz

    def ToSystem ( self ):
        """Return a system."""
        system = None
        if self.QPARSED:
            try:
                system              = System.FromAtoms ( self.atomicNumbers, bonds = getattr ( self, "bonds", None ) )
                system.label        = self.label
                system.coordinates3 = self.xyz
                if hasattr ( self, "crystalClass" ) and hasattr ( self, "symmetryParameters" ):
                    definitions = self.crystalClass.GetUniqueSymmetryParameters ( self.symmetryParameters )
                    definitions["crystalClass"] = self.crystalClass
                    system.DefineSymmetry ( **definitions )
            except:
                pass
        return system

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def MOL2File_ToAtomNames ( filename, log = logFile ):
    """Helper function that reads atom names from a MOL2 file."""
    infile = MOL2FileReader ( filename )
    infile.Parse ( log = log )
    charges = infile.ToAtomNames ( )
    return charges

def MOL2File_ToCharges ( filename, log = logFile ):
    """Helper function that reads charges from a MOL2 file."""
    infile = MOL2FileReader ( filename )
    infile.Parse ( log = log )
    charges = infile.ToCharges ( )
    return charges

def MOL2File_ToCoordinates3 ( filename, log = logFile ):
    """Helper function that reads coordinates from a MOL2 file."""
    infile = MOL2FileReader ( filename )
    infile.Parse ( log = log )
    coordinates3 = infile.ToCoordinates3 ( )
    return coordinates3

def MOL2File_ToSystem ( filename, log = logFile ):
    """Helper function that reads a system from a MOL2 file."""
    infile = MOL2FileReader ( filename )
    infile.Parse ( log = log )
    system = infile.ToSystem ( )
    return system

# . Importer definitions.
_Importer.AddHandler ( { Coordinates3 : MOL2File_ToCoordinates3 ,
                         System       : MOL2File_ToSystem       } , [ "mol2", "MOL2" ], "Tripos MOL2", defaultFunction = MOL2File_ToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
