#-------------------------------------------------------------------------------
# . File      : MOLFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
#===================================================================================================================================
# . Classes and functions to read MOL and SDF files.
#===================================================================================================================================

from pCore        import Coordinates3, logFile, LogFileActive, TextFileReader
from ExportImport import _Importer
from pMolecule    import Atom, PeriodicTable, System, \
                         AromaticSingleBond, DoubleBond, SingleBond, TripleBond, UndefinedBond

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Bond type definitions.
_MOLBONDTYPES = { 1 : SingleBond ( ), \
                  2 : DoubleBond ( ), \
                  3 : TripleBond ( ), \
                  4 : AromaticSingleBond ( ) }

# . Charge definitions.
_CHARGECODES = { 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3 }

#===================================================================================================================================
# . MOL file reader class.
#===================================================================================================================================
class MOLFileReader ( TextFileReader ):
    """MOLFileReader is the class for MOL or SDF files that are to be read."""

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            self.molrecords = []
            # . Open the file.
            self.Open ( )
            # . Parse the data.
            try:
                # . Parse all entries.
                QCONTINUE = True
                while QCONTINUE:
                    # . Header block.
                    # . Molecule name.
                    label = self.GetLine ( QWARNING = False )
                    # . Skip lines 1 and 2 which are not parsed.
                    self.GetLine ( )
                    self.GetLine ( )
                    # . Connection table.
                    # . Counts line - only atom and bond numbers are parsed.
                    ( natoms, nbonds ) = self.GetFixedFormatTokens ( ( 0, 3, int, 0 ), ( 3, 6, int, 0 ) )
                    # . Initialize atom and coordinates data structures.
                    atomicNumbers = []
                    bonds         = []
                    charges       = []
                    coordinates3  = Coordinates3.WithExtent ( natoms )
                    mchg          = {}
                    miso          = {}
                    mrad          = {}
                    coordinates3.Set ( 0.0 )
                    # . Read the atom lines.
                    for n in range ( natoms ):
                        ( x, y, z, atomicNumber, charge ) = self.GetFixedFormatTokens ( ( 0, 10, float, 0.0 ), ( 10, 20, float, 0.0 ), ( 20, 30, float, 0.0 ), ( 31, 34, PeriodicTable.AtomicNumber, -1 ), ( 36, 39, int, 0 ) )
                        atomicNumbers.append ( atomicNumber )
                        charges.append ( _CHARGECODES.get ( charge, 0 ) )
                        coordinates3[n,0] = x
                        coordinates3[n,1] = y
                        coordinates3[n,2] = z
                    # . Read the bond lines.
                    for n in range ( nbonds ):
                        ( atom1, atom2, code ) = self.GetFixedFormatTokens ( ( 0, 3, int, 0 ), ( 3, 6, int, 0 ), ( 6, 9, int, 0 ) )
                        if ( atom1 <= 0 ) or ( atom1 > natoms ) or ( atom2 <= 0 ) or ( atom2 > natoms ): self.Warning ( "Bond atom indices out of range: {:d}, {:d}.".format ( atom1, atom2 ), True )
                        bonds.append ( ( atom1 - 1, atom2 - 1, _MOLBONDTYPES.get ( code, UndefinedBond ( ) ) ) )
                    # . Properties lines.
                    while True:
                        line = self.GetLine ( )
                        if   line.startswith ( "M  CHG" ): self.ParsePropertiesLine ( line, mchg )
                        elif line.startswith ( "M  ISO" ): self.ParsePropertiesLine ( line, miso )
                        elif line.startswith ( "M  RAD" ): self.ParsePropertiesLine ( line, mrad )
                        elif line.startswith ( "M  END" ): break
                    # . Store the data.
                    self.molrecords.append ( ( label, atomicNumbers, bonds, charges, coordinates3, mchg, miso, mrad ) )
                    # . Check for an SDF terminator.
                    QCONTINUE = False
                    while True:
                        line = self.GetLine ( QWARNING = False )
                        if line == "$$$$":
                            QCONTINUE = True
                            break
            except EOFError:
                pass
	    # . Close the file.
            self.WarningStop ( )
	    self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.QPARSED = True

    def ParsePropertiesLine ( self, line, mdict ):
        """Parse a CHG, ISO or RAD properties line."""
        try:    n = int ( line[8:9] )
        except: n = 0
        for i in range ( n ):
            o = i * 8 + 10
            try:
                a = int ( line[o  :o+3] ) - 1
                v = int ( line[o+4:o+7] )
                if a >= 0: mdict[a] = v
            except:
                pass

    def ToAtoms ( self, atomicNumbers, charges, mchg, miso, mrad ):
        """Set some atom properties."""
        atoms = []
        if self.QPARSED:
            for ( i, atomicNumber ) in enumerate ( atomicNumbers ):
                atoms.append ( Atom ( atomicNumber = atomicNumber ) )
            if len ( mchg ) > 0:
                for ( i, v ) in mchg.iteritems ( ):
                    atoms[i].formalCharge = v
        return atoms

    def ToCoordinates3 ( self, index = 0 ):
        """Return a coordinates3 object."""
        if self.QPARSED:
            if index in range ( len ( self.molrecords ) ):
                data = self.molrecords[index]
                return data[4]
            else:
                raise IndexError ( "MOL record index {:d} not in range [0,{:d}].".format ( index, len ( self.molrecords ) - 1 ) )
        else:
            return None

    def ToSystem ( self, index = 0 ):
        """Return a system."""
        if self.QPARSED:
            if index in range ( len ( self.molrecords ) ):
                # . Get data.
                ( label, atomicNumbers, bonds, charges, coordinates3, mchg, miso, mrad ) = self.molrecords[index]
                # . Make atoms.
                atoms = self.ToAtoms ( atomicNumbers, charges, mchg, miso, mrad )
                # . Make system.
                system              = System.FromAtoms ( atoms, bonds = bonds )
                system.label        = label
                system.coordinates3 = coordinates3
                return system
            else:
                raise IndexError ( "MOL record index {:d} not in range [0,{:d}].".format ( index, len ( self.molrecords ) - 1 ) )
        else:
            return None

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def MOLFile_ToCoordinates3 ( filename, index = 0, log = logFile ):
    """Helper function that reads coordinates from a MOL file."""
    infile = MOLFileReader ( filename )
    infile.Parse ( log = log )
    coordinates3 = infile.ToCoordinates3 ( index = index )
    return coordinates3

def MOLFile_ToSystem ( filename, index = 0, log = logFile ):
    """Helper function that reads a system from a MOL file."""
    infile = MOLFileReader ( filename )
    infile.Parse ( log = log )
    system = infile.ToSystem ( index = index )
    return system

def SDFFile_ToSystems ( filename, log = logFile ):
    """Helper function that reads all the systems from an SDF file."""
    infile = MOLFileReader ( filename )
    infile.Parse ( log = log )
    systems = []
    for i in range ( len ( infile.molrecords ) ):
        systems.append ( infile.ToSystem ( index = i ) )
    return systems

# . Importer definitions.
_Importer.AddHandler ( { Coordinates3 : MOLFile_ToCoordinates3 ,
                         System       : MOLFile_ToSystem       } , [ "mol", "MOL" ], "MDL MOL", defaultFunction = MOLFile_ToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
