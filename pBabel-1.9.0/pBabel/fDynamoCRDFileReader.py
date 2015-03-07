#-------------------------------------------------------------------------------
# . File      : fDynamoCRDFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Read data from an fDynamo CRD file."""

from pCore        import Clone, Coordinates3, logFile, LogFileActive, TextFileReader
from ExportImport import _Importer
from pMolecule    import Sequence, System

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class fDynamoCRDFileReader ( TextFileReader ):
    """Class for reading fDynamo CRD files."""

    def GetLine ( self, QWARNING = True ):
        """Get a non-empty line removed of comments."""
        try:
            QFOUND = False
            while not QFOUND:
                line  = next ( self.file )
                index = line.find ( "!" )
                if index >= 0: line = line[:index]
                line  = line.strip ( )
                self.nlines += 1
                QFOUND = ( len ( line ) > 0 )
            return line
        except:
            if QWARNING: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Read data.
            self.atoms      = []
            self.residues   = []
            self.subsystems = []
	    # . Open the file.
	    self.Open ( )
	    # . Parse all the lines.
            try:
                ( natoms, nresidues, nsubsystems ) = self.__GetCounters ( )
                # . Subsystems.
                for isub in range ( nsubsystems ):
                    nres = self.__SubsystemHeader ( )
                    for ires in range ( nres ):
                        natm = self.__ResidueHeader ( )
                        for iatm in range ( natm ):
                            self.__AtomRecord ( )
                # . Do a final logic check.
                if ( len ( self.atoms ) != natoms ) or ( len ( self.residues ) != nresidues ) or ( len ( self.subsystems ) != nsubsystems ):
                    self.Warning ( "Counter mismatch after parsing: {:d}/{:d}, {:d}/{:d}, {:d}/{:d}.".format ( len ( self.atoms      ) , natoms      , \
                                                                                                               len ( self.residues   ) , nresidues   , \
                                                                                                               len ( self.subsystems ) , nsubsystems ), True )
            except EOFError:
                pass
	    # . Close the file.
            self.WarningStop ( )
	    self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.QPARSED = True

    def Summary ( self, log = logFile ):
        """Print a summary of the stored data."""
        if self.QPARSED and LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "CRD File Summary" )
            summary.Entry ( "Number of Atoms",      "{:d}".format ( len ( self.atoms      ) ) )
            summary.Entry ( "Number of Residues",   "{:d}".format ( len ( self.residues   ) ) )
            summary.Entry ( "Number of Subsystems", "{:d}".format ( len ( self.subsystems ) ) )
            summary.Stop ( )

    def ToCoordinates3 ( self ):
        """Return the coordinates."""
        if self.QPARSED:
            coordinates3 = Coordinates3.WithExtent ( len ( self.atoms ) )
            for ( i, ( aname, atomicNumber, x, y, z ) ) in enumerate ( self.atoms ):
                coordinates3[i,0] = x
                coordinates3[i,1] = y
                coordinates3[i,2] = z
            return coordinates3
        else:
            return None

    def ToSequence ( self ):
        """Return a sequence."""
        sequence = None
        if self.QPARSED:
            majorSeparator = Sequence.defaultAttributes["labelSeparator"]
            minorSeparator = Sequence.defaultAttributes["fieldSeparator"]
            atomPaths = []
            for ( sname, sresidues ) in self.subsystems:
                for ires in range ( sresidues ):
                    ( rname, ratoms ) = residues.pop ( 0 )
                    pathHead = sname + majorSeparator + rname + minorSeparator + "{:d}".format ( ires + 1 ) + majorSeparator
                    for iatom in range ( ratoms ):
                        atom    = atoms.pop ( 0 )
                        atomPaths.append ( pathHead + atom )
            sequence = Sequence.FromAtomPaths ( atomPaths )
        return sequence

    def ToSystem ( self ):
        """Return a system from the CRD data but without sequence information."""
        if self.QPARSED:
            atoms = []
            for ( aname, atomicNumber, x, y, z ) in self.atoms:
                atoms.append ( atomicNumber )
            system = System.FromAtoms ( atoms )
            system.coordinates3 = self.ToCoordinates3 ( )
            return system
        else:
            return None

    # . Private methods.
    def __AtomRecord ( self ):
        """Extract atom data from the next record."""
        tokens = self.GetTokens ( converters = [ int, None, int, float, float, float ] )
        if ( len ( tokens ) >= 6 ):
            self.atoms.append ( ( tokens[1], tokens[2], tokens[3], tokens[4], tokens[5] ) )
        else:
            self.Warning ( "Invalid atom data line." )
            self.atoms.append ( ( "", -1, 0.0, 0.0, 0.0 ) )

    def __GetCounters ( self ):
        """Get the counters from the first line of the file."""
        tokens = self.GetTokens ( converters = [ int, int, int ] )
        if ( len ( tokens ) >= 3 ) and ( tokens[0] >= 0 ) and ( tokens[1] >= 0 ) and ( tokens[2] >= 0 ):
            return ( tokens[0], tokens[1], tokens[2] ) # . natoms, nresidues, nsubsystems.
        else:
            self.Warning ( "Invalid natoms, nresidues or nsubsystems counter.", True )
            return ( 0, 0, 0 )

    def __ResidueHeader ( self ):
        """Parse a residue header."""
        # . Initialization.
        name = ""
        natm = ""
        # . Residue line.
        tokens = self.GetTokens ( converters = [ None, int, None ] )
        if ( len ( tokens ) >= 3 ) and ( tokens[0].upper ( ) == "RESIDUE" ): name = tokens[2]
        else: self.Warning ( "Invalid residue header line.", True )
        # . Atom counter line.
        tokens = self.GetTokens ( converters = [ int ] )
        if ( len ( tokens ) >= 1 ) and ( tokens[0] >= 0 ): natm = tokens[0]
        else: self.Warning ( "Invalid residue atom counter.", True )
        self.residues.append ( ( name, natm ) )
        return natm

    def __SubsystemHeader ( self ):
        """Parse a subsystem header."""
        # . Initialization.
        name = ""
        nres = ""
        # . Subsystem line.
        tokens = self.GetTokens ( converters = [ None, int, None ] )
        if ( len ( tokens ) >= 3 ) and ( tokens[0].upper ( ) == "SUBSYSTEM" ): name = tokens[2]
        else: self.Warning ( "Invalid subsystem header line.", True )
        # . Residue counter line.
        tokens = self.GetTokens ( converters = [ int ] )
        if ( len ( tokens ) >= 1 ) and ( tokens[0] >= 0 ): nres = tokens[0]
        else: self.Warning ( "Invalid subsystem residue counter.", True )
        self.subsystems.append ( ( name, nres ) )
        return nres

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def fDynamoCRDFile_ToCoordinates3 ( fileName, log = logFile ):
    """Helper function that returns a set of coordinates from an fDynamo CRD file."""
#    atomTranslation      = kwargs.get ( "atomTranslation",       {} ) # . For future use (maybe).
#    residueTranslation   = kwargs.get ( "residueTranslation",    {} )
#    subsystemTranslation = kwargs.get ( "subsystemTranslation ", {} )
    infile = fDynamoCRDFileReader ( fileName )
    infile.Parse   ( log = log )
    infile.Summary ( log = log )
    return infile.ToCoordinates3 ( )

def fDynamoCRDFile_ToSequence ( fileName, log = logFile ):
    """Helper function that returns a sequence from an fDynamo CRD file."""
#    atomTranslation      = kwargs.get ( "atomTranslation",       {} ) # . For future use (maybe).
#    residueTranslation   = kwargs.get ( "residueTranslation",    {} )
#    subsystemTranslation = kwargs.get ( "subsystemTranslation ", {} )
    infile = fDynamoCRDFileReader ( fileName )
    infile.Parse   ( log = log )
    infile.Summary ( log = log )
    return infile.ToSequence ( )

def fDynamoCRDFile_ToSystem ( fileName, log = logFile ):
    """Helper function that returns a system from an fDynamo CRD file."""
#    atomTranslation      = kwargs.get ( "atomTranslation",       {} ) # . For future use (maybe).
#    residueTranslation   = kwargs.get ( "residueTranslation",    {} )
#    subsystemTranslation = kwargs.get ( "subsystemTranslation ", {} )
    infile = fDynamoCRDFileReader ( fileName )
    infile.Parse   ( log = log )
    infile.Summary ( log = log )
    return infile.ToSystem ( )

# . Importer definitions.
_Importer.AddHandler ( { Coordinates3 : fDynamoCRDFile_ToCoordinates3 ,
                         System       : fDynamoCRDFile_ToSystem       } , [ "fcrd", "FCRD" ], "Fortran Dynamo Coordinates", defaultFunction = fDynamoCRDFile_ToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

