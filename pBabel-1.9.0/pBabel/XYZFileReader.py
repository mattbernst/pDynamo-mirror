#-------------------------------------------------------------------------------
# . File      : XYZFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for reading XYZ files."""

from pCore                    import Coordinates3, logFile, LogFileActive, TextFileReader
from ExportImport             import _Importer
from pMolecule                import PeriodicTable, System
from SystemGeometryTrajectory import SystemGeometryTrajectory

#===================================================================================================================================
# . XYZ file reader class.
#===================================================================================================================================
class XYZFileReader ( TextFileReader ):
    """XYZFileReader is the class for XYZ files that are to be read."""

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Number of atoms.
                items  = self.GetTokens ( converters = [ int ] )
                natoms = items[0]
                # . Title line.
                self.title = self.GetLine ( )
                # . XYZ lines.
                self.atomicNumbers = []
                self.coordinates3  = Coordinates3.WithExtent ( natoms )
                for i in range ( natoms ):
                    items = self.GetTokens ( converters = [ PeriodicTable.AtomicNumber, float, float, float ] )
                    self.atomicNumbers.append ( items[0] )
                    self.coordinates3[i,0] = items[1]
                    self.coordinates3[i,1] = items[2]
                    self.coordinates3[i,2] = items[3]
            except EOFError:
                pass
	    # . Close the file.
            self.WarningStop ( )
	    self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.QPARSED = True

    def ToCoordinates3 ( self ):
        """Return a coordinates3 object."""
        if self.QPARSED and hasattr ( self, "coordinates3" ): return self.coordinates3
        else:                                                 return None

    def ToSystem ( self ):
        """Return a system."""
        system = None
        if self.QPARSED and hasattr ( self, "atomicNumbers" ):
            system              = System.FromAtoms ( self.atomicNumbers )
            system.label        = self.title
            system.coordinates3 = self.ToCoordinates3 ( )
        return system

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def XYZFile_ToCoordinates3 ( filename, log = logFile ):
    """Helper function that reads the coordinates from a XYZ file."""
    infile = XYZFileReader ( filename )
    infile.Parse ( )
    coordinates3 = infile.ToCoordinates3 ( )
    return coordinates3

def XYZFile_ToSystem ( filename, log = logFile ):
    """Helper function that reads a system from a XYZ file."""
    infile = XYZFileReader ( filename )
    infile.Parse ( )
    system = infile.ToSystem ( )
    return system

def XYZFiles_ToSystemGeometryTrajectory ( inPaths, outPath, system ):
    """Convert XYZ files to a SystemGeometryTrajectory.

    Files are added in the order that they are supplied.
    """
    # . Define the output trajectory.
    outTrajectory = SystemGeometryTrajectory  ( outPath, system, mode = "w" )
    # . Loop over paths.
    for inPath in inPaths:
        system.coordinates3 = XYZFile_ToCoordinates3 ( inPath )
        outTrajectory.WriteOwnerData ( )
    outTrajectory.Close ( )

# . Importer definitions.
_Importer.AddHandler ( { Coordinates3 : XYZFile_ToCoordinates3 ,
                         System       : XYZFile_ToSystem       } , [ "xyz" ], "XYZ", defaultFunction = XYZFile_ToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
