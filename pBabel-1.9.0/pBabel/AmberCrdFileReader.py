#-------------------------------------------------------------------------------
# . File      : AmberCrdFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for reading Amber crd files."""

from pCore        import Coordinates3, logFile, LogFileActive, TextFileReader
from ExportImport import _Importer
from pMolecule    import SymmetryParameters

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . The width of the natoms field (second line).
_NATOMSWIDTH = 6

# . The number and the widths of XYZ elements on a line.
_XYZNUMBER =  6
_XYZWIDTH  = 12

# . The symmetry width.
_SYMMETRYNUMBER =  3
_SYMMETRYWIDTH  = 12

#===================================================================================================================================
# . AmberCrd file reader class.
#===================================================================================================================================
class AmberCrdFileReader ( TextFileReader ):
    """AmberCrdFileReader is the class for Amber crd files that are to be read."""

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Keyword line.
                self.title = self.GetLine ( )
	        # . Number of atoms.
                items    = self.GetTokens ( converters = ( int, ) )
	        natoms   = items[0]
	        # . The coordinate data.
                items    = self.GetFixedFormatArray ( 3 * natoms, _XYZNUMBER, _XYZWIDTH, converter = float, default = 0.0 )
                self.xyz = Coordinates3.WithExtent ( natoms )
                for n in range ( natoms ):
                    for i in range ( 3 ): self.xyz[n,i] = items[3*n+i]
                # . Symmetry data - optional.
                items = self.GetFixedFormatArray ( 3, _SYMMETRYNUMBER, _SYMMETRYWIDTH, converter = float, default = 0.0, QWARNING = False )
                self.symmetryParameters = SymmetryParameters ( )
                self.symmetryParameters.SetCrystalParameters ( items[0], items[1], items[2], 90.0, 90.0, 90.0 )
            except EOFError:
                pass
	    # . Close the file.
            self.WarningStop ( )
	    self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.QPARSED = True

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ) and self.QPARSED:
            summary = log.GetSummary ( )
            summary.Start ( "Amber Crd File Summary" )
            summary.Entry ( "Number of Atoms", "{:d}".format ( self.xyz.rows ) )
            summary.Stop ( )

    def ToCoordinates3 ( self ):
        """Return the coordinates."""
        if self.QPARSED: return self.xyz
        else:            return None

    def ToSymmetryParameters ( self ):
        """Return the symmetry parameters."""
        if self.QPARSED: return self.symmetryParameters
        else:            return None

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def AmberCrdFile_ToCoordinates3 ( filename, log = logFile ):
    """Helper function that reads coordinates from an Amber crd file."""
    infile = AmberCrdFileReader ( filename )
    infile.Parse   ( log = log )
    infile.Summary ( log = log )
    return infile.ToCoordinates3 ( )

def AmberCrdFile_ToSymmetryParameters ( filename, log = logFile ):
    """Helper function that reads symmetryparamters from an Amber crd file."""
    infile = AmberCrdFileReader ( filename )
    infile.Parse   ( log = log )
    infile.Summary ( log = log )
    return infile.ToSymmetryParameters ( )

# . Importer definitions.
_Importer.AddHandler ( { Coordinates3 : AmberCrdFile_ToCoordinates3 } , [ "crd", "CRD" ], "Amber Coordinates", defaultFunction = AmberCrdFile_ToCoordinates3 )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
