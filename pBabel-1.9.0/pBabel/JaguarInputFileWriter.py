#-------------------------------------------------------------------------------
# . File      : JaguarInputFileWriter.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for writing Jaguar input files."""

from pCore        import Coordinates3, TextFileWriter, TextFileWriterError
from ExportImport import _Exporter
from pMolecule    import PeriodicTable, System

#===================================================================================================================================
# . Jaguar input file writer class.
#===================================================================================================================================
class JaguarInputFileWriter ( TextFileWriter ):
    """JaguarInputFileWriter is the class for Jaguar input files that are to be written."""

    def WriteFrame ( self, system, label = None, xyz = None ):
        """Write a single frame."""
        # . Check the data.
        if xyz == None: xyz = system.coordinates3
        if ( xyz is None ) or ( not isinstance ( xyz, Coordinates3 ) ) or ( xyz.rows != len ( system.atoms ) ): raise TextFileWriterError ( "Invalid or missing data to write to Jaguar input file." )
        # . Check the label.
        if label is None: label = system.label
        if label is None: label = "Jaguar input file."
        # . Check the electronic state.
        if system.electronicState is None:
            charge = 0
            multip = 1
        else:
            charge = system.electronicState.charge
            multip = system.electronicState.multiplicity
        # . Header.
        self.file.write ( label + "\n\n" )
        # . General section.
        self.file.write ( "&gen\n" )
        self.file.write ( "molchg={:d}\n".format ( charge ) )
        self.file.write ( "multip={:d}\n".format ( multip ) )
        self.file.write ( "&\n" )
        # . Coordinates.
        self.file.write ( "&zmat\n" )
        numbers = system.atoms.GetItemAttributes ( "atomicNumber" )
        for ( i, n ) in enumerate ( numbers ):
            symbol = PeriodicTable.Symbol ( n, index = i+1 )
            self.file.write ( "{:<5s}{:25.15f}{:25.15f}{:25.15f}\n".format ( symbol, xyz[i,0], xyz[i,1], xyz[i,2] ) )
        self.file.write ( "&\n" )

    def WriteSingleSystem ( self, system, label = None, xyz = None ):
        """Write a complete file."""
        self.Open  ( )
        self.WriteFrame ( system, label = label, xyz = xyz )
        self.Close ( )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def JaguarInputFile_FromSystem ( filename, system, label = None, xyz = None ):
    """Helper function that writes a system to a Jaguar input file."""
    outfile = JaguarInputFileWriter ( filename )
    outfile.WriteSingleSystem ( system, label = label, xyz = xyz )

# . Exporter definitions.
_Exporter.AddHandler ( { System : JaguarInputFile_FromSystem } , [ "jagin", "jin", "JAGIN", "JIN" ], "Jaguar Input" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
