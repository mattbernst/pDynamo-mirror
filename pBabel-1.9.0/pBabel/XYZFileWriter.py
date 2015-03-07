#-------------------------------------------------------------------------------
# . File      : XYZFileWriter.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for writing XYZ files."""

import os

from pCore                    import Coordinates3, TextFileWriter, TextFileWriterError
from ExportImport             import _Exporter
from pMolecule                import PeriodicTable, System
from SystemGeometryTrajectory import SystemGeometryTrajectory

# . Maximum label length (needed to avoid errors in some programs - e.g. VMD).
_MaximumLabelLength = 80

#===================================================================================================================================
# . XYZ file writer class.
#===================================================================================================================================
class XYZFileWriter ( TextFileWriter ):
    """XYZFileWriter is the class for XYZ files that are to be written."""

    def WriteFrame ( self, system, label = None, xyz = None ):
        """Write a single frame."""
        # . Check the data.
        if xyz == None: xyz = system.coordinates3
        if ( xyz is None ) or ( not isinstance ( xyz, Coordinates3 ) ) or ( xyz.rows != len ( system.atoms ) ): raise TextFileWriterError ( "Invalid or missing data to write to XYZ file." )
        # . Check the label.
        if label is None: label = system.label
        # . Write the frame.
        self.file.write ( "{:6d}\n".format ( len ( system.atoms ) ) )
        if label is None:
            self.file.write ( "\n" )
        else:
            if len ( label ) > _MaximumLabelLength: label = label[0:_MaximumLabelLength]
            self.file.write ( label + "\n" )
        numbers = system.atoms.GetItemAttributes ( "atomicNumber" )
        for i in range ( xyz.rows ):
            symbol = PeriodicTable.Symbol ( numbers[i] )
            self.file.write ( "{:<5s}{:25.15f}{:25.15f}{:25.15f}\n".format ( symbol, xyz[i,0], xyz[i,1], xyz[i,2] ) )

    def WriteSingleSystem ( self, system, label = None, xyz = None ):
        """Write a complete file."""
        self.Open  ( )
        self.WriteFrame ( system, label = label, xyz = xyz )
        self.Close ( )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def XYZFile_FromSystem ( filename, system, label = None, xyz = None ):
    """Helper function that writes a system to a XYZ file."""
    outfile = XYZFileWriter ( filename )
    outfile.WriteSingleSystem ( system, label = label, xyz = xyz )

def XYZFiles_FromSystemGeometryTrajectory ( outPath, inPath, system ):
    """Convert a SystemGeometryTrajectory to XYZ files."""
    # . Define the paths.
    inTrajectory = SystemGeometryTrajectory  ( inPath, system, mode = "r" )
    if not os.path.exists ( outPath ): os.mkdir ( outPath )
    # . Find format.
    n      = len ( "{:d}".format ( len ( inTrajectory ) ) )
    format = "frame{:0" + repr ( n ) + "d}.xyz"
    # . Loop over frames.
    i = 0
    while inTrajectory.RestoreOwnerData ( ):
        XYZFile_FromSystem ( os.path.join ( outPath, format.format ( i+1 ) ), system )
        i += 1
    inTrajectory.Close ( )

# . Exporter definitions.
_Exporter.AddHandler ( { System : XYZFile_FromSystem } , [ "xyz" ], "XYZ" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
