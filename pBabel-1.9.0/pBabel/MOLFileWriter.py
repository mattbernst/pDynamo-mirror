#-------------------------------------------------------------------------------
# . File      : MOLFileWriter.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for writing MOL files."""

from pCore        import Coordinates3, TextFileWriter, TextFileWriterError
from ExportImport import _Exporter
from pMolecule    import PeriodicTable, System, \
                         AromaticSingleBond, DoubleBond, SingleBond, TripleBond

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The maximum line length.
_MAXIMUMLINELENGTH = 80

# . Bond type definitions.
_MOLBONDTYPES = { SingleBond ( )          : 1, \
                  DoubleBond ( )          : 2, \
                  TripleBond ( )          : 3, \
                  AromaticSingleBond ( )  : 4  }

#===================================================================================================================================
# . MOL file writer class.
#===================================================================================================================================
class MOLFileWriter ( TextFileWriter ):
    """MOLFileWriter is the class for MOL files that are to be written."""

    def WriteFrame ( self, system, label = None, xyz = None ):
        """Write a single frame."""
        # . Check the data.
        if xyz == None: xyz = system.coordinates3
        if ( xyz is None ) or ( not isinstance ( xyz, Coordinates3 ) ) or ( xyz.rows != len ( system.atoms ) ): raise TextFileWriterError ( "Invalid or missing data to write to MOL file." )
        # . Check the label.
        if label is None: label = system.label
        # . Write the header block.
        if label is None: self.file.write ( "\n\n\n" )
        else:             self.file.write ( label.strip ( )[0:min ( len ( label ), _MAXIMUMLINELENGTH )] + "\n\n\n" )
        # . Counts line.
        self.file.write ( "{:3d}{:3d}".format ( len ( system.atoms ), len ( system.connectivity.bonds ) ) + 8 * "  0" + "999 V2000\n" )
        # . Atom block.
        for ( i, atom ) in enumerate ( system.atoms ):
            self.file.write ( "{:10.4f}{:10.4f}{:10.4f} {:<2s}".format ( xyz[i,0], xyz[i,1], xyz[i,2], PeriodicTable.Symbol ( atom.atomicNumber ) ) + 5 * "  0" + "\n" )
        # . Bond block.
        for ( i, bond ) in enumerate ( system.connectivity.bonds ):
            self.file.write ( "{:3d}{:3d}{:3d}".format ( bond.i + 1, bond.j + 1, _MOLBONDTYPES[bond.type] ) + 3 * "  0" + "\n" )
        # . Properties block.
        for ( i, atom ) in enumerate ( system.atoms ):
            formalCharge = getattr ( atom, "formalCharge", 0 )
            if formalCharge != 0: self.file.write ( "M  CHG  1{:4d}{:4d}\n".format ( i+1, formalCharge ) )
        self.file.write ( "M  END\n" )

    def WriteSingleSystem ( self, system, label = None, xyz = None ):
        """Write a complete file."""
        self.Open  ( )
        self.WriteFrame ( system, label = label, xyz = xyz )
        self.Close ( )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def MOLFile_FromSystem ( filename, system, label = None, xyz = None ):
    """Helper function that writes a system to a MOL file."""
    outfile = MOLFileWriter ( filename )
    outfile.WriteSingleSystem ( system, label = label, xyz = xyz )

# . Exporter definitions.
_Exporter.AddHandler ( { System : MOLFile_FromSystem } , [ "mol", "MOL" ], "MDL MOL" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
