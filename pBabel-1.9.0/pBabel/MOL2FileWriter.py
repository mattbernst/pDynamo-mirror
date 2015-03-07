#-------------------------------------------------------------------------------
# . File      : MOL2FileWriter.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for writing MOL2 files."""

from pCore        import Coordinates3, TextFileWriter, TextFileWriterError
from ExportImport import _Exporter
from pMolecule    import PeriodicTable, System, \
                         AromaticDoubleBond, AromaticSingleBond, DoubleBond, SingleBond, TripleBond, UndefinedBond

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Atom type.
_ATOMTYPE = "Any"

# . The charge type.
_CHARGETYPE = "NO_CHARGES"

# . The maximum line length.
_MAXIMUMLINELENGTH = 80

# . Bond type definitions.
_MOL2BONDTYPES = { SingleBond         ( ) : "1" , \
                   DoubleBond         ( ) : "2" , \
                   TripleBond         ( ) : "3" , \
                   AromaticSingleBond ( ) : "ar", \
                   AromaticDoubleBond ( ) : "ar", \
                   UndefinedBond      ( ) : "un"  }

# . The molecule type.
_MOLECULETYPE = "SMALL"

# . Section header.
_RTISTRING  = "@<TRIPOS>"

#===================================================================================================================================
# . MOL2 file writer class.
#===================================================================================================================================
class MOL2FileWriter ( TextFileWriter ):
    """MOL2FileWriter is the class for MOL2 files that are to be written."""

    def WriteFrame ( self, system, label = None, xyz = None ):
        """Write a single frame."""
        # . Get data.
        natoms = len ( system.atoms )
        nbonds = len ( system.connectivity.bonds )
        # . Check the data.
        if xyz == None: xyz = system.coordinates3
        if ( xyz is None ) or ( not isinstance ( xyz, Coordinates3 ) ) or ( xyz.rows != natoms ): raise TextFileWriterError ( "Invalid or missing data to write to MOL2 file." )
        # . Check the label.
        if label is None: label = system.label
        # . Molecule block.
        self.file.write ( _RTISTRING + "MOLECULE\n" )
        if label is None: self.file.write ( "\n" )
        else:             self.file.write ( label.strip ( )[0:min ( len ( label ), _MAXIMUMLINELENGTH )] + "\n" )
        self.file.write ( "{:d} {:d}\n".format ( natoms, nbonds ) )
        self.file.write ( _MOLECULETYPE + "\n" )
        self.file.write ( _CHARGETYPE   + "\n" )
        self.file.write ( "\n" )
        # . Atom block.
        self.file.write ( _RTISTRING + "ATOM\n" )
        for ( i, atom ) in enumerate ( system.atoms ):
            self.file.write ( "{:6d} {:<6s} {:10.4f} {:10.4f} {:10.4f}   {:s}\n".format ( i+1,  PeriodicTable.Symbol ( atom.atomicNumber, index = i ), xyz[i,0], xyz[i,1], xyz[i,2], _ATOMTYPE ) )
        self.file.write ( "\n" )
        # . Bond block.
        self.file.write ( _RTISTRING + "BOND\n" )
        for ( i, bond ) in enumerate ( system.connectivity.bonds ):
            self.file.write ( "{:6d} {:6d} {:6d} {:<s}\n".format ( i+1, bond.i + 1, bond.j + 1, _MOL2BONDTYPES[bond.type] ) )
        self.file.write ( "\n" )
        # . Symmetry block.
        if system.symmetryParameters is not None:
            p = system.symmetryParameters
            self.file.write ( _RTISTRING + "CRYSIN\n" )
            self.file.write ( "{:.3f} {:.3f} {:.3f} {:.1f} {:.1f} {:.1f}  1  1\n".format ( p.a, p.b, p.c, p.alpha, p.beta, p.gamma ) )

    def WriteSingleSystem ( self, system, label = None, xyz = None ):
        """Write a complete file."""
        self.Open  ( )
        self.WriteFrame ( system, label = label, xyz = xyz )
        self.Close ( )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def MOL2File_FromSystem ( filename, system, label = None, xyz = None ):
    """Helper function that writes a system to a MOL2 file."""
    outfile = MOL2FileWriter ( filename )
    outfile.WriteSingleSystem ( system, label = label, xyz = xyz )

# . Exporter definitions.
_Exporter.AddHandler ( { System : MOL2File_FromSystem } , [ "mol2", "MOL2" ], "Tripos MOL2" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
