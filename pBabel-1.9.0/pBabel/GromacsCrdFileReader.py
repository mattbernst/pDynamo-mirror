#-------------------------------------------------------------------------------
# . File      : GromacsCrdFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
#===========================================================================================
#
# . File      : GromacsCrdFileReader.py
#
# . Author    : Guilherme M. Arantes (University of Sao Paulo, Brazil, 2011)
#
# . Based on and to be used with the pDynamo library, copyright CEA, CNRS, Martin J. Field 
#
# . Web       : http://www.pdynamo.org
#
#===========================================================================================
"""Read data from a Gromacs .gro coordinate file."""

import sys

from pCore        import Coordinates3, logFile, LogFileActive, TextFileReader
from ExportImport import _Importer
from pMolecule    import CrystalClassCubic, CrystalClassTriclinic, CrystalClassOrthorhombic, SymmetryParameters

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GromacsCrdFileReader ( TextFileReader ):
    """Class for reading Gromacs .gro coordinate files."""

    def GetAtomLineFormat ( self ):
        """Get the format of the atom lines.

        This is:
	RESNO RES   TYPE ATOMNO   X     Y     Z  
	I5    A4 2X A4   I5       F8.3  F8.3  F8.3 
        """
        a = 4 ; f = 8 ; i = 5 ; x = 2
        p = 0
        format =      [ ( p, p+i, int  , 0   ) ] ; p +=  i    # . Res. number.
        format.append ( ( p, p+a, None , ""  ) ) ; p += (a+x) # . Res. name.
        format.append ( ( p, p+a, None , ""  ) ) ; p += (a+i) # . Atom name.
        format.append ( ( p, p+f, float, 0.0 ) ) ; p +=  f    # . X.
        format.append ( ( p, p+f, float, 0.0 ) ) ; p +=  f    # . Y.
        format.append ( ( p, p+f, float, 0.0 ) )              # . Z.
        return format

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Get the atom line format.
            atomlineformat = self.GetAtomLineFormat ( )
	    # . Open the file.
	    self.Open ( )
	    # . Parse all the lines.
            try:
                # . Keyword line.
                self.title = self.GetLine ( )
                # . Number of atoms.
                items    = self.GetTokens ( converters = ( int, ) )
                natoms   = items[0]
                # . The coordinate data.
                self.xyz = Coordinates3 ( natoms )
                for n in range ( natoms ):
                    tokens = self.GetFixedFormatTokens ( *atomlineformat )
                    for i in range ( 3 ): self.xyz[n,i] = float ( tokens[i+3]*10.0 )
                # . Symmetry data.
                self.symmetryItems = self.GetTokens ( )
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
            summary.Start ( "Gromacs .gro File Summary" )
            summary.Entry ( "Number of Atom Lines", "{:d}".format ( self.xyz.rows ) )
            summary.Stop ( )

    def ToCoordinates3 ( self ):
        """Return the coordinates."""
        if self.QPARSED: return self.xyz
        else:            return None

    def ToSymmetryParameters ( self ):
        """Return the symmetry parameters."""
        # . Will assign only cubic, orthorhombic, dodecahedron or octahedron boxes. 
	# . Will fail for general triclinic, but this could be hard-coded once sizes/angles are known.
        if self.QPARSED: 
	    # . triclinic box
	    items = self.symmetryItems
            if len ( items ) == 9: 
	        specialTriclinic = items[8] == items[7] or "-" + items[7] == items[5]
	        if not specialTriclinic: self.Warning ( "Invalid general triclinic box symmetry.", True )
	        # . Dodecahedron 
	        if   items[0] == items[1]:
                    alpha = 60.0
	            beta  = 60.0
	            gamma = 90.0
	        # . Octahedron
	        else: 
                    alpha = 70.53
	            beta  = 109.47
	            gamma = 70.53
	        a = float ( items[0] ) * 10.0
	        b = a
	        c = a
            # . Cubic or orthorhombic box
	    elif len ( items ) == 3  :
                alpha = 90.0
	        beta  = 90.0
	        gamma = 90.0
	        items = [ float ( items[i] ) * 10.0 for i in range(3) ]
	        a = items[0]
	        b = items[1]
	        c = items[2]
	    else: self.Warning ( "Invalid or unrecognized box symmetry.", True )
            self.symmetryParameters = SymmetryParameters ( )
            self.symmetryParameters.SetCrystalParameters (  a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma )
	    return self.symmetryParameters
        else:            return None

    def ToSymmetry ( self, system ):
        """Assign symmetry to a system."""
        # . Check for symmetry.
        if not self.QPARSED: return None
	else: 
            # . Unpack data.
            alpha = self.symmetryParameters.alpha 
            beta  = self.symmetryParameters.beta
            gamma = self.symmetryParameters.gamma
            a     = self.symmetryParameters.a
            b     = self.symmetryParameters.b
            c     = self.symmetryParameters.c
            # . Get the crystal class.
	    if alpha > 80.0:
                if   ( a == b ) and ( a == c ): crystalClass = CrystalClassCubic        ( )
                else:                           crystalClass = CrystalClassOrthorhombic ( )
	    else           :                    crystalClass = CrystalClassTriclinic    ( )
            # . Define symmetry.
            system.DefineSymmetry ( crystalClass = crystalClass,  a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma )

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def GromacsCrdFile_Process ( filename, log = logFile, system = None ):
    """Helper function that returns a set of coordinates from a Gromacs .gro file."""
    infile = GromacsCrdFileReader ( filename )
    infile.Parse   ( log = log )
    infile.Summary ( log = log )
    infile.ToSymmetryParameters ( )
    infile.ToSymmetry ( system ) 
    return infile.ToCoordinates3 ( )

def GromacsCrdFile_ToCoordinates3 ( filename, log = logFile ):
    """Helper function that returns a set of coordinates from a Gromacs .gro file."""
    infile = GromacsCrdFileReader ( filename )
    infile.Parse   ( log = log )
    infile.Summary ( log = log )
    return infile.ToCoordinates3 ( )

def GromacsCrdFile_ToSymmetry ( filename, log = logFile ):
    """Helper function that reads symmetryparamters from a Gromacs .gro file."""
    infile = GromacsCrdFileReader ( filename )
    infile.Parse   ( log = log )
    infile.Summary ( log = log )
    return infile.ToSymmetryParameters ( )

# . Importer definitions.
_Importer.AddHandler ( { Coordinates3 : GromacsCrdFile_ToCoordinates3 } , [ "gro", "GRO" ], "Gromacs Coordinates", defaultFunction = GromacsCrdFile_ToCoordinates3 )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

