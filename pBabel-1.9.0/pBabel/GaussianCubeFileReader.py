#-------------------------------------------------------------------------------
# . File      : GaussianCubeFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for reading Gaussian cube files.

See the preamble to GaussianCubeFileReader for more information about the format.
"""

import math

from pCore        import Coordinates3, logFile, LogFileActive, RealNDArray, RegularGrid_FromDimensionData, TextFileReader, UNITS_LENGTH_BOHRS_TO_ANGSTROMS
from ExportImport import _Importer
from pMolecule    import System

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GaussianCubeFileReader ( TextFileReader ):
    """Class for reading Gaussian cube files."""

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            # . Parse the data.
            try:
                # . Header.
                self.label = self.GetLine ( )
                # . More description - ignored.
                self.GetLine ( )
                # . First definition line.
                ( natoms, ox, oy, oz ) = self.GetFixedFormatTokens ( ( 0, 5, int, 0 ), ( 5, 17, float, 0.0 ), ( 17, 29, float, 0.0 ), ( 29, 41, float, 0.0 ) )
                hasOrbitals = ( natoms < 0 )
                if hasOrbitals: natoms = abs ( natoms )
                # . Grid axis definitions - number of points and increment vector.
                griddata = []
                ngrid    = 1
                for ( i, o ) in enumerate ( ( ox, oy, oz ) ):
                    items = self.GetFixedFormatTokens ( ( 0, 5, int, 0 ), ( 5, 17, float, 0.0 ), ( 17, 29, float, 0.0 ), ( 29, 41, float, 0.0 ) )
                    n     = items.pop ( 0 )
                    h     = items.pop ( i ) * UNITS_LENGTH_BOHRS_TO_ANGSTROMS
                    o    *= UNITS_LENGTH_BOHRS_TO_ANGSTROMS
                    if ( items[0] != 0.0 ) or ( items[1] != 0.0 ): self.Warning ( "Unable to handle grids whose increment vectors are not aligned along the Cartesian axes.", False )
                    griddata.append ( { "bins" : n, "binSize" : h, "lower" : o - 0.5 * h } )
                    ngrid *= n
                # . Get the grid.
                self.grid = RegularGrid_FromDimensionData ( griddata )
                # . Atom data.
                self.atomicNumbers = []
                self.coordinates3  = Coordinates3.WithExtent ( natoms )
                for i in range ( natoms ):
                    ( n, q, x, y, z ) = self.GetFixedFormatTokens ( ( 0, 5, int, 0 ), ( 5, 17, float, 0.0 ), ( 17, 29, float, 0.0 ), ( 29, 41, float, 0.0 ), ( 41, 53, float, 0.0 ) )
                    self.atomicNumbers.append ( n )
                    self.coordinates3[i,0] = x
                    self.coordinates3[i,1] = y
                    self.coordinates3[i,2] = z
                self.coordinates3.Scale ( UNITS_LENGTH_BOHRS_TO_ANGSTROMS )
                # . Orbital data.
                if hasOrbitals:
                    indices = self.GetTokens ( )
                    for ( i, index ) in enumerate ( indices ): indices[i] = int ( index )
                    norbitals = indices.pop ( 0 )
                    if norbitals != len ( indices ): self.Warning ( "The number of orbitals does not match the number of orbital indices.", True )
                    self.orbitalindices = indices
                else:
                    self.orbitalindices = None
                # . Field data.
                # . Unfortunately each row is written out separately so reading in has to be done by row as well.
                keys = []
                if hasOrbitals:
                    for index in self.orbitalindices:
                        keys.append ( "Orbital {:d}".format ( index ) )
                else:
                    keys.append ( "Density" )
                self.fielddata = {}
                for key in keys:
                    self.fielddata[key] = RealNDArray.WithExtents ( *self.grid.shape )
                nfields = len ( keys )
                if nfields == 1:
                    fielddata = self.fielddata[keys[0]] # []
                else:
                    fielddata = []
                    for key in keys:
                        fielddata.append ( self.fielddata[key] )
                nitems  = griddata[2]["bins"] * nfields
                for ix in range ( griddata[0]["bins"] ):
                    for iy in range ( griddata[1]["bins"] ):
                        newdata = self.GetFixedFormatArray ( nitems, 6, 13, converter = float, default = 0.0 )
                        if nfields == 1:
#                            fielddata.extend ( newdata )
                            for iz in range ( nitems ):
                                fielddata[ix,iy,iz] = newdata[iz]
                        else:
                            for ifield in range ( nfields ):
                                for ( iz, i ) in enumerate ( range ( ifield, nitems, nfields ) ):
                                    fielddata[ifield][ix,iy,iz] = newdata[i]
#                            for ifield in range ( nfields ):
#                                fielddata[ifield].extend ( newdata[ifield::nfields] )
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
        if self.QPARSED: return self.coordinates3
        else:            return None

    def ToSystem ( self ):
        """Return a system."""
        if self.QPARSED:
            system = System.FromAtoms ( self.atomicNumbers )
            system.label        = self.label
            system.coordinates3 = self.ToCoordinates3 ( )
            return system
        else:
            return None

    def ToVolumetricData ( self ):
        """Return the volumetric data."""
        # . Temporarily the grid will go here too.
        data = {}
        if self.QPARSED:
            data = self.fielddata
            data["Grid"] = self.grid
        return data

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def GaussianCubeFile_ToCoordinates3 ( path, log = logFile ):
    """Helper function that reads coordinates3 from a Gaussian cube file."""
    infile = GaussianCubeFileReader ( path )
    infile.Parse ( log = log )
    system = infile.ToCoordinates3 ( )
    return system

def GaussianCubeFile_ToSystem ( path, log = logFile, volumetricdata = None ):
    """Helper function that reads a system from a Gaussian cube file."""
    infile = GaussianCubeFileReader ( path )
    infile.Parse ( log = log )
    system = infile.ToSystem ( )
    if volumetricdata is not None: volumetricdata.update ( infile.ToVolumetricData ( ) )
    return system

# . Importer definitions.
_Importer.AddHandler ( { Coordinates3 : GaussianCubeFile_ToCoordinates3 ,
                         System       : GaussianCubeFile_ToSystem       } , [ "cub", "CUB", "cube", "CUBE" ], "Gaussian Cube File", defaultFunction = GaussianCubeFile_ToSystem )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
