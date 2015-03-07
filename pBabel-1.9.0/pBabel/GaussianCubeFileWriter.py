#-------------------------------------------------------------------------------
# . File      : GaussianCubeFileWriter.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for writing Gaussian cube files.

The following is taken from VMD 1.8:

Write a formatted cubefile very similar to those created by the gaussian program or the cubegen utility.
The format is as follows (last checked against gaussian 98):

Line   Format      Contents
===============================================================
 1     a           Title.
 2     a           Description of the quantity stored in the cubefile.
 3     i5,3f12.6   #atoms, x-,y-,z-coordinates of origin, X0, Y0, Z0
 4     i5,3f12.6   #gridPoints, increment vector, N1, X1, Y1, Z1 - slowest running direction.
 5     i5,3f12.6   #gridPoints, increment vector, N2, X2, Y2, Z2
 6     i5,3f12.6   #gridPoints, increment vector, N3, X3, Y3, Z3 - fastest running direction.

 #atoms lines of atom coordinates:
 ...   i5,4f12.6   Atom number, charge, x-,y-,z-coordinate.

 if #atoms is negative:
 1 line 10i5       Number of orbitals and their numbers.

 rest: 6e13.5      Cube data (with z increment moving fastest,
                   then y and then x) unless there is more than
                   one orbital in which case the orbital number
                   increments fastest. Each row is written out
                   separately so there will be partly empty lines
                   unless N3 is a multiple of 6.

 Other formats that appear to be possible are 6f12.6 and 3d23.16.

All coordinates are in atomic units.

The coordinates of the grid point (I1,I2,I3) are (Fortran convention):

X0+(I1-1)*X1+(I2-1)*X2+(I3-1)*X3
Y0+(I1-1)*Y1+(I2-1)*Y2+(I3-1)*Y3
Z0+(I1-1)*Z1+(I2-1)*Z2+(I3-1)*Z3

"""

import math

from pCore import Coordinates3, logFile, LogFileActive, Real1DArray, TextFileWriter, UNITS_LENGTH_BOHRS_TO_ANGSTROMS, Vector3

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DEFAULT_GRIDSPACING  = 0.2 # . Atomic units.
_DEFAULT_RADIUSFACTOR = 3.0
_DEFAULT_ZERODENSITY  = 1.0e-10

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class RegularCubicGrid3 ( object ):
    """Defines a regular cubic grid."""

    defaultattributes = { "coordinates3" : None ,
                          "gridspacing"  : 0.0  ,
                          "npoints"      : None ,
                          "origin"       : None }

    def __init__ ( self, gridspacing, coordinates3, radii ):
        """Constructor.

        |gridspacing| defines the grid spacing. It is fixed.
        |coordinates3| are the coordinates of the points to be enclosed.
        |radii| are the ranges around each point that must be included in the grid.

        The grid takes the same units as those of coordinates3 and radii.
        """
        self.__dict__.update ( self.__class__.defaultattributes )
        self.gridspacing = gridspacing
        self.npoints     = []
        ( self.origin, extents ) = coordinates3.EnclosingOrthorhombicBox ( radii = radii )
        # . Determine the number of points and readjust the origin and extents accordingly.
        for i in range ( len ( extents ) ):
            n = int ( math.ceil ( extents[i] / gridspacing ) )
            l = float ( n ) * gridspacing
            self.origin[i] -= 0.5 * ( l - extents[i] )
            extents[i] = l
            self.npoints.append ( n )

    def GetGridPointCoordinates ( self ):
        """Return the grid point coordinates."""
        if self.coordinates3 is None:
            coordinates3 = Coordinates3.WithExtent ( self.NumberOfPoints ( ) )
            n = 0
            for ix in range ( self.npoints[0] ):
                for iy in range ( self.npoints[1] ) :
                    for iz in range ( self.npoints[2] ) :
                        point = self.Point ( ix, iy, iz )
                        coordinates3[n,0] = point[0]
                        coordinates3[n,1] = point[1]
                        coordinates3[n,2] = point[2]
                        n += 1
            self.coordinates3 = coordinates3
        return self.coordinates3

    def NumberOfPoints ( self ):
        """Return the number of points on the grid."""
        return ( self.npoints[0] * self.npoints[1] * self.npoints[2] )

    def Point ( self, ix, iy, iz ):
        """Return the coordinates of a point."""
        point   = Vector3.Null ( )
        point[0] = self.origin[0] + float ( ix ) * self.gridspacing
        point[1] = self.origin[1] + float ( iy ) * self.gridspacing
        point[2] = self.origin[2] + float ( iz ) * self.gridspacing
        return point

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "Regular 3D Cubic Grid" )
            summary.Entry ( "Number of Points"   , "{:d}"  .format (         self.NumberOfPoints ( ) ) )
            summary.Entry ( "Grid Spacing"       , "{:.2f}".format (         self.gridspacing        ) )
            summary.Entry ( "Origin x"           , "{:.2f}".format (         self.origin [0]         ) )
            summary.Entry ( "Extent x"           , "{:.2f}".format ( float ( self.npoints[0] ) * self.gridspacing ) )
            summary.Entry ( "Origin y"           , "{:.2f}".format (         self.origin [1]         ) )
            summary.Entry ( "Extent y"           , "{:.2f}".format ( float ( self.npoints[1] ) * self.gridspacing ) )
            summary.Entry ( "Origin z"           , "{:.2f}".format (         self.origin [2]         ) )
            summary.Entry ( "Extent z"           , "{:.2f}".format ( float ( self.npoints[2] ) * self.gridspacing ) )
            summary.Entry ( "Number of x Points" , "{:d}"  .format (         self.npoints[0]         ) )
            summary.Entry ( "Number of y Points" , "{:d}"  .format (         self.npoints[1]         ) )
            summary.Entry ( "Number of z Points" , "{:d}"  .format (         self.npoints[2]         ) )
            summary.Stop ( )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GaussianCubeFileWriter ( TextFileWriter ):
    """Class for writing Gaussian cube files."""

    defaultattributes = { "grid"           : None  ,
                          "qcCoordinates3" : None  ,
                          "qcSelection"    : None  ,
                          "system"         : None  }
    defaultattributes.update ( TextFileWriter.defaultattributes )

    def __init__ ( self, path, system ):
        """Constructor."""
        super ( GaussianCubeFileWriter, self ).__init__ ( path )
        self.system      = system
        self.qcSelection = system.energyModel.qcAtoms.GetFullSelection ( )

    def DefineGrid ( self, gridspacing = _DEFAULT_GRIDSPACING, radiusfactor = _DEFAULT_RADIUSFACTOR ):
        """Define the grid."""
        if ( self.system is not None ) and ( self.system.coordinates3 is not None ) and ( self.system.energyModel is not None ) and ( self.system.energyModel.qcModel is not None ):
            self.qcCoordinates3 = self.system.energyModel.qcAtoms.GetCoordinates3 ( self.system.coordinates3, toBohrs = True )
            allRadii            = self.system.atoms.GetItemAttributes ( "vdwRadius" )
            if len ( self.qcSelection ) == len ( self.system.atoms ):
                radii = allRadii
            else:
                radii = Real1DArray ( len ( self.qcSelection ) )
                for ( i, s ) in enumerate ( self.qcSelection ):
                    radii[i] = allRadii[s]
            radii.Scale ( radiusfactor / UNITS_LENGTH_BOHRS_TO_ANGSTROMS )
            self.grid = RegularCubicGrid3 ( gridspacing, self.qcCoordinates3, radii )

    def WriteSystemData ( self, datatype, orbitals = None, useDensityP = True, spinDensities = False ):
        """Get the system data."""
        # . Grid points.
        gridPoints = self.grid.GetGridPointCoordinates ( )
        # . Generate the data.
        if   datatype == "Density"  : data = self.system.energyModel.qcModel.GridPointDensities  ( self.system.configuration, gridPoints, spinDensities = spinDensities )
        elif datatype == "Orbitals" : data = self.system.energyModel.qcModel.GridPointOrbitals   ( self.system.configuration, gridPoints, orbitals = orbitals, useDensityP = useDensityP )
        elif datatype == "Potential": data = self.system.energyModel.qcModel.GridPointPotentials ( self.system.configuration, gridPoints )
        else: raise ValueError ( "Unknown grid data type: " + datatype + "." )
        # . Get the multiplicity of the data.
        if datatype == "Orbitals": ndatum = len ( orbitals )
        else:                      ndatum = 1
        # . Write the data - explicit loop required here as each row is written out separately.
        m = 0
        n = 0
        for ix in range ( self.grid.npoints[0] ):
            for iy in range ( self.grid.npoints[1] ) :
                for iz in range ( self.grid.npoints[2] ) :
                    for i in range ( ndatum ):
                        self.file.write ( "{:13.5e}".format ( data[m] ) )
                        m += 1
                        n += 1
                        if ( n > 5 ):
                            self.file.write ( "\n" )
                            n = 0
                if n != 0 :
                    self.file.write ( "\n" )
                    n = 0
        # . Finish up.
        self.Close ( )

    def WriteHeader ( self, QORBITALS = False ):
        """Write the header to the file."""
        self.Open ( )
        if self.system.label is None: self.file.write ( "Gaussian Cube File.\n" )
        else:                         self.file.write ( "Gaussian Cube File for " + self.system.label + ".\n" )
        self.file.write ( "X: outer loop, Y: middle loop, Z: inner loop.\n" )
        if QORBITALS: self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( - len ( self.system.energyModel.qcAtoms.QCAtomSelection ( ) ), self.grid.origin[0], self.grid.origin[1], self.grid.origin[2] ) )
        else:         self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format (   len ( self.system.energyModel.qcAtoms.QCAtomSelection ( ) ), self.grid.origin[0], self.grid.origin[1], self.grid.origin[2] ) )
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( self.grid.npoints[0], self.grid.gridspacing, 0.0, 0.0 ) )
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( self.grid.npoints[1], 0.0, self.grid.gridspacing, 0.0 ) )
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( self.grid.npoints[2], 0.0, 0.0, self.grid.gridspacing ) )
        numbers = self.system.atoms.GetItemAttributes ( "atomicNumber" )
        for ( i, s ) in enumerate ( self.qcSelection ):
            self.file.write ( "{:<5d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n".format ( numbers[s], float ( numbers[s] ), self.qcCoordinates3[i,0], self.qcCoordinates3[i,1], self.qcCoordinates3[i,2] ) )

    def WriteOrbitalHeader ( self, orbitals ):
        """Write the orbital header."""
        if orbitals is None: orbitals = range ( 0, 1 )
        if len ( orbitals ) > 0:
            self.file.write ( "{:5d}".format ( len ( orbitals ) ) )
            for i in orbitals: self.file.write ( "{:5d}".format ( i ) )
            self.file.write ( "\n" )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def GaussianCubeFile_FromSystemDensity ( path, system, gridspacing = _DEFAULT_GRIDSPACING, spinDensities = False, radiusfactor = _DEFAULT_RADIUSFACTOR, log = logFile ):
    """Helper function that writes a system's electronic density to a Gaussian cube file."""
    outfile = GaussianCubeFileWriter ( path, system )
    outfile.DefineGrid ( gridspacing, radiusfactor )
    outfile.grid.Summary    ( log = log )
    outfile.WriteHeader     ( )
    outfile.WriteSystemData ( "Density", spinDensities = spinDensities )

def GaussianCubeFile_FromSystemOrbitals ( path, system, gridspacing = _DEFAULT_GRIDSPACING, orbitals = None, useDensityP = True, radiusfactor = _DEFAULT_RADIUSFACTOR, log = logFile ):
    """Helper function that writes a system's orbitals to a Gaussian cube file."""
    outfile = GaussianCubeFileWriter ( path, system )
    outfile.DefineGrid ( gridspacing, radiusfactor )
    outfile.grid.Summary       ( log = log )
    outfile.WriteHeader        ( QORBITALS = True )
    outfile.WriteOrbitalHeader ( orbitals )
    outfile.WriteSystemData    ( "Orbitals", orbitals = orbitals, useDensityP = useDensityP )

def GaussianCubeFile_FromSystemPotential ( path, system, gridspacing = _DEFAULT_GRIDSPACING, radiusfactor = _DEFAULT_RADIUSFACTOR, log = logFile ):
    """Helper function that writes a system's electrostatic potential to a Gaussian cube file."""
    outfile = GaussianCubeFileWriter ( path, system )
    outfile.DefineGrid ( gridspacing, radiusfactor )
    outfile.grid.Summary    ( log = log )
    outfile.WriteHeader     ( )
    outfile.WriteSystemData ( "Potential" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
