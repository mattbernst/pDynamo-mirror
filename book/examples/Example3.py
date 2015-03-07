"""Example 3."""

from Definitions import *

# . Read in a system.
molecule = XYZFile_ToSystem ( os.path.join ( xyzPath, "bala_c7eq.xyz" ) )
molecule.Summary ( )

# . Define a table for the results.
table = logFile.GetTable ( columns = [ 15, 15 ] )
table.Start ( )
table.Title ( "Bond Analysis" )
table.Heading ( "Safety Factor" )
table.Heading ( "Bonds Found"   )

# . Loop over the buffer sizes.
for i in range ( 21 ):
    safety = 0.1 * float ( i )
    molecule.BondsFromCoordinates3 ( safety = safety )
    table.Entry ( "{:4.1f}".format ( safety ) )
    table.Entry ( "{:d}".format ( len ( molecule.connectivity.bonds ) ) )

# . Finish up the table.
table.Stop ( )

# . Generate the bonds with the default safety factor.
molecule.BondsFromCoordinates3 ( safety = 0.5 )
molecule.Summary ( )

# . Print the bonds.
table = logFile.GetTable ( columns = 4 * [ 5, 5, 10 ] )
table.Start ( )
table.Title ( "Bond Lengths (Angstroms)" )
for bond in molecule.connectivity.bonds:
    table.Entry ( "{:d}".format ( bond.i ) )
    table.Entry ( "{:d}".format ( bond.j ) )
    table.Entry ( "{:6.3f}".format ( molecule.coordinates3.Distance ( bond.i, bond.j ) ) )
table.Stop ( )

# . Print the angles.
table = logFile.GetTable ( columns = 3 * [ 5, 5, 5, 10 ] )
table.Start ( )
table.Title ( "Angles (Degrees)" )
for angle in molecule.connectivity.angles:
    table.Entry ( "{:d}".format ( angle.i ) )
    table.Entry ( "{:d}".format ( angle.j ) )
    table.Entry ( "{:d}".format ( angle.k ) )
    table.Entry ( "{:6.1f}".format ( molecule.coordinates3.Angle ( angle.i, angle.j, angle.k ) ) )
table.Stop ( )

# . Print the dihedrals.
table = logFile.GetTable ( columns = 4 * [ 5, 5, 5, 5, 10 ] )
table.Start ( )
table.Title ( "Dihedrals (Degrees)" )
for dihedral in molecule.connectivity.dihedrals:
    table.Entry ( "{:d}".format ( dihedral.i ) )
    table.Entry ( "{:d}".format ( dihedral.j ) )
    table.Entry ( "{:d}".format ( dihedral.k ) )
    table.Entry ( "{:d}".format ( dihedral.l ) )
    table.Entry ( "{:6.1f}".format ( molecule.coordinates3.Dihedral ( dihedral.i, dihedral.j, dihedral.k, dihedral.l ) ) )
table.Stop ( )
