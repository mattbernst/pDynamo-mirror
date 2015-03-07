"""Example 17."""

from Definitions import *

# . Read the molecule definition.
molecule = MOLFile_ToSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )
molecule.Summary ( )

# . Define the trajectory.
trajectory = SystemGeometryTrajectory ( os.path.join ( scratchPath, "bala_c7eq.trj" ), molecule, mode = "r" )

# . Loop over the frames in the trajectory.
phi = []
psi = []
while trajectory.RestoreOwnerData ( ):
    phi.append ( molecule.coordinates3.Dihedral ( 4, 6,  8, 14 ) )
    psi.append ( molecule.coordinates3.Dihedral ( 6, 8, 14, 16 ) )

# . Set up the statistics calculation.
phiStatistics = Statistics ( phi )
psiStatistics = Statistics ( psi )

# . Output the results.
table = logFile.GetTable ( columns = [ 20, 20, 20 ] )
table.Start   ( )
table.Title   ( "Phi/Psi Angles" )
table.Heading ( "Frame" )
table.Heading ( "Phi"   )
table.Heading ( "Psi"   )
for ( i, ( h, s ) ) in enumerate ( zip ( phi, psi ) ):
    table.Entry ( "{:d}"  .format ( i ) )
    table.Entry ( "{:.2f}".format ( h ) )
    table.Entry ( "{:.2f}".format ( s ) )
table.Entry ( "Mean:",               alignment = "l" )
table.Entry ( "{:.2f}".format ( phiStatistics.mean ) )
table.Entry ( "{:.2f}".format ( psiStatistics.mean ) )
table.Entry ( "Standard Deviation:", alignment = "l" )
table.Entry ( "{:.2f}".format ( phiStatistics.standardDeviation ) )
table.Entry ( "{:.2f}".format ( psiStatistics.standardDeviation ) )
table.Stop ( )
