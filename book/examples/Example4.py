"""Example 4."""

from Definitions import *

# . Define the list of structures.
xyzFiles = [ "bala_alpha.xyz", "bala_c5.xyz", "bala_c7ax.xyz", "bala_c7eq.xyz" ]

# . Define a molecule.
xyzFile  = xyzFiles.pop ( )
molecule = XYZFile_ToSystem ( os.path.join ( xyzPath, xyzFile ) )
molecule.Summary ( )

# . Translate the system to its center of mass.
masses = molecule.atoms.GetItemAttributes ( "mass" )
molecule.coordinates3.TranslateToCenter ( weights = masses )

# . Calculate and print the inertia matrix before reorientation.
inertia = molecule.coordinates3.InertiaMatrix ( weights = masses )
inertia.Print ( title = "Inertia Matrix Before Reorientation" )

# . Transform to principal axes.
molecule.coordinates3.ToPrincipalAxes ( weights = masses )

# . Calculate and print the inertia matrix after reorientation.
inertia = molecule.coordinates3.InertiaMatrix ( weights = masses )
inertia.Print ( title = "Inertia Matrix After Reorientation" )

# . Define a table for the results.
table = logFile.GetTable ( columns = [ 20, 10, 10 ] )
table.Start ( )
table.Title ( "RMS Coordinate Deviations" )
table.Heading ( "Structure" )
table.Heading ( "Before"    )
table.Heading ( "After"     )

# . Loop over the remaining structures.
for xyzFile in xyzFiles:
    coordinates3 = XYZFile_ToCoordinates3 ( os.path.join ( xyzPath, xyzFile ) )
    rms0 = coordinates3.RMSDeviation ( molecule.coordinates3, weights = masses )
    coordinates3.Superimpose ( molecule.coordinates3, weights = masses )
    rms1 = coordinates3.RMSDeviation ( molecule.coordinates3, weights = masses )
    table.Entry ( xyzFile[0:-4], alignment = "l" )
    table.Entry ( "{:.2f}".format ( rms0 ) )
    table.Entry ( "{:.2f}".format ( rms1 ) )

# . Finish up the table.
table.Stop ( )
