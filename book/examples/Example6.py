"""Example 6."""

from Definitions import *

# . Generate the molecule.
molecule = XYZFile_ToSystem ( os.path.join ( xyzPath, "bala_c7eq.xyz" ) )
molecule.DefineQCModel ( QCModelMNDO ( ) )
molecule.Summary ( )

# . Create an objective function for the molecule.
of = SystemGeometryObjectiveFunction.FromSystem ( molecule )

# . Test the gradients.
of.TestGradients ( )
