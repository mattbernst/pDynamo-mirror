"""Example 21."""

from Definitions import *

# . Read the system definition.
solvent = Unpickle ( os.path.join ( scratchPath, "water216_cubicBox.pkl" ) )
solvent.Summary ( )

# . Select all oxygens.
indices = []
for ( i, atom ) in enumerate ( solvent.atoms ):
    if atom.atomicNumber == 8: indices.append ( i )
oxygens = Selection.FromIterable ( indices )

# . Analyse the trajectory data.
trajectory = SystemGeometryTrajectory ( os.path.join ( scratchPath, "water216_cubicBox.trj" ), solvent, mode = "r" )

# . Self-diffusion function.
SelfDiffusionFunction ( trajectory, oxygens )

# . Radial distribution function.
RadialDistributionFunction ( trajectory, oxygens )
