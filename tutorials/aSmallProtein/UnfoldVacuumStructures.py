"""Constrained minimization to unfold vacuum structures."""

from Definitions import *

# . PDB files.
_PDBPaths = ( "1UAO", "2E4E" )

# . Define parameters for the minimization and soft constraint.
_ForceConstant  = 100.0
_SCName         = "scMD"
_TargetDistance =  20.0

# . Loop over structures.
for pdbPath in _PDBPaths:

    # . Retrieve system.
    system = Unpickle ( os.path.join ( outPath, pdbPath + ".pkl" ) )
    system.coordinates3 = XYZFile_ToCoordinates3 ( os.path.join ( outPath, pdbPath + "_folded.xyz" ) )
    system.Summary ( )
    system.Energy ( )

    # . Get atom indices.
    nTer = system.sequence.AtomIndex ( "A:GLY.1:N"  )
    cTer = system.sequence.AtomIndex ( "A:GLY.10:C" )

    # . Get the starting constraint distance (rounded to the nearest tenth of an Angstrom).
    distance0 = system.coordinates3.Distance ( nTer, cTer )

    # . Set up the constraint.
    constraints = SoftConstraintContainer ( )
    system.DefineSoftConstraints ( constraints )
    scModel    = SoftConstraintEnergyModelHarmonic ( _TargetDistance, _ForceConstant )
    constraint = SoftConstraintDistance ( nTer, cTer, scModel )
    constraints[_SCName] = constraint

    # . Optimize.
    system.Energy ( doGradients = True )
    ConjugateGradientMinimize_SystemGeometry ( system,                      \
                                               logFrequency         =  100, \
                                               maximumIterations    = 5000, \
                                               rmsGradientTolerance =  0.5  )

    # . Get the energy and constraint distance.
    system.DefineSoftConstraints ( None )
    energy   = system.Energy ( doGradients = True )
    distance1 = system.coordinates3.Distance ( nTer, cTer )

    # . Print the starting and stopping distances.
    print ( "\n\nStarting and stopping distances = {:.2f} {:.2f}.\n".format ( distance0, distance1 ) )

    # . Save structures.
    XYZFile_FromSystem ( os.path.join ( outPath, pdbPath + "_unfolded.xyz" ), system )
