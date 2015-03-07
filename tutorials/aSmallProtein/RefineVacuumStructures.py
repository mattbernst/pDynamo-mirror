"""Refinement of vacuum structures."""

from Definitions import *

# . PDB files.
_PDBPaths = ( "1UAO", "2E4E" )

# . Parameters for the minimization.
_Iterations   = 1000
_Logfrequency =  100
_Tolerance    =  1.0

# . Loop over structures.
for pdbPath in _PDBPaths:

    # . Retrieve the system.
    system = Unpickle ( os.path.join ( outPath, pdbPath + ".pkl" ) )
    system.Summary ( )
    system.Energy ( )

    # . Get selection for the hydrogens and heavy atoms.
    hydrogens = Selection.FromIterable ( AtomSelection.FromAtomAttributes ( system, "atomicNumber", 1 ) )
    heavies   = hydrogens.Complement ( upperBound = len ( system.atoms ) )

    # . Minimize the coordinates of the hydrogens.
    system.DefineFixedAtoms ( hydrogens )
    ConjugateGradientMinimize_SystemGeometry ( system,                               \
                                               maximumIterations    = _Iterations,   \
                                               logFrequency         = _Logfrequency, \
                                               rmsGradientTolerance = _Tolerance     )
    system.DefineFixedAtoms ( None )

    # . Loop over minimizations.
    for forceConstant in ( 1000.0, 500.0, 250.0, 100.0, 50.0 ):

        # . Harmonically constrain heavy atoms.
        tethers = None
        if ( forceConstant is not None ):
            reference          = Clone ( system.coordinates3 )
            tetherEnergyModel  = SoftConstraintEnergyModelHarmonic ( 0.0, forceConstant )
            tethers            = SoftConstraintContainer ( )
            tethers["tethers"] = SoftConstraintMultipleTether ( heavies, reference, tetherEnergyModel )
        system.DefineSoftConstraints ( tethers )

        # . Minimize.
        ConjugateGradientMinimize_SystemGeometry ( system,                               \
                                                   maximumIterations    = _Iterations,   \
                                                   logFrequency         = _Logfrequency, \
                                                   rmsGradientTolerance = _Tolerance     )
        system.Energy ( doGradients = True )

    # . Save the refined coordinates.
    XYZFile_FromSystem ( os.path.join ( outPath, pdbPath + "_folded.xyz" ), system )

