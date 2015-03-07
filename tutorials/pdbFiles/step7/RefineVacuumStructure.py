"""Refine the vacuum structure for PKA."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions      import outPath
from pCore            import Clone, Selection, Pickle, Unpickle
from pMolecule        import AtomSelection, NBModelABFS, SoftConstraintContainer, SoftConstraintEnergyModelHarmonic, SoftConstraintMultipleTether
from pMoleculeScripts import ConjugateGradientMinimize_SystemGeometry

# . Parameters.
_Iterations   = 1000
_Logfrequency =  100
_Tolerance    =  1.0

# . Get the system with an appropriate NB model.
system = Unpickle ( os.path.join ( outPath, "step6.pkl" ) )
system.DefineNBModel ( NBModelABFS ( ) )
system.Summary ( )

# . Get a selection with the indices of the atoms whose coordinates were built.
built = AtomSelection.FromAtomAttributes ( system, "atomicNumber", 1 )
built.Add ( system.sequence.AtomIndex ( "A:LYS.8:CB"  ) )
built.Add ( system.sequence.AtomIndex ( "I:SER.17:OG" ) )
built = Selection.FromIterable ( built )

# . Minimize the coordinates of these atoms.
system.DefineFixedAtoms ( built.Complement ( upperBound = len ( system.atoms ) ) )
ConjugateGradientMinimize_SystemGeometry ( system,                               \
                                           maximumIterations    = _Iterations,   \
                                           logFrequency         = _Logfrequency, \
                                           rmsGradientTolerance = _Tolerance     )
system.DefineFixedAtoms ( None )

# . Get a selection corresponding to heavy atoms.
hydrogens = Selection.FromIterable ( AtomSelection.FromAtomAttributes ( system, "atomicNumber", 1 ) )
heavies   = hydrogens.Complement ( upperBound = len ( system.atoms ) )

# . Loop over minimizations.
for forceConstant in ( 1000.0, 500.0, 250.0, 100.0, 50.0, 10.0, 4.0, 1.0, 0.25, None ):

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

# . Save the system.
system.configuration.Clear ( )
Pickle ( os.path.join ( outPath, "step7.pkl" ), system )
