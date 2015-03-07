"""Solvate the system using a water box."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions      import outPath                                                                                                                         
from pBabel           import PDBFile_FromSystem                                                                                                                   
from pCore            import Clone, NormalDeviateGenerator, Pickle, RandomNumberGenerator, Selection, Unpickle
from pMolecule        import AtomSelection, MMModelOPLS, NBModelABFS, SoftConstraintContainer, SoftConstraintEnergyModelHarmonic, SoftConstraintMultipleTether
from pMoleculeScripts import ConjugateGradientMinimize_SystemGeometry, LangevinDynamics_SystemGeometry, SolvateSystemBySuperposition                                                                                                         

# . Parameters.
# . Refine options.
_Optimize = True
_Refine   = True
_Steps    = 2000

# . Define the MM and NB models.
nbModel = NBModelABFS ( )

# . Get the solute system.
solute = Unpickle ( os.path.join ( outPath, "step8_a.pkl" ) )
solute.Summary ( )

# . Get the solvent system.
solvent = Unpickle ( os.path.join ( outPath, "waterBox.pkl" ) )
solvent.DefineNBModel ( nbModel )
solvent.Summary ( )

# . Create the solvated system.
solution       = SolvateSystemBySuperposition ( solute, solvent, reorientSolute = False )
solution.label = "Solvated PKA Simulation System - Chains A and I Only"
solution.DefineNBModel  ( nbModel )
solution.Summary ( )

# . Refinement.
if _Refine:

    # . Check energies.
    solute.Energy ( )
    for system in ( solvent, solution ):
        system.DefineNBModel ( nbModel )
        system.Energy ( )

    # . Get selections for the protein atoms in the A and I chains.
    protein   = AtomSelection.FromAtomPattern ( solution, "A:*:*" ) | AtomSelection.FromAtomPattern ( solution, "I:*:*" )
    hydrogens = AtomSelection.FromAtomAttributes ( solution, "atomicNumber", 1 )
    heavies   = Selection.FromIterable ( protein - hydrogens )
    protein   = Selection.FromIterable ( protein )

    # . Fix all protein atoms.
    solution.DefineFixedAtoms ( protein )
    solution.Summary ( )

    # . Optimization.
    if _Optimize:
        ConjugateGradientMinimize_SystemGeometry ( solution                    ,
                                                   maximumIterations    =  200 ,
                                                   logFrequency         =   10 ,
                                                   rmsGradientTolerance = 10.0 )

    # . Define a normal deviate generator in a given state.
    normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 517152 ) )

    # . Dynamics loop - fixed atoms then tether constraints.
    for forceConstant in ( None, 500.0, 100.0, 20.0, 4.0 ):

        # . Harmonically constrain protein heavy atoms.
        tethers = None
        if ( forceConstant is not None ):
            reference          = Clone ( solution.coordinates3 )
            tetherEnergyModel  = SoftConstraintEnergyModelHarmonic ( 0.0, forceConstant )
            tethers            = SoftConstraintContainer ( )
            tethers["tethers"] = SoftConstraintMultipleTether ( heavies, reference, tetherEnergyModel )
        solution.DefineSoftConstraints ( tethers )

        # . Dynamics.
        LangevinDynamics_SystemGeometry ( solution                        ,
                                          collisionFrequency     =   25.0 ,
                                          logFrequency           =    100 ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  = _Steps ,
                                          temperature            =  300.0 ,
                                          timeStep               =  0.001 )

        # . Unfix atoms and remove tethers.
        solution.DefineFixedAtoms      ( None )
        solution.DefineSoftConstraints ( None )

# . Save the system.
Pickle ( os.path.join ( outPath, "step8_b.pkl" ), solution )

# . Print PDB file.
PDBFile_FromSystem ( os.path.join ( outPath, "step8_b.pdb" ), solution )
