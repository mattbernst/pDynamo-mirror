"""Add counterions to the systems and solvate."""

from Definitions import *

# . PDB files and the number of cations to add.
_PDBPaths = (  ( "1UAO", 2 ), ( "2E4E", 1 ) )

# . Structures.
_Structures = ( "folded", "unfolded" )

# . Box sizes.
_XBox = 28.0
_YBox = 28.0
_ZBox = 28.0

# . Reorient option.
_Reorient = False

# . Define the solvent MM and NB models.
mmModel = MMModelOPLS ( "bookSmallExamples" )
nbModel = NBModelABFS ( )

# . Refine options.
_NSteps0  = 2500
_NSteps1  =  500
_Optimize = True
_Refine   = True

# . The type of cation to add.
_PositiveIon = "sodium"
cation = MOLFile_ToSystem ( os.path.join ( dataPath, _PositiveIon + ".mol" ) )
cation.DefineMMModel ( mmModel )
cation.Summary ( )

# . Get the solvent system.
solvent = Unpickle ( os.path.join ( outPath, "waterBox.pkl" ) )
solvent.DefineNBModel ( nbModel )
solvent.Summary ( )

# . Loop over structures.
for ( pdbPath, _NPositive ) in _PDBPaths:

    # . Retrieve the system.
    system = Unpickle ( os.path.join ( outPath, pdbPath + ".pkl" ) )
    system.Summary ( )
    masses = system.atoms.GetItemAttributes ( "mass" )

    # . Loop over xyz files.
    for structure in _Structures:

        # . Get the coordinates.
        system.coordinates3 = XYZFile_ToCoordinates3 ( os.path.join ( outPath, pdbPath + "_" + structure + ".xyz" ) )
        system.Energy ( )

        # . Reorient the system if necessary (see the results of GetSolvationInformation.py).
        if _Reorient: system.coordinates3.ToPrincipalAxes ( weights = masses )

        # . Add the counterions.
        solute = AddCounterIons ( system, 0, None, _NPositive, cation, ( _XBox, _YBox, _ZBox ) )
        solute.Summary ( )

        # . Create the solvated system.
        solution       = SolvateSystemBySuperposition ( solute, solvent, reorientSolute = False )
        solution.label = "Solvated " + pdbPath
        solution.DefineNBModel  ( nbModel )
        solution.Summary ( )

        # . Refinement.
        if _Refine:

            # . Check energies.
            solute.Energy ( )
            for temporary in ( solvent, solution ):
                temporary.DefineNBModel ( nbModel )
                temporary.Energy ( )

            # . Fix solute atoms - protein only.
            fixedAtoms = Selection.FromIterable ( AtomSelection.FromAtomPattern ( solution, "A:*:*" ) )
            solution.DefineFixedAtoms ( fixedAtoms )
            solution.Summary ( )

            # . Optimization.
            if _Optimize:
                ConjugateGradientMinimize_SystemGeometry ( solution                    ,
                                                           maximumIterations    =  200 ,
                                                           logFrequency         =    1 ,
                                                           rmsGradientTolerance = 10.0 )

            # . Define a normal deviate generator in a given state.
            normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 517151 ) )

            # . Do a Langevin dynamics calculation to refine the solvent and counterion coordinates.
            LangevinDynamics_SystemGeometry ( solution                          ,
                                              collisionFrequency     =     25.0 ,
                                              logFrequency           =      100 ,
                                              normalDeviateGenerator = normalDeviateGenerator ,
                                              steps                  = _NSteps0 ,
                                              temperature            =    300.0,
                                              timeStep               =    0.001 )

            # . Unfix atoms.
            solution.DefineFixedAtoms ( None )

            # . Do a Langevin dynamics calculation for all atoms.
            LangevinDynamics_SystemGeometry ( solution                          ,
                                              collisionFrequency     =     25.0 ,
                                              logFrequency           =      100 ,
                                              normalDeviateGenerator = normalDeviateGenerator ,
                                              steps                  = _NSteps1 ,
                                              temperature            =    300.0 ,
                                              timeStep               =    0.001 )

        # . Save the system.
        solution.configuration.Clear ( )
        Pickle ( os.path.join ( outPath, pdbPath + "_" + structure + "_solvated.pkl" ), solution )

        # . Save the PDB file for visualization.
        PDBFile_FromSystem ( os.path.join ( outPath, pdbPath + "_" + structure + "_solvated.pdb" ), solution )
