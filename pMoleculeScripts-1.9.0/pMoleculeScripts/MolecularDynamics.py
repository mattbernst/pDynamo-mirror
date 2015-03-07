#-------------------------------------------------------------------------------
# . File      : MolecularDynamics.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Helper functions for performing molecular dynamics simulations."""

from pCore     import logFile, LangevinVelocityVerletIntegrator, LeapFrogIntegrator, VelocityVerletIntegrator
from pMolecule import CrystalClassCubic, SystemGeometryObjectiveFunction

#===================================================================================================================================
# . Helper function for setting up a molecular dynamics simulation.
#===================================================================================================================================
def _SetUpSimulation ( system, defaultOptions, defaultsToChange, inputOptions ):
    """Generic function to set up an optimization."""
    # . Get default options - overridden with more suitable values if necessary.
    options = dict ( defaultOptions )
    options.update ( { "logFrequency" : 1, "steps" : 1000, "timeStep" : 0.001 } )
    options.update ( defaultsToChange )
    # . Update with the input options.
    options.update ( inputOptions )
    # . Change steps to maximumIterations.
    options["maximumIterations"] = options.pop ( "steps" )
    # . Get some non-optimizer options.
    log                       = options.pop ( "log"                       , logFile )
    removeRotationTranslation = options.pop ( "removeRotationTranslation" , True    )
    trajectories              = options.pop ( "trajectories"              , []      )
    # . Set up the objective function.
    of = SystemGeometryObjectiveFunction.FromSystem ( system )
    of.DefineWeights ( )
    if removeRotationTranslation: of.RemoveRotationTranslation ( )
    for ( trajectory, saveFrequency ) in trajectories:
        of.DefineTrajectory ( trajectory, saveFrequency )
    # . Finish up.
    return ( of, options, log )

#===================================================================================================================================
# . Langevin molecular dynamics.
#===================================================================================================================================
# . Keyword arguments are those from LangevinVelocityVerletIntegrator.defaultAttributes along with "log", "removeRotationTranslation" and "trajectories".
def LangevinDynamics_SystemGeometry ( system, **keywordArguments ):
    """Molecular dynamics using the velocity Verlet algorithm."""
    ( of, options, log ) = _SetUpSimulation ( system, LangevinVelocityVerletIntegrator.defaultAttributes, {}, keywordArguments )
    of.VelocitiesAssign ( options["temperature"], normalDeviateGenerator = options.get ( "normalDeviateGenerator", None ) )
    optimizer = LangevinVelocityVerletIntegrator ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Molecular dynamics using the leapfrog algorithm.
#===================================================================================================================================
# . Keyword arguments are those from LeapFrogIntegrator.defaultAttributes along with "log", "normalDeviateGenerator", "removeRotationTranslation" and "trajectories".
def LeapFrogDynamics_SystemGeometry ( system, **keywordArguments ):
    """Molecular dynamics using the leap-frog algorithm."""
    ( of, options, log ) = _SetUpSimulation ( system, LeapFrogIntegrator.defaultAttributes, {}, keywordArguments )
    # . The system has symmetry.
    if hasattr ( system, "symmetry" ) and ( system.symmetry is not None ):
        # . For the moment restrict calculation to systems with cubic symmetry.
        if options.get ( "pressureControl", False ) and not isinstance ( system.symmetry.crystalClass, CrystalClassCubic ):
            raise ValueError ( "Pressure coupling only works currently for systems with cubic symmetry." )
    # . Turn off pressure coupling for systems with no symmetry.
    else:
        options["pressureControl"] = False
    of.VelocitiesAssign ( options["temperature"], normalDeviateGenerator = options.pop ( "normalDeviateGenerator", None ) )
    optimizer = LeapFrogIntegrator ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Molecular dynamics using the velocity Verlet algorithm.
#===================================================================================================================================
# . Keyword arguments are those from VelocityVerletIntegrator.defaultAttributes along with "log", "normalDeviateGenerator", "removeRotationTranslation" and "trajectories".
def VelocityVerletDynamics_SystemGeometry ( system, **keywordArguments ):
    """Molecular dynamics using the velocity Verlet algorithm."""
    ( of, options, log ) = _SetUpSimulation ( system, VelocityVerletIntegrator.defaultAttributes, {}, keywordArguments )
    of.VelocitiesAssign ( options["temperatureStart"], normalDeviateGenerator = options.pop ( "normalDeviateGenerator", None ) )
    optimizer = VelocityVerletIntegrator ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
