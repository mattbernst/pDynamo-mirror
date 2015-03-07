#-------------------------------------------------------------------------------
# . File      : GeometryOptimization.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Helper functions for geometry optimizations."""

from pCore     import BakerOptimizer, ConjugateGradientMinimizer, FIREMinimizer, LBFGSMinimizer, QuasiNewtonMinimizer, SteepestDescentMinimizer, logFile
from pMolecule import SystemGeometryObjectiveFunction

#===================================================================================================================================
# . Helper function for setting up an optimization.
#===================================================================================================================================
def _SetUpOptimization ( system, defaultOptions, defaultsToChange, inputOptions, removeRotationTranslation ):
    """Generic function to set up an optimization."""
    # . Get default options - overridden with more suitable values if necessary.
    options = dict ( defaultOptions )
    options.update ( { "logFrequency" : 1, "maximumIterations" : 50, "rmsGradientTolerance" : 1.5 } )
    options.update ( defaultsToChange )
    # . Update with the input options.
    options.update ( inputOptions )
    # . Get some non-optimizer options.
    log                        = options.pop ( "log"                       , logFile )
    optimizeSymmetryParameters = options.pop ( "optimizeSymmetryParameters", False   )
    trajectories               = options.pop ( "trajectories"              , []      )
    # . Set up the objective function.
    of = SystemGeometryObjectiveFunction.FromSystem ( system )
    if removeRotationTranslation : of.RemoveRotationTranslation ( )
    if optimizeSymmetryParameters: of.IncludeSymmetryParameters ( )
    for ( trajectory, saveFrequency ) in trajectories:
        of.DefineTrajectory ( trajectory, saveFrequency )
    # . Finish up.
    return ( of, options, log )

#===================================================================================================================================
# . Baker saddle optimization.
#===================================================================================================================================
# . Keyword arguments are those from BakerOptimizer.defaultAttributes along with "findMinimum", "log", "optimizeSymmetryParameters" and "trajectories".
def BakerSaddleOptimize_SystemGeometry ( system, **keywordArguments ):
    """Baker stationary point geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, BakerOptimizer.defaultAttributes, {}, keywordArguments, True )
    # . Minimization.
    if options.pop ( "findMinimum", False ):
        options["followMode"           ] = -1
        options["numberOfNegativeModes"] =  0
        if "hessianUpdatingOption" not in keywordArguments: options["hessianUpdatingOption"] = "BFGS"
    # . Saddle point.
    else:
        options["followMode"           ] = max ( options["followMode"], 0 )
        options["numberOfNegativeModes"] = 1
        if "hessianUpdatingOption" not in keywordArguments: options["hessianUpdatingOption"] = "BOFILL"
    optimizer = BakerOptimizer ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Conjugate gradient minimization.
#===================================================================================================================================
# . Keyword arguments are those from ConjugateGradientMinimizer.defaultAttributes along with "log", "optimizeSymmetryParameters" and "trajectories".
# . It seems to be necessary to limit the initialStep here so as to prevent problems with SCF convergence when there are QC atoms.
def ConjugateGradientMinimize_SystemGeometry ( system, **keywordArguments ):
    """Conjugate gradient geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, ConjugateGradientMinimizer.defaultAttributes, { "initialStep" : 0.1 }, keywordArguments, True )
    optimizer = ConjugateGradientMinimizer ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . FIRE minimization.
#===================================================================================================================================
# . Keyword arguments are those from FIREMinimizer.defaultAttributes along with "log", "optimizeSymmetryParameters" and "trajectories".
def FIREMinimize_SystemGeometry ( system, **keywordArguments ):
    """FIRE geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, FIREMinimizer.defaultAttributes, { "maximumTimeStep" : 0.01, "timeStep" : 0.001 }, keywordArguments, True )
    optimizer = FIREMinimizer ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . LBFGS minimization.
#===================================================================================================================================
# . Keyword arguments are those from LBFGSMinimizer.defaultAttributes along with "log", "optimizeSymmetryParameters" and "trajectories".
def LBFGSMinimize_SystemGeometry ( system, **keywordArguments ):
    """LBFGS geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, LBFGSMinimizer.defaultAttributes, {}, keywordArguments, True )
    optimizer = LBFGSMinimizer ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Quasi-Newton minimization.
#===================================================================================================================================
# . Keyword arguments are those from QuasiNewtonMinimizer.defaultAttributes along with "hessian", "log", "optimizeSymmetryParameters" and "trajectories".
def QuasiNewtonMinimize_SystemGeometry ( system, **keywordArguments ):
    """Quasi-Newton geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, QuasiNewtonMinimizer.defaultAttributes, {}, keywordArguments, True )
    hessian   = options.pop ( "hessian", None )
    optimizer = QuasiNewtonMinimizer ( **options )
    optimizer.Summary ( log = log )
    if hessian is None:
        variables = of.VariablesAllocate ( )
        of.VariablesGet ( variables )
        hessian   = of.StartingHessian ( variables )
    of.SetStartingHessian ( hessian )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Steepest descent minimization.
#===================================================================================================================================
# . Keyword arguments are those from SteepestDescentMinimizer.defaultAttributes along with "log", "optimizeSymmetryParameters" and "trajectories".
def SteepestDescentMinimize_SystemGeometry ( system, **keywordArguments ):
    """Steepest descent geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, SteepestDescentMinimizer.defaultAttributes, {}, keywordArguments, False )
    optimizer = SteepestDescentMinimizer ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
