#-------------------------------------------------------------------------------
# . File      : GrowingStringPath.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Path generation using the growing string method."""

from pBabel    import SystemGeometryTrajectory
from pCore     import Clone, LBFGSMinimizer, logFile, LogFileActive, Real1DArray
from pMolecule import SystemGeometryObjectiveFunction

#===================================================================================================================================
# . Notes:
#
# . The number of images is nimages and are numbered from 0 to nimages-1. The first (0) and last images (nimages-1) are fixed so
#   only those in the range [1,nimages-2] move.
#
# . The overall characteristics of the optimization seem to be similar to that of the original although some of the rms gradient
#   tolerance and stepSize formulae may need playing with. The module, however, should provide a foundation for future work.
#
#===================================================================================================================================

# . Errors.
class ArgumentError ( Exception ): pass

#===================================================================================================================================
# . "Growing string" method for generating an initial path.
#===================================================================================================================================
# . Basic minimization options.
_MinimizerOptions = { "logFrequency"         :   1 ,
                      "maximumIterations"    :  50 ,
                      "rmsGradientTolerance" : 1.5 }

# . Change path to trajectory object (as for all other scripts)?

def GrowingStringInitialPath ( system, npoints, point0, pointn, path, **keywordArguments ):
    """Use a \"growing string\" method for generating an initial path."""

    # . Get some options.
    log      = keywordArguments.pop ( "log", logFile )
    sequence = keywordArguments.pop ( "sequence", "Forward" ).lower ( )

    # . Get the number of points.
    npoints = max ( npoints, 2 )

    # . Check for fixed atoms.
    QFIXED = ( system.hardConstraints is not None ) and ( system.hardConstraints.NumberOfFixedAtoms ( ) > 0 )

    # . Save the current system coordinates.
    temporary = Clone ( system.coordinates3 )

    # . Create the trajectory and write the first frame.
    trajectory = SystemGeometryTrajectory ( path, system, mode = "w" )
    trajectory.WriteFrame ( point0, frame = 0 )

    # . Get the first and last energies.
    energies = {}
    point0.CopyTo ( system.coordinates3 )
    energies[0] = system.Energy ( log = log )
    pointn.CopyTo ( system.coordinates3 )
    if not QFIXED: system.coordinates3.Superimpose ( point0 )
    energies[npoints-1] = system.Energy ( log = log )

    # . Save the last point.
    trajectory.WriteFrame ( system.coordinates3, frame = npoints - 1 )

    # . Generate intermediate points.
    if npoints > 2:

        # . Set up the optimizer.
        algorithm = keywordArguments.pop ( "optimizer", LBFGSMinimizer )
        options   = dict ( _MinimizerOptions )
        options.update ( keywordArguments )
        optimizer = algorithm ( **options )
        optimizer.Summary ( log = log )

        # . Set up the objective function.
        of = SystemGeometryObjectiveFunction.FromSystem ( system )

        # . Get variables for the first and last points and the tangent constraint.
        first   = Real1DArray.WithExtent ( of.nvariables )
        last    = Real1DArray.WithExtent ( of.nvariables )
        tangent = Real1DArray.WithExtent ( of.nvariables )

        # . Save the last point (given it is already in system.coordinates3).
        of.VariablesGet ( last )

        # . Save the first point.
        point0.CopyTo ( system.coordinates3 )
        of.VariablesGet ( first )

        # . Generate the sequence over which iterations are to occur.
        sequenceData = []
        if   sequence == "alternate" :
            for i in range ( ( npoints - 2 ) // 2 ):
                sequenceData.append ( ( first, last, i + 1, npoints - 2 - ( 2 * i ) ) )
                sequenceData.append ( ( last, first, npoints - 2 - i, npoints - 2 - ( 2 * i + 1 ) ) )
            if ( npoints - 2 ) % 2 != 0: sequenceData.append ( ( first, last, ( npoints - 2 ) // 2 + 1, 1 ) )
        elif sequence == "backward"  :
            for i in range ( npoints - 2 ): sequenceData.append ( ( last, first, npoints - 2 - i, npoints - 2 - i ) )
        elif sequence == "forward"   :
            for i in range ( npoints - 2 ): sequenceData.append ( ( first, last, i + 1, npoints - 2 - i ) )
        else: raise ArgumentError ( "Unknown sequence option - " + sequence + "." )

        # . Loop over points.
        for ( start, stop, index, remainder ) in sequenceData:

            # . Get the new point to optimize and the tangent constraint.
            stop.CopyTo ( tangent )
            tangent.AddScaledArray ( -1.0, start )
            start.AddScaledArray ( 1.0 / float ( remainder + 1 ), tangent )
            tangent.Normalize ( )

            # . Set up the object function.
            of.RemoveRotationTranslation ( reference = point0 )
            of.AddLinearConstraint ( tangent )
            of.VariablesPut ( start )

            # . Minimize and save the point.
            optimizer.Iterate ( of, log = log )
            of.VariablesGet ( start )
            energies[index] = system.Energy ( log = log )
            trajectory.WriteFrame ( system.coordinates3, frame = index )

    # . Finish up.
    system.coordinates3 = temporary
    trajectory.Close ( )

    # . Output.
    if LogFileActive ( log ):
        table = log.GetTable ( columns = [ 30, 20 ] )
        table.Start   ( )
        table.Title   ( "Growing String Structure Energies" )
        table.Heading ( "Structure" )
        table.Heading ( "Energy"    )
        keys = energies.keys ( )
        keys.sort ( )
        for key in keys:
            table.Entry ( repr ( key ) )
            table.Entry ( "{:.3f}".format ( energies[key] ) )
        table.Stop ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
