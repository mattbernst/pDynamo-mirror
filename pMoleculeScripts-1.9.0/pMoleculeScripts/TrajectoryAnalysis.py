#-------------------------------------------------------------------------------
# . File      : TrajectoryAnalysis.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Helper functions for analyzing trajectories."""

import math

from pCore     import Clone, Coordinates3, logFile, LogFileActive, Real1DArray, Real2DArray, Selection, SymmetricMatrix
from pMolecule import CrystalClassCubic

# . The trajectories to AveragePositions, CoordinateFluctuations and CovarianceMatrix should have had their rotation/translational motion removed.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Safety factor for estimating upper bound for RDF calculation.
_DEFAULTSAFETY = 1.0

#===================================================================================================================================
# . Helper functions for argument processing.
#===================================================================================================================================
def _CheckForUnknownOptions ( options ):
    """Check to see if there are unknown options."""
    unknowns = sorted ( options.keys ( ) )
    if len ( unknowns ) > 0: raise ValueError ( "Invalid keyword arguments: " + ", ".join ( unknowns ) + "." )

def _GetSystemAndTrajectoriesFromArguments ( *arguments ):
    """Get a system and its trajectories from an argument list."""
    if isinstance ( arguments[0], ( list, tuple ) ): trajectories =   arguments[0]
    else:                                            trajectories = [ arguments[0] ]
    system = trajectories[0].owner
    return ( system, trajectories )

#===================================================================================================================================
# . Calculate average positions.
#===================================================================================================================================
#AveragePositions ( trajectory/trajectories, selection = None )
def AveragePositions ( *arguments, **keywordArguments ):
    """Calculate the average positions for selected particles."""

    # . Initialization.
    positions = None

    # . Process arguments.
    ( system, trajectories ) = _GetSystemAndTrajectoriesFromArguments ( *arguments )

    # . Get the selection (or all particles otherwise).
    selection = keywordArguments.pop ( "selection", None )
    if selection is None: selection = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    _CheckForUnknownOptions ( keywordArguments )

    # . Get the size of the problem.
    n = len ( selection )
    if n > 0:

        # . Allocate space.
        positions = Coordinates3.WithExtent ( n )
        positions.Set ( 0.0 )

        # . Loop over trajectory frames.
        numberFrames = 0
        for trajectory in trajectories:
            trajectory.ReadHeader ( )
            while trajectory.RestoreOwnerData ( ):
	        frame = system.coordinates3
                for ( p, i ) in enumerate ( selection ):
                    positions[p,0] += frame[i,0]
                    positions[p,1] += frame[i,1]
                    positions[p,2] += frame[i,2]
            trajectory.ReadFooter ( )
            trajectory.Close ( )
            numberFrames += len ( trajectory )

        # . Scale.
        if numberFrames > 0: positions.Scale ( 1.0 / float ( numberFrames ) )

    return positions

#===================================================================================================================================
# . Coordinate fluctuations.
#===================================================================================================================================
#CoordinateFluctuations ( trajectory/trajectories, anisotropic = False, asBFactors = False, selection = None )
def CoordinateFluctuations ( *arguments, **keywordArguments ):
    """Calculate the coordinate fluctuations for selected particles."""

    # . Initialization.
    fluctuations = None

    # . Get the trajectory and associated system data.
    ( system, trajectories ) = _GetSystemAndTrajectoriesFromArguments ( *arguments )

    # . Get the selection and the size of the problem.
    selection = keywordArguments.pop ( "selection", None )
    if selection is None: selection = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    n         = len ( selection )

    # . Get the average positions.
    averagePositions = keywordArguments.pop ( "averagePositions", None )
    if averagePositions is None:
        averagePositions = AveragePositions ( trajectories, selection = selection )

    # . Various other options.
    anisotropic = keywordArguments.pop ( "anisotropic", False )
    asBFactors  = keywordArguments.pop ( "asBFactors" , False )
    _CheckForUnknownOptions ( keywordArguments )

    # . Continue processing.
    if ( n > 0 ) and ( averagePositions is not None ):

        # . Allocate space.
        displacement = Coordinates3.WithExtent ( n )
        if anisotropic: fluctuations = Real2DArray.WithExtents ( n, 6 )
        else:           fluctuations = Real1DArray.WithExtent  ( n    )
        displacement.Set ( 0.0 )
        fluctuations.Set ( 0.0 )

        # . Loop over trajectory frames.
        numberFrames = 0
        for trajectory in trajectories:
            trajectory.ReadHeader ( )
            while trajectory.RestoreOwnerData ( ):
	        frame = system.coordinates3
                if anisotropic:
                    for ( p, i ) in enumerate ( selection ):
                        dx = frame[i,0] - averagePositions[p,0]
                        dy = frame[i,1] - averagePositions[p,1]
                        dz = frame[i,2] - averagePositions[p,2]
                        fluctuations[p,0] += dx * dx
                        fluctuations[p,1] += dy * dx
                        fluctuations[p,2] += dy * dy
                        fluctuations[p,3] += dz * dx
                        fluctuations[p,4] += dz * dy
                        fluctuations[p,5] += dz * dz
                else:
                    for ( p, i ) in enumerate ( selection ):
                        fluctuations[p] += ( ( frame[i,0] - averagePositions[p,0] )**2 + \
                                             ( frame[i,1] - averagePositions[p,1] )**2 + \
                                             ( frame[i,2] - averagePositions[p,2] )**2 )
            trajectory.ReadFooter ( )
            trajectory.Close ( )
            numberFrames += len ( trajectory )

        # . Scale.
        if numberFrames > 0: fluctuations.Scale ( 1.0 / float ( numberFrames ) )

        # . Convert to B-factors if necessary.
        if asBFactors:
            conversionFactor = 8.0 * math.pi**2
            if not anisotropic: conversionFactor /= 3.0
            fluctuations.Scale ( conversionFactor )

    # . Finish up.
    return fluctuations

#===================================================================================================================================
# . Calculate the covariance matrix.
#===================================================================================================================================
#CovarianceMatrix ( trajectory/trajectories, selection = None )
def CovarianceMatrix ( *arguments, **keywordArguments ):
    """Calculate the covariance matrix for selected particles."""

    # . Initialization.
    covariance = None

    # . Get the trajectory and associated system data.
    ( system, trajectories ) = _GetSystemAndTrajectoriesFromArguments ( *arguments )

    # . Get the selection and the size of the problem.
    selection = keywordArguments.pop ( "selection", None )
    if selection is None: selection = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    n         = 3 * len ( selection )

    # . Get the average positions.
    averagePositions = keywordArguments.pop ( "averagePositions", None )
    if averagePositions is None:
        averagePositions = AveragePositions ( trajectories, selection = selection )
    _CheckForUnknownOptions ( keywordArguments )

    # . Continue processing.
    if ( n > 0 ) and ( averagePositions is not None ):

        # . Allocate space.
        covariance   = SymmetricMatrix.WithExtent ( n ) ; covariance.Set   ( 0.0 )
        displacement = Real1DArray.WithExtent     ( n ) ; displacement.Set ( 0.0 )

        # . Loop over trajectory frames.
        numberFrames = 0
        for trajectory in trajectories:
            trajectory.ReadHeader ( )
            while trajectory.RestoreOwnerData ( ):
	        frame = system.coordinates3
                for ( p, i ) in enumerate ( selection ):
                    displacement[3*p  ] = frame[i,0] - averagePositions[p,0]
                    displacement[3*p+1] = frame[i,1] - averagePositions[p,1]
                    displacement[3*p+2] = frame[i,2] - averagePositions[p,2]
                for i in range ( n ):
                    dI = displacement[i]
                    for j in range ( i + 1 ):
                        covariance[i,j] += ( dI * displacement[j] )
            trajectory.ReadFooter ( )
            trajectory.Close ( )
            numberFrames += len ( trajectory )

        # . Scale.
        if numberFrames > 0: covariance.Scale ( 1.0 / float ( numberFrames ) )

    return covariance

#===================================================================================================================================
# . Calculate a radial distribution function.
#===================================================================================================================================
#RadialDistributionFunction ( trajectory/trajectories, bins = 100, log = logFile, maximumR = None, selection1 = None, selection2 = None )
def RadialDistributionFunction ( *arguments, **keywordArguments ):
    """Calculate a radial distribution function."""

    # . Get the trajectory and associated system data.
    ( system, trajectories ) = _GetSystemAndTrajectoriesFromArguments ( *arguments )

    # . Keyword arguments.
    bins       = keywordArguments.pop ( "bins"       ,  100           )
    log        = keywordArguments.pop ( "log"        ,  logFile       )
    safety     = keywordArguments.pop ( "safety"     , _DEFAULTSAFETY )
    selection1 = keywordArguments.pop ( "selection1" , None           )
    selection2 = keywordArguments.pop ( "selection2" , None           )
    upper      = keywordArguments.pop ( "maximumR"   , None           )
    if selection1 is None: selection1 = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    if selection2 is None: selection2 = selection1
    _CheckForUnknownOptions ( keywordArguments )

    # . Get the numbers of particles.
    np1 = len ( selection1 )
    np2 = len ( selection2 )

    # . Check for a cubic system.
    try:
        cc = system.symmetry.crystalClass
        if not isinstance ( cc, CrystalClassCubic ): raise
    except:
        raise ValueError ( "System does not have cubic symmetry." )

    # . Estimate an upper bound.
    if upper is None:
        try:    upper = ( system.symmetryParameters.a - safety ) / 2.0
        except: raise ValueError ( "Please supply an upper bound for the RDF calculation." )

    # . Initialization.
    lower     = 0.0
    distances = []
    histogram = [ 0 for i in range ( bins ) ]
    rdf       = []

    # . Calculate the width of each bin.
    width = ( upper - lower ) / float ( bins )
    vacc  = 0.0

    # . Loop over trajectory frames.
    numberFrames = 0
    for trajectory in trajectories:
        trajectory.ReadHeader ( )
        while trajectory.RestoreOwnerData ( ):

            # . Get the frame.
	    frame         = system.coordinates3
            numberFrames += 1

            # . Accumulate the volume.
            a     = system.symmetryParameters.a
            vacc += a**3

            # . Check a.
            if upper > 0.5 * a: raise ValueError ( "Box does not satisfy minimum image convention." )

            # . Bin the distances.
	    # . This loop involves duplicate work for self rdfs.
            for i in selection1:
	        for j in selection2:
		    if i != j:
                        dr = frame.Displacement ( i, j )
                        for c in range ( 3 ):
                            x      = a * round ( dr[c] / a )
                            dr[c] -= x
                        r = dr.Norm2 ( )
                        if ( r >= lower ) and ( r < upper ):
		            b = int ( ( r - lower ) / width )
			    histogram[b] += 1

        # . Finish up.
        trajectory.ReadFooter ( )
        trajectory.Close ( )

    # . Calculate g(r).
    fact = 4.0 * math.pi * float ( np1 * np2 * numberFrames**2 ) / ( 3.0 * vacc )
    for i in range ( bins ):
        rlower = lower  + float ( i ) * width
        rupper = rlower + width
	distances.append ( rlower + 0.5 * width )
        rdf.append ( float ( histogram[i] ) / ( fact * ( rupper**3 - rlower**3 ) ) )

    # . Output the results.
    if LogFileActive ( log ):
        table = log.GetTable ( columns = [ 20, 20 ] )
        table.Start ( )
        table.Title ( "Radial Distribution Function" )
        table.Heading ( "Distance" )
	table.Heading ( "g(r)"     )
        for ( r, g ) in zip ( distances, rdf ):
	    table.Entry ( "{:20.4f}".format ( r ) )
	    table.Entry ( "{:20.4f}".format ( g ) )
        table.Stop ( )

    # . Return results.
    return ( distances, rdf )

#===================================================================================================================================
# . Remove rotational and translational motion from a trajectory.
#===================================================================================================================================
def RemoveRotationTranslation ( inTrajectory, outTrajectory, system, reference3 = None, useWeights = True ):
    """Remove rotational and translational motion from a trajectory."""

    # . Reference structure.
    # . Use the first trajectory structure if none specified.
    if reference3 is system.coordinates3: reference3 = Clone ( system.coordinates3 )

    # . Get weights.
    if useWeights: masses = system.atoms.GetItemAttributes ( "mass" )
    else:          masses = None

    # . Loop over trajectory frames.
    inTrajectory.ReadHeader   ( )
    outTrajectory.WriteHeader ( )
    while inTrajectory.RestoreOwnerData ( ):
        if reference3 is None:
            reference3 = Clone ( system.coordinates3 )
        else:
            system.coordinates3.Superimpose ( reference3, weights = masses )
        outTrajectory.WriteOwnerData ( )
    inTrajectory.ReadFooter   ( )
    outTrajectory.WriteFooter ( )
    inTrajectory.Close  ( )
    outTrajectory.Close ( )

#===================================================================================================================================
# . Calculate the self diffusion function.
#===================================================================================================================================
#SelfDiffusionFunction ( trajectory/trajectories, log = logFile, maximumTime = None, selection = None )
def SelfDiffusionFunction ( *arguments, **keywordArguments ):
    """Calculate the self diffusion function."""

    # . Get the trajectory and associated system data.
    ( system, trajectories ) = _GetSystemAndTrajectoriesFromArguments ( *arguments )

    # . Keyword arguments.
    log         = keywordArguments.pop ( "log"         , logFile )
    maximumTime = keywordArguments.pop ( "maximumTime" , None    )
    selection   = keywordArguments.pop ( "selection"   , None    )
    timeStep    = keywordArguments.pop ( "timeStep"    , 1.0     )
    if selection is None: selection = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    _CheckForUnknownOptions ( keywordArguments )

    # . Initialization.
    np     = len ( selection ) # . The number of particles.
    dSelf  = []
    frames = []
    times  = []

    # . Loop over trajectory frames.
    for trajectory in trajectories:
        trajectory.ReadHeader ( )
        while trajectory.RestoreOwnerData ( ):
            frames.append ( system.coordinates3.Prune ( selection ) )
        trajectory.ReadFooter ( )
        trajectory.Close ( )

    # . Get the number of frames and tStop.
    numberFrames = len ( frames )
    if maximumTime is None:
        tStop = numberFrames - 1
    else:
        tStop = int ( maximumTime / timeStep )
        tStop = min ( numberFrames - 1, tStop )

    # . Calculate the function - slow version.
    # . Loop over time differences to calculate.
    for t in range ( tStop + 1 ):

        # . Initialization.
        tmax  = numberFrames - t
	total = 0.0

        # . Loop over allowed increments.
	for t0 in range ( tmax ):
	    # . Calculate differences.
	    fa = frames[t0+t]
	    fb = frames[t0]
	    for i in range ( np ):
	        total += ( fa[i,0] - fb[i,0] )**2 + ( fa[i,1] - fb[i,1] )**2 + ( fa[i,2] - fb[i,2] )**2

        # . Calculate dSelf.
	dSelf.append ( total / float ( 3 * np * tmax ) )
        times.append ( float ( t ) * timeStep )

    # . Output the results.
    if LogFileActive ( log ):
        table = log.GetTable ( columns = [ 20, 20 ] )
        table.Start ( )
        table.Title ( "Self-Diffusion Function" )
        table.Heading ( "Time"  )
	table.Heading ( "Dself" )
        for ( t, d ) in zip ( times, dSelf ):
	    table.Entry ( "{:20.4f}".format ( t ) )
	    table.Entry ( "{:20.4f}".format ( d ) )
        table.Stop ( )

    # . Return results.
    return ( times, dSelf )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
