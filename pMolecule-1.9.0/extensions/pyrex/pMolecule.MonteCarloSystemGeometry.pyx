#-------------------------------------------------------------------------------
# . File      : pMolecule.MonteCarloSystemGeometry.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Functions for the Monte Carlo simulation of molecular systems."""

from pCore.cDefinitions               cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3               cimport CCoordinates3, Coordinates3
from pCore.SelectionContainer         cimport SelectionContainer, CSelectionContainer
from pCore.Status                     cimport Status, Status_Success
from pMolecule.MonteCarloState        cimport MonteCarloState_AdjustMoveSizes, MonteCarloState_Allocate, MonteCarloState_MoveIsolate, MonteCarloState_MoveVolume,   \
                                              MonteCarloState_StatisticsBlockAccumulate, MonteCarloState_StatisticsBlockStart, MonteCarloState_StatisticsBlockStop, \
                                              MonteCarloState_StatisticsStart, MonteCarloState_StatisticsStop, CMonteCarloState
from pMolecule.NBModelMonteCarlo      cimport NBModelMonteCarlo, NBModelMonteCarlo_MMMMEnergyFull
from pMolecule.NBModelMonteCarloState cimport NBModelMonteCarloState
from pMolecule.SymmetryParameters     cimport SymmetryParameters, SymmetryParameters_Volume

from pCore       import CONSTANT_MOLAR_GAS, logFile, LogFileActive, RandomNumberGenerator, UNITS_PRESSURE_ATMOSPHERES_TO_KILOJOULES_PER_MOLE                     
from EnergyTerms import EnergyTerms                                                                               
from System      import System                                                                                    

#===============================================================================
# . Parameters.
#===============================================================================
# . The number of random numbers required per Monte Carlo move.
_NRANDOMISOLATE = 7
_NRANDOMVOLUME  = 2

# . The simulation options.
_SIMULATIONOPTIONS = { "acceptanceRatio"       :     0.4 ,
                       "adjustFrequency"       :    1000 ,
                       "blocks"                :      10 ,
                       "log"                   : logFile ,
                       "logFrequency"          :       1 ,
                       "moves"                 :   10000 ,
                       "pressure"              :     1.0 ,
                       "randomNumberGenerator" :    None ,
                       "rotation"              :    15.0 ,
                       "temperature"           :   300.0 ,
                       "trajectories"          :    None ,
                       "translation"           :    0.15 ,
                       "volumeChange"          :   400.0 ,
                       "volumeFrequency"       :     500 }

# . The simulation option names.
_SIMULATIONOPTIONNAMES = { "Acceptance Ratio"   : "acceptanceRatio" ,
                           "Adjust Frequency"   : "adjustFrequency" ,
                           "Blocks"             : "blocks"          ,
                           "Input RNG"          : "QRNG"            ,
                           "Log Frequency"      : "logFrequency"    ,
                           "Max. Rotation"      : "rotation"        ,
                           "Max. Translation"   : "translation"     ,
                           "Max. Volume Change" : "volumeChange"    ,
                           "Moves"              : "moves"           ,
                           "Pressure"           : "pressure"        ,
                           "Temperature"        : "temperature"     ,
                           "Trajectories"       : "QTRAJECTORIES"   ,
                           "Volume Frequency"   : "volumeFrequency" }

#===============================================================================
# . Private functions.
#===============================================================================
def _CheckSystem ( system, QCHECKSTATE = False ):
    """Check to see if a system is set up for a Monte Carlo calculation."""
    QOK = isinstance ( system, System )
    if QOK:
        # . Basic checks.
        QOK = hasattr ( system, "configuration" ) and hasattr ( system, "energyModel" ) and hasattr ( system.energyModel, "nbModel" ) and isinstance ( getattr ( system.energyModel, "nbModel" ), NBModelMonteCarlo )
        # . Check the state.
        if QCHECKSTATE and QOK:
            # . Set up the state if it does not exist.
            if ( not hasattr ( system.configuration, "nbState" ) ) or ( not isinstance ( getattr ( system.configuration, "nbState" ), NBModelMonteCarloState ) ): system.Energy ( log = None )
    if not QOK: raise ValueError ( "System not set up for Monte Carlo calculation." )

#===============================================================================
# . Calculate the interaction energy between one isolate and the remainder.
#===============================================================================
def MonteCarlo_IsolateInteractionEnergy ( system, isolate, log = None ):
    """Calculate the interaction energy between one isolate and all the others."""
    _CheckSystem ( system, QCHECKSTATE = True )
    # . Set up the system.
    system.configuration.nbState.Initialize ( system.configuration )
    # . Calculate the energy.
    energyTerms = EnergyTerms ( )
    energyTerms.Extend ( system.energyModel.nbModel.IsolateInteractionEnergy ( isolate, system.configuration.nbState ) )
    energyTerms.Total ( )
    # . Do some printing.
    energyTerms.Summary ( log = log )
    # . Finish up.
    return energyTerms.potentialEnergy

#===============================================================================
# . Scale the interaction parameters for an isolate.
#===============================================================================
def MonteCarlo_ScaleIsolateInteractionParameters ( system, isolate, chargeScale = 1.0, epsilonScale = 1.0, sigmaScale = 1.0, log = None ):
    """Scale the interaction parameters for an isolate."""
    _CheckSystem ( system, QCHECKSTATE = True )
    # . Set the scaling parameters in nbState.
    system.configuration.nbState.ScaleIsolateInteractionParameters ( isolate, chargeScale, epsilonScale, sigmaScale, log = log )

#===============================================================================
# . The simulation function.
#===============================================================================
def MonteCarlo_SystemGeometry ( system, **keywordArguments ):
    """Metropolis Monte Carlo simulation of a system."""

    # . Declarations.
    cdef Coordinates3           coordinates3
    cdef NBModelMonteCarlo      nbModel
    cdef NBModelMonteCarloState nbState
    cdef SelectionContainer     isolates
    cdef SymmetryParameters     symmetryParameters
    cdef Integer                    blocks, i, iblock, imove, moves, tmove
    cdef CMonteCarloState  *mcstate
    cdef Status            status

    # . Get the options.
    options  = dict ( _SIMULATIONOPTIONS )
    unknowns = set ( )
    for ( key, value ) in keywordArguments.items ( ):
        if key in options: options[key] = value
        else: unknowns.add ( key )
    if len ( unknowns ) > 0: raise ValueError ( "Invalid keyword arguments: " + ", ".join ( unknowns ) + "." )

    # . Get various counters.
    adjustfrequency = options["adjustFrequency"]
    blocks          = options["blocks"]
    moves           = options["moves" ]
    volumefrequency = options["volumeFrequency"]

    # . Check for printing.
    log          = options["log"]
    logFrequency = options["logFrequency"]
    QPRINT       = LogFileActive ( log ) and ( logFrequency > 0 ) and ( logFrequency <= blocks )
    if not QPRINT: options["logFrequency"] = 0

    # . Check for a random number generator.
    randomNumberGenerator = options["randomNumberGenerator"]
    options["QRNG"]       = ( randomNumberGenerator is not None )
    if randomNumberGenerator is None: randomNumberGenerator = RandomNumberGenerator.WithRandomSeed ( )

    # . Check for trajectories.
    logging      = []
    trajectories = options["trajectories"]
    if trajectories is not None:
        for ( trajectory, savefrequency ) in trajectories:
            if ( trajectory is not None ) and ( savefrequency > 0 ) and ( savefrequency <= blocks * moves ):
                logging.append ( ( savefrequency, trajectory ) )
                trajectory.WriteHeader ( pressure = options["pressure"], temperature = options["temperature"] )
    QLOGGING = ( len ( logging ) > 0 )
    options["QTRAJECTORIES"] = QLOGGING

    # . Print a summary.
    if QPRINT:
        # . Collect the options for printing.
        keys = _SIMULATIONOPTIONNAMES.keys ( )
        keys.sort ( )
        # . Do the summary.
        summary = log.GetSummary ( valueWidth = 16 )
        summary.Start ( "Monte Carlo Simulation Options" )
        for key in keys:
            value = options[_SIMULATIONOPTIONNAMES[key]]
            if   isinstance ( value, basestring ):
                valuestring = value
            elif isinstance ( value, bool ):
                if value: valuestring = "True"
                else:     valuestring = "False"
            elif isinstance ( value, float      ):
                valuestring = "{:.6g}".format ( value )
            else:
                valuestring = "{!r}".format ( value )
            summary.Entry ( key, valuestring )
        summary.Stop ( )

    # . Make sure the system is set up correctly.
    # . Basic checks.
    _CheckSystem ( system )
    if ( system.hardConstraints is not None ) and ( system.hardConstraints.NumberOfFixedAtoms ( ) > 0 ): raise ValueError ( "Monte Carlo simulations do not currently work with fixed atoms." )
    # . An energy (an easy way of making sure everything is OK).
    system.Energy ( log = None )

    # . Check that there are moves.
    if ( blocks * moves ) > 0:

        # . Set up the Monte Carlo state.
        # . Allocation.
        mcstate = MonteCarloState_Allocate ( len ( system.atoms ), max ( _NRANDOMISOLATE, _NRANDOMVOLUME ) )

        # . Options.
        mcstate.acceptanceratio = options["acceptanceRatio"]
        mcstate.blocks          = blocks
        mcstate.moves           = moves
        mcstate.rmax            = options["rotation"]
        mcstate.tmax            = options["translation"]
        mcstate.vmax            = options["volumeChange"]

        # . Aliases.
        coordinates3       = system.coordinates3          ; mcstate.coordinates3       = coordinates3.cObject
        nbModel            = system.energyModel.nbModel   ; mcstate.nbModel            = nbModel.cObject
        nbState            = system.configuration.nbState ; mcstate.nbState            = nbState.cObject
        isolates           = system.connectivity.isolates ; mcstate.isolates           = isolates.cObject
        symmetryParameters = system.symmetryParameters    ; mcstate.symmetryParameters = symmetryParameters.cObject

        # . Constants.
        mcstate.beta     = 1.0e+03 / ( CONSTANT_MOLAR_GAS * options["temperature"] )
        mcstate.pressure = UNITS_PRESSURE_ATMOSPHERES_TO_KILOJOULES_PER_MOLE * options["pressure"]
        mcstate.tfact    = float ( mcstate.isolates.nitems ) / mcstate.beta

        # . The initial energy and volume.
        mcstate.ecurrent = NBModelMonteCarlo_MMMMEnergyFull ( mcstate.nbModel, mcstate.nbState )
        mcstate.volume   = SymmetryParameters_Volume ( mcstate.symmetryParameters )

        # . Initialize statistics.
        MonteCarloState_StatisticsStart ( mcstate )

        if QPRINT:
            table = log.GetTable ( columns = [ 8, 16, 16, 16, 16, 16, 16, 16 ] )
            table.Start   ( )
            table.Title   ( "Monte Carlo Run Statistics" )
            table.Heading ( "Block"   )
            table.Heading ( "Accept." )
            table.Heading ( "<E>"     )
            table.Heading ( "<dE^2>"  )
            table.Heading ( "<H>"     )
            table.Heading ( "<dH^2>"  )
            table.Heading ( "<V>"     )
            table.Heading ( "<dV^2>"  )

        # . Initialize the total move counter.
        tmove = 0

        # . Loop over the blocks.
        for iblock from 0 <= iblock < blocks:

            # . Initialize the block statistics.
            MonteCarloState_StatisticsBlockStart ( mcstate )

            # . Loop over the moves per block.
            for imove from 0 <= imove < moves:

                # . Increment the total move counter.
                tmove = tmove + 1

                # . Check for a volume move.
                QVOLUME = ( volumefrequency > 0 ) and ( ( tmove % volumefrequency ) == 0 )

                # . Perform a move.
                if QVOLUME:
                    for i from 0 <= i < _NRANDOMVOLUME:
                        mcstate.random[i] = randomNumberGenerator.NextReal ( )
                    status = MonteCarloState_MoveVolume  ( mcstate )
                else:
                    for i from 0 <= i < _NRANDOMISOLATE:
                        mcstate.random[i] = randomNumberGenerator.NextReal ( )
                    status = MonteCarloState_MoveIsolate ( mcstate )
                if status != Status_Success: raise ValueError ( "Monte Carlo move error." )

                # . Accumulate statistics for the configuration.
                MonteCarloState_StatisticsBlockAccumulate ( mcstate )

                # . Adjust move sizes.
                if ( adjustfrequency > 0 ) and ( ( tmove % adjustfrequency ) == 0 ):
                    MonteCarloState_AdjustMoveSizes ( mcstate )

                # . Save the configuration.
                if QLOGGING:
                    for ( savefrequency, trajectory ) in logging:
                        if ( tmove % savefrequency ) == 0: trajectory.WriteOwnerData ( )

            # . Do the statistics for the block.
            MonteCarloState_StatisticsBlockStop ( mcstate )
            if QPRINT and ( ( iblock % logFrequency ) == 0 ):
                table.Entry ( "{:d}".format ( iblock ) )
                table.Entry ( "{:16.4f}".format ( float ( moves - mcstate.nreject ) / float ( moves ) ) )
                table.Entry ( "{:16.4f}".format ( mcstate.eav  ) )
                table.Entry ( "{:16.4f}".format ( mcstate.eav2 ) )
                table.Entry ( "{:16.4f}".format ( mcstate.hav  ) )
                table.Entry ( "{:16.4f}".format ( mcstate.hav2 ) )
                table.Entry ( "{:16.4f}".format ( mcstate.vav  ) )
                table.Entry ( "{:16.4f}".format ( mcstate.vav2 ) )

        # . Finish up.
        # . Do the statistics for the run.
        MonteCarloState_StatisticsStop ( mcstate )
        if QPRINT:
            if ( blocks > 1 ):
                # . Complete run.
                table.Entry ( "Run" )
                table.Entry ( "{:16.4f}".format ( float ( blocks * moves - mcstate.nrejectt ) / float ( blocks * moves ) ) )
                table.Entry ( "{:16.4f}".format ( mcstate.etot , ) )
                table.Entry ( "{:16.4f}".format ( mcstate.etot2 ) )
                table.Entry ( "{:16.4f}".format ( mcstate.htot  ) )
                table.Entry ( "{:16.4f}".format ( mcstate.htot2 ) )
                table.Entry ( "{:16.4f}".format ( mcstate.vtot  ) )
                table.Entry ( "{:16.4f}".format ( mcstate.vtot2 ) )
                # . Block statistics.
                table.Entry ( "Block" )
                table.Entry ( "-" )
                table.Entry ( "{:16.4f}".format ( mcstate.etotb  ) )
                table.Entry ( "{:16.4f}".format ( mcstate.etotb2 ) )
                table.Entry ( "{:16.4f}".format ( mcstate.htotb  ) )
                table.Entry ( "{:16.4f}".format ( mcstate.htotb2 ) )
                table.Entry ( "{:16.4f}".format ( mcstate.vtotb  ) )
                table.Entry ( "{:16.4f}".format ( mcstate.vtotb2 ) )
            table.Stop ( )

        # . Write out the final move sizes.
        if QPRINT and ( adjustfrequency > 0 ):
            summary = log.GetSummary ( valueWidth = 16 )
            summary.Start ( "Adjusted Move Sizes" )
            summary.Entry ( "Max. Rotation",    "{:16.6g}".format ( mcstate.rmax ) )
            summary.Entry ( "Max. Translation", "{:16.6g}".format ( mcstate.tmax ) )
            summary.Entry ( "Max. Vol. Move",   "{:16.6g}".format ( mcstate.vmax ) )
            summary.Stop ( )

    # . Deactivate logging.
    if QLOGGING:
        for ( savefrequency, trajectory ) in logging:
            trajectory.WriteFooter ( )
            trajectory.Close ( )
        QLOGGING = False
