#-------------------------------------------------------------------------------
# . File      : CoordinateBuilding.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Helper functions for building coordinates.

These functions only use connectivity information. An energy model of a particular type is not required.
"""

import math

from pCore     import logFile, LogFileActive, RandomNumberGenerator, Selection, Vector3
from pMolecule import PeriodicTable, System

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Hydrogen building.
# . Default heavy atom - hydrogen atom bond distance (in Angstroms).
_DEFAULTBONDDISTANCE     = 1.0

# . Default angles for different coordinations (in degrees).
_LINEARANGLE             = 180.0
_TETRAHEDRALANGLE        = math.degrees ( 2.0 * math.asin ( math.sqrt ( 2.0 ) / math.sqrt ( 3.0 ) ) )
_TRIGONALPLANARANGLE     = 120.0

# . The default angles for a given coordination.
_COORDINATIONANGLES      = { 1: 0.0, 2: _LINEARANGLE, 3: _TRIGONALPLANARANGLE, 4: _TETRAHEDRALANGLE }

# . The default plane angle for a given coordination.
_COORDINATIONPLANEANGLES = { 3: 180.0, 4: 180.0 - 0.5 * _TETRAHEDRALANGLE }

# . Possible undefined coordinate value and tolerance.
#_UNDEFINEDCOORDINATE = 9999.0
#_UNDEFINEDTOLERANCE  = 0.1

#===================================================================================================================================
# . Build hydrogen coordinates using connectivity information only.
#===================================================================================================================================
def BuildHydrogenCoordinates3FromConnectivity ( system, log = logFile, randomNumberGenerator = None ):
    """Build hydrogen coordinates.

    The coordinates are built using connectivity information only (bonds
    and angles) which means that no account is taken of non-connectivity
    information (such as hydrogen-bonding). These interactions will have
    to be optimized separately using energy minimization or dynamics.

    Note that bonds to other hydrogens are ignored and unbound hydrogens
    or hydrogens linked to more than one heavy atom will not be built."""

    # . Check for a system object.
    if isinstance ( system, System ) and ( system.connectivity is not None ) and ( system.connectivity.HasFullConnectivity ( ) ):

        # . Check whether there are undefined coordinates.
        coordinates3     = system.coordinates3
        numberUndefined0 = coordinates3.numberUndefined
        if numberUndefined0 > 0:

            # . Initialization.
            bonds     = system.connectivity.bonds
            direction = Vector3.Null ( )
            if randomNumberGenerator is None:
                randomNumberGenerator = RandomNumberGenerator.WithRandomSeed ( )

            # . Loop over heavy atoms with defined coordinates.
            for ( c, atom ) in enumerate ( system.atoms ):
                if ( atom.atomicNumber != 1 ) and ( c not in coordinates3.undefined ):

                    # . Initialization.
                    built       = []
                    builth      = []
                    tobuild     = []
                    unbuildable = []

                    # . Loop over the connected atoms.
                    others = bonds.GetConnectedAtoms ( c )
                    for i in others:
                        other  = system.atoms[i]
                        QBUILT = ( i not in coordinates3.undefined )
                        # . Hydrogens.
                        if other.atomicNumber == 1:
                            if   QBUILT                                    : builth.append      ( i )
                            elif len ( bonds.GetConnectedAtoms ( i ) ) == 1: tobuild.append     ( i )
                            else:                                            unbuildable.append ( i )
                        # . Other atoms.
                        else:
                            if QBUILT: built.append       ( i )
                            else:      unbuildable.append ( i )

                    # . Skip this atom if the number of connections is greater than four, there are no hydrogens to build or there are unbuildable atoms.
                    if ( len ( others ) > 4 ) or ( len ( tobuild ) == 0 ) or ( len ( unbuildable ) > 0 ): continue

                    # . Order the lists and put built hydrogens after built heavy atoms as it is reasoned that heavy atom coordinates will be more reliable.
                    built.sort   ( )
                    builth.sort  ( )
                    tobuild.sort ( )
                    built += builth

                    # . Get coordination data for the center.
                    nconnections = len ( built ) + len ( tobuild )
                    bondlength   = PeriodicTable.Element ( atom.atomicNumber ).GetSingleBondDistance ( 1 )
                    if bondlength is None: bondlength = _DEFAULTBONDDISTANCE
                    angle        = PeriodicTable.Element ( atom.atomicNumber ).GetCoordinationAngle ( nconnections )
                    if angle is None: angle = _COORDINATIONANGLES.get      ( nconnections, 0.0 )
                    planeangle              = _COORDINATIONPLANEANGLES.get ( nconnections, 0.0 )

                    # . Build the hydrogens.
                    while len ( tobuild ) > 0:

                        # . Get the hydrogen index.
                        h = tobuild.pop  ( 0 )

                        # . Build according to the number of built connected atoms.
                        nbuilt = len ( built )

                        # . Get a random normalized vector.
                        if ( nbuilt == 0 ) or ( nbuilt == 1 ):
                            for i in range ( 3 ): direction[i] = 2.0 * ( randomNumberGenerator.NextReal ( ) - 0.5 )
                            direction.Normalize ( )

                        # . Put the hydrogen in a random direction from the center.
                        # . Works for all cases.
                        if   nbuilt == 0: coordinates3.BuildPointFromDistance ( h, c, bondlength, direction )

                        # . Put the hydrogen at the correct angle from the center and built atom but in a random plane.
                        # . Works for all cases given correct choice of angle.
                        elif nbuilt == 1:
                            coordinates3.BuildPointFromDistanceAngle ( h, c, built[0], bondlength, angle, direction )

                        # . Put the hydrogen away from the other built points at an appropriate angle from their plane.
                        # . The sign of the plane angle is arbitrary.
                        # . Works for cases 3, 4, 5 (square pyramidal), 6 with correct choice of planeangle.
                        elif nbuilt == 2:
                            coordinates3.BuildPointFromDistancePlaneAngle ( h, c, built[0], built[1], bondlength, planeangle )

                        # . Put the hydrogen using a tetrahedral tripod.
                        # . Only works for tetrahedral coordination.
                        elif nbuilt == 3:
                            coordinates3.BuildPointFromDistanceTetrahedralTripod ( h, c, built[0], built[1], built[2], bondlength )

                        # . Cannot handle valencies greater than 4 for the moment.
                        else: break

                        # . The hydrogen has been built.
                        built.append ( h )
                        coordinates3.FlagCoordinateAsDefined ( h )

            # . Output a summary.
            if LogFileActive ( log ):
                numberToBuild = coordinates3.numberUndefined
                numberBuilt   = ( numberUndefined0 - numberToBuild )
                if   numberBuilt <= 0: log.Paragraph ( "Coordinates for no hydrogens were built." )
                elif numberBuilt == 1: log.Paragraph ( "Coordinates for one hydrogen were built." )
                else:                  log.Paragraph ( "Coordinates for {:d} hydrogens were built.".format ( numberBuilt ) )

#===================================================================================================================================
# . Identify undefined coordinates.
#===================================================================================================================================
def IdentifyUndefinedCoordinates3 ( system, log = logFile, printHeavies = True, printHydrogens = False, usePDBNotation = True ):
    """Identify any undefined coordinates."""
    # . Initialization.
    numberOfHeavies   = 0
    numberOfHydrogens = 0

    # . See if there are undefined coordinates.
    if system.coordinates3.numberUndefined > 0:

        # . Get the indices of the undefined coordinates.
        undefined = system.coordinates3.undefined

        # . Find numbers.
        heavies   = []
        hydrogens = []
        for i in undefined:
            n = system.atoms[i].atomicNumber
            if n == 1: hydrogens.append ( i )
            else:      heavies.append   ( i )
        numberOfHeavies   = len ( heavies   )
        numberOfHydrogens = len ( hydrogens )

        # . Set some options.
        usePDBNotation = usePDBNotation and ( system.sequence is not None )

        # . Print the results.
        if LogFileActive ( log ):

            # . Basic summary.
            summary = log.GetSummary ( )
            summary.Start ( "Undefined Coordinates3 Summary" )
            summary.Entry ( "Undefined Heavy Atoms", "{:d}".format ( numberOfHeavies   ) )
            summary.Entry ( "Undefined Hydrogens",   "{:d}".format ( numberOfHydrogens ) )
            summary.Stop ( )

            # . Explicit listings.
            for ( doPrinting, indices, tag ) in ( ( printHeavies, heavies, "Heavy" ), ( printHydrogens, hydrogens, "Hydrogen" ) ):
                n = len ( indices )
                if doPrinting and ( n > 0 ):
                    if usePDBNotation: table = log.GetTable ( columns = min (  5, n ) * [ 20 ] )
                    else:              table = log.GetTable ( columns = min ( 10, n ) * [ 10 ] )
                    table.Start ( )
                    table.Title ( "Undefined " + tag + " Atoms" )
                    for i in indices:
                        table.Entry ( system.atoms[i].path )
                    table.Stop ( )
    # . Return.
    return ( numberOfHeavies, numberOfHydrogens )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
