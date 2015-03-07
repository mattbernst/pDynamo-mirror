#-------------------------------------------------------------------------------
# . File      : CrystalUtilities.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Utilities for dealing with crystals."""

import math

from pCore     import Clone, CrossPairList_FromDoubleCoordinates3, logFile, LogFileActive
from pMolecule import CrystalClassTriclinic, System

from Merge     import MergeByAtom

#===============================================================================
# . Analyze the transformations of a crystal.
#===============================================================================
def CrystalAnalyzeTransformations ( system, log = logFile ):
    """Analyze the transformations for a crystal.

    Transformations must be either proper or improper rotations and have inverses.

    An inverse check needs to be added.
    """

    # . Basic checks.
    if LogFileActive ( log ) and isinstance ( system, System ) and ( system.symmetry is not None ):

        # . Get the lattice matrices for the system.
        sp       = system.configuration.symmetryParameters
        M        = sp.M
        inverseM = sp.inverseM

        # . Loop over the transformations.
        for ( i, t3 ) in enumerate ( system.symmetry.transformations ):

            # . Output the transformation.
            newt3 = Clone ( t3 )
	    newt3.Orthogonalize ( M, inverseM )
            newt3.Print ( log = log, title = "Transformation {:d}".format ( i ) )

            # . Check for a rotation of some sort.
            if   newt3.rotation.IsProperRotation   ( ): log.Paragraph ( "Transformation is a proper rotation."    )
            elif newt3.rotation.IsImproperRotation ( ): log.Paragraph ( "Transformation is an improper rotation." )
            else: log.Paragraph ( "Transformation is neither a proper nor an improper rotation." )

#===============================================================================
# . Center the coordinates of a system by putting them inside the primary
# . image (i.e. their fractional coordinates are in the range [0,1]).
# . Centering can be done in three ways:
# 1. By atom.
# 2. By isolate.
# 3. By mmisolate but keeping all QC and fixed atom isolates as one.
# . The last method is the default as this is the one that will be needed
# . most often for simulations. 1 and 2 may be useful for analysis but will
# . often disrupt any energy calculations.
#
# . Centering of isolates is done with respect to the center of geometry of
# . the isolate.
#
# . If |selection| is present, only selected atoms or isolates with selected
# . atoms will be centered.
#
# . The return values of the function are the centered coordinates and the
# . translations effected (simply the difference between the new and old
# . coordinates).
#===============================================================================
_CENTERINGOPTIONS = [ "byatom", "byisolate", "bymmisolate" ]

def _GetMMIsolates ( system ):
    """Get the MM isolates for a system.

    Fixed and QC atoms are kept together within the same isolate so that their atoms' relative positions are maintained.
    """
    # . Initialization.
    mmisolates = None
    # . Basic checks.
    if isinstance ( system, System ) and ( system.energyModel is not None ) and ( system.energyModel.mmModel is not None ) and ( system.energyModel.exclusions is not None ):
        # . Get the basic set of isolates from the exclusions.
        mmisolates = system.energyModel.exclusions.GetIsolates ( upperBound = len ( system.atoms ) )
        # . Merge fixed atom and QC atom isolates.
        if system.hardConstraints     is not None: mmisolates.MergeIsolates ( system.hardConstraints.fixedAtoms               )
        if system.energyModel.qcAtoms is not None: mmisolates.MergeIsolates ( system.energyModel.qcAtoms.GetFullSelection ( ) )
    return mmisolates

def CrystalCenterCoordinates ( system, log = logFile, mode = "bymmisolate", selection = None ):
    """Center the coordinates of a system in the primary image."""

    # . Initialization.
    coordinates3 = None
    translations = None

    # . Basic checks.
    if isinstance ( system, System ) and ( system.symmetry is not None ) and ( system.coordinates3 is not None ) and ( system.symmetryParameters is not None ):

        # . Check the centering mode.
        option = mode.lower ( )
        if ( option not in _CENTERINGOPTIONS ): raise ValueError ( "Unknown centering mode: " + option + "." )

        # . Get a set of cloned coordinates.
        coordinates3 = Clone ( system.coordinates3 )

        # . By atom.
        if ( option == "byatom" ):
            system.symmetryParameters.CenterCoordinatesByAtom ( coordinates3, selection = selection )
        # . By isolate.
        elif ( option == "byisolate" ):
            if ( system.isolates is None ): raise ValueError ( "Isolates for the system are not defined." )
            else: system.symmetryParameters.CenterCoordinatesByIsolate ( coordinates3, system.isolates, selection = selection )
        # . By MM isolate.
        elif ( option == "bymmisolate" ):
            mmisolates = _GetMMIsolates ( system )
            system.symmetryParameters.CenterCoordinatesByIsolate ( coordinates3, mmisolates, selection = selection )

        # . Get the translations.
        translations = Clone ( coordinates3 )
        translations.AddScaledMatrix ( -1.0, system.coordinates3 )

    return ( coordinates3, translations )

#===============================================================================
# . Non-P1 to P1 with an arbitrary number of cells in each direction.
#===============================================================================
def CrystalExpandToP1 ( system, arange = range ( 1 ), brange = range ( 1 ), crange = range ( 1 ), imposeTriclinic = False, log = logFile ):
    """Create a system with P1 symmetry from one with arbitrary crystal symmetry.

    |system| is the input system.

    |arange| is the range of unit cells to build in the a direction.
    |brange| is the range of unit cells to build in the a direction.
    |crange| is the range of unit cells to build in the a direction.

    The return value is the new system but note that:
    |None|   is returned if |system| does not have crystal symmetry.
    |system| is returned if it is already P1.

    Although the system is converted to P1, the original crystal class is retained.

    Other options to include are:
    |QTRICLINIC| for selecting a triclinic crystal class for the output system.
    """
    result = None

    # . Get the number of cells.
    ncells = len ( arange ) * len ( brange ) * len ( crange )

    # . Basic checks.
    if isinstance ( system, System ) and ( system.symmetry is not None ):

        # . The system is already P1 and only one cell is required.
        if system.symmetry.transformations.IsIdentity ( ) and ( ncells == 1 ):
            result = system

        # . Build the new system.
        else:

            # . Get an alias for the system's symmetry parameters.
            sp = system.configuration.symmetryParameters

            # . Define the lattice matrix M and its inverse.
            M        = sp.M
            inverseM = sp.inverseM

            # . Create the images of the original system by building over each cell in turn.
            images = []
            for a in arange:
                for b in brange:
                    for c in crange:
                        for t3 in system.symmetry.transformations:
                            # . Create the transformation.
                            newt3 = Clone ( t3 )
                            newt3.translation[0] += a
                            newt3.translation[1] += b
                            newt3.translation[2] += c
	                    newt3.Orthogonalize ( M, inverseM )
                            # . Create the image.
                            image = Clone ( system )
                            image.coordinates3.Transform ( newt3 )
	                    images.append ( image )

            # . Get the new box lengths.
            a = len ( arange ) * sp.a
            b = len ( brange ) * sp.b
            c = len ( crange ) * sp.c

            # . Create the new system with the correct symmetry.
            if imposeTriclinic: crystalClass = CrystalClassTriclinic ( )
            else:               crystalClass = system.symmetry.crystalClass
            result = MergeByAtom ( images )
            result.DefineSymmetry ( crystalClass = crystalClass, a = a, b = b, c = c, alpha = sp.alpha, beta = sp.beta, gamma = sp.gamma )
            if system.label is None: result.label = "P1 Crystal System"
            else:                    result.label = system.label + " - P1 Crystal Symmetry"
            result.Summary ( log = log )

    return result

#===============================================================================
# . Get lists of possible bonds between the primary and secondary images.
# . It may be advisable to center the coordinates first (or use a large
# . range for the image search).
#===============================================================================
def CrystalGetImageBondPairs ( system, arange = range ( -1, 2 ), brange = range ( -1, 2 ), crange = range ( -1, 2 ), radii = None, safety = 0.45 ):
    """Create pairlists of interactions between the primary and secondary images."""

    # . Initialization.
    results = None

    # . Basic checks.
    if isinstance ( system, System ) and ( system.symmetry is not None ) and hasattr ( system.symmetry, "transformations" ):

        # . Need to have centering option!

        # . Get an alias for the system's symmetry parameters.
        sp = system.configuration.symmetryParameters

        # . Define the lattice matrix M and its inverse.
        M        = sp.M
        inverseM = sp.inverseM

        # . Create the images of the original system by building over each cell in turn.
        results = []
        for a in arange:
            for b in brange:
                for c in crange:
                    for ( t, t3 ) in enumerate ( system.symmetry.transformations ):
                        # . Create the transformation.
                        newt3 = Clone ( t3 )
                        newt3.translation[0] += a
                        newt3.translation[1] += b
                        newt3.translation[2] += c
	                newt3.Orthogonalize ( M, inverseM )
                        # . Skip the identity transformation.
                        if not newt3.IsIdentity ( ):
                            # . Create the image coordinates.
                            icoordinates3 = Clone ( system.coordinates3 )
                            icoordinates3.Transform ( newt3 )
                            # . Create the pair list.
                            pairlist = CrossPairList_FromDoubleCoordinates3 ( system.coordinates3, icoordinates3, radii1 = radii, radii2 = radii, safety = safety )
                            # . Calculate and save the pairs.
                            if ( pairlist is not None ) and ( len ( pairlist ) > 0 ):
                                pairs = []
                                for ( i, j ) in pairlist:
                                    dx = system.coordinates3[i,0] - icoordinates3[j,0]
                                    dy = system.coordinates3[i,1] - icoordinates3[j,1]
                                    dz = system.coordinates3[i,2] - icoordinates3[j,2]
                                    d  = math.sqrt ( dx * dx + dy * dy + dz * dz )
                                    pairs.append ( ( i, j, d ) )
                                results.append ( ( t, a, b, c, pairs ) )
    return results

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
