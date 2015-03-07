#-------------------------------------------------------------------------------
# . File      : Aromaticity.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Functions for determining aromaticity."""

#===================================================================================================================================
#                           Aromaticity Determination
#===================================================================================================================================
#
# . The following rules are based upon those used in the SMILES definition.
#
# . An aromatic ring (or combination of rings) obeys the following rules:
#
# . 1. All atoms in the ring or ring system must be aromatic.
#
# . 2. The number of electrons donated to the ring system must obey Huckel's
#      rule - i.e. 4n + 2.
#
# . An aromatic atom obeys the following rules:
#
# . 1. It must be one of the possible aromatic elements (C, N, O, P, S, As, Se).
#
# . 2. It must be in a ring and have no more than 3 connections.
#
# . 3. The number of electrons an atom can donate depends upon the element type:
#
#        C: The third connection can be a single bond to any other atom (I),
#           although double bonds to non-carbon, non-ring atoms are allowed
#           too (as long as the total valence is 4) (II).
#
#        N, P, As: These atoms can be 3-coordinate and 3-valence (I) or
#                  2-coordinate or 3-coordinate and 4-valence in which case
#                  the double bond must be to an extra-cyclic oxygen (II).
#
#        O: Oxygen can only ever have two connections (I).
#
#        S, Se: These atoms can be 2-coordinate or 3-coordinate and
#               4-valence with the third connection a double bond to an
#               extra-cyclic oxygen (I).
#
#        Donated electrons and extra valencies:
#
#        Element            Case     Donated Electrons      Extra Valence
#
#        C                    I              1                    1
#                             II             0                    0
#        N, P, As             I              2                    0
#                             II             1                    1
#        O                    I              2                    0
#        S, Se                I              2                    0
#
#        The number of donated electrons is augmented or diminished by the
#        charge on the atom.
#
# . At the moment no attempt is made to aromaticize/kekularize the bonds but
# . this should be added eventually.
#
#===================================================================================================================================

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Elements which can be aromatic.
_AROMATICELEMENTS = ( 6, 7, 8, 15, 16, 33, 34 )

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def DetermineAromaticAtoms ( connectivity ):
    """Determine which atoms of a connectivity are aromatic."""

    # . This is overly simple and clearly needs to be made more robust.
    # . The minimum amount of information is used - particularly that concerning bond types.

    # . Initialization.
    atoms    = connectivity.atoms
    bonds    = connectivity.bonds
    ringSets = connectivity.ringSets

    # . Initialize the isAromatic flags for all atoms.
    for atom in atoms: atom.isAromatic = False

    # . Check that there are rings.
    if ( ringSets is not None ) and ( len ( ringSets ) > 0 ):

        # . Initialization.
        aromaticValencies = {} # . Not needed.

        # . Loop over each ringSet separately.
        for ringSet in ringSets:

            # . Get the ring atoms.
            ringatoms = set ( [ x for ring in ringSet for x in ring ] )

            # . Determine whether an atom can be aromatic and, if so, how many electrons it can donate.
            electrons = {}
            possible  = set ( )
            for atom in ringatoms:
                atomicNumber = atoms[atom].atomicNumber
                connections  = atoms[atom].connections
                # . Basic checks.
                if ( atomicNumber in _AROMATICELEMENTS ) and ( connections in [ 2, 3 ] ):
                    # . Initialization.
                    QEXTRACYCLICDOUBLE         = False
                    QEXTRACYCLICDOUBLETOOXYGEN = False
                    QEXTRACYCLICSINGLE         = False
                    QNOEXTRACYCLIC             = True
                    # . Get the (single) extracyclic atom bound to the atom if there is one.
                    for bindex in bonds.GetConnections ( atom ):
                        bond  = bonds[bindex]
                        other = bond.Other ( atom )
                        if other not in ringatoms:
                            QNOEXTRACYCLIC = False
                            if   bond.type.bondOrder == 1:
                                QEXTRACYCLICSINGLE         = True
                            elif bond.type.bondOrder == 2:
                                QEXTRACYCLICDOUBLE         = True
                                QEXTRACYCLICDOUBLETOOXYGEN = ( atoms[other].atomicNumber == 8 )
                            break
                    # . C.
                    if atomicNumber == 6:
                        if QEXTRACYCLICDOUBLE:
                            electrons[atom] = 0
                            possible.add ( atom )
                        elif ( QNOEXTRACYCLIC or QEXTRACYCLICSINGLE ):
                            electrons[atom] = 1
                            possible.add ( atom )
                    # . N, P, As.
                    elif atomicNumber in ( 7, 15, 33 ):
                        if ( connections == 2 ) or QEXTRACYCLICDOUBLETOOXYGEN:
                            electrons[atom] = 1
                            possible.add ( atom )
                        elif ( connections == 3 ) and ( QNOEXTRACYCLIC or QEXTRACYCLICSINGLE ):
                            electrons[atom] = 2
                            possible.add ( atom )
                    # . O, S, Se with 2 connections or S, Se with 3 connections in which the third is an extracyclic double bond to oxygen.
                    elif ( ( atomicNumber in ( 8, 16, 34 ) ) and ( connections == 2 ) ) or ( ( atomicNumber in ( 16, 34 ) ) and ( QEXTRACYCLICDOUBLETOOXYGEN ) ):
                        electrons[atom] = 2
                        possible.add ( atom )

            # . Modify the donated electrons by the atom's formal charge.
            for atom in electrons: electrons[atom] -= getattr ( atoms[atom], "formalCharge", 0 )

            # . This needs to be replaced by a loop over combinations of rings.
            # . Loop over each ring separately, check that each atom is a possible aromatic and then apply Huckel's rule.
            for ring in ringSet:
                nelectrons = 0
                QPOSSIBLE  = True
                for atom in ring:
                    if atom in possible:
                        nelectrons += electrons[atom]
                    else:
                        QPOSSIBLE = False
                        break
                # . Possible aromatic ring so apply Huckel's rule.
                if QPOSSIBLE:
                    n = ( nelectrons - 2 ) // 4
                    if ( 4 * n + 2 ) == nelectrons:
                        # . Flag the atoms as aromatic and set the aromatic valence.
                        # . The regular valence is unchanged as this will be altered due to bond changes.
                        for atom in ring:
                            atoms[atom].isAromatic = True
                            if electrons[atom] == 1: aromaticValencies[atom] = 1

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
