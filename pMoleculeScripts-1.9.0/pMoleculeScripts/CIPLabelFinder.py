#-------------------------------------------------------------------------------
# . File      : CIPLabelFinder.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Script to find CIP labels for a molecule.

A concise, clear nice summary may be found in:

"Preferred IUPAC Names", Chapter 9, IUPAC (http://www.iupac.org)

The last CIP paper published by one of the eponymous authors:

V. Prelog and G. Helmchen
"Basic Principles of the CIP-System and Proposals for a Revision."
Angew. Chem. Int. Ed. Engl. 1982, 21, 567-583.

A paper indicating some of the complications of applying them:

P. Mata, A. M. Lobo, C. Marshall and A. P. Johnson
"Implementation of the Cahn-Ingold-Prelog System for Stereochemical
Perception in the LHASA Program."
J. Chem. Inf. Comp. Sci. 1994, 34, 491-504.

The rules for assignment of priority are complicated and only those
that make use of connectivity information are implemented here. This
means that priorities are determined using atomic numbers, atomic
masses and bond connectivities converted to hierarchical digraphs.

Omitted are complications involving multiple Kekule structures for
a system (e.g. fractional dummy atoms) and all priority comparisons
involving stereocenters. The latter include things such as the fact
that like descriptor pairs (R/R, S/S) take precedence over (R/S, S/R).

It would be good to try and differentiate between centers that are
clearly achiral (e.g. with > 1 bound Hs) and those that cannot be
determined.

Tetrahedral centers - lowest priority substituent at back with
                      a right-handed coordinate system:

   R (rectus)   - clockwise (highest to lowest).
   S (sinister) - anticlockwise.

Double bond centers - two highest priority substituents:

   E (Entgegen) - trans.
   Z (Zusammen) - cis.

Amino acids in proteins are S at their Calpha atoms, except for
cysteine which is R. Isoleucine/threonine - sidechains?

"""

import copy, math, operator

from pCore   import LogFileActive, logFile
from pMolecule import DoubleBond

#===============================================================================
# . Public functions.
#===============================================================================
def CIPLabelFinder ( molecule, log = logFile ):
    """Determine CIP labels for a molecule."""

    # . The molecule must have full connectivity.
    if molecule.connectivity.HasFullConnectivity ( ):

        # . Make the connection representation of the bonds.
        molecule.connectivity.bonds.MakeConnections ( )

        # . The molecule should be kekularized.
        #molecule.Kekularize ( )

        # . Treat the tetrahedral centers.
        ( tcenters, rtcenters, stcenters, utcenters ) = _TetrahedralCenters ( molecule )

        # . Treat the double bond centers.
        ( dcenters, edcenters, zdcenters, udcenters ) = _DoubleBondCenters  ( molecule )

        # . Output results.
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "Cahn-Ingold-Prelog Labels" )
            summary.Entry ( "Tetrahedral Centers"  , "{:d}".format ( len ( tcenters  ) ) )
            if len ( tcenters ) > 0:
               summary.Entry ( "Undefined Centers" , "{:d}".format ( len ( utcenters ) ) )
               summary.Entry ( "R Centers"         , "{:d}".format ( len ( rtcenters ) ) )
               summary.Entry ( "S Centers"         , "{:d}".format ( len ( stcenters ) ) )
            summary.Entry ( "Double Bond Centers"  , "{:d}".format ( len ( dcenters  ) ) )
            if len ( dcenters ) > 0:
               summary.Entry ( "Undefined Centers" , "{:d}".format ( len ( udcenters ) ) )
               summary.Entry ( "E Centers"         , "{:d}".format ( len ( edcenters ) ) )
               summary.Entry ( "Z Centers"         , "{:d}".format ( len ( zdcenters ) ) )
            summary.Stop ( )

        # . Finish up.
        return ( ( tcenters, rtcenters, stcenters, utcenters ), ( dcenters, edcenters, zdcenters, udcenters ) )

    # . Do nothing.
    else:
        return None

#===============================================================================
# . Private functions.
#===============================================================================
def _AssignCIPPriorities ( branches, highpriority, lowpriority ):
    """Order a set of atoms with respect to their CIP priorities."""
    # . Sort the branches in order of increasing priority.
    branches.sort ( key = operator.attrgetter ( "levelData" ) )

    # . Check for branches of lower and higher priority.
    while True:
        if len ( branches ) > 1:
            if branches[0].levelData < branches[1].levelData:
                branch = branches.pop ( 0 )
                lowpriority.insert ( 0, branch.level1 )
            else: break
        else: break
    while True:
        if len ( branches ) > 1:
            if branches[-1].levelData > branches[-2].levelData:
                branch = branches.pop ( -1 )
                highpriority.append ( branch.level1 )
            else: break
        else: break

    # . Check for an odd branch.
    if len ( branches ) == 1:
        branch = branches.pop ( )
        highpriority.append ( branch.level1 )

    # . There are no more branches to assign.
    if len ( branches ) == 0:
        return highpriority + lowpriority
    # . The next level is needed.
    else:
        # . Construct the next level for each branch.
        newbranches = []
        for branch in branches:
            branch.IncrementLevel ( )
            if len ( branch ) > 0: newbranches.append ( branch )
        if len ( newbranches ) > 0: return _AssignCIPPriorities ( newbranches, highpriority, lowpriority )
        else:                       return None

def _CIPPriorityOrderer ( molecule, rootatom, branchatoms ):
    """Order a set of atoms with respect to their CIP priorities."""
    # . Initialization.
    branches     = []
    highpriority = []
    lowpriority  = []
    # . Generate the branches.
    for level1 in branchatoms:
        branches.append ( _CIPBranch ( molecule.connectivity, rootatom, level1 ) )
     # . Assign priorities.
    return _AssignCIPPriorities ( branches, highpriority, lowpriority )

def _DoubleBondCenters ( molecule ):
    """Treat the double bond centers."""
    # . Find all double bonds between two sp2 carbons.
    dcenters  = []
    edcenters = []
    zdcenters = []
    udcenters = []
    for ( b, bond ) in enumerate ( molecule.connectivity.bonds ):
        if bond.type is DoubleBond ( ):
            i = bond.i
            j = bond.j
            atomi = molecule.atoms[i]
            atomj = molecule.atoms[j]
            if ( atomi.atomicNumber == 6 ) and ( atomi.connections == 3 ) and \
               ( atomj.atomicNumber == 6 ) and ( atomj.connections == 3 ):
               # . Save the bond.
               dcenters.append ( b )
               # . Get the connected atoms.
               iatoms = molecule.connectivity.bonds.GetConnectedAtoms ( i )
               iatoms.remove ( j )
               jatoms = molecule.connectivity.bonds.GetConnectedAtoms ( j )
               jatoms.remove ( i )
               # . Assign priorities.
               oiatoms = _CIPPriorityOrderer ( molecule, i, iatoms )
               ojatoms = _CIPPriorityOrderer ( molecule, j, jatoms )
               if ( oiatoms is not None ) and ( ojatoms is not None ):
                   label = _DoubleBondCenterLabel ( molecule.coordinates3, i, oiatoms[0], j, ojatoms[0] )
                   if label == "E": edcenters.append ( b )
                   else:            zdcenters.append ( b )
               # . Cannot assign priorities to the atoms.
               else: udcenters.append ( b )
    return ( dcenters, edcenters, zdcenters, udcenters )

def _DoubleBondCenterLabel ( coordinates3, d1, s1, d2, s2 ):
    """Determine the double bond label - E or Z.

    d1 and d2 are the double bond atoms and s1 and s2 are
    the highest priority substitutents on each atom.
    """
    v1 = coordinates3.Displacement ( s1, d1 )
    v2 = coordinates3.Displacement ( s2, d2 )
    if v1.Dot ( v2 ) < 0.0: return "E"
    else:                   return "Z"

def _SignedVolume ( a, b, c ):
    """Get the signed volume of three vectors."""
    sv = a[0] * ( b[1] * c[2] - b[2] * c[1] ) + \
         a[1] * ( b[2] * c[0] - b[0] * c[2] ) + \
         a[2] * ( b[0] * c[1] - b[1] * c[0] )
    return sv

def _TetrahedralCenters ( molecule ):
    """Treat the tetrahedral centers."""
    # . Find all tetrahedral carbons.
    tcenters  = []
    rtcenters = []
    stcenters = []
    utcenters = []
    for ( i, atom ) in enumerate ( molecule.atoms ):
        if ( atom.atomicNumber == 6 ) and ( atom.connections == 4 ):
            # . Save the atom.
            tcenters.append ( i )
            # . Order the connected atoms by their priorities.
            iatoms = molecule.connectivity.bonds.GetConnectedAtoms ( i )
            oatoms = _CIPPriorityOrderer ( molecule, i, iatoms )
            if oatoms is not None:
                label = _TetrahedralCenterLabel ( molecule.coordinates3, i, oatoms[0], oatoms[1], oatoms[2] )
                if label == "R": rtcenters.append ( i )
                else:            stcenters.append ( i )
            # . Cannot assign priorities to the atoms.
            else: utcenters.append ( i )
    return ( tcenters, rtcenters, stcenters, utcenters )

def _TetrahedralCenterLabel ( coordinates3, t, s1, s2, s3 ):
    """Determine the tetrahedral bond label - R or S.

    t is the center and s1, s2 and s3 are the highest
    priority substituents (1 > 2 > 3).

    Assigning a signed volume of a particular sign to "R" or "S"
    relies on the definition of the coordinate axes. Here the
    standard x (towards), y (right), z (up) system is used.
    """
    a = coordinates3.Displacement ( s1, t )
    b = coordinates3.Displacement ( s2, t )
    c = coordinates3.Displacement ( s3, t )
    sv = _SignedVolume ( a, b, c )
    if sv < 0.0: return "R"
    else:        return "S"

#===============================================================================
# . Class.
#===============================================================================
# . Comparison is done on the levelData.

class _CIPBranch ( object ):
    """Class to hold a hierarchical digraph originating from an atom bound to a stereocenter."""

    def __init__ ( self, connectivity, level0, level1 ):
        """Constructor.

        |connectivity| is the connectivity.
        |level0|       is the stereocenter.
        |level1|       is the atom from which the branch originates.
        """
        self.connectivity = connectivity
        self.level0       = level0
        self.level1       = level1
        self.paths        = [ [ level0, level1 ] ]
        self.levelData    = [ ( self.connectivity.atoms[level1].atomicNumber, self.connectivity.atoms[level1].mass ) ]

    def __len__ ( self ):
        """Return the number of paths in the branch."""
        return len ( self.paths )

    def IncrementLevel ( self ):
        """Go to the next level."""
        # . Get the new paths.
        oldpaths   = self.paths
        self.paths = []
        for oldpath in oldpaths:
            # . Get the top atom and its source.
            oldsource = oldpath[-2]
            oldtop    = oldpath[-1]
            # . A dummy atom terminates a path.
            if isinstance ( oldtop, _CIPDummyAtom ):
                pass
            # . A regular atom.
            else:
                # . Loop over the bonds for the atom.
                for b in self.connectivity.bonds.GetConnections ( oldtop ):
                    bond      = self.connectivity.bonds[b]
                    otheratom = bond.Other ( oldtop )
                    # . Multiple bonds add dummy atom paths no matter whether the bond is old or new.
                    for i in range ( bond.type.bondOrder - 1 ):
                        newpath = copy.copy ( oldpath )
                        newpath.append ( _CIPDummyAtom ( self.connectivity.atoms[otheratom].atomicNumber, self.connectivity.atoms[otheratom].mass ) )
                        self.paths.append ( newpath )
                    # . This is a new bond.
                    if otheratom != oldsource:
                        # . A ring adds a dummy atom path.
                        if otheratom in oldpath:
                            newpath = copy.copy ( oldpath )
                            newpath.append ( _CIPDummyAtom ( self.connectivity.atoms[otheratom].atomicNumber, self.connectivity.atoms[otheratom].mass ) )
                            self.paths.append ( newpath )
                        # . Add a normal path.
                        else:
                            newpath = copy.copy ( oldpath )
                            newpath.append ( otheratom )
                            self.paths.append ( newpath )
        # . Get the sorted atomic number and mass data for the level.
        self.levelData = []
        for path in self.paths:
            leveln = path[-1]
            if isinstance ( leveln, _CIPDummyAtom ): topatom = leveln
            else:                                    topatom = self.connectivity.atoms[leveln]
            self.levelData.append ( ( topatom.atomicNumber, topatom.mass ) )
        self.levelData.sort    ( )
        self.levelData.reverse ( )

class _CIPDummyAtom ( object ):
    """Class to represent a CIP dummy atom."""

    def __init__ ( self, atomicNumber, mass ):
        """Constructor."""
        self.atomicNumber = atomicNumber
        self.mass         = mass

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
