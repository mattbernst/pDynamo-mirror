#-------------------------------------------------------------------------------
# . File      : Dihedral.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes for handling dihedrals in a connectivity."""

import copy, operator

from pCore import Clone, logFile, LogFileActive, SingleObjectContainer

#===============================================================================
# . Class.
#===============================================================================
class Dihedral ( object ):
    """A class to represent a dihedral."""

    defaultattributes = { "i" : -1, \
                          "j" : -1, \
                          "k" : -1, \
                          "l" : -1  }

    def __init__ ( self, i = -1, j = -1, k = -1, l = -1 ):
        """Constructor."""
        self.i = i
        self.j = j
        self.k = k
        self.l = l

    def Order ( self ):
        """Ordering."""
        if self.j < self.k:
            self.i, self.l = self.l, self.i
            self.j, self.k = self.k, self.j

    @property
    def sortKey ( self ): return ( self.j, self.k, self.i, self.l )

#===============================================================================
# . Class.
#===============================================================================
class DihedralContainer ( SingleObjectContainer ):
    """A container class for dihedrals."""

    def __init__ ( self, *arguments ):
        """Constructor."""
        super ( DihedralContainer, self ).__init__ ( )
        self.QSORTED = False
        if len ( arguments ) > 0:
            # . Test for an iterable sequence.
            try:    items = iter ( arguments[0] )
            except: TypeError ( "Argument initializer must be an iterable sequence." )
            # . Create the dihedrals.
            dihedrals = []
            for ( i, j, k, l ) in items:
                dihedrals.append ( Dihedral ( i = i, j = j, k = k, l = l ) )
            self.items = dihedrals
            self.Sort ( )

    @classmethod
    def FromBondContainer ( selfClass, bonds ):
        """Constructor from a bond container."""
        bonds.MakeConnections ( )
        items = []
        for bond in bonds:
            j = bond.i
            k = bond.j
            jconnections = bonds.GetConnectedAtoms ( j )
            kconnections = bonds.GetConnectedAtoms ( k )
            for i in jconnections:
                if ( i != k ):
                    for l in kconnections:
                        if ( l != i ) and ( l != j ): items.append ( Dihedral ( i = i, j = j, k = k, l = l ) )
        self         = selfClass ( )
        self.items   = items
        self.QSORTED = True
        return self

    def ItemClass ( self ): return Dihedral

    def ItemName ( self ): return "Dihedral"

    def Sort ( self ):
        """Sorting."""
        if not self.QSORTED:
            items = self.items
            for item in items: item.Order ( )
            items.sort ( key = operator.attrgetter ( "sortKey" ) )
            for i in range ( len ( self ), 1, -1 ):
                if items[i-1].sortKey == items[i-2].sortKey: del items[i-1]
            self.QSORTED = True

    def UpperBound ( self ):
        """Return the upperbound of the container (the highest index + 1)."""
        self.Sort ( )
        if len ( self ) > 0: return max ( self[-1].i, self[-1].j, self[-1].l ) + 1
        else:                return 0

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
