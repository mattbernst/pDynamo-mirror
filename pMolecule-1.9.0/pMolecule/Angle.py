#-------------------------------------------------------------------------------
# . File      : Angle.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes for handling angles in a connectivity."""

import copy, operator

from pCore import Clone, logFile, LogFileActive, SingleObjectContainer

#===============================================================================
# . Class.
#===============================================================================
class Angle ( object ):
    """A class to represent an angle."""

    def __init__ ( self, i = -1, j = -1, k = -1 ):
        """Constructor."""
        self.i = i
        self.j = j
        self.k = k

    def Order ( self ):
        """Ordering."""
        if self.i < self.k: self.i, self.k = self.k, self.i

    @property
    def sortKey ( self ): return ( self.j, self.i, self.k )

#===============================================================================
# . Class.
#===============================================================================
class AngleContainer ( SingleObjectContainer ):
    """A container class for angles."""

    def __init__ ( self, *arguments ):
        """Constructor."""
        super ( AngleContainer, self ).__init__ ( )
        self.QSORTED = False
        if len ( arguments ) > 0:
            # . Test for an iterable sequence.
            try:    items = iter ( arguments[0] )
            except: TypeError ( "Argument initializer must be an iterable sequence." )
            # . Create the angles.
            angles = []
            for ( i, j, k ) in items:
                angles.append ( Angle ( i = i, j = j, k = k ) )
            self.items = angles
            self.Sort ( )

    @classmethod
    def FromBondContainer ( selfClass, bonds ):
        """Constructor from a bond container."""
        self = selfClass ( )
        if len ( bonds ) > 0:
            bonds.MakeConnections ( )
            items = []
            for j in range ( bonds.UpperBound ( ) ):
                connections  = bonds.GetConnectedAtoms ( j )
                nconnections = len ( connections )
                for i in range ( 1, nconnections ):
                    for k in range ( 0, i ):
                        items.append ( Angle ( i = connections[i], j = j, k = connections[k] ) )
            self.items   = items
            self.QSORTED = True
        return self

    def ItemClass ( self ): return Angle

    def ItemName ( self ): return "Angle"

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
        if len ( self ) > 0: return max ( self[-1].i, self[-1].j ) + 1
        else:                return 0

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
