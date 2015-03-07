#-------------------------------------------------------------------------------
# . File      : pMolecule.LebedevLaikovGrid.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Functions for Lebedev-Laikov grids."""

from pCore.Coordinates3 cimport Coordinates3
from pCore.Memory       cimport Memory_Allocate_Array_Real, Memory_Deallocate_Real
from pCore.Real1DArray  cimport Real1DArray

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def LebedevLaikovGrid_GetGridPoints ( gridAngularMomentum ):
    """Get grid points given an angular momentum."""
    cdef Integer      i, npts
    cdef Real        *w, *x, *y, *z
    cdef Coordinates3 gridPoints
    cdef Real1DArray  weights
    gridPoints = None
    weights    = None
    npts       = LebedevLaikov_Number_Of_Points ( gridAngularMomentum )
    if npts > 0:
        # . Allocate space.
        x = Memory_Allocate_Array_Real ( npts )
        y = Memory_Allocate_Array_Real ( npts )
        z = Memory_Allocate_Array_Real ( npts )
        w = Memory_Allocate_Array_Real ( npts )
        # . Get the grid points and weights.
        LebedevLaikov_Points ( npts, x, y, z, w )
        # . Fill final arrays.
        gridPoints = Coordinates3.WithExtent ( npts )
        weights    = Real1DArray.WithExtent  ( npts )
        for i from 0 <= i < npts:
            gridPoints[i,0] = x[i]
            gridPoints[i,1] = y[i]
            gridPoints[i,2] = z[i]
            weights[i]      = w[i]
        # . Deallocate space.
        Memory_Deallocate_Real ( &x )
        Memory_Deallocate_Real ( &y )
        Memory_Deallocate_Real ( &z )
        Memory_Deallocate_Real ( &w )
    return ( gridPoints, weights )
