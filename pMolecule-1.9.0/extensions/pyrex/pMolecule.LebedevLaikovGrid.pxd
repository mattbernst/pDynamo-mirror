#-------------------------------------------------------------------------------
# . File      : pMolecule.LebedevLaikovGrid.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3 cimport Coordinates3
from pCore.Memory       cimport Memory_Allocate_Array_Integer
from pCore.Selection    cimport Selection, CSelection
from pCore.Status       cimport Status

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "Lebedev.h":

    cdef Integer LebedevLaikov_Angular_Momentum_Value ( Integer npts   )
    cdef Integer LebedevLaikov_Number_Of_Points       ( Integer lvalue )
    cdef Integer LebedevLaikov_Points                 ( Integer N, Real *X, Real *Y, Real *Z, Real *W )
