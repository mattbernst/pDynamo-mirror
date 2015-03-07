#-------------------------------------------------------------------------------
# . File      : pCore.RandomDeviateGenerator.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions          cimport Boolean, Cardinal, Integer, Real
from pCore.RandomNumberGenerator cimport CRandomNumberGenerator, RandomNumberGenerator
from pCore.Real1DArray           cimport CReal1DArray, Real1DArray, Real1DArray_Length, Real1DArray_SetItem

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RandomNumberDistribution.h":

    # . Functions.
    cdef Real RandomNumberDistribution_GaussianBoxMueller ( CRandomNumberGenerator *rng, Real mu, Real sigma )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NormalDeviateGenerator:

    cdef public RandomNumberGenerator randomNumberGenerator
    cdef public object                mu
    cdef public object                sigma
