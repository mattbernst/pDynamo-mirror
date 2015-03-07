#-------------------------------------------------------------------------------
# . File      : pCore.RandomNumberGenerator.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, Cardinal, Integer, Real
from pCore.Real1DArray  cimport CReal1DArray, Real1DArray, Real1DArray_Length, Real1DArray_SetItem

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RandomNumberGenerator.h":

    # . The random number generator type type.
    ctypedef struct CRandomNumberGeneratorType "RandomNumberGeneratorType":
        pass

    # . The random number generator type.
    ctypedef struct CRandomNumberGenerator "RandomNumberGenerator":
        pass

    # . Functions.
    cdef CRandomNumberGenerator *RandomNumberGenerator_Allocate     ( CRandomNumberGeneratorType *type )
    cdef CRandomNumberGenerator *RandomNumberGenerator_Clone        ( CRandomNumberGenerator  *self )
    cdef void                    RandomNumberGenerator_Deallocate   ( CRandomNumberGenerator **self )
    cdef Cardinal                RandomNumberGenerator_NextCardinal ( CRandomNumberGenerator  *self )
    cdef Real                    RandomNumberGenerator_NextReal     ( CRandomNumberGenerator  *self )
    cdef Real                    RandomNumberGenerator_NextRealOpen ( CRandomNumberGenerator  *self )
    cdef void                    RandomNumberGenerator_SetSeed      ( CRandomNumberGenerator  *self, Cardinal seed )

    # . Type declarations.
    cdef CRandomNumberGeneratorType *RandomNumberGeneratorType_MersenneTwister

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RandomNumberGenerator:

    cdef CRandomNumberGenerator *cObject
    cdef public object           initialSeed
    cdef public object           isOwner
    cdef public object           owner
