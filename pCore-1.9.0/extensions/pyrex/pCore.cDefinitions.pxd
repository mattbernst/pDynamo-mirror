#-------------------------------------------------------------------------------
# . File      : pCore.cDefinitions.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Definitions.h":

    ctypedef short unsigned int Integer16
    ctypedef       unsigned int Integer32

    ctypedef size_t CSize

cdef extern from "Boolean.h":

    ctypedef enum Boolean:
        CFalse "False" = 0
        CTrue  "True"  = 1

cdef extern from "Cardinal.h":

    ctypedef unsigned int Cardinal

cdef extern from "Integer.h":

    ctypedef int Integer

cdef extern from "Real.h":

    ctypedef double Real
