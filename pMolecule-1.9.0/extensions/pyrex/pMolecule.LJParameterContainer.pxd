#-------------------------------------------------------------------------------
# . File      : pMolecule.LJParameterContainer.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "LJParameterContainer.h":

    ctypedef struct CLJParameterContainer "LJParameterContainer":
        Integer  ntypes
        Real    *epsilon
        Real    *sigma
        Real    *tableA
        Real    *tableB
        Integer *tableindex

    cdef CLJParameterContainer *LJParameterContainer_Allocate              ( Integer ntypes )
    cdef CLJParameterContainer *LJParameterContainer_Clone                 ( CLJParameterContainer  *self )
    cdef void                   LJParameterContainer_Deallocate            ( CLJParameterContainer **self )
    cdef void                   LJParameterContainer_MakeEpsilonSigmaAMBER ( CLJParameterContainer  *self )
    cdef void                   LJParameterContainer_MakeEpsilonSigmaOPLS  ( CLJParameterContainer  *self )
    cdef void                   LJParameterContainer_MakeTableAMBER        ( CLJParameterContainer  *self )
    cdef void                   LJParameterContainer_MakeTableOPLS         ( CLJParameterContainer  *self )
    cdef CLJParameterContainer *LJParameterContainer_MergeEpsilonSigma     ( CLJParameterContainer  *self, CLJParameterContainer *other )
    cdef void                   LJParameterContainer_Scale                 ( CLJParameterContainer  *self, Real scale )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class LJParameterContainer:

    cdef CLJParameterContainer *cObject
    cdef public object          analyticForm
    cdef public object          isOwner
    cdef public object          label
    cdef public object          parameterKeys
