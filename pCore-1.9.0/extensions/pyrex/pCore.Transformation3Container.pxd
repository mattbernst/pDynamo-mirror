#-------------------------------------------------------------------------------
# . File      : pCore.Transformation3Container.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions    cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Transformation3 cimport CTransformation3, Transformation3, Transformation3_IsEqual,  Transformation3_IsIdentity

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "Transformation3Container.h":

    ctypedef struct CTransformation3Container "Transformation3Container":
        Boolean  isOwner
        Integer  identity
        Integer  nitems
        Integer *inverses
        CTransformation3 **items

    cdef CTransformation3Container *Transformation3Container_Allocate     ( Integer nitems )
    cdef void                       Transformation3Container_Deallocate   ( CTransformation3Container **self )
    cdef void                       Transformation3Container_FindIdentity ( CTransformation3Container  *self )
    cdef void                       Transformation3Container_FindInverses ( CTransformation3Container  *self )

#===============================================================================
# . Class.
#===============================================================================
cdef class Transformation3Container:

    cdef CTransformation3Container *cObject
    cdef public object              isOwner
    cdef public object              items
