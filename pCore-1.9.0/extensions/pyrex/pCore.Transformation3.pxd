#-------------------------------------------------------------------------------
# . File      : pCore.Transformation3.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Matrix33     cimport CMatrix33, Matrix33
from pCore.Vector3      cimport CVector3, Vector3

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "Transformation3.h":

    ctypedef struct CTransformation3 "Transformation3":
        CMatrix33 *rotation
        CVector3  *translation

    cdef CTransformation3 *Transformation3_Allocate      ( )
    cdef CTransformation3 *Transformation3_Clone         ( CTransformation3  *self )
    cdef void              Transformation3_Copy          ( CTransformation3  *self, CTransformation3 *other )
    cdef void              Transformation3_Deallocate    ( CTransformation3 **self )
    cdef Boolean           Transformation3_IsEqual       ( CTransformation3  *self, CTransformation3 *other )
    cdef Boolean           Transformation3_IsIdentity    ( CTransformation3  *self )
    cdef void              Transformation3_Orthogonalize ( CTransformation3  *self, CMatrix33 *A, CMatrix33 *B )

#===============================================================================
# . Class.
#===============================================================================
cdef class Transformation3:

    cdef CTransformation3 *cObject
    cdef public object     isOwner
    cdef public object     rotation
    cdef public object     translation
