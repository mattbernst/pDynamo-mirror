#-------------------------------------------------------------------------------
# . File      : pCore.Vector3.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Status       cimport Status, Status_Continue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Vector3.h":

    ctypedef struct CVector3 "Vector3":
        Boolean  isOwner
        Boolean  isView
        Integer  length
        Integer  offset
        Integer  size
        Integer  stride
        Real    *data

    cdef Real      Vector3_AbsoluteMaximum      ( CVector3  *self )
    cdef Integer   Vector3_AbsoluteMaximumIndex ( CVector3  *self )
    cdef void      Vector3_AddScalar            ( CVector3  *self, Real alpha )
    cdef void      Vector3_AddScaledVector      ( CVector3  *self, Real alpha, CVector3 *other, Status *status )
    cdef CVector3 *Vector3_Allocate             ( )
    cdef CVector3 *Vector3_Clone                ( CVector3  *self )
    cdef void      Vector3_CopyTo               ( CVector3  *self, CVector3 *other, Status *status )
    cdef void      Vector3_CrossProduct         ( CVector3  *self, CVector3 *other )
    cdef void      Vector3_Deallocate           ( CVector3 **self )
    cdef void      Vector3_Divide               ( CVector3  *self, CVector3 *other, Status *status )
    cdef Real      Vector3_Dot                  ( CVector3  *self, CVector3 *other, Status *status )
    cdef Integer   Vector3_Length               ( CVector3  *self )
    cdef Real      Vector3_GetItem              ( CVector3  *self, Integer i, Status *status )
    cdef Real      Vector3_Maximum              ( CVector3  *self )
    cdef Real      Vector3_Minimum              ( CVector3  *self )
    cdef void      Vector3_Multiply             ( CVector3  *self, CVector3 *other, Status *status )
    cdef Real      Vector3_Norm2                ( CVector3  *self )
    cdef void      Vector3_Normalize            ( CVector3  *self, Real *tolerance, Status *status )
    cdef Real      Vector3_RootMeanSquare       ( CVector3  *self )
    cdef void      Vector3_Scale                ( CVector3  *self, Real alpha )
    cdef void      Vector3_Set                  ( CVector3  *self, Real alpha )
    cdef void      Vector3_SetItem              ( CVector3  *self, Integer i, Real alpha, Status *status )
    cdef Real      Vector3_Sum                  ( CVector3  *self )
    cdef void      Vector3_ViewOfRaw            ( CVector3  *self, Integer offset, Integer length, Integer stride, Real *data, Integer size, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Vector3:

    cdef CVector3     *cObject
    cdef public object isOwner
    cdef public object owner
