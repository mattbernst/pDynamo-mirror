#-------------------------------------------------------------------------------
# . File      : pCore.Matrix33.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Status       cimport Status, Status_Continue
from pCore.Vector3      cimport CVector3, Vector3, Vector3_GetItem

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Matrix33.h":

    ctypedef struct CMatrix33 "Matrix33":
        pass

    cdef CMatrix33 *Matrix33_Allocate                 ( )
    cdef void       Matrix33_ApplyToVector3           ( CMatrix33  *self, CVector3 *vector3 )
    cdef CMatrix33 *Matrix33_Clone                    ( CMatrix33  *self )
    cdef void       Matrix33_CopyTo                   ( CMatrix33  *self, CMatrix33 *other, Status *status )
    cdef void       Matrix33_Deallocate               ( CMatrix33 **self )
    cdef Real       Matrix33_GetItem                  ( CMatrix33  *self, Integer i1, Integer i2, Status *status )
    cdef void       Matrix33_Invert                   ( CMatrix33  *self, CMatrix33 *other )
    cdef Boolean    Matrix33_IsEqual                  ( CMatrix33  *self, CMatrix33 *other )
    cdef Boolean    Matrix33_IsIdentity               ( CMatrix33  *self )
    cdef Boolean    Matrix33_IsImproperRotation       ( CMatrix33  *self )
    cdef Boolean    Matrix33_IsOrthogonal             ( CMatrix33  *self )
    cdef Boolean    Matrix33_IsProperRotation         ( CMatrix33  *self )
    cdef Integer    Matrix33_Length                   ( CMatrix33  *self, Integer dimension )
    cdef void       Matrix33_PostMultiplyBy           ( CMatrix33  *self, CMatrix33 *other )
    cdef void       Matrix33_PreMultiplyBy            ( CMatrix33  *self, CMatrix33 *other )
    cdef void       CMatrix33_Reflection             "Matrix33_Reflection"             ( CMatrix33 **self, CVector3 *normal )
    cdef Status     CMatrix33_RotationAboutAxis      "Matrix33_RotationAboutAxis"      ( CMatrix33 **self, Real angle, Real x, Real y, Real z )
    cdef Status     CMatrix33_RotationFromQuaternion "Matrix33_RotationFromQuaternion" ( CMatrix33 **self, Real q0, Real q1, Real q2, Real q3 )
    cdef void       Matrix33_Scale                    ( CMatrix33  *self, Real value )
    cdef void       Matrix33_Set                      ( CMatrix33  *self, Real value )
    cdef void       Matrix33_SetItem                  ( CMatrix33  *self, Integer i1, Integer i2, Real value, Status *status )
    cdef void       Matrix33_Transpose                ( CMatrix33  *self, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Matrix33:

    cdef CMatrix33     *cObject
    cdef public object  isOwner
    cdef public object  owner
