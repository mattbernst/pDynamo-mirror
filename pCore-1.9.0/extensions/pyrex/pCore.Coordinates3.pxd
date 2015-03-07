#-------------------------------------------------------------------------------
# . File      : pCore.Coordinates3.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions    cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Matrix33        cimport CMatrix33, Matrix33
from pCore.Real1DArray     cimport CReal1DArray, Real1DArray, Real1DArray_ViewOfRaw
from pCore.Real2DArray     cimport CReal2DArray, Real2DArray, Real2DArray_VectorMultiply, Real2DArray_ViewOfRaw
from pCore.RegularGrid     cimport CRegularGrid, RegularGrid
from pCore.Selection       cimport CSelection, Selection
from pCore.Status          cimport Status, Status_Continue, Status_Success
from pCore.SymmetricMatrix cimport CSymmetricMatrix, SymmetricMatrix
from pCore.Transformation3 cimport CTransformation3, Transformation3
from pCore.Vector3         cimport CVector3, Vector3, Vector3_ViewOfRaw

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "Coordinates3.h":

    ctypedef struct CCoordinates3 "Coordinates3":
        Boolean  isOwner
        Boolean  isView
        Integer  length
        Integer  length0
        Integer  length1
        Integer  offset
        Integer  size
        Integer  stride0
        Integer  stride1
        Real    *data

    # . Real2DArray procedures.
    cdef void              Coordinates3_1DSlice                                 ( CCoordinates3  *self, Integer start0, Integer stop0, Integer stride0, Integer start1, Integer stop1, Integer stride1, CReal1DArray  *slice, Status *status )
    cdef void              Coordinates3_2DSlice                                 ( CCoordinates3  *self, Integer start0, Integer stop0, Integer stride0, Integer start1, Integer stop1, Integer stride1, CReal2DArray  *slice, Status *status )
    cdef Real              Coordinates3_AbsoluteMaximum                         ( CCoordinates3  *self )
    cdef void              Coordinates3_AddScaledArray                          ( CCoordinates3  *self, Real value, CCoordinates3 *other, Status *status )
    cdef Status            Coordinates3_CopyFromArray                           ( CCoordinates3  *self, CReal1DArray *vector, CSelection *selection )
    cdef void              Coordinates3_CopyTo                                  ( CCoordinates3  *self, CCoordinates3 *other, Status *status )
    cdef Status            Coordinates3_CopyToArray                             ( CCoordinates3  *self, CReal1DArray *vector, CSelection *selection )
    cdef void              Coordinates3_Deallocate                              ( CCoordinates3 **self )
    cdef Real              Coordinates3_GetItem                                 ( CCoordinates3  *self, Integer i, Integer j, Status *status )
    cdef Integer           Coordinates3_Length                                  ( CCoordinates3  *self, Integer dimension )
    cdef CCoordinates3    *Coordinates3_Prune                                   ( CCoordinates3  *self, CSelection *selection )
    cdef Real              Coordinates3_RootMeanSquare                          ( CCoordinates3  *self )
    cdef void              Coordinates3_SetItem                                 ( CCoordinates3  *self, Integer i, Integer j, Real value, Status *status )
    cdef void              Coordinates3_Scale                                   ( CCoordinates3  *self, Real value )
    cdef void              Coordinates3_Set                                     ( CCoordinates3  *self, Real value )
    cdef void              Coordinates3_SetByRow                                ( CCoordinates3  *self, CSelection *selection, Real alpha )
    cdef void              Coordinates3_SliceVector3                            ( CCoordinates3  *self, Integer start0, Integer stop0, Integer stride0, Integer start1, Integer stop1, Integer stride1, CVector3      *slice, Status *status )
    cdef void              Coordinates3_SliceCoordinates3                       ( CCoordinates3  *self, Integer start0, Integer stop0, Integer stride0, Integer start1, Integer stop1, Integer stride1, CCoordinates3 *slice, Status *status )
    cdef void              Coordinates3_ViewOfRaw                               ( CCoordinates3  *self, Integer offset, Integer length0, Integer stride0, Integer length1, Integer stride1, Real *data, Integer size, Status *status )

    # . Coordinate specific procedures.
    cdef CCoordinates3    *Coordinates3_Allocate                                ( Integer extent )
    cdef Real              Coordinates3_Angle                                   ( CCoordinates3  *self, Integer i, Integer j, Integer k )
    cdef Status            Coordinates3_BuildPointFromDistance                  ( CCoordinates3  *self, Integer i, Integer j, Real r, CVector3 *direction )
    cdef Status            Coordinates3_BuildPointFromDistanceAngle             ( CCoordinates3  *self, Integer i, Integer j, Integer k, Real r, Real theta, CVector3 *direction )
    cdef Status            Coordinates3_BuildPointFromDistanceAngleDihedral     ( CCoordinates3  *self, Integer i, Integer j, Integer k, Integer l, Real r, Real theta, Real phi )
    cdef Status            Coordinates3_BuildPointFromDistancePlaneAngle        ( CCoordinates3  *self, Integer i, Integer j, Integer k, Integer l, Real r, Real planeangle )
    cdef Status            Coordinates3_BuildPointFromDistanceTetrahedralTripod ( CCoordinates3  *self, Integer i, Integer j, Integer k, Integer l, Integer m, Real r )
    cdef Status            Coordinates3_Center                                  ( CCoordinates3  *self, CSelection *selection, CReal1DArray *weights, CVector3 **center )
    cdef Real              Coordinates3_Dihedral                                ( CCoordinates3  *self, Integer i, Integer j, Integer k, Integer l )
    cdef Real              Coordinates3_Distance                                ( CCoordinates3  *self, Integer i, Integer j )
    cdef void              Coordinates3_EnclosingOrthorhombicBox                ( CCoordinates3  *self, CSelection *selection, CReal1DArray *radii, CVector3 *origin, CVector3 *extents )
    cdef void              Coordinates3_Gather                                  ( CCoordinates3  *self, CCoordinates3 *other, CSelection *selection )
    cdef void              Coordinates3_GatherAddScaledMatrix                   ( CCoordinates3  *self, Real alpha, CCoordinates3 *other, CSelection *selection )
    cdef Status            CCoordinates3_FromRegularGrid "Coordinates3_FromRegularGrid" ( CCoordinates3 **self, CRegularGrid *grid, CSelection *selection )
    cdef Status            Coordinates3_IdentifyOccupiedGridPoints              ( CCoordinates3  *self, CRegularGrid *grid, CReal1DArray *radii, Boolean QMIDPOINTOVERLAP, CSelection **occupied )
    cdef CSymmetricMatrix *Coordinates3_InertiaMatrix                           ( CCoordinates3  *self, CSelection *selection, CReal1DArray *weights )
    cdef void              Coordinates3_MomentsOfInertia                        ( CCoordinates3  *self, CSelection *selection, CReal1DArray *weights, CVector3 *moments, CMatrix33 *axes )
    cdef Real              Coordinates3_RadiusOfGyration                        ( CCoordinates3  *self, CSelection *selection, CReal1DArray *weights )
    cdef Real              Coordinates3_RMSDeviation                            ( CCoordinates3  *self, CCoordinates3 *other, CSelection *selection, CReal1DArray *weights )
    cdef void              Coordinates3_Rotate                                  ( CCoordinates3  *self, CMatrix33 *rotation, CSelection *selection )
    cdef CReal2DArray     *Coordinates3_RotationTranslationVectors              ( CCoordinates3  *self, CReal1DArray *weights, Boolean QRx, Boolean QRy, Boolean QRz, Boolean QTx, Boolean QTy, Boolean QTz, Integer dimension )
    cdef void              Coordinates3_ScaleRows                               ( CCoordinates3  *self, CReal1DArray *rowScalingFactors, Status *status )
    cdef void              Coordinates3_Scatter                                 ( CCoordinates3  *self, CCoordinates3 *other, CSelection *selection )
    cdef void              Coordinates3_ScatterAddScaledMatrix                  ( CCoordinates3  *self, Real alpha, CCoordinates3 *other, CSelection *selection )
    cdef void              Coordinates3_Superimpose                             ( CCoordinates3  *self, CCoordinates3 *other, CSelection *selection, CReal1DArray *weights, CMatrix33 *rotation, CVector3 *translation )
    cdef void              Coordinates3_ToPrincipalAxes                         ( CCoordinates3  *self, CSelection *selection, CReal1DArray *weights )
    cdef void              Coordinates3_Transform                               ( CCoordinates3  *self, CTransformation3 *transformation, CSelection *selection )
    cdef void              Coordinates3_Translate                               ( CCoordinates3  *self, CVector3 *translation, CSelection *selection )
    cdef void              Coordinates3_TranslateToCenter                       ( CCoordinates3  *self, CSelection *selection, CReal1DArray *weights )

#===============================================================================
# . Class.
#===============================================================================
cdef class Coordinates3:

    cdef CCoordinates3 *cObject
    cdef public object  isOwner
    cdef public object  owner
    cdef public object  undefined

