#-------------------------------------------------------------------------------
# . File      : pCore.SymmetricMatrix.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Real1DArray  cimport CReal1DArray, Real1DArray
from pCore.Real2DArray  cimport CReal2DArray, Real2DArray
from pCore.Status       cimport Status, Status_Continue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "SymmetricMatrix.h":

    ctypedef struct CSymmetricMatrix "SymmetricMatrix":
        Integer  dimension
        Integer  size
        Integer  dimensionP
        Integer  sizeP
        Real    *data

    ctypedef enum SYMMETRICMATRIXUPDATING_OPTION:
        SYMMETRICMATRIXUPDATING_BFGS   = 0
        SYMMETRICMATRIXUPDATING_BOFILL = 1
        SYMMETRICMATRIXUPDATING_MS     = 2
        SYMMETRICMATRIXUPDATING_POWELL = 3

    cdef CSymmetricMatrix *SymmetricMatrix_Allocate               ( Integer n )
    cdef Status            SymmetricMatrix_AddScaledMatrix        ( CSymmetricMatrix  *self, Real alpha, CSymmetricMatrix *other )
    cdef CSymmetricMatrix *SymmetricMatrix_Clone                  ( CSymmetricMatrix  *self )
    cdef Status            SymmetricMatrix_CopyTo                 ( CSymmetricMatrix  *self, CSymmetricMatrix *other )
    cdef void              SymmetricMatrix_Deallocate             ( CSymmetricMatrix **self )
    cdef void              SymmetricMatrix_Diagonalize            ( CSymmetricMatrix  *self, CReal1DArray *eigenvalues, CReal2DArray *eigenvectors, Status *status )
    cdef Integer           SymmetricMatrix_Dimension              ( CSymmetricMatrix  *self )
    cdef Status            SymmetricMatrix_Fill                   ( CSymmetricMatrix  *self, CSymmetricMatrix *other )
    cdef Real              SymmetricMatrix_Get_Component          ( CSymmetricMatrix  *self, Integer i1, Integer i2 )
    cdef Real              SymmetricMatrix_GetItem                ( CSymmetricMatrix  *self, Integer i, Integer j, Status *status )
    cdef Real              SymmetricMatrix_MaximumAbsoluteValue   ( CSymmetricMatrix  *self )
    cdef void              SymmetricMatrix_ProjectOut             ( CSymmetricMatrix **self, CReal2DArray *vectors, Status *status )
    cdef void              SymmetricMatrix_Raise                  ( CSymmetricMatrix  *self, CReal2DArray *vectors, Real value, Status *status )
    cdef void              SymmetricMatrix_Scale                  ( CSymmetricMatrix  *self, Real value )
    cdef void              SymmetricMatrix_Scale_DB               ( CSymmetricMatrix  *self, Integer start, Integer n, Real scale )
    cdef void              SymmetricMatrix_Scale_OB               ( CSymmetricMatrix  *self, Integer istart, Integer ni, Integer jstart, Integer nj, Real scale )
    cdef void              SymmetricMatrix_Set                    ( CSymmetricMatrix  *self, Real value )
    cdef void              SymmetricMatrix_SetItem                ( CSymmetricMatrix  *self, Integer i, Integer j, Real value, Status *status )
    cdef void              SymmetricMatrix_Set_Column_From_Vector ( CSymmetricMatrix  *self, Integer column, Integer mstart, CReal1DArray *vector, Integer vstart, Integer n )
    cdef void              SymmetricMatrix_Set_Component          ( CSymmetricMatrix  *self, Integer i1, Integer i2, Real value )
    cdef void              SymmetricMatrix_Set_Diagonal           ( CSymmetricMatrix  *self, Integer start, Integer n, Real value )
    cdef void              SymmetricMatrix_Set_DB_From_SM_DB      ( CSymmetricMatrix  *matrix1, Integer start1,
                                                                    CSymmetricMatrix  *matrix2, Integer start2, Integer ncolumns )
    cdef void              SymmetricMatrix_Set_Zero               ( CSymmetricMatrix  *self )
    cdef Integer           SymmetricMatrix_Size                   ( CSymmetricMatrix  *self )
    cdef void              SymmetricMatrix_Transfer               ( CSymmetricMatrix  *matrix1, CSymmetricMatrix *matrix2 )
    cdef Status            SymmetricMatrix_Update                 ( CSymmetricMatrix  *self, CReal1DArray *dx, CReal1DArray *dg, SYMMETRICMATRIXUPDATING_OPTION option, Real *tolerance )
    cdef void              SymmetricMatrix_VectorMultiply         ( CSymmetricMatrix  *self, CReal1DArray *other, CReal1DArray *result, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetricMatrix:

    cdef CSymmetricMatrix *cObject
    cdef public object     isOwner
    cdef public object     owner
