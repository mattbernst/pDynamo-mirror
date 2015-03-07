/*------------------------------------------------------------------------------
! . File      : SymmetricMatrix.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SYMMETRICMATRIX
# define _SYMMETRICMATRIX

# include "Definitions.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The symmetric matrix type. */
/* . dimensionP and sizeP are the allocated dimension and size whereas
! . dimension and size are the actual dimension and size (less than or
! . equal to dimensionP and sizeP). */
typedef struct {
   Integer dimension,  size  ;
   Integer dimensionP, sizeP ;
   Real   *data ;
} SymmetricMatrix ;

/* . Updating options. */
typedef enum {
    SYMMETRICMATRIXUPDATING_BFGS   = 0,
    SYMMETRICMATRIXUPDATING_BOFILL = 1,
    SYMMETRICMATRIXUPDATING_MS     = 2,
    SYMMETRICMATRIXUPDATING_POWELL = 3
} SYMMETRICMATRIXUPDATING_OPTION ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define SymmetricMatrix_Data( self ) ( &((self)->data[0]) )

/* . Index to an item on the diagonal. */
# define SymmetricMatrix_DiagonalIndex(i) ( i * ( i + 1 ) ) / 2

/* . Access an item (i >= j). */
# define SymmetricMatrix_Item( self, i, j ) ( (self)->data[( (i) * ( i + 1 ) ) / 2 + j] )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Status           SymmetricMatrix_AddScaledMatrix               (       SymmetricMatrix  *self, const Real alpha, const SymmetricMatrix *other ) ;
extern SymmetricMatrix *SymmetricMatrix_Allocate                      ( const Integer n ) ;
extern SymmetricMatrix *SymmetricMatrix_AllocateN                     ( const Integer n, Status *status ) ;
extern void             SymmetricMatrix_AnticommutatorSS              (       SymmetricMatrix  *self, const SymmetricMatrix *a, const SymmetricMatrix *b, Status *status ) ;
extern SymmetricMatrix *SymmetricMatrix_Clone                         ( const SymmetricMatrix  *self ) ;
extern Status           SymmetricMatrix_CopyTo                        ( const SymmetricMatrix  *self, SymmetricMatrix *other ) ;
extern void             SymmetricMatrix_CopyFromReal2DArray           (       SymmetricMatrix  *self, const Real2DArray *other, Status *status ) ;
extern void             SymmetricMatrix_CopyToReal2DArray             ( const SymmetricMatrix  *self,       Real2DArray *other, Status *status ) ;
extern void             SymmetricMatrix_Deallocate                    (       SymmetricMatrix **self ) ;
extern void             SymmetricMatrix_Decrement                     (       SymmetricMatrix  *self, SymmetricMatrix *other ) ;
extern void             SymmetricMatrix_Diagonalize                   (       SymmetricMatrix  *self, Real1DArray *eigenValues, Real2DArray *eigenVectors, Status *status ) ;
extern void             SymmetricMatrix_DiagonalizePartial            (       SymmetricMatrix  *self, Integer lower, Integer upper, Real1DArray *eigenValues, Real2DArray *eigenVectors, Status *status ) ;
extern Integer          SymmetricMatrix_Dimension                     ( const SymmetricMatrix  *self ) ;
extern void             SymmetricMatrix_GetColumn                     ( const SymmetricMatrix  *self, const Integer n, Real1DArray *column, Status *status ) ;
extern Real             SymmetricMatrix_GetItem                       ( const SymmetricMatrix  *self, const Integer i, const Integer j, Status *status ) ;
extern Real             SymmetricMatrix_Get_Component                 (       SymmetricMatrix  *self, const Integer i, const Integer j ) ;
extern void             SymmetricMatrix_Increment                     (       SymmetricMatrix  *self, SymmetricMatrix *other ) ;
extern void             SymmetricMatrix_IncrementComponent            (       SymmetricMatrix  *self, const Integer i, const Integer j, const Real value ) ;
extern void             SymmetricMatrix_IncrementDiagonal             (       SymmetricMatrix  *self, const Real value ) ;
extern void             SymmetricMatrix_IncrementDiagonalFromArray    ( const SymmetricMatrix  *self, const Real1DArray *diagonal, Status *status ) ;
extern void             SymmetricMatrix_IncrementItem                 (       SymmetricMatrix  *self, const Integer i, const Integer j, const Real value, Status *status ) ;
extern void             SymmetricMatrix_Increment_DB                  (       SymmetricMatrix  *self, const Integer start, const Integer n, const Real *dblock ) ;
extern void             SymmetricMatrix_Increment_DBlockM             (       SymmetricMatrix  *self, const Integer istart, const Integer ni, const Real *dblock ) ;
extern void             SymmetricMatrix_Increment_OB                  (       SymmetricMatrix  *self, const Integer istart, const Integer ni,
                                                                                                                     const Integer jstart, const Integer nj, const Real *oblock ) ;
extern void             SymmetricMatrix_Increment_OBlockM             (       SymmetricMatrix  *self, const Integer istart, const Integer ni, const Integer jstart, const Integer nj, const Real *oblock ) ;
extern void             SymmetricMatrix_IndexedCopyToReal2DArray      ( const SymmetricMatrix  *self, const Integer1DArray *indices, Real2DArray *target, Status *status ) ;
extern void             SymmetricMatrix_Invert                        ( const SymmetricMatrix  *self, const Real *tolerance, SymmetricMatrix *inverse, Status *status ) ;
extern Status           SymmetricMatrix_InverseSquareRoot             (       SymmetricMatrix  *self, Real1DArray *eigenValues, const Real2DArray *eigenVectors ) ;
extern Boolean          SymmetricMatrix_IsCompact                     ( const SymmetricMatrix  *self ) ;
extern Boolean          SymmetricMatrix_IsDiagonal                    ( const SymmetricMatrix  *self, const Real tolerance ) ;
extern Boolean          SymmetricMatrix_IsUniform                     ( const SymmetricMatrix  *self ) ;
extern void             SymmetricMatrix_MakeFromEigensystem           (       SymmetricMatrix  *self            ,
                                                                        const Boolean           zeroMatrix      ,
                                                                        const Integer           numberOfVectors ,
                                                                        const Real1DArray      *eigenvalues     ,
                                                                        const Real2DArray      *eigenvectors    ,
                                                                              Status           *status          ) ;
extern Real             SymmetricMatrix_MaximumAbsoluteValue          ( const SymmetricMatrix  *self ) ;
extern Status           SymmetricMatrix_Multiply2                     ( const SymmetricMatrix  *self, const SymmetricMatrix *other, Real2DArray *result ) ;
extern Real             SymmetricMatrix_Multiply2_Trace               ( const SymmetricMatrix  *self, const SymmetricMatrix *other ) ;
extern Real             SymmetricMatrix_MultiplyASBS_Trace            ( const SymmetricMatrix  *a, const SymmetricMatrix *b, const SymmetricMatrix *s, Status *status ) ;
extern Real2DArray     *SymmetricMatrix_OrthogonalizingTransformation (       SymmetricMatrix  *self, const Real *diagonalTolerance, const Real *eigenValueTolerance, Status *status ) ;
extern Status           SymmetricMatrix_PostMatrixMultiply            ( const SymmetricMatrix  *self, const Real2DArray *matrix, const Boolean useTranspose, Real2DArray *result ) ;
extern void             SymmetricMatrix_PreMatrixMultiply             ( const SymmetricMatrix  *self, const Real2DArray *matrix, const Boolean useTranspose, Real2DArray *result, Status *status ) ;
extern void             SymmetricMatrix_Print                         (       SymmetricMatrix  *self ) ;
extern void             SymmetricMatrix_ProjectOut                    (       SymmetricMatrix **self, const Real2DArray *vectors, Status *status ) ;
extern SymmetricMatrix *SymmetricMatrix_ProjectionMatrix              ( const Real2DArray *vectors ) ;
extern void             SymmetricMatrix_Rank1Update                   (       SymmetricMatrix  *self, const Real alpha, const Real1DArray *vector, Status *status ) ;
extern void             SymmetricMatrix_Raise                         (       SymmetricMatrix  *self, const Real2DArray *vectors, const Real value, Status *status ) ;
extern void             SymmetricMatrix_Scale                         (       SymmetricMatrix  *self, const Real value ) ;
extern void             SymmetricMatrix_Scale_Diagonal                (       SymmetricMatrix  *self, const Integer start, const Integer n, const Real scale ) ;
extern void             SymmetricMatrix_Scale_DB                      (       SymmetricMatrix  *self, const Integer start, const Integer n, const Real scale ) ;
extern void             SymmetricMatrix_Scale_OB                      (       SymmetricMatrix  *self, const Integer istart, const Integer ni,
                                                                                                                     const Integer jstart, const Integer nj, const Real scale ) ;
extern void             SymmetricMatrix_Scale_OD                      (       SymmetricMatrix  *self, const Real value ) ;
extern void             SymmetricMatrix_Set                           (       SymmetricMatrix  *self, const Real value ) ;
extern void             SymmetricMatrix_SetItem                       (       SymmetricMatrix  *self, const Integer i, const Integer j, const Real value, Status *status ) ;
extern void             SymmetricMatrix_Set_Column_From_Vector        (       SymmetricMatrix  *self,  const Integer column, const Integer mstart,
                                                                                                             const Real1DArray *vector, const Integer vstart, const Integer n ) ;
extern void             SymmetricMatrix_Set_Component                 (       SymmetricMatrix  *self, const Integer i, const Integer j, const Real value ) ;
extern void             SymmetricMatrix_Set_DBlockM                   (       SymmetricMatrix  *self, const Integer istart, const Integer ni, const Real *dblock ) ;
extern void             SymmetricMatrix_Set_DB_From_SM_DB             (       SymmetricMatrix  *self , const Integer start1,
                                                                              SymmetricMatrix  *other, const Integer start2,
                                                                                                       const Integer ncolumns ) ;
extern void             SymmetricMatrix_Set_Diagonal                  (       SymmetricMatrix  *self, const Integer start, const Integer n, const Real value ) ;
extern void             SymmetricMatrix_Set_Diagonal_V                (       SymmetricMatrix  *self, const Integer start, const Integer n, const Real *values ) ;
extern void             SymmetricMatrix_Set_OBlockM                   (       SymmetricMatrix  *self, const Integer istart, const Integer ni, const Integer jstart, const Integer nj, const Real *oblock ) ;
extern void             SymmetricMatrix_Set_OB_From_SM_OB             (       SymmetricMatrix  *self , const Integer istart1, const Integer jstart1,
                                                                              SymmetricMatrix  *other, const Integer istart2, const Integer jstart2,
                                                                                                            const Integer ni, const Integer nj ) ;
extern void             SymmetricMatrix_Set_Zero                      (       SymmetricMatrix  *self ) ;
extern Integer          SymmetricMatrix_Size                          ( const SymmetricMatrix  *self ) ;
extern void             SymmetricMatrix_SolveLinearEquations          (       SymmetricMatrix  *self, Real1DArray *rhs, Status *status ) ;
extern Real             SymmetricMatrix_Sparsity                      ( const SymmetricMatrix  *self, const Real tolerance ) ;
extern Status           SymmetricMatrix_SquareRoot                    (       SymmetricMatrix  *self, Real1DArray *eigenValues, const Real2DArray *eigenVectors ) ;
extern Status           SymmetricMatrix_SymmetricMultiply2            ( const SymmetricMatrix  *self, const SymmetricMatrix *other, SymmetricMatrix *result ) ;
extern void             SymmetricMatrix_SymmetricTransform            ( const SymmetricMatrix  *a, const SymmetricMatrix *b, SymmetricMatrix *bab, Status *status ) ;
extern Real             SymmetricMatrix_Trace                         ( const SymmetricMatrix  *self ) ;
extern Real             SymmetricMatrix_Trace2                        ( const SymmetricMatrix  *self, const SymmetricMatrix *other, Status *status ) ;
extern void             SymmetricMatrix_Transfer                      (       SymmetricMatrix  *matrix1, SymmetricMatrix *matrix2 ) ;
extern Status           SymmetricMatrix_Transform                     ( const SymmetricMatrix  *self, const Real2DArray *matrix, const Boolean useTranspose, SymmetricMatrix *result ) ;
extern Real             SymmetricMatrix_TransformByMatrixColumns      (       SymmetricMatrix  *self, Real2DArray *matrix, const Integer i, const Integer l ) ;
extern void             SymmetricMatrix_Transform_In_Place            (       SymmetricMatrix  *self, Real2DArray *matrix ) ;
extern void             SymmetricMatrix_Unweight                      (       SymmetricMatrix  *self ) ;
extern Status           SymmetricMatrix_Update                        (       SymmetricMatrix  *self, const Real1DArray *dx, const Real1DArray *dg, const SYMMETRICMATRIXUPDATING_OPTION option, const Real *tolerance ) ;
extern void             SymmetricMatrix_VectorMultiply                ( const SymmetricMatrix  *self, const Real1DArray *other, Real1DArray *result, Status *status ) ;
extern void             SymmetricMatrix_ViewOfRaw                     (       SymmetricMatrix  *self, const Integer length, const Integer stride, Real *data, const Integer rawstride, Status *status ) ;

# endif
