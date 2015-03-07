/*------------------------------------------------------------------------------
! . File      : Real2DArray.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _REAL2DARRAY
# define _REAL2DARRAY

# include "Boolean.h"
# include "cblas.h"
# include "Integer.h"
# include "Integer1DArray.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Selection.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The array type. */
typedef struct {
    Boolean  isOwner ;
    Boolean  isView  ;
    Integer  length  ;
    Integer  length0 ;
    Integer  length1 ;
    Integer  offset  ;
    Integer  size    ;
    Integer  stride0 ;
    Integer  stride1 ;
    Real    *data    ;
} Real2DArray ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The number of columns. */
# define Real2DArray_Columns( self ) ( (self)->length1 )

/* . A pointer to the start of the data. */
# define Real2DArray_Data( self ) ( &((self)->data[(self)->offset]) )

/* . An item. */
# define Real2DArray_Item( self, i, j ) ( (self)->data[(self)->offset+(i)*(self)->stride0+(j)*(self)->stride1] )

/* . An item by index. */
# define Real2DArray_ItemByIndex( self, ij ) ( (self)->data[ij] )

/* . An item index. */
# define Real2DArray_ItemIndex( self, i, j ) ( (self)->offset+(i)*(self)->stride0+(j)*(self)->stride1 )

/* . A pointer to an item. */
# define Real2DArray_ItemPointer( self, i, j ) ( &((self)->data[(self)->offset+(i)*(self)->stride0+(j)*(self)->stride1]) )

/* . A pointer to a row. */
# define Real2DArray_RowPointer( self, i ) Real2DArray_ItemPointer( self, i, 0 )

/* . The number of rows. */
# define Real2DArray_Rows( self ) ( (self)->length0 )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void         Real2DArray_1DSlice                  ( const Real2DArray  *self, const Integer start0, const Integer stop0, const Integer stride0, const Integer start1, const Integer stop1, const Integer stride1, Real1DArray *slice, Status *status ) ;
extern Real         Real2DArray_AbsoluteMaximum          ( const Real2DArray  *self ) ;
extern void         Real2DArray_AddScaledArray           (       Real2DArray  *self, const Real value, const Real2DArray *other, Status *status ) ;
extern Real2DArray *Real2DArray_Allocate                 ( const Integer length0, const Integer length1, Status *status ) ;
extern Real2DArray *Real2DArray_Clone                    ( const Real2DArray  *self, Status *status ) ;
extern void         Real2DArray_ColumnDotProducts        ( const Boolean      initialize ,
                                                           const Real2DArray *a          ,
                                                           const Real2DArray *b          ,
                                                                 Real1DArray *c          ) ;
extern void         Real2DArray_ColumnSlice              ( const Real2DArray  *self, const Integer column, Real1DArray *slice, Status *status ) ;
extern Status       Real2DArray_CopyFromArrayByRow       (       Real2DArray  *self, const Real1DArray *vector, const Selection *selection ) ;
extern void         Real2DArray_CopyTo                   ( const Real2DArray  *self, Real2DArray *other, Status *status ) ;
extern Status       Real2DArray_CopyToArrayByRow         ( const Real2DArray  *self, Real1DArray *vector, const Selection *selection ) ;
extern void         Real2DArray_Deallocate               (       Real2DArray **self ) ;
extern Real         Real2DArray_Determinant              (       Real2DArray  *self, Status *status ) ;
extern void         Real2DArray_Exp                      (       Real2DArray  *self ) ;
extern Real         Real2DArray_FrobeniusNorm            ( const Real2DArray  *self ) ;
extern Real         Real2DArray_GetItem                  ( const Real2DArray  *self, const Integer i, const Integer j, Status *status ) ;
extern Integer      Real2DArray_GramSchmidtOrthogonalize ( const Real2DArray  *self, const Integer *maximumIterations, const Integer *numberConstant, const Real *tolerance, Status *status ) ;
extern Boolean      Real2DArray_IsCompact                ( const Real2DArray  *self, const Integer dimension, Status *status ) ;
extern Boolean      Real2DArray_IsOrthogonal             ( const Real2DArray  *self, const Real *tolerance, Real *deviation, Status *status ) ;
extern Boolean      Real2DArray_IsSymmetric              ( const Real2DArray  *self, const Real *tolerance, Real *deviation ) ;
extern Boolean      Real2DArray_IsUniform                ( const Real2DArray  *self ) ;
extern Integer      Real2DArray_Length                   ( const Real2DArray  *self, const Integer dimension ) ;
extern void         Real2DArray_LUPFactorizationInPlace  (       Real2DArray  *self, Integer1DArray *pivots, Status *status ) ;
extern void         Real2DArray_MatrixMultiply           ( const Boolean aTranspose, const Boolean bTranspose, const Real alpha, const Real2DArray *a, const Real2DArray *b, const Real beta, Real2DArray *c, Status *status ) ;
extern void         Real2DArray_Multiply                 (       Real2DArray  *self, const Real2DArray *other, Status *status ) ;
extern void         Real2DArray_Print                    ( const Real2DArray  *self ) ;
extern void         Real2DArray_ProjectOutOf1DArray      ( const Real2DArray  *self, Real1DArray *vector, Status *status ) ;
extern Real2DArray *Real2DArray_PruneByRow               ( const Real2DArray  *self, const Selection *selection ) ;
extern void         Real2DArray_Resize                   (       Real2DArray  *self, const Integer length0, const Real *initializer, Status *status ) ;
extern Real         Real2DArray_RootMeanSquare           ( const Real2DArray  *self ) ;
extern void         Real2DArray_RowCopyTo                ( const Real2DArray  *self, const Integer m, Real2DArray *other, const Integer n, Status *status ) ;
extern void         Real2DArray_RowSlice                 ( const Real2DArray  *self, const Integer row, Real1DArray *slice, Status *status ) ;
extern void         Real2DArray_Scale                    (       Real2DArray  *self, Real value ) ;
extern void         Real2DArray_Set                      (       Real2DArray  *self, Real value ) ;
extern void         Real2DArray_SetByRow                 (       Real2DArray  *self, const Selection *selection, const Real alpha ) ;
extern void         Real2DArray_SetItem                  (       Real2DArray  *self, const Integer i, const Integer j, const Real value, Status *status ) ;
extern void         Real2DArray_Slice                    ( const Real2DArray  *self, const Integer start0, const Integer stop0, const Integer stride0, const Integer start1, const Integer stop1, const Integer stride1, Real2DArray *slice, Status *status ) ;
extern Real         Real2DArray_Trace                    ( const Real2DArray  *self, const Real2DArray *other, Status *status ) ;
extern void         Real2DArray_Transpose                (       Real2DArray  *self, Status *status ) ;
extern void         Real2DArray_VectorMultiply           ( const Boolean aTranspose, const Real alpha, const Real2DArray *a, const Real1DArray *x, const Real beta, Real1DArray *y, Status *status ) ;
extern void         Real2DArray_ViewOfRaw                (       Real2DArray  *self, const Integer offset, const Integer length0, const Integer stride0, const Integer length1, const Integer stride1, Real *data, const Integer size, Status *status ) ;

# endif
