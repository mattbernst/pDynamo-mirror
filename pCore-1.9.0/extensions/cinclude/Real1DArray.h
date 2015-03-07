/*------------------------------------------------------------------------------
! . File      : Real1DArray.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _REAL1DARRAY
# define _REAL1DARRAY

# include "Boolean.h"
# include "Integer.h"
# include "Integer1DArray.h"
# include "Real.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The array type. */
typedef struct {
    Boolean  isOwner ;
    Boolean  isView  ;
    Integer  length  ;
    Integer  offset  ;
    Integer  size    ;
    Integer  stride  ;
    Real    *data    ;
} Real1DArray ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define Real1DArray_Data( self ) ( &((self)->data[(self)->offset]) )

/* . An item. */
# define Real1DArray_Item( self, i ) ( (self)->data[(self)->offset+(i)*(self)->stride] )

/* . A pointer to an item. */
# define Real1DArray_ItemPointer( self, i ) ( &((self)->data[(self)->offset+(i)*(self)->stride]) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real         Real1DArray_AbsoluteMaximum      ( const Real1DArray  *self ) ;
extern Integer      Real1DArray_AbsoluteMaximumIndex ( const Real1DArray  *self ) ;
extern void         Real1DArray_Add                  (       Real1DArray  *self, const Real1DArray *other, Status *status ) ;
extern void         Real1DArray_AddScalar            (       Real1DArray  *self, const Real value ) ;
extern void         Real1DArray_AddScaledArray       (       Real1DArray  *self, const Real value, const Real1DArray *other, Status *status ) ;
extern Real1DArray *Real1DArray_Allocate             ( const Integer length, Status *status ) ;
extern Real1DArray *Real1DArray_Clone                ( const Real1DArray  *self, Status *status ) ;
extern void         Real1DArray_CopyTo               ( const Real1DArray  *self, Real1DArray *other, Status *status ) ;
extern void         Real1DArray_Deallocate           (       Real1DArray **self ) ;
extern void         Real1DArray_Divide               (       Real1DArray  *self, const Real1DArray *other, Status *status ) ;
extern Real         Real1DArray_Dot                  ( const Real1DArray  *self, const Real1DArray *other, Status *status ) ;
extern void         Real1DArray_Exp                  (       Real1DArray  *self ) ;
extern Real         Real1DArray_GetItem              ( const Real1DArray  *self, const Integer i, Status *status ) ;
extern Boolean      Real1DArray_IsCompact            ( const Real1DArray  *self ) ;
extern Boolean      Real1DArray_IsUniform            ( const Real1DArray  *self ) ;
extern Integer      Real1DArray_Length               ( const Real1DArray  *self ) ;
extern void         Real1DArray_Ln                   (       Real1DArray  *self ) ;
extern Real         Real1DArray_Maximum              ( const Real1DArray  *self ) ;
extern Real         Real1DArray_Minimum              ( const Real1DArray  *self ) ;
extern void         Real1DArray_Multiply             (       Real1DArray  *self, const Real1DArray *other, Status *status ) ;
extern void         Real1DArray_Normalize            (       Real1DArray  *self, const Real *tolerance, Status *status ) ;
extern Real         Real1DArray_Norm2                ( const Real1DArray  *self ) ;
extern void         Real1DArray_Permute              (       Real1DArray  *self, Integer1DArray *permutation, Status *status ) ;
extern void         Real1DArray_Print                ( const Real1DArray  *self ) ;
extern void         Real1DArray_Range                ( const Real1DArray  *self, Real *lower, Real *upper ) ;
extern void         Real1DArray_Reciprocate          (       Real1DArray  *self ) ;
extern void         Real1DArray_Resize               (       Real1DArray  *self, const Integer length, const Real *initializer, Status *status ) ;
extern void         Real1DArray_Reverse              (       Real1DArray  *self ) ;
extern Real         Real1DArray_RootMeanSquare       ( const Real1DArray  *self ) ;
extern void         Real1DArray_Scale                (       Real1DArray  *self, const Real scale ) ;
extern void         Real1DArray_Set                  (       Real1DArray  *self, Real value ) ;
extern void         Real1DArray_SetItem              (       Real1DArray  *self, const Integer i, const Real value, Status *status ) ;
extern void         Real1DArray_Slice                ( const Real1DArray  *self, const Integer start, const Integer stop, const Integer stride, Real1DArray *slice, Status *status ) ;
extern void         Real1DArray_Sort                 (       Real1DArray  *self ) ;
extern void         Real1DArray_SortIndex            ( const Real1DArray  *self, Integer1DArray *indices, Status *status ) ;
extern void         Real1DArray_SortUnique           (       Real1DArray  *self ) ;
extern Real         Real1DArray_Sparsity             ( const Real1DArray  *self, const Real tolerance ) ;
extern Real         Real1DArray_Sum                  ( const Real1DArray  *self ) ;
extern void         Real1DArray_ViewOfRaw            (       Real1DArray  *self, const Integer offset, const Integer length, const Integer stride, Real *data, const Integer size, Status *status ) ;

# endif
