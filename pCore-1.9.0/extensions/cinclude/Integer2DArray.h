/*------------------------------------------------------------------------------
! . File      : Integer2DArray.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _INTEGER2DARRAY
# define _INTEGER2DARRAY

# include "Boolean.h"
# include "Integer.h"
# include "Integer1DArray.h"
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
    Integer *data    ;
} Integer2DArray ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define Integer2DArray_Data( self ) ( &((self)->data[(self)->offset]) )

/* . An item. */
# define Integer2DArray_Item( self, i, j ) ( (self)->data[(self)->offset+(i)*(self)->stride0+(j)*(self)->stride1] )

/* . An item index. */
# define Integer2DArray_ItemIndex( self, i, j ) ( (self)->offset+(i)*(self)->stride0+(j)*(self)->stride1 )

/* . A pointer to an item. */
# define Integer2DArray_ItemPointer( self, i, j ) ( &((self)->data[(self)->offset+(i)*(self)->stride0+(j)*(self)->stride1]) )

/* . A pointer to a row. */
# define Integer2DArray_RowPointer( self, i ) Integer2DArray_ItemPointer( self, i, 0 )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void            Integer2DArray_1DSlice     ( const Integer2DArray  *self, const Integer start0, const Integer stop0, const Integer stride0, const Integer start1, const Integer stop1, const Integer stride1, Integer1DArray *slice, Status *status ) ;
extern Integer2DArray *Integer2DArray_Allocate    ( const Integer length0, const Integer length1, Status *status ) ;
extern Integer2DArray *Integer2DArray_Clone       ( const Integer2DArray  *self, Status *status ) ;
extern void            Integer2DArray_ColumnSlice ( const Integer2DArray  *self, const Integer column, Integer1DArray *slice, Status *status ) ;
extern void            Integer2DArray_CopyTo      ( const Integer2DArray  *self, Integer2DArray *other, Status *status ) ;
extern void            Integer2DArray_Deallocate  (       Integer2DArray **self ) ;
extern Integer         Integer2DArray_GetItem     ( const Integer2DArray  *self, const Integer i, const Integer j, Status *status ) ;
extern Integer         Integer2DArray_Length      ( const Integer2DArray  *self, const Integer dimension ) ;
extern void            Integer2DArray_Resize      (       Integer2DArray  *self, const Integer length0, const Integer *initializer, Status *status ) ;
extern void            Integer2DArray_RowSlice    ( const Integer2DArray  *self, const Integer row, Integer1DArray *slice, Status *status ) ;
extern void            Integer2DArray_Set         (       Integer2DArray  *self, Integer value ) ;
extern void            Integer2DArray_SetItem     (       Integer2DArray  *self, const Integer i, const Integer j, const Integer value, Status *status ) ;
extern void            Integer2DArray_Slice       ( const Integer2DArray  *self, const Integer start0, const Integer stop0, const Integer stride0, const Integer start1, const Integer stop1, const Integer stride1, Integer2DArray *slice, Status *status ) ;
extern void            Integer2DArray_ViewOfRaw   (       Integer2DArray  *self, const Integer offset, const Integer length0, const Integer stride0, const Integer length1, const Integer stride1, Integer *data, const Integer size, Status *status ) ;

# endif
