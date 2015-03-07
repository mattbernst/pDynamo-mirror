/*------------------------------------------------------------------------------
! . File      : Integer1DArray.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _INTEGER1DARRAY
# define _INTEGER1DARRAY

# include "Boolean.h"
# include "Integer.h"
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
    Integer *data    ;
} Integer1DArray ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define Integer1DArray_Data( self ) ( &((self)->data[(self)->offset]) )

/* . An item. */
# define Integer1DArray_Item( self, i ) ( (self)->data[(self)->offset+(i)*(self)->stride] )

/* . A pointer to an item. */
# define Integer1DArray_ItemPointer( self, i ) ( &((self)->data[(self)->offset+(i)*(self)->stride]) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void            Integer1DArray_AddScaledArray     (       Integer1DArray  *self, const Integer value, const Integer1DArray *other, Status *status ) ;
extern Integer1DArray *Integer1DArray_Allocate           ( const Integer length, Status *status ) ;
extern Integer1DArray *Integer1DArray_Clone              ( const Integer1DArray  *self, Status *status ) ;
extern void            Integer1DArray_CopyTo             ( const Integer1DArray  *self, Integer1DArray *other, Status *status ) ;
extern void            Integer1DArray_Deallocate         (       Integer1DArray **self ) ;
extern Integer         Integer1DArray_Dot                ( const Integer1DArray  *self, const Integer1DArray *other, Status *status ) ;
extern Integer         Integer1DArray_GetItem            ( const Integer1DArray  *self, const Integer i, Status *status ) ;
extern void            Integer1DArray_LeftCircularShift  (       Integer1DArray  *self ) ;
extern Integer         Integer1DArray_Length             ( const Integer1DArray  *self ) ;
extern Integer         Integer1DArray_Maximum            ( const Integer1DArray  *self ) ;
extern void            Integer1DArray_Print              ( const Integer1DArray  *self ) ;
extern void            Integer1DArray_Resize             (       Integer1DArray  *self, const Integer length, const Integer *initializer, Status *status ) ;
extern void            Integer1DArray_Reverse            (       Integer1DArray  *self ) ;
extern void            Integer1DArray_RightCircularShift (       Integer1DArray  *self ) ;
extern void            Integer1DArray_Set                (       Integer1DArray  *self, Integer value ) ;
extern void            Integer1DArray_SetItem            (       Integer1DArray  *self, const Integer i, const Integer value, Status *status ) ;
extern void            Integer1DArray_Slice              ( const Integer1DArray  *self, const Integer start, const Integer stop, const Integer stride, Integer1DArray *slice, Status *status ) ;
extern void            Integer1DArray_Sort               (       Integer1DArray  *self ) ;
extern void            Integer1DArray_SortIndex          ( const Integer1DArray  *self, Integer1DArray *indices, Status *status ) ;
extern void            Integer1DArray_SortUnique         (       Integer1DArray  *self ) ;
extern Integer         Integer1DArray_Sum                ( const Integer1DArray  *self ) ;
extern void            Integer1DArray_ViewOfRaw          (       Integer1DArray  *self, const Integer offset, const Integer length, const Integer stride, Integer *data, const Integer size, Status *status ) ;

# endif
