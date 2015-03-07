/*------------------------------------------------------------------------------
! . File      : IntegerNDArray.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _INTEGERNDARRAY
# define _INTEGERNDARRAY

# include "Boolean.h"
# include "Integer.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The N-D array type. */
typedef struct {
    Boolean  isOwner     ;
    Boolean  isView      ;
    Integer  length      ;
    Integer  ndimensions ;
    Integer  offset      ;
    Integer  size        ;
    Integer *lengths     ;
    Integer *strides     ;
    Integer *data        ;
} IntegerNDArray ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
! . These are both general (any array dimension) and specific.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define IntegerNDArray_Data( self ) ( &((self)->data[(self)->offset]) )

/* . Items. */
/* . These can be used on arrays of higher dimension as all higher arrays indices will have the value 0. */
# define IntegerNDArray_Item1D( self, i )          ( (self)->data[(self)->offset+(i)*(self)->strides[0]] )
# define IntegerNDArray_Item2D( self, i, j )       ( (self)->data[(self)->offset+(i)*(self)->strides[0]+(j)*(self)->strides[1]] )
# define IntegerNDArray_Item3D( self, i, j, k )    ( (self)->data[(self)->offset+(i)*(self)->strides[0]+(j)*(self)->strides[1]+(k)*(self)->strides[2]] )
# define IntegerNDArray_Item4D( self, i, j, k, l ) ( (self)->data[(self)->offset+(i)*(self)->strides[0]+(j)*(self)->strides[1]+(k)*(self)->strides[2]+(l)*(self)->strides[3]] )

/* . Pointers to items. */
# define IntegerNDArray_ItemPointer1D( self, i )          ( &((self)->data[(self)->offset+(i)*(self)->strides[0]]) )
# define IntegerNDArray_ItemPointer2D( self, i, j )       ( &((self)->data[(self)->offset+(i)*(self)->strides[0]+(j)*(self)->strides[1]]) )
# define IntegerNDArray_ItemPointer3D( self, i, j, k )    ( &((self)->data[(self)->offset+(i)*(self)->strides[0]+(j)*(self)->strides[1]+(k)*(self)->strides[2]]) )
# define IntegerNDArray_ItemPointer4D( self, i, j, k, l ) ( &((self)->data[(self)->offset+(i)*(self)->strides[0]+(j)*(self)->strides[1]+(k)*(self)->strides[2]+(l)*(self)->strides[3]]) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern IntegerNDArray *IntegerNDArray_Allocate           ( const Integer ndimensions, const Integer *lengths, Status *status ) ;
extern IntegerNDArray *IntegerNDArray_Clone              ( const IntegerNDArray  *self, Status *status ) ;
extern void            IntegerNDArray_CopyTo             ( const IntegerNDArray  *self, IntegerNDArray *other, Status *status ) ;
extern void            IntegerNDArray_Deallocate         (       IntegerNDArray **self ) ;
extern Integer         IntegerNDArray_Length             ( const IntegerNDArray  *self, const Integer dimension ) ;
extern Integer         IntegerNDArray_NumberOfDimensions ( const IntegerNDArray  *self ) ;
extern void            IntegerNDArray_Resize             (       IntegerNDArray  *self, const Integer length0, const Integer *initializer, Status *status ) ;
extern void            IntegerNDArray_Set                (       IntegerNDArray  *self, Integer value ) ;

# endif
