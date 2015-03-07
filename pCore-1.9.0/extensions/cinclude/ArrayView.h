/*------------------------------------------------------------------------------
! . File      : ArrayView.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _ARRAYVIEW
# define _ARRAYVIEW

# include "Boolean.h"
# include "Integer.h"
# include "SliceOld.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The array view type. */
typedef struct {
    Integer  offset  ;
    Integer  rank    ;
    Integer  size    ;
    Integer *extents ;
    Integer *strides ;
} ArrayView ;

/* . The array view item iterator type. */
typedef struct {
    Integer    current ;
    Integer   *indices ;
    ArrayView *target  ;
} ArrayViewItemIterator ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Items. */
/* . These can be used on arrays of higher dimension as all higher arrays indices will have the value 0. */
# define ArrayView_ItemIndex1D( self, i )          ( (self)->offset+(i)*(self)->strides[0] )
# define ArrayView_ItemIndex2D( self, i, j )       ( (self)->offset+(i)*(self)->strides[0]+(j)*(self)->strides[1] )
# define ArrayView_ItemIndex3D( self, i, j, k )    ( (self)->offset+(i)*(self)->strides[0]+(j)*(self)->strides[1]+(k)*(self)->strides[2] )
# define ArrayView_ItemIndex4D( self, i, j, k, l ) ( (self)->offset+(i)*(self)->strides[0]+(j)*(self)->strides[1]+(k)*(self)->strides[2]+(l)*(self)->strides[3] )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern ArrayView *ArrayView_Allocate            ( const Integer rank, Status *status ) ;
extern ArrayView *ArrayView_AllocateRaw         ( Status *status ) ;
extern ArrayView *ArrayView_AllocateWithExtents ( const Integer rank, const Integer *extents, Status *status ) ;
extern Boolean    ArrayView_CheckCapacity       ( const ArrayView  *self, const Integer capacity ) ;
extern ArrayView *ArrayView_Clone               ( const ArrayView  *self, Status *status ) ;
extern void       ArrayView_Deallocate          (       ArrayView **self ) ;
extern Integer    ArrayView_Extent              ( const ArrayView  *self, const Integer dimension, Status *status ) ;
extern Boolean    ArrayView_IsCompact           ( const ArrayView  *self ) ;
extern Boolean    ArrayView_IsConformable       ( const ArrayView  *self, const ArrayView *other ) ;
extern Boolean    ArrayView_IsUniform           ( const ArrayView  *self ) ;
extern Integer    ArrayView_ItemIndex           ( const ArrayView  *self, const Integer rank, const Integer *indices ) ;
extern Integer    ArrayView_Rank                ( const ArrayView  *self ) ;
extern Boolean    ArrayView_Reshape1D           ( const ArrayView  *self, Integer *offset, Integer *extent, Integer *stride, Status *status ) ;
extern Integer    ArrayView_Size                ( const ArrayView  *self ) ;
extern Status     ArrayView_Slice               ( const ArrayView  *self, const MultiSliceX *multiSlice, ArrayView **view ) ;

/* . Iterator functions. */
extern ArrayViewItemIterator *ArrayViewItemIterator_Allocate        ( ArrayView *target, Status *status ) ;
extern void                   ArrayViewItemIterator_Deallocate      ( ArrayViewItemIterator **self ) ;
extern void                   ArrayViewItemIterator_Initialize      ( ArrayViewItemIterator  *self ) ;
extern Integer                ArrayViewItemIterator_Next            ( ArrayViewItemIterator  *self ) ;
extern Integer                ArrayViewItemIterator_NextWithIndices ( ArrayViewItemIterator  *self, Integer *indices ) ;

# endif
