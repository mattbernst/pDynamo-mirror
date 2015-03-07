/*------------------------------------------------------------------------------
! . File      : RealNDArray.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _REALNDARRAY
# define _REALNDARRAY

# include "ArrayView.h"
# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    Boolean    isOwner  ;
    Boolean    isView   ;
    Integer    capacity ;
    ArrayView *view     ;
    Real      *data     ;
} RealNDArray ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
! . These are both general (any array dimension) and specific.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define RealNDArray_Data( self ) ( &((self)->data[(self)->view->offset]) )

/* . Items. */
/* . These can be used on arrays of higher dimension as all higher arrays indices will have the value 0. */
# define RealNDArray_Item1D( self, i          ) ( (self)->data[ ArrayView_ItemIndex1D ( (self)->view, i          ) ] )
# define RealNDArray_Item2D( self, i, j       ) ( (self)->data[ ArrayView_ItemIndex2D ( (self)->view, i, j       ) ] )
# define RealNDArray_Item3D( self, i, j, k    ) ( (self)->data[ ArrayView_ItemIndex3D ( (self)->view, i, j, k    ) ] )
# define RealNDArray_Item4D( self, i, j, k, l ) ( (self)->data[ ArrayView_ItemIndex4D ( (self)->view, i, j, k, l ) ] )

/* . Pointers to items. */
# define RealNDArray_ItemPointer1D( self, i          ) ( &( RealNDArray_Item1D( self, i          ) ) )
# define RealNDArray_ItemPointer2D( self, i, j       ) ( &( RealNDArray_Item2D( self, i, j       ) ) )
# define RealNDArray_ItemPointer3D( self, i, j, k    ) ( &( RealNDArray_Item3D( self, i, j, k    ) ) )
# define RealNDArray_ItemPointer4D( self, i, j, k, l ) ( &( RealNDArray_Item4D( self, i, j, k, l ) ) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern RealNDArray *RealNDArray_Allocate    ( const Integer rank, const Integer *extents, Status *status ) ;
extern RealNDArray *RealNDArray_AllocateRaw ( Status *status ) ;
extern RealNDArray *RealNDArray_Clone       ( const RealNDArray  *self, const Boolean doShallow, Status *status ) ;
extern void         RealNDArray_CopyTo      ( const RealNDArray  *self, RealNDArray *other, Status *status ) ;
extern void         RealNDArray_Deallocate  (       RealNDArray **self ) ;
extern Integer      RealNDArray_Extent      ( const RealNDArray  *self, const Integer dimension, Status *status ) ;
extern Real         RealNDArray_GetItem     (       RealNDArray  *self, const Integer *indices, Status *status ) ;
extern Boolean      RealNDArray_IsCompact   ( const RealNDArray  *self ) ;
extern Boolean      RealNDArray_IsUniform   ( const RealNDArray  *self ) ;
extern Integer      RealNDArray_Rank        ( const RealNDArray  *self ) ;
extern void         RealNDArray_Reshape1D   ( const RealNDArray  *self, Real1DArray *view, Status *status ) ;
extern void         RealNDArray_Scale       (       RealNDArray  *self, Real value ) ;
extern void         RealNDArray_Set         (       RealNDArray  *self, Real value ) ;
extern void         RealNDArray_SetItem     (       RealNDArray  *self, const Integer *indices, const Real value, Status *status ) ;
extern Integer      RealNDArray_Size        ( const RealNDArray  *self ) ;
extern void         RealNDArray_Slice1D     ( const RealNDArray  *self, const ArrayView *view, Real1DArray *slice, Status *status ) ;
extern void         RealNDArray_Slice2D     ( const RealNDArray  *self, const ArrayView *view, Real2DArray *slice, Status *status ) ;
extern RealNDArray *RealNDArray_SliceND     ( const RealNDArray  *self,       ArrayView *view, Status *status ) ;
extern void         RealNDArray_TailSlice2D ( const RealNDArray  *self, Integer *indices, Real2DArray *slice, Status *status ) ;
extern RealNDArray *RealNDArray_ViewOfRaw   ( ArrayView *view, Real *data, const Integer capacity, Status *status ) ;

# endif
