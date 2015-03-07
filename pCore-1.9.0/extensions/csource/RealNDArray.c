/*------------------------------------------------------------------------------
! . File      : RealNDArray.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . N-D arrays.
!=================================================================================================================================*/

# include "Macros.h"
# include "Memory.h"
# include "RealNDArray.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RealNDArray *RealNDArray_Allocate ( const Integer rank, const Integer *extents, Status *status )
{
    ArrayView   *view = NULL ;
    RealNDArray *self = NULL ;
    /* . Basic allocation. */
    self = RealNDArray_AllocateRaw ( status ) ;
    view = ArrayView_AllocateWithExtents ( rank, extents, status ) ;
    if ( ( self == NULL ) || ( view == NULL ) )
    {
        ArrayView_Deallocate   ( &view ) ;
        RealNDArray_Deallocate ( &self ) ;
    }
    else
    {
        self->view = view ;
        /* . Data allocation. */
        if ( view->size > 0 )
        {
            MEMORY_ALLOCATEARRAY ( self->data, view->size, Real ) ;
            if ( self->data == NULL )
            {
                RealNDArray_Deallocate ( &self ) ;
                Status_Set ( status, Status_MemoryAllocationFailure ) ;
            }
            else self->capacity = view->size ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Raw allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RealNDArray *RealNDArray_AllocateRaw ( Status *status )
{
    RealNDArray *self = NULL ;
    MEMORY_ALLOCATE ( self, RealNDArray ) ;
    if ( self != NULL )
    {
        self->capacity = 0 ;
        self->data     = NULL  ;
        self->isOwner  = True  ;
        self->isView   = False ;
        self->view     = NULL  ;
    }
    else Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
RealNDArray *RealNDArray_Clone ( const RealNDArray *self, const Boolean doShallow, Status *status )
{
    RealNDArray *clone = NULL ;
    if ( self != NULL )
    {
        /* . Shallow cloning does not clone non-owner data. */
        if ( doShallow && ( ! self->isOwner ) )
        {
            ArrayView *view ;
            clone = RealNDArray_AllocateRaw ( status ) ;
            view  = ArrayView_Clone ( self->view, status ) ;
            if ( ( clone == NULL ) || ( view == NULL ) )
            {
                ArrayView_Deallocate   ( &view  ) ;
                RealNDArray_Deallocate ( &clone ) ;
            }
            else
            {
                clone->capacity = self->capacity ;
                clone->data     = self->data     ;
                clone->isOwner  = False          ;
                clone->isView   = self->isView   ;
                clone->view     = view           ;
            }
        }
        /* . Other cases. */
        else
        {
            clone = RealNDArray_Allocate ( self->view->rank, self->view->extents, status ) ;
            RealNDArray_CopyTo ( self, clone, status ) ;
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealNDArray_CopyTo ( const RealNDArray *self, RealNDArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( ArrayView_IsConformable ( self->view, other->view ) )
        {
            auto Integer m, n ;
            auto ArrayViewItemIterator *from, *to ;
            from = ArrayViewItemIterator_Allocate ( self->view , NULL ) ;
            to   = ArrayViewItemIterator_Allocate ( other->view, NULL ) ;
            ArrayViewItemIterator_Initialize ( from ) ;
            ArrayViewItemIterator_Initialize ( to   ) ;
            while ( ( m = ArrayViewItemIterator_Next ( from ) ) >= 0 )
            {
                n = ArrayViewItemIterator_Next ( to ) ;
                other->data[n] = self->data[m] ;
            }
            ArrayViewItemIterator_Deallocate ( &from ) ;
            ArrayViewItemIterator_Deallocate ( &to   ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealNDArray_Deallocate ( RealNDArray **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        if ( (*self)->isOwner ) MEMORY_DEALLOCATE ( (*self)->data ) ;
        MEMORY_DEALLOCATE ( (*self)->view ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Extent along a dimension.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RealNDArray_Extent ( const RealNDArray *self, const Integer dimension, Status *status )
{
    Integer extent = 0 ;
    if ( self != NULL ) extent = ArrayView_Extent ( self->view, dimension, status ) ;
    return extent ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real RealNDArray_GetItem ( RealNDArray *self, const Integer *indices, Status *status )
{
    Real value = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer n ;
        n = ArrayView_ItemIndex ( self->view, self->view->rank, indices ) ;
        if ( n < 0 ) Status_Set ( status, Status_IndexOutOfRange ) ;
        else value = self->data[n] ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for compactness (minimal spacing of items).
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean RealNDArray_IsCompact ( const RealNDArray *self )
{
    Boolean isCompact = False ;
    if ( self != NULL ) isCompact = ArrayView_IsCompact ( self->view ) ;
    return isCompact ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for uniformity (equal spacing of items).
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean RealNDArray_IsUniform ( const RealNDArray *self )
{
    Boolean isUniform = False ;
    if ( self != NULL ) isUniform = ArrayView_IsUniform ( self->view ) ;
    return isUniform ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Rank.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RealNDArray_Rank ( const RealNDArray *self )
{
    Integer rank = 0 ;
    if ( self != NULL ) rank = ArrayView_Rank ( self->view ) ;
    return rank ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a 1-D reshaped view of the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealNDArray_Reshape1D ( const RealNDArray *self, Real1DArray *view, Status *status )
{
    if ( ( self != NULL ) && ( view != NULL ) )
    {
        auto Integer extent, offset, stride ;
        if ( ArrayView_Reshape1D ( self->view, &offset, &extent, &stride, status ) )
        {
            /* . View. */
            view->length  = extent ;
            view->offset  = offset ;
            view->stride  = stride ;
            /* . Data. */
            view->isOwner = False          ;
            view->isView  = True           ;
            view->size    = self->capacity ;
            view->data    = self->data     ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealNDArray_Scale ( RealNDArray *self, Real value )
{
    if ( ( self != NULL ) && ( self->view->size > 0 ) )
    {
        auto Integer n ;
        auto ArrayViewItemIterator *iterator ;
        iterator = ArrayViewItemIterator_Allocate ( self->view, NULL ) ;
        ArrayViewItemIterator_Initialize ( iterator ) ;
        while ( ( n = ArrayViewItemIterator_Next ( iterator ) ) >= 0 ) self->data[n] *= value ;
        ArrayViewItemIterator_Deallocate ( &iterator ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealNDArray_Set ( RealNDArray *self, Real value )
{
    if ( ( self != NULL ) && ( self->view->size > 0 ) )
    {
        auto Integer n ;
        auto ArrayViewItemIterator *iterator ;
        iterator = ArrayViewItemIterator_Allocate ( self->view, NULL ) ;
        ArrayViewItemIterator_Initialize ( iterator ) ;
        while ( ( n = ArrayViewItemIterator_Next ( iterator ) ) >= 0 ) self->data[n] = value ;
        ArrayViewItemIterator_Deallocate ( &iterator ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealNDArray_SetItem ( RealNDArray *self, const Integer *indices, const Real value, Status *status )
{
    if ( self != NULL )
    {
        Integer n ;
        n = ArrayView_ItemIndex ( self->view, self->view->rank, indices ) ;
        if ( n < 0 ) Status_Set ( status, Status_IndexOutOfRange ) ;
        else self->data[n] = value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Size.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RealNDArray_Size ( const RealNDArray *self )
{
    Integer size = 0 ;
    if ( self != NULL ) size = ArrayView_Size ( self->view ) ;
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Slice - 1-D.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealNDArray_Slice1D ( const RealNDArray *self, const ArrayView *view, Real1DArray *slice, Status *status )
{
    if ( ( self != NULL ) && ( view != NULL ) && ( slice != NULL ) )
    {
        if ( view->rank == 1 ) Real1DArray_ViewOfRaw ( slice, view->offset, view->extents[0], view->strides[0],
                                                                         self->data, self->capacity, status ) ;
        else Status_Set ( status, Status_InvalidRank ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Slice - 2-D.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealNDArray_Slice2D ( const RealNDArray *self, const ArrayView *view, Real2DArray *slice, Status *status )
{
    if ( ( self != NULL ) && ( view != NULL ) && ( slice != NULL ) )
    {
        if ( view->rank == 2 ) Real2DArray_ViewOfRaw ( slice, view->offset, view->extents[0], view->strides[0],
                                                                            view->extents[1], view->strides[1], 
                                                                         self->data, self->capacity, status ) ;
        else Status_Set ( status, Status_InvalidRank ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Slice - N-D.
!---------------------------------------------------------------------------------------------------------------------------------*/
RealNDArray *RealNDArray_SliceND ( const RealNDArray *self, ArrayView *view, Status *status )
{
    RealNDArray *slice = NULL ;
    if ( ( self != NULL ) && ( view != NULL ) ) slice = RealNDArray_ViewOfRaw ( view, self->data, self->capacity, status ) ;
    return slice ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a 2-D slice of the last two dimensions of given index.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealNDArray_TailSlice2D ( const RealNDArray *self, Integer *indices, Real2DArray *slice, Status *status )
{
    if ( ( self != NULL ) && ( indices != NULL ) && ( slice != NULL ) )
    {
        auto ArrayView *view ;
        auto Integer    n ;

        /* . Check the indices. */
        view = self->view ;
        n    = ArrayView_ItemIndex ( view, view->rank-2, indices ) ;

        /* . OK. */
        if ( n >= 0 )
        {
            auto Integer a, b ;
            a = view->rank - 2 ;
            b = view->rank - 1 ;
            slice->data    = self->data ;
            slice->isOwner = False ;
            slice->isView  = True  ;
            slice->length  = view->extents[a] * view->extents[b] ;
            slice->length0 = view->extents[a] ;
            slice->length1 = view->extents[b] ;
            slice->offset  = n ;
            slice->size    = self->capacity ;
            slice->stride0 = view->strides[a] ;
            slice->stride1 = view->strides[b] ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create an array view of raw data.
!---------------------------------------------------------------------------------------------------------------------------------*/
RealNDArray *RealNDArray_ViewOfRaw ( ArrayView *view, Real *data, const Integer capacity, Status *status )
{
    RealNDArray *self = NULL ;
    if ( ( view != NULL ) && ( data != NULL ) )
    {
        if ( ArrayView_CheckCapacity ( view, capacity ) )
        {
            self = RealNDArray_AllocateRaw ( status ) ;
            if ( self != NULL )
            {
                self->capacity = capacity ;
                self->data     = data     ;
                self->isOwner  = False    ;
                self->isView   = True     ;
                self->view     = view     ;
            }
            else Status_Set ( status, Status_MemoryAllocationFailure ) ;
        }
        else Status_Set ( status, Status_ArrayOutOfBounds ) ;
    }
    return self ;
}
