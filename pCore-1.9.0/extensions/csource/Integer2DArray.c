/*------------------------------------------------------------------------------
! . File      : Integer2DArray.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . 2-D arrays.
!=================================================================================================================================*/

# include "Integer2DArray.h"
# include "Macros.h"
# include "Memory.h"
# include "SliceOld.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a 1D slice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_1DSlice ( const Integer2DArray *self, const Integer start0, const Integer stop0, const Integer stride0, const Integer start1, const Integer stop1, const Integer stride1, Integer1DArray *slice, Status *status )
{
    if ( ( self != NULL ) && ( slice != NULL ) )
    {
        auto Integer length0, length1 ;
        auto Status  localstatus ;
        length0 = SliceIndices_GetLength ( start0, stop0, stride0, self->length0, &localstatus ) ;
        length1 = SliceIndices_GetLength ( start1, stop1, stride1, self->length1, &localstatus ) ;
        if ( ! Status_OK ( &localstatus ) ) Status_Set ( status, localstatus ) ;
        else if ( ( length0 != 1 ) && ( length1 != 1 ) ) Status_Set ( status, Status_InvalidSliceIndices ) ;
        else
        {
            slice->data    = self->data ;
            slice->isOwner = False ;
            slice->isView  = True  ;
            slice->length  = length0 * length1 ;
            slice->offset  = Integer2DArray_ItemIndex ( self, start0, start1 ) ;
            slice->size    = self->size ;
            if ( length0 == 1 ) slice->stride = self->stride1 * stride1 ;
            else                slice->stride = self->stride0 * stride0 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer2DArray *Integer2DArray_Allocate ( const Integer length0, const Integer length1, Status *status )
{
    Integer2DArray *self = NULL ;
    MEMORY_ALLOCATE ( self, Integer2DArray ) ;
    if ( self != NULL )
    {
        self->data    = NULL  ;
        self->isOwner = True  ;
        self->isView  = False ;
        self->length0 = Maximum ( length0, 0 ) ;
        self->length1 = Maximum ( length1, 0 ) ;
        self->length  = self->length0 * self->length1 ;
        self->offset  = 0 ;
        self->size    = self->length  ;
        self->stride0 = self->length1 ;
        self->stride1 = 1 ;
        if ( self->length > 0 )
        {
            MEMORY_ALLOCATEARRAY ( self->data, self->length, Integer ) ;
            if ( self->data == NULL ) Integer2DArray_Deallocate ( &self ) ;
        }
        else if ( self->length < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
! . Data is copied and kept in the most compact representation possible.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer2DArray *Integer2DArray_Clone ( const Integer2DArray *self, Status *status )
{
    Integer2DArray *clone = NULL ;
    if ( self != NULL )
    {
        clone = Integer2DArray_Allocate ( self->length0, self->length1, status ) ;
        Integer2DArray_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a slice of a column.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_ColumnSlice ( const Integer2DArray  *self, const Integer column, Integer1DArray *slice, Status *status )
{
    if ( ( self != NULL ) && ( slice != NULL ) )
    {
        if ( ( column < 0 ) || ( column >= self->length1 ) ) Status_Set ( status, Status_IndexOutOfRange ) ;
        else
        {
            slice->data    = self->data ;
            slice->isOwner = False ;
            slice->isView  = True  ;
            slice->length  = self->length0 ;
            slice->offset  = self->offset + column * self->stride1 ;
            slice->size    = self->size    ;
            slice->stride  = self->stride0 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
! . As much data as possible is copied.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_CopyTo ( const Integer2DArray *self, Integer2DArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer n0, n1 ;
        if ( self->length0 == other->length0 ) n0 = self->length0 ;
        else
        {
            Status_Set ( status, Status_ArrayLengthMismatch ) ;
            n0 = Minimum ( self->length0, other->length0 ) ;
        }
        if ( self->length1 == other->length1 ) n1 = self->length1 ;
        else
        {
            Status_Set ( status, Status_ArrayLengthMismatch ) ;
            n1 = Minimum ( self->length1, other->length1 ) ;
        }
        if ( ( n0 > 0 ) && ( n1 > 0 ) )
        {
            auto Integer i, j ;
            for ( i = 0 ; i < n0 ; i++ )
            {
                for ( j = 0 ; j < n1 ; j++ ) Integer2DArray_Item ( other, i, j ) = Integer2DArray_Item ( self, i, j ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_Deallocate ( Integer2DArray **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        if ( (*self)->isOwner ) MEMORY_DEALLOCATE ( (*self)->data ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Integer2DArray_GetItem ( const Integer2DArray *self, const Integer i, const Integer j, Status *status )
{
    Integer value = 0 ;
    if ( self != NULL )
    {
        if ( self->length > 0 )
        {
            if ( ( i >= 0 ) && ( i < self->length0 ) && ( j >= 0 ) && ( j < self->length1 ) ) value = Integer2DArray_Item ( self, i, j ) ;
            else Status_Set ( status, Status_IndexOutOfRange ) ;
        }
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Integer2DArray_Length ( const Integer2DArray *self, const Integer dimension )
{
    Integer length = 0 ;
    if ( self != NULL )
    {
        switch ( dimension )
        {
            case 0  : length = self->length0 ; break ;
            case 1  : length = self->length1 ; break ;
            default : length = self->length  ; break ;
        }
    }
    return length ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
! . This is only done for the slowest changing dimension.
! . This can only be done if the array owns the data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_Resize ( Integer2DArray *self, const Integer length0, const Integer *initializer, Status *status )
{
    if ( self != NULL )
    {
        if ( self->isOwner )
        {
            /* . Only do something if the new and old lengths differ. */
            if ( length0 != self->length0 )
            {
                if ( length0 > 0 )
                {
                    auto Integer n ;
                    auto Integer   *new = NULL ;
                    /* . The minimal number of items needed is, in fact, offset + length0 * stride0 - ( stride0 - 1 ).
                    ! . However, the longer, stride-consistent, value is left here for the moment just in case. */
                    n = self->offset + length0 * self->stride0 ;
                    MEMORY_REALLOCATEARRAY ( new, self->data, n, Integer ) ;
                    /* . Failure so leave as is. */
                    if ( new == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
                    /* . OK. */
                    else
                    {
                        auto Integer oldlength0 ;
                        /* . Save the old length. */
                        oldlength0 = self->length0 ;
                        /* . Reset data values. */
                        self->data    = new     ;
                        self->length0 = length0 ;
                        self->length  = self->length0 * self->length1 ;
                        self->size    = n       ;
                        /* . Initialize. */
                        if ( ( length0 > oldlength0 ) && ( initializer != NULL ) )
                        {
                            auto Integer i, j ;
                            for ( i = oldlength0 ; i < length0 ; i++ )
                            {
                                for ( j = 0 ; j < self->length1 ; j++ ) Integer2DArray_Item ( self, i, j ) = (*initializer) ;
                            }
                        }
                    }
                }
                else
                {
                    MEMORY_DEALLOCATE ( self->data ) ;
                    self->data    = NULL  ;
                    self->isOwner = True  ;
                    self->isView  = False ;
                    self->length  = 0     ;
                    self->length0 = 0     ;
                    self->size    = 0     ;
                    if ( length0 < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
                }
            }
        }
        else Status_Set ( status, Status_InvalidArrayOperation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a slice of a row.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_RowSlice ( const Integer2DArray  *self, const Integer row, Integer1DArray *slice, Status *status )
{
    if ( ( self != NULL ) && ( slice != NULL ) )
    {
        if ( ( row < 0 ) || ( row >= self->length0 ) ) Status_Set ( status, Status_IndexOutOfRange ) ;
        else
        {
            slice->data    = self->data ;
            slice->isOwner = False ;
            slice->isView  = True  ;
            slice->length  = self->length1 ;
            slice->offset  = self->offset + row * self->stride0 ;
            slice->size    = self->size    ;
            slice->stride  = self->stride1 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_Set ( Integer2DArray *self, Integer value )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i, j ;
        for ( i = 0 ; i < self->length0 ; i++ )
        {
            for ( j = 0 ; j < self->length1 ; j++ ) Integer2DArray_Item ( self, i, j ) = value ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_SetItem ( Integer2DArray *self, const Integer i, const Integer j, const Integer value, Status *status )
{
    if ( self != NULL )
    {
        if ( self->length > 0 )
        {
            if ( ( i >= 0 ) && ( i < self->length0 ) && ( j >= 0 ) && ( j < self->length1 ) ) Integer2DArray_Item ( self, i, j ) = value ;
            else Status_Set ( status, Status_IndexOutOfRange ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a 2D slice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_Slice ( const Integer2DArray *self, const Integer start0, const Integer stop0, const Integer stride0, const Integer start1, const Integer stop1, const Integer stride1, Integer2DArray *slice, Status *status )
{
    if ( ( self != NULL ) && ( slice != NULL ) )
    {
        auto Integer length0, length1 ;
        auto Status  localstatus ;
        length0 = SliceIndices_GetLength ( start0, stop0, stride0, self->length0, &localstatus ) ;
        length1 = SliceIndices_GetLength ( start1, stop1, stride1, self->length1, &localstatus ) ;
        if ( ! Status_OK ( &localstatus ) ) Status_Set ( status, localstatus ) ;
        else
        {
            slice->data    = self->data ;
            slice->isOwner = False ;
            slice->isView  = True  ;
            slice->length  = length0 * length1 ;
            slice->length0 = length0 ;
            slice->length1 = length1 ;
            slice->offset  = Integer2DArray_ItemIndex ( self, start0, start1 ) ;
            slice->size    = self->size ;
            slice->stride0 = self->stride0 * stride0 ;
            slice->stride1 = self->stride1 * stride1 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create an array view of a raw array.
! . The input offset, length and stride values are final.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer2DArray_ViewOfRaw ( Integer2DArray *self, const Integer offset, const Integer length0, const Integer stride0, const Integer length1, const Integer stride1, Integer *data, const Integer size, Status *status )
{
    if ( ( self != NULL ) && ( data != NULL ) )
    {
        auto Integer n ;
        /* . Find space required. */
        if ( length0 * length1 == 0 ) n = 0 ;
        else
        {
            n = offset + 1 ;
            if ( stride0 > 0 ) n += ( length0 - 1 ) * stride0 ;
            if ( stride1 > 0 ) n += ( length1 - 1 ) * stride1 ;
        }
        if ( n > size ) Status_Set ( status, Status_ArrayOverFlow ) ;
        else
        {
            self->data    = data    ;
            self->isOwner = False   ;
            self->isView  = True    ;
            self->length  = length0 * length1 ;
            self->length0 = length0 ;
            self->length1 = length1 ;
            self->offset  = offset  ;
            self->size    = size    ;
            self->stride0 = stride0 ;
            self->stride1 = stride1 ;
        }
    }
}
