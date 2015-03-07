/*------------------------------------------------------------------------------
! . File      : Integer1DArray.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . 1-D arrays.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "Integer1DArray.h"
# include "Macros.h"
# include "Memory.h"
# include "SliceOld.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer Item_Compare ( const void *vterm1, const void *vterm2 ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add a scaled array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_AddScaledArray ( Integer1DArray *self, const Integer value, const Integer1DArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && ( value != 0 ) )
    {
        auto Integer n ;
        if ( self->length == other->length ) n = self->length ;
        else
        {
            Status_Set ( status, Status_ArrayLengthMismatch ) ;
            n = Minimum ( self->length, other->length ) ;
        }
        if ( n > 0 )
        {
            auto Integer i ;
            for ( i = 0 ; i < n ; i++ ) Integer1DArray_Item ( self, i ) += value * Integer1DArray_Item ( other, i ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer1DArray *Integer1DArray_Allocate ( const Integer length, Status *status )
{
    Integer1DArray *self = NULL ;
    MEMORY_ALLOCATE ( self, Integer1DArray ) ;
    if ( self != NULL )
    {
        self->data    = NULL   ;
        self->isOwner = True   ;
        self->isView  = False  ;
        self->length  = Maximum ( length, 0 ) ;
        self->offset  = 0      ;
        self->size    = Maximum ( length, 0 ) ;
        self->stride  = 1      ;
        if ( length > 0 )
        {
            MEMORY_ALLOCATEARRAY ( self->data, length, Integer ) ;
            if ( self->data == NULL ) Integer1DArray_Deallocate ( &self ) ;
        }
        else if ( length < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
! . Data is copied and kept in the most compact representation possible.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer1DArray *Integer1DArray_Clone ( const Integer1DArray *self, Status *status )
{
    Integer1DArray *clone = NULL ;
    if ( self != NULL )
    {
        clone = Integer1DArray_Allocate ( self->length, status ) ;
        Integer1DArray_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
! . As much data as possible is copied.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_CopyTo ( const Integer1DArray *self, Integer1DArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer n ;
        if ( self->length == other->length ) n = self->length ;
        else
        {
            Status_Set ( status, Status_ArrayLengthMismatch ) ;
            n = Minimum ( self->length, other->length ) ;
        }
        if ( n > 0 )
        {
            auto Integer i ;
            for ( i = 0 ; i < n ; i++ ) Integer1DArray_Item ( other, i ) = Integer1DArray_Item ( self, i ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_Deallocate ( Integer1DArray **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        if ( (*self)->isOwner ) MEMORY_DEALLOCATE ( (*self)->data ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dot product.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Integer1DArray_Dot ( const Integer1DArray *self, const Integer1DArray *other, Status *status )
{
    Integer dot = 0 ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer n ;
        if ( self->length == other->length ) n = self->length ;
        else
        {
            Status_Set ( status, Status_ArrayLengthMismatch ) ;
            n = Minimum ( self->length, other->length ) ;
        }
        if ( n > 0 )
        {
            auto Integer i ;
            for ( i = 0 ; i < n ; i++ ) dot += Integer1DArray_Item ( other, i ) * Integer1DArray_Item ( self, i ) ;
        }
    }
    return dot ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Integer1DArray_GetItem ( const Integer1DArray *self, const Integer i, Status *status )
{
    Integer value = 0.0e+00 ;
    if ( self != NULL )
    {
        if ( self->length > 0 )
        {
            if ( ( i >= 0 ) && ( i < self->length ) ) value = Integer1DArray_Item ( self, i ) ;
            else Status_SafeSet ( status, Status_IndexOutOfRange ) ;
        }
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Left circular shift.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_LeftCircularShift ( Integer1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        auto Integer first, i ;
        first = Integer1DArray_Item ( self, 0 ) ;
        for ( i = 1 ; i < self->length ; i++ ) Integer1DArray_Item ( self, i-1 ) = Integer1DArray_Item ( self, i ) ;
        Integer1DArray_Item ( self, self->length-1 ) = first ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Integer1DArray_Length ( const Integer1DArray *self )
{
    if ( self == NULL ) return 0 ;
    else                return self->length ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the maximum value in the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Integer1DArray_Maximum ( const Integer1DArray *self )
{
    Integer value = 0 ; /* . Should really be largest negative number representable by Integer. */
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        value = Integer1DArray_Item ( self, 0 ) ;
        for ( i = 0 ; i < self->length ; i++ ) value = Maximum ( value, Integer1DArray_Item ( self, i ) ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_Print ( const Integer1DArray *self )
{
    if ( self == NULL ) printf ( "Null integer 1-D array.\n" ) ;
    else
    {
        auto Integer i ;
        for ( i = 0 ; i < self->length ; i++ ) printf ( "%6d", Integer1DArray_Item ( self, i ) ) ;
        printf ( "\n" ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
! . This can only be done if the array owns the data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_Resize ( Integer1DArray *self, const Integer length, const Integer *initializer, Status *status )
{
    if ( self != NULL )
    {
        if ( self->isOwner )
        {
            /* . Only do something if the new and old lengths differ. */
            if ( length != self->length )
            {
                if ( length > 0 )
                {
                    auto Integer n, *new = NULL ;
                    /* . The minimal number of items needed is, in fact, offset + length * stride - ( stride - 1 ).
                    ! . However, the longer, stride-consistent, value is left here for the moment just in case. */
                    n = self->offset + length * self->stride ;
                    MEMORY_REALLOCATEARRAY ( new, self->data, n, Integer ) ;
                    /* . Failure so leave as is. */
                    if ( new == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
                    /* . OK. */
                    else
                    {
                        auto Integer oldlength ;
                        /* . Save the old length. */
                        oldlength = self->length ;
                        /* . Reset data values. */
                        self->data   = new    ;
                        self->length = length ;
                        self->size   = n      ;
                        /* . Initialize. */
                        if ( ( length > oldlength ) && ( initializer != NULL ) )
                        {
                            auto Integer i ;
                            for ( i = oldlength ; i < length ; i++ ) Integer1DArray_Item ( self, i ) = (*initializer) ;
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
                    self->offset  = 0     ;
                    self->size    = 0     ;
                    self->stride  = 1     ;
                    if ( length < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
                }
            }
        }
        else Status_Set ( status, Status_InvalidArrayOperation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Reverse - in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_Reverse ( Integer1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        auto Integer i, j ;
        auto Integer t ;
        for ( i = 0, j = self->length - 1 ; i < ( self->length / 2 ) ; i++, j-- )
        {
            t = Integer1DArray_Item ( self, i ) ;
            Integer1DArray_Item ( self, i ) = Integer1DArray_Item ( self, j ) ;
            Integer1DArray_Item ( self, j ) = t ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Right circular shift.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_RightCircularShift ( Integer1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        auto Integer i, last ;
        last = Integer1DArray_Item ( self, self->length-1 ) ;
        for ( i = self->length-1 ; i > 0 ; i-- ) Integer1DArray_Item ( self, i ) = Integer1DArray_Item ( self, i-1 ) ;
        Integer1DArray_Item ( self, 0 ) = last ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_Set ( Integer1DArray *self, Integer value )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->length ; i++ ) Integer1DArray_Item ( self, i ) = value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_SetItem ( Integer1DArray *self, const Integer i, const Integer value, Status *status )
{
    if ( self != NULL )
    {
        if ( self->length > 0 )
        {
            if ( ( i >= 0 ) && ( i < self->length ) ) Integer1DArray_Item ( self, i ) = value ;
            else Status_SafeSet ( status, Status_IndexOutOfRange ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Slicing - 1D.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_Slice ( const Integer1DArray  *self, const Integer start, const Integer stop, const Integer stride, Integer1DArray *slice, Status *status )
{
    if ( ( self != NULL ) && ( slice != NULL ) )
    {
        auto Integer length ;
        auto Status  localstatus ;
        length = SliceIndices_GetLength ( start, stop, stride, self->length, &localstatus ) ;
        if ( ! Status_OK ( &localstatus ) ) Status_Set ( status, localstatus ) ;
        else
        {
            slice->data    = self->data ;
            slice->isOwner = False ;
            slice->isView  = True  ;
            slice->length  = length ;
            slice->offset  = self->offset + start * self->stride ;
            slice->size    = self->size ;
            slice->stride  = self->stride * stride ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting - ascending in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_Sort ( Integer1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 1 ) ) qsort ( ( void * ) Integer1DArray_Data ( self ), ( size_t ) self->length, self->stride * sizeof ( Integer ), ( void * ) Item_Compare ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting - return an array of sorted indices. The array itself remains unchanged.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_SortIndex ( const Integer1DArray *self, Integer1DArray *indices, Status *status )
{
    if ( ( self != NULL ) && ( indices != NULL ) )
    {
        /* . Length mismatch. */
        if ( indices->length != self->length ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
        /* . Continue as long as there is space for all indices. */
        if ( ( indices->length >= self->length ) && ( self->length ) )
        {
            auto Integer j, k, ki, l, m, n, tmp ;
            n = self->length ;
            /* . Initialize indices. */
            for ( k = 0 ; k < n               ; k++ ) Integer1DArray_Item ( indices, k ) =  k ;
            for ( k = n ; k < indices->length ; k++ ) Integer1DArray_Item ( indices, k ) = -1 ;
            /* . Use heapsort. */
            m = n - 1 ;
            k = m / 2;
            k ++ ;
            do
            {
                k -- ;
                ki = Integer1DArray_Item ( indices, k ) ;
                l  = k ;
                while ( l <= ( m / 2 ) )
                {
                    j = 2 * l ;
                    if ( j < m && ( Integer1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) < Integer1DArray_Item ( self, Integer1DArray_Item ( indices, j+1 ) ) ) ) j++ ;
                    if ( Integer1DArray_Item ( self, ki ) >= Integer1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) ) break ;
                    Integer1DArray_Item ( indices, l ) = Integer1DArray_Item ( indices, j );
                    l = j;
                }
                Integer1DArray_Item ( indices, l ) = ki;
            }
            while ( k > 0 ) ;
            while ( m > 0 )
            {
                tmp = Integer1DArray_Item ( indices, 0 ) ; Integer1DArray_Item ( indices, 0 ) = Integer1DArray_Item ( indices, m ) ; Integer1DArray_Item ( indices, m ) = tmp ;
                m -- ;
                ki = Integer1DArray_Item ( indices, 0 ) ;
                l  = k ;
                while ( l <= ( m / 2 ) )
                {
                    j = 2 * l ;
                    if ( j < m && ( Integer1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) < Integer1DArray_Item ( self, Integer1DArray_Item ( indices, j+1 ) ) ) ) j++ ;
                    if ( Integer1DArray_Item ( self, ki ) >= Integer1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) ) break ;
                    Integer1DArray_Item ( indices, l ) = Integer1DArray_Item ( indices, j );
                    l = j;
                }
                Integer1DArray_Item ( indices, l ) = ki;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting of unique items only (this means duplicates are removed).
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_SortUnique ( Integer1DArray *self )
{
    Integer1DArray_Sort ( self ) ;
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        auto Integer i, n ;
        for ( i = 1, n = 1 ; i < self->length ; i++ )
        {
            if ( Integer1DArray_Item ( self, i-1 ) != Integer1DArray_Item ( self, i ) )
            {
                if ( n < i ) Integer1DArray_Item ( self, n ) = Integer1DArray_Item ( self, i ) ;
                n++ ;
            }
        }
        self->length = n ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Summing.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Integer1DArray_Sum ( const Integer1DArray *self )
{
    Integer sum = 0 ;
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->length ; i++ ) sum += Integer1DArray_Item ( self, i ) ;
    }
    return sum ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create an array view of a raw array with a given size.
! . The input offset, length and stride are the final values.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Integer1DArray_ViewOfRaw ( Integer1DArray *self, const Integer offset, const Integer length, const Integer stride, Integer *data, const Integer size, Status *status )
{
    if ( ( self != NULL ) && ( data != NULL ) )
    {
        auto Integer n ;
        /* . Find space required. */
        if ( length == 0 ) n = 0 ;
        else
        {
            n = offset + 1 ;
            if ( stride > 0 ) n += ( length - 1 ) * stride ;
        }
        if ( n > size ) Status_Set ( status, Status_ArrayOverFlow ) ;
        else
        {
            self->data    = data   ;
            self->isOwner = False  ;
            self->isView  = True   ;
            self->length  = length ;
            self->offset  = offset ;
            self->size    = size   ;
            self->stride  = stride ;
        }
    }
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
static Integer Item_Compare ( const void *vterm1, const void *vterm2 )
{
    Integer i ;
    Integer *term1, *term2 ;
    term1 = ( Integer * ) vterm1 ;
    term2 = ( Integer * ) vterm2 ;
         if ( term1[0] < term2[0] ) i = -1 ;
    else if ( term1[0] > term2[0] ) i =  1 ;
    else i = 0 ;
    return i ;
}
