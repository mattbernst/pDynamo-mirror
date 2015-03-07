/*------------------------------------------------------------------------------
! . File      : Boolean1DArray.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . 1-D arrays.
!=================================================================================================================================*/

# include <stdlib.h>

# include "Boolean1DArray.h"
# include "Macros.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer Item_Compare ( const void *vterm1, const void *vterm2 ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean1DArray *Boolean1DArray_Allocate ( const Integer length, Status *status )
{
    Boolean1DArray *self = NULL ;
    MEMORY_ALLOCATE ( self, Boolean1DArray ) ;
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
            MEMORY_ALLOCATEARRAY ( self->data, length, Boolean ) ;
            if ( self->data == NULL ) Boolean1DArray_Deallocate ( &self ) ;
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
Boolean1DArray *Boolean1DArray_Clone ( const Boolean1DArray *self, Status *status )
{
    Boolean1DArray *clone = NULL ;
    if ( self != NULL )
    {
        clone = Boolean1DArray_Allocate ( self->length, status ) ;
        Boolean1DArray_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
! . As much data as possible is copied.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Boolean1DArray_CopyTo ( const Boolean1DArray *self, Boolean1DArray *other, Status *status )
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
            for ( i = 0 ; i < n ; i++ ) Boolean1DArray_Item ( other, i ) = Boolean1DArray_Item ( self, i ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Boolean1DArray_Deallocate ( Boolean1DArray **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        if ( (*self)->isOwner ) MEMORY_DEALLOCATE ( (*self)->data ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Boolean1DArray_Length ( const Boolean1DArray *self )
{
    if ( self == NULL ) return 0 ;
    else                return self->length ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
! . This can only be done if the array owns the data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Boolean1DArray_Resize ( Boolean1DArray *self, const Integer length, const Boolean *initializer, Status *status )
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
                    auto Boolean *new = NULL ;
                    auto Integer  n ;
                    /* . The minimal number of items needed is, in fact, offset + length * stride - ( stride - 1 ).
                    ! . However, the longer, stride-consistent, value is left here for the moment just in case. */
                    n = self->offset + length * self->stride ;
                    MEMORY_REALLOCATEARRAY ( new, self->data, n, Boolean ) ;
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
                            for ( i = oldlength ; i < length ; i++ ) Boolean1DArray_Item ( self, i ) = (*initializer) ;
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
void Boolean1DArray_Reverse ( Boolean1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        auto Boolean t ;
        auto Integer i, j ;
        for ( i = 0, j = self->length - 1 ; i < ( self->length / 2 ) ; i++, j-- )
        {
            t = Boolean1DArray_Item ( self, i ) ;
            Boolean1DArray_Item ( self, i ) = Boolean1DArray_Item ( self, j ) ;
            Boolean1DArray_Item ( self, j ) = t ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Boolean1DArray_Set ( Boolean1DArray *self, Boolean value )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->length ; i++ ) Boolean1DArray_Item ( self, i ) = value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting - ascending in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Boolean1DArray_Sort ( Boolean1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 1 ) ) qsort ( ( void * ) Boolean1DArray_Data ( self ), ( size_t ) self->length, self->stride * sizeof ( Boolean ), ( void * ) Item_Compare ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting - return an array of sorted indices. The array itself remains unchanged.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Boolean1DArray_SortIndex ( const Boolean1DArray *self, Integer1DArray *indices, Status *status )
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
                    if ( j < m && ( Boolean1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) < Boolean1DArray_Item ( self, Integer1DArray_Item ( indices, j+1 ) ) ) ) j++ ;
                    if ( Boolean1DArray_Item ( self, ki ) >= Boolean1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) ) break ;
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
                    if ( j < m && ( Boolean1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) < Boolean1DArray_Item ( self, Integer1DArray_Item ( indices, j+1 ) ) ) ) j++ ;
                    if ( Boolean1DArray_Item ( self, ki ) >= Boolean1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) ) break ;
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
void Boolean1DArray_SortUnique ( Boolean1DArray *self )
{
    Boolean1DArray_Sort ( self ) ;
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        auto Integer i, n ;
        for ( i = 1, n = 1 ; i < self->length ; i++ )
        {
            if ( Boolean1DArray_Item ( self, i-1 ) != Boolean1DArray_Item ( self, i ) )
            {
                if ( n < i ) Boolean1DArray_Item ( self, n ) = Boolean1DArray_Item ( self, i ) ;
                n++ ;
            }
        }
        self->length = n ;
    }
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
static Integer Item_Compare ( const void *vterm1, const void *vterm2 )
{
    Integer i ;
    Boolean *term1, *term2 ;
    term1 = ( Boolean * ) vterm1 ;
    term2 = ( Boolean * ) vterm2 ;
         if ( term1[0] < term2[0] ) i = -1 ;
    else if ( term1[0] > term2[0] ) i =  1 ;
    else i = 0 ;
    return i ;
}
