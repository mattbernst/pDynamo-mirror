/*------------------------------------------------------------------------------
! . File      : IntegerNDArray.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . N-D arrays.
!=================================================================================================================================*/

# include "IntegerNDArray.h"
# include "Macros.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
IntegerNDArray *IntegerNDArray_Allocate ( const Integer ndimensions, const Integer *lengths, Status *status )
{
    IntegerNDArray *self = NULL ;
    /* . Check the dimensions. */
    if ( ndimensions >= 0 )
    {
        MEMORY_ALLOCATE ( self, IntegerNDArray ) ;
        if ( self != NULL )
        {
            /* . Initialization. */
            self->data        = NULL  ;
            self->lengths     = NULL  ;
            self->strides     = NULL  ;
            self->isOwner     = True  ;
            self->isView      = False ;
            self->length      = 0 ;
            self->ndimensions = ndimensions ;
            self->offset      = 0 ;
            self->size        = 0 ;
            /* . Dimensions. */
            if ( ndimensions > 0 )
            {
                MEMORY_ALLOCATEARRAY ( self->lengths, ndimensions, Integer ) ;
                MEMORY_ALLOCATEARRAY ( self->strides, ndimensions, Integer ) ;
                if ( ( self->lengths == NULL ) || ( self->strides == NULL ) ) IntegerNDArray_Deallocate ( &self ) ;
                else
                {
                    auto Integer d, length, n ;
                    /* . Initialize the dimensions - loop over in reverse order. */
                    for ( d = ndimensions - 1, n = 1 ; d >= 0 ; d-- )
                    {
                        length = lengths[d] ;
                        if ( length < 0 )
                        {
                            length = 0 ;
                            Status_Set ( status, Status_NegativeArrayLength ) ;
                        }
                        self->lengths[d] = length ;
                        self->strides[d] = n ;
                        n *= length ;
                    }
                    self->length = n ;
                    self->size   = n ;
                    /* . Allocate the data. */
                    if ( n > 0 )
                    {
                        MEMORY_ALLOCATEARRAY ( self->data, n, Integer ) ;
                        if ( self->data == NULL ) IntegerNDArray_Deallocate ( &self ) ;
                    }
                }
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    else Status_Set ( status, Status_InvalidDimension ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
! . Data is copied and kept in the most compact representation possible.
!---------------------------------------------------------------------------------------------------------------------------------*/
IntegerNDArray *IntegerNDArray_Clone ( const IntegerNDArray *self, Status *status )
{
    IntegerNDArray *clone = NULL ;
    if ( self != NULL )
    {
        clone = IntegerNDArray_Allocate ( self->ndimensions, self->lengths, status ) ;
        IntegerNDArray_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
! . As much data as possible is copied as long as the array dimensions are the same.
!---------------------------------------------------------------------------------------------------------------------------------*/
void IntegerNDArray_CopyTo ( const IntegerNDArray *self, IntegerNDArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( self->ndimensions != other->ndimensions )
        {
            auto Integer d, length, ndimensions = self->ndimensions ;
            auto Integer lengths[self->ndimensions] ;
            /* . Find the lengths to copy. */
            for ( d = 0, length = 1 ; d < ndimensions ; d++ )
            {
                if ( self->lengths[d] == other->lengths[d] ) lengths[d] = self->lengths[d] ;
                else
                {
                    Status_Set ( status, Status_ArrayLengthMismatch ) ;
                    lengths[d] = Minimum ( self->lengths[d], other->lengths[d] ) ;
                }
                length *= lengths[d] ;
            }
            /* . Copy. */
            if ( length > 0 )
            {
                auto Integer i, index, nother, nself ;
                auto Integer indices[self->ndimensions], *otherstrides = other->strides, *selfstrides = self->strides ;
                for ( d = 0 ; d < ndimensions ; d++ ) indices[d] = 0 ;
                for ( i = 0, nother = other->offset, nself = self->offset ; i < length ; i++ )
                {
                    other->data[nother] = self->data[nself] ;
                    for ( d = ndimensions-1 ; d >= 0 ; d-- )
                    {
                       index = indices[d] + 1 ;
                       if ( index >= lengths[d] )
                       {
                            indices[d] = 0 ;
                            nother -= ( lengths[d] - 1 ) * otherstrides[d] ;
                            nself  -= ( lengths[d] - 1 ) * selfstrides[d] ;
                       }
                       else
                       {
                            indices[d] = index ;
                            nother += otherstrides[d] ;
                            nself  += selfstrides[d] ;
                            break ;
                        }
                    }
                }
            }
        }
        else Status_Set ( status, Status_ArrayDimensionMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void IntegerNDArray_Deallocate ( IntegerNDArray **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        if ( (*self)->isOwner ) MEMORY_DEALLOCATE ( (*self)->data ) ;
        MEMORY_DEALLOCATE ( (*self)->lengths ) ;
        MEMORY_DEALLOCATE ( (*self)->strides ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer IntegerNDArray_Length ( const IntegerNDArray *self, const Integer dimension )
{
    Integer length = 0 ;
    if ( self != NULL )
    {
        if ( ( dimension < 0 ) || ( dimension >= self->ndimensions ) ) length = self->length ;
        else                                                           length = self->lengths[dimension] ;
    }
    return length ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Number of dimensions.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer IntegerNDArray_NumberOfDimensions ( const IntegerNDArray *self )
{
    if ( self == NULL ) return 0 ;
    else                return self->ndimensions ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
! . This is only done for the slowest changing dimension.
! . This can only be done if the array owns the data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void IntegerNDArray_Resize ( IntegerNDArray *self, const Integer length0, const Integer *initializer, Status *status )
{
    if ( self != NULL )
    {
        if ( self->isOwner )
        {
            /* . Only do something if the new and old lengths differ. */
            if ( ( self->ndimensions > 0 ) && ( length0 != self->lengths[0] ) )
            {
                if ( length0 > 0 )
                {
                    auto Integer n, *new = NULL ;
                    /* . The minimal number of items needed is, in fact, offset + length0 * stride0 - ( stride0 - 1 ).
                    ! . However, the longer, stride-consistent, value is left here for the moment just in case. */
                    n = self->offset + length0 * self->strides[0] ;
                    MEMORY_REALLOCATEARRAY ( new, self->data, n, Integer ) ;
                    /* . Failure so leave as is. */
                    if ( new == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
                    /* . OK. */
                    else
                    {
                        auto Integer oldlength, oldlength0 ;
                        /* . Save the old lengths. */
                        oldlength  = self->length    ;
                        oldlength0 = self->lengths[0] ;
                        /* . Reset data values. */
                        self->data       = new        ;
                        self->length     = ( oldlength / oldlength0 ) * length0 ;
                        self->lengths[0] = length0    ;
                        self->size       = n          ;
                        /* . Initialize. */
                        if ( ( length0 > oldlength0 ) && ( initializer != NULL ) )
                        {
                            /* . Better to use a slice here when have worked out how to do it! */
/* . Start. */
                            auto Integer d, i, index, n, ndimensions = self->ndimensions ;
                            auto Integer indices[self->ndimensions], lengths[self->ndimensions], *strides = self->strides ;
                            for ( d = 0 ; d < ndimensions ; d++ )
                            {
                                indices[d] = 0 ;
                                lengths[d] = self->lengths[d] ;
                            }
                            lengths[0] -= oldlength0 ;
                            for ( i = 0, n = self->offset + oldlength0 * strides[0] ; i < ( self->length - oldlength ) ; i++ )
                            {
                                self->data[n] = (*initializer) ;
                                for ( d = ndimensions-1 ; d >= 0 ; d-- )
                                {
                                   index = indices[d] + 1 ;
                                   if ( index >= lengths[d] ) { indices[d] = 0 ; n -= ( lengths[d] - 1 ) * strides[d] ; }
                                   else { indices[d] = index ; n += strides[d] ; break ; }
                                }
                            }
/* . Stop. */
                        }
                    }
                }
                else
                {
                    MEMORY_DEALLOCATE ( self->data ) ;
                    self->data       = NULL  ;
                    self->isOwner    = True  ;
                    self->isView     = False ;
                    self->length     = 0     ;
                    self->lengths[0] = 0     ;
                    self->size       = 0     ;
                    if ( length0 < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
                }
            }
        }
        else Status_Set ( status, Status_InvalidArrayOperation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void IntegerNDArray_Set ( IntegerNDArray *self, Integer value )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer d, i, index, n, ndimensions = self->ndimensions ;
        auto Integer indices[self->ndimensions], *lengths = self->lengths, *strides = self->strides ;
        for ( d = 0 ; d < ndimensions ; d++ ) indices[d] = 0 ;
        for ( i = 0, n = self->offset ; i < self->length ; i++ )
        {
            self->data[n] = value ;
            for ( d = ndimensions-1 ; d >= 0 ; d-- )
            {
               index = indices[d] + 1 ;
               if ( index >= lengths[d] ) { indices[d] = 0 ; n -= ( lengths[d] - 1 ) * strides[d] ; }
               else { indices[d] = index ; n += strides[d] ; break ; }
            }
        }
    }
}
