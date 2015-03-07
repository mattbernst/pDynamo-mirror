/*------------------------------------------------------------------------------
! . File      : Real1DArray.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . 1-D arrays.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "cblas.h"
# include "Macros.h"
# include "Memory.h"
# include "Real1DArray.h"
# include "SliceOld.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer Item_Compare ( const void *vterm1, const void *vterm2 ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . The abolute maximum value.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real1DArray_AbsoluteMaximum ( const Real1DArray *self )
{
    Real value = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        i = Real1DArray_AbsoluteMaximumIndex ( self ) ;
        if ( i >= 0 ) value = fabs ( Real1DArray_Item ( self, i ) ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The index of the abolute maximum value.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Real1DArray_AbsoluteMaximumIndex ( const Real1DArray *self )
{
    Integer index = -1 ;
    if ( ( self != NULL ) && ( self->length > 0 ) ) index = cblas_idamax ( self->length, Real1DArray_Data ( self ), self->stride ) ;
    return index ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Element-wise addition.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Add ( Real1DArray *self, const Real1DArray *other, Status *status )
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
            for ( i = 0 ; i < n ; i++ ) Real1DArray_Item ( self, i ) += Real1DArray_Item ( other, i ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add a scalar.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_AddScalar ( Real1DArray *self, const Real value )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->length ; i++ ) Real1DArray_Item ( self, i ) += value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add a scaled array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_AddScaledArray ( Real1DArray  *self, const Real value, const Real1DArray *other, Status *status )
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
        if ( n > 0 ) cblas_daxpy ( n, value, Real1DArray_Data ( other ), other->stride, Real1DArray_Data ( self ), self->stride ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real1DArray *Real1DArray_Allocate ( const Integer length, Status *status )
{
    Real1DArray *self = NULL ;
    MEMORY_ALLOCATE ( self, Real1DArray ) ;
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
            MEMORY_ALLOCATEARRAY ( self->data, length, Real ) ;
            if ( self->data == NULL ) Real1DArray_Deallocate ( &self ) ;
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
Real1DArray *Real1DArray_Clone ( const Real1DArray *self, Status *status )
{
    Real1DArray *clone = NULL ;
    if ( self != NULL )
    {
        clone = Real1DArray_Allocate ( self->length, status ) ;
        Real1DArray_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
! . As much data as possible is copied.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_CopyTo ( const Real1DArray *self, Real1DArray *other, Status *status )
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
        if ( n > 0 ) cblas_dcopy ( n, Real1DArray_Data ( self ), self->stride, Real1DArray_Data ( other ), other->stride ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Deallocate ( Real1DArray **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        if ( (*self)->isOwner ) MEMORY_DEALLOCATE ( (*self)->data ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Element-wise divistion.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Divide ( Real1DArray *self, const Real1DArray *other, Status *status )
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
            for ( i = 0 ; i < n ; i++ ) Real1DArray_Item ( self, i ) /= Real1DArray_Item ( other, i ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dot product.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real1DArray_Dot ( const Real1DArray *self, const Real1DArray *other, Status *status )
{
    Real dot = 0.0e+00 ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer n ;
        if ( self->length == other->length ) n = self->length ;
        else
        {
            Status_Set ( status, Status_ArrayLengthMismatch ) ;
            n = Minimum ( self->length, other->length ) ;
        }
        if ( n > 0 ) dot = cblas_ddot ( n, Real1DArray_Data ( self ), self->stride, Real1DArray_Data ( other ), other->stride ) ;
    }
    return dot ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Element-wise exponentiation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Exp ( Real1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        auto Real    v ;
        for ( i = 0 ; i < self->length ; i++ )
        {
            v = Real1DArray_Item ( self, i ) ;
            if ( v > Real_MaximumExponent ) v = Real_Huge ;
            else                            v = exp ( v ) ;
            Real1DArray_Item ( self, i ) = v ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real1DArray_GetItem ( const Real1DArray *self, const Integer i, Status *status )
{
    Real value = 0.0e+00 ;
    if ( self != NULL )
    {
        if ( self->length > 0 )
        {
            if ( ( i >= 0 ) && ( i < self->length ) ) value = Real1DArray_Item ( self, i ) ;
            else Status_Set ( status, Status_IndexOutOfRange ) ;
        }
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for compactness.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean Real1DArray_IsCompact ( const Real1DArray *self )
{
    Boolean isCompact = True ;
    if ( self != NULL ) isCompact = ( self->stride == 1 ) ;
    return isCompact ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for uniformness - always true.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean Real1DArray_IsUniform ( const Real1DArray *self ) { return True ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Real1DArray_Length ( const Real1DArray *self )
{
    if ( self == NULL ) return 0 ;
    else                return self->length ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Element-wise natural logarithms.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Ln ( Real1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        auto Real    v ;
        for ( i = 0 ; i < self->length ; i++ )
        {
            v = Real1DArray_Item ( self, i ) ;
            if ( v <= 0.0e+00 ) v = - Real_Huge ;
            else                v = log ( v ) ;
            Real1DArray_Item ( self, i ) = v ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the maximum value in the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real1DArray_Maximum ( const Real1DArray *self )
{
    Real value = 0.0e+00 ; /* . Should really be largest negative number representable by Real. */
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        value = Real1DArray_Item ( self, 0 ) ;
        for ( i = 0 ; i < self->length ; i++ ) value = Maximum ( value, Real1DArray_Item ( self, i ) ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the minimum value in the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real1DArray_Minimum ( const Real1DArray *self )
{
    Real value = 0.0e+00 ; /* . Should really be largest positive number representable by Real. */
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        value = Real1DArray_Item ( self, 0 ) ;
        for ( i = 0 ; i < self->length ; i++ ) value = Minimum ( value, Real1DArray_Item ( self, i ) ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Element-wise multiplication.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Multiply ( Real1DArray *self, const Real1DArray *other, Status *status )
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
            for ( i = 0 ; i < n ; i++ ) Real1DArray_Item ( self, i ) *= Real1DArray_Item ( other, i ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Normalize ( Real1DArray *self, const Real *nullnormvalue, Status *status )
{
    if ( self != NULL )
    {
        auto Real delta, norm2, safeMinimum ;
        /* . Get the 2-norm. */
        norm2 = Real1DArray_Norm2 ( self ) ;
        /* . Get the null-norm value for normalization. */
        safeMinimum = Real_SafeMinimum ( ) ;
        if ( nullnormvalue == NULL ) delta = safeMinimum ;
        else                         delta = Maximum ( safeMinimum, fabs ( (*nullnormvalue) ) ) ;
        /* . Normalize or set to zero. */
        if ( norm2 > delta )
        {
            Real1DArray_Scale ( self, 1.0e+00 / norm2 ) ;
            Status_Set ( status, Status_Success ) ;
        }
        else
        {
            Real1DArray_Set ( self, 0.0e+00 ) ;
            Status_Set ( status, Status_NullVector ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the 2-norm of the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real1DArray_Norm2 ( const Real1DArray *self )
{
    Real norm2 = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->data != NULL ) ) norm2 = cblas_dnrm2 ( self->length, Real1DArray_Data ( self ), self->stride ) ;
    return norm2 ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Permute the items in an array in place given a permutation vector.
! . The permutation vector is assumed to be valid.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Also inverse permutation? */

void Real1DArray_Permute ( Real1DArray *self, Integer1DArray *permutation, Status *status )
{
    if ( ( self != NULL ) && ( permutation != NULL ) )
    {
        if ( self->length != permutation->length ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
        else
        {
            auto Integer i, j, k, n ;
            auto Real    a, b ;

            /* . Initialization. */
            n = self->length ;
            for ( i = 0 ; i < n ; i++ )
            {
                j = Integer1DArray_Item ( permutation, i ) ;
                if ( i != j ) Integer1DArray_Item ( permutation, i ) = -(j+1) ; /* . Add one to deal with case j = 0. */
            }

            /* . Permutation. */
            for ( i = 0 ; i < n ; i++ )
            {
                k = Integer1DArray_Item ( permutation, i ) ;
                if ( k < 0 )
                {
                    j = i ;
                    a = Real1DArray_Item ( self, i ) ;
                    while ( k < 0 )
                    {
                        k = -(k+1) ;
                        Integer1DArray_Item ( permutation, j ) = k ;
                        b = Real1DArray_Item ( self, k ) ;
                        Real1DArray_Item ( self, k ) = a ;
                        j = k ;
                        a = b ;
                        k = Integer1DArray_Item ( permutation, j ) ;
                    }
                    if ( j != i ) { Status_Set ( status, Status_InvalidPermutation ) ; break ; }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Print ( const Real1DArray *self )
{
    if ( self == NULL ) printf ( "Null real 1-D array.\n" ) ;
    else
    {
        auto Integer i ;
        for ( i = 0 ; i < self->length ; i++ ) printf ( "%15.10f", Real1DArray_Item ( self, i ) ) ;
        printf ( "\n" ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the range spanned by the values in the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Range ( const Real1DArray *self, Real *lower, Real *upper )
{
    if ( lower != NULL ) (*lower) = Real1DArray_Minimum ( self ) ;
    if ( upper != NULL ) (*upper) = Real1DArray_Maximum ( self ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Element-wise reciprocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Reciprocate ( Real1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        auto Real    v ;
        for ( i = 0 ; i < self->length ; i++ )
        {
            v = Real1DArray_Item ( self, i ) ;
            if ( v == 0.0e+00 ) v = Real_Huge   ;
            else                v = 1.0e+00 / v ;
            Real1DArray_Item ( self, i ) = v ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
! . This can only be done if the array owns the data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Resize ( Real1DArray *self, const Integer length, const Real *initializer, Status *status )
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
                    auto Integer n ;
                    auto Real   *new = NULL ;
                    /* . The minimal number of items needed is, in fact, offset + length * stride - ( stride - 1 ).
                    ! . However, the longer, stride-consistent, value is left here for the moment just in case. */
                    n = self->offset + length * self->stride ;
                    MEMORY_REALLOCATEARRAY ( new, self->data, n, Real ) ;
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
                            for ( i = oldlength ; i < length ; i++ ) Real1DArray_Item ( self, i ) = (*initializer) ;
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
void Real1DArray_Reverse ( Real1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        auto Integer i, j ;
        auto Real t ;
        for ( i = 0, j = self->length - 1 ; i < ( self->length / 2 ) ; i++, j-- )
        {
            t = Real1DArray_Item ( self, i ) ;
            Real1DArray_Item ( self, i ) = Real1DArray_Item ( self, j ) ;
            Real1DArray_Item ( self, j ) = t ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the root mean square of the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real1DArray_RootMeanSquare ( const Real1DArray *self )
{
    Real rms = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->length > 0 ) ) rms = Real1DArray_Norm2 ( self ) / sqrt ( ( Real ) self->length ) ;
    return rms ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scaling.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Scale ( Real1DArray *self, const Real scale )
{
    if ( ( self != NULL ) && ( self->data != NULL ) ) cblas_dscal ( self->length, scale, Real1DArray_Data ( self ), self->stride ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Set ( Real1DArray *self, Real value )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->length ; i++ ) Real1DArray_Item ( self, i ) = value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_SetItem ( Real1DArray *self, const Integer i, const Real value, Status *status )
{
    if ( self != NULL )
    {
        if ( self->length > 0 )
        {
            if ( ( i >= 0 ) && ( i < self->length ) ) Real1DArray_Item ( self, i ) = value ;
            else Status_Set ( status, Status_IndexOutOfRange ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Slicing - 1D.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_Slice ( const Real1DArray  *self, const Integer start, const Integer stop, const Integer stride, Real1DArray *slice, Status *status )
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
void Real1DArray_Sort ( Real1DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 1 ) ) qsort ( ( void * ) Real1DArray_Data ( self ), ( size_t ) self->length, self->stride * sizeof ( Real ), ( void * ) Item_Compare ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting - return an array of sorted indices. The array itself remains unchanged.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_SortIndex ( const Real1DArray *self, Integer1DArray *indices, Status *status )
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
                    if ( j < m && ( Real1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) < Real1DArray_Item ( self, Integer1DArray_Item ( indices, j+1 ) ) ) ) j++ ;
                    if ( Real1DArray_Item ( self, ki ) >= Real1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) ) break ;
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
                    if ( j < m && ( Real1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) < Real1DArray_Item ( self, Integer1DArray_Item ( indices, j+1 ) ) ) ) j++ ;
                    if ( Real1DArray_Item ( self, ki ) >= Real1DArray_Item ( self, Integer1DArray_Item ( indices, j ) ) ) break ;
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
void Real1DArray_SortUnique ( Real1DArray *self )
{
    Real1DArray_Sort ( self ) ;
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        auto Integer i, n ;
        for ( i = 1, n = 1 ; i < self->length ; i++ )
        {
            if ( Real1DArray_Item ( self, i-1 ) != Real1DArray_Item ( self, i ) )
            {
                if ( n < i ) Real1DArray_Item ( self, n ) = Real1DArray_Item ( self, i ) ;
                n++ ;
            }
        }
        self->length = n ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sparsity (as a percentage).
! . Sparsity is defined as the number of items with a magnitude less than a given tolerance.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real1DArray_Sparsity ( const Real1DArray *self, const Real tolerance )
{
    Real sparsity = 1.0e+00 ;
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i, n ;
        for ( i = n = 0 ; i < self->length ; i++ )
        {
            if ( fabs ( Real1DArray_Item ( self, i ) ) <= tolerance ) n += 1 ;
        }
        sparsity = ( ( Real ) n ) / ( ( Real ) self->length ) ;
    }
    return 100.0e+00 * sparsity ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sum of array items.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real1DArray_Sum ( const Real1DArray *self )
{
    Real sum = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->length ; i++ ) sum += Real1DArray_Item ( self, i ) ;
    }
    return sum ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create an array view of a raw array with a given size.
! . The input offset, length and stride are the final values.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real1DArray_ViewOfRaw ( Real1DArray *self, const Integer offset, const Integer length, const Integer stride, Real *data, const Integer size, Status *status )
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
    Real *term1, *term2 ;
    term1 = ( Real * ) vterm1 ;
    term2 = ( Real * ) vterm2 ;
         if ( term1[0] < term2[0] ) i = -1 ;
    else if ( term1[0] > term2[0] ) i =  1 ;
    else i = 0 ;
    return i ;
}
