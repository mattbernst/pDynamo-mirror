/*------------------------------------------------------------------------------
! . File      : DoubleSymmetricMatrix.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Double symmetric arrays.
!
! . These are symmetric four-dimensional arrays. The following items are equivalent: ijkl, ijlk, jikl, jilk, klij, klji, lkij, lkji.
!
! . Packing is similar to SymmetricMatrix: i >= j, k >= l and (ij) >= (kl).
!
!=================================================================================================================================*/

# include <stdio.h>

# include "cblas.h"
# include "Memory.h"
# include "DoubleSymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DoubleSymmetricMatrix *DoubleSymmetricMatrix_Allocate ( const Integer length, Status *status )
{
    DoubleSymmetricMatrix *self = NULL ;
    MEMORY_ALLOCATE ( self, DoubleSymmetricMatrix ) ;
    if ( self != NULL )
    {
        auto Integer n ;
        n = Maximum ( length, 0 ) ;
        self->data     = NULL  ;
        self->isOwner  = True  ;
        self->length0  = n ;
        self->length01 = ( n * ( n + 1 ) ) / 2 ;
        self->length   = self->length01 * self->length01 ;
        self->size     = self->length ;
        if ( self->length > 0 )
        {
            MEMORY_ALLOCATEARRAY ( self->data, self->length, Real ) ;
            if ( self->data == NULL ) DoubleSymmetricMatrix_Deallocate ( &self ) ;
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
DoubleSymmetricMatrix *DoubleSymmetricMatrix_Clone ( const DoubleSymmetricMatrix *self, Status *status )
{
    DoubleSymmetricMatrix *clone = NULL ;
    if ( self != NULL )
    {
        clone = DoubleSymmetricMatrix_Allocate ( self->length0, status ) ;
        DoubleSymmetricMatrix_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_CopyTo ( const DoubleSymmetricMatrix *self, DoubleSymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( self->length != other->length ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
        else cblas_dcopy ( self->length, self->data, 1, other->data, 1 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_Deallocate ( DoubleSymmetricMatrix **self )
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
Real DoubleSymmetricMatrix_GetItem ( const DoubleSymmetricMatrix *self, const Integer i, const Integer j, const Integer k, const Integer l, Status *status )
{
    Real value = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer ijkl ;
        ijkl = DoubleSymmetricMatrix_Index ( i, j, k, l ) ;
        if ( ( ijkl >= 0 ) && ( ijkl < self->length ) ) value = self->data[ijkl] ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_IncrementItem ( const DoubleSymmetricMatrix *self, const Integer i, const Integer j, const Integer k, const Integer l, const Real value, Status *status )
{
    if ( self != NULL )
    {
        auto Integer ijkl ;
        ijkl = DoubleSymmetricMatrix_Index ( i, j, k, l ) ;
        if ( ( ijkl >= 0 ) && ( ijkl < self->length ) ) self->data[ijkl] += value ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an index into the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define IJIndex(i,j) ( ( i * ( i + 1 ) ) / 2 + j )
Integer DoubleSymmetricMatrix_Index ( const Integer i, const Integer j, const Integer k, const Integer l )
{
    auto Integer pqrs, p, pq, q, r, rs, s ;
    p  = Maximum ( i, j ) ; q  = Minimum ( i, j ) ;
    r  = Maximum ( k, l ) ; s  = Minimum ( k, l ) ;
    pq = IJIndex ( p, q ) ; rs = IJIndex ( r, s ) ;
    if ( pq >= rs ) pqrs = IJIndex ( pq, rs ) ;
    else            pqrs = IJIndex ( rs, pq ) ;
    return pqrs ;
}
# undef IJIndex

/*----------------------------------------------------------------------------------------------------------------------------------
! . Length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer DoubleSymmetricMatrix_Length ( const DoubleSymmetricMatrix *self )
{
    Integer length = 0 ;
    if ( self != NULL ) length = self->length0 ;
    return length ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_Print ( const DoubleSymmetricMatrix *self )
{
    if ( self == NULL ) printf ( "Null double symmetric matrix.\n" ) ;
    else
    {
        auto Integer i, ij, ijkl, j, k, kl, l, upper ;
        for ( i = ij = ijkl = 0 ; i < self->length0 ; i++ )
        {
            for ( j = 0 ; j <= i ; ij++, j++ )
            {
                for ( k = kl = 0 ; k <= i ; k++ )
                {
                    if ( k == i ) upper = j ;
                    else          upper = k ;
                    for ( l = 0 ; l <= upper ; ijkl++, kl++, l++ )
                    {
                        printf ( "%5d %5d %5d %5d %5d %5d %5d %12.6f\n", ijkl, i, j, k, l, ij, kl, self->data[ijkl] ) ;
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_Set ( DoubleSymmetricMatrix *self, Real value )
{
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->length ; i++ ) self->data[i] = value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_SetItem ( const DoubleSymmetricMatrix *self, const Integer i, const Integer j, const Integer k, const Integer l, const Real value, Status *status )
{
    if ( self != NULL )
    {
        auto Integer ijkl ;
        ijkl = DoubleSymmetricMatrix_Index ( i, j, k, l ) ;
        if ( ( ijkl >= 0 ) && ( ijkl < self->length ) ) self->data[ijkl] = value ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unweight the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_Unweight ( DoubleSymmetricMatrix *self )
{
    if ( ( self != NULL ) && ( self->length0 > 0 ) )
    {
        auto Integer i, ij, ijkl, j, k, kl, l, upper ;
        auto Real    w ;
        for ( i = ij = ijkl = 0 ; i < self->length0 ; i++ )
        {
            for ( j = 0 ; j <= i ; ij++, j++ )
            {
                for ( k = kl = 0 ; k <= i ; k++ )
                {
                    if ( k == i ) upper = j ;
                    else          upper = k ;
                    for ( l = 0 ; l <= upper ; ijkl++, kl++, l++ )
                    {
                        w = 0.125e+00 ;
                        if ( i == j ) w *= 2.0e+00 ;
                        if ( k == l ) w *= 2.0e+00 ;
                        if ( ij == kl ) w *= 2.0e+00 ;
                        self->data[ijkl] *= w ;
                    }
                }
            }
        }
    }
}
