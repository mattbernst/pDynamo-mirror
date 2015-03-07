/*------------------------------------------------------------------------------
! . File      : AntisymmetricMatrix.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Real antisymmetric matrices.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Notes:
!
!   In future it may be useful to have two strides, i.e. stride0 and stride1, which give the difference between row and column
!   elements, respectively. This would complicate the code but could be employed to have ASM diagonal block slices of existing
!   ASMs. Items in this scheme would be accessed as: (i-1)*s0 + ((i-1)*(i-2))*s1/2 + j*s1 (i > j). 1D views would only then be
!   possible if s0 = s1.
!
!---------------------------------------------------------------------------------------------------------------------------------*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "AntisymmetricMatrix.h"
# include "Memory.h"
# include "SliceOld.h"

/*
extern void                 AntisymmetricMatrix_AnticommutatorAS         ( const AntisymmetricMatrix  *self, const SymmetricMatrix *a, AntisymmetricMatrix *result, Status *status ) ;
*/

/*
AS Comm. self * other + other * self.
S Comm.  self * other - other * self.

From2Symm. a * b - b * a (and, optionally, a * b + b * a).
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . The abolute maximum value.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real AntisymmetricMatrix_AbsoluteMaximum ( const AntisymmetricMatrix *self )
{
    return Real1DArray_AbsoluteMaximum ( &(self->view) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add a scaled matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_AddScaledMatrix ( AntisymmetricMatrix *self, const Real value, const AntisymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) ) Real1DArray_AddScaledArray ( &(self->view), value, &(other->view), status ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
AntisymmetricMatrix *AntisymmetricMatrix_Allocate ( const Integer length, Status *status )
{
    AntisymmetricMatrix *self = NULL ;
    MEMORY_ALLOCATE ( self, AntisymmetricMatrix ) ;
    if ( self != NULL )
    {
        auto Integer n ;
        n = Maximum ( length, 0 ) ;
        /* . Basic data. */
        self->data    = NULL  ;
        self->isOwner = True  ;
        self->isView  = False ;
        self->length  = n ;
        self->offset  = 0 ;
        self->size    = ( n * ( n - 1 ) ) / 2 ;
        self->stride  = 1 ;
        /*. View. */
        self->view.isOwner = False ;
        self->view.isView  = True  ;
        self->view.length  = self->size ;
        self->view.offset  = 0 ;
        self->view.size    = self->size ;
        self->view.stride  = 1 ;
        self->view.data    = NULL ;
        if ( length > 0 )
        {
            MEMORY_ALLOCATEARRAY ( self->data, self->size, Real ) ;
            if ( self->data == NULL ) AntisymmetricMatrix_Deallocate ( &self ) ;
            else self->view.data = self->data ;
        }
        else if ( length < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Anticommutator of antisymmetric and symmetric matrices A * S + S * A.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_AnticommutatorAS ( const AntisymmetricMatrix *self, const SymmetricMatrix *a, AntisymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( result != NULL ) )
    {
        if ( ( self->length == a->dimension ) && ( a->dimension == result->length ) )
        {
            auto Integer i, j, k ;
            auto Real    sum ;
            for ( i = 1 ; i < self->length ; i++ )
            {
                for ( k = 0 ; k < i ; k++ )
                {
                    sum = 0.0e+00 ;
                    /* . Care needed with signs when swapping indices. */
                    for ( j = 0         ; j < k            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, k, j ) - SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, k, j ) ) ;
                    for ( j = ( k + 1 ) ; j < i            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, j, k ) + SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, j, k ) ) ;
                    for ( j = ( i + 1 ) ; j < self->length ; j++ ) sum += ( - AntisymmetricMatrix_Item ( self, j, i ) * SymmetricMatrix_Item ( a, j, k ) + SymmetricMatrix_Item ( a, j, i ) * AntisymmetricMatrix_Item ( self, j, k ) ) ;
                    sum += AntisymmetricMatrix_Item ( self, i, k ) * ( SymmetricMatrix_Item ( a, i, i ) + SymmetricMatrix_Item ( a, k, k ) ) ;
                    AntisymmetricMatrix_Item ( result, i, k ) = sum ;
                }
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of antisymmetric and symmetric matrices A * S - S * A.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorAS ( const AntisymmetricMatrix *self, const SymmetricMatrix *a, SymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( result != NULL ) )
    {
        if ( ( self->length == a->dimension ) && ( a->dimension == result->dimension ) )
        {
            auto Integer i, j, k ;
            auto Real    sum ;
            for ( i = 0 ; i < self->length ; i++ )
            {
                for ( k = 0 ; k < i ; k++ )
                {
                    sum = 0.0e+00 ;
                    for ( j = 0         ; j < k            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, k, j ) + SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, k, j ) ) ;
                    for ( j = ( k + 1 ) ; j < i            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, j, k ) - SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, j, k ) ) ;
                    for ( j = ( i + 1 ) ; j < self->length ; j++ ) sum += ( - AntisymmetricMatrix_Item ( self, j, i ) * SymmetricMatrix_Item ( a, j, k ) - SymmetricMatrix_Item ( a, j, i ) * AntisymmetricMatrix_Item ( self, j, k ) ) ;
                    sum += AntisymmetricMatrix_Item ( self, i, k ) * ( SymmetricMatrix_Item ( a, k, k ) - SymmetricMatrix_Item ( a, i, i ) );
                    SymmetricMatrix_Item ( result, i, k ) = sum ;
                }
                sum = 0.0e+00 ;
                for ( j = 0         ; j < i            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, i, j ) + SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, i, j ) ) ;
                for ( j = ( i + 1 ) ; j < self->length ; j++ ) sum += ( - AntisymmetricMatrix_Item ( self, j, i ) * SymmetricMatrix_Item ( a, j, i ) - SymmetricMatrix_Item ( a, j, i ) * AntisymmetricMatrix_Item ( self, j, i ) ) ;
                SymmetricMatrix_Item ( result, i, i ) = sum ;
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of two symmetric matrices A * B - B * A. Faster version that requires much more memory.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorSS_Fast (       AntisymmetricMatrix *self   ,
                                             const SymmetricMatrix     *a      ,
                                             const SymmetricMatrix     *b      ,
                                                   Real2DArray         *mA     ,
                                                   Real2DArray         *mB     ,
                                                   Real2DArray         *mC     ,
                                                   Status              *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) && ( mA != NULL ) && ( mB != NULL ) && ( mC != NULL ) )
    {
        if ( ( self->length == a->dimension ) && ( a->dimension == b->dimension ) )
        {
            SymmetricMatrix_CopyToReal2DArray ( a, mA, status ) ;
            SymmetricMatrix_CopyToReal2DArray ( b, mB, status ) ;
            Real2DArray_MatrixMultiply ( False, False, 1.0e+00, mA, mB, 0.0e+00, mC, status ) ;
            AntisymmetricMatrix_CopyFromReal2DArray ( self, mC, False, status ) ;
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of two symmetric matrices A * B - B * A. Reference (slow) version.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorSS_Reference ( AntisymmetricMatrix *self, const SymmetricMatrix *a, const SymmetricMatrix *b, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) )
    {
# ifdef USEOPENMP
	# pragma omp parallel /*num_threads(MAXIMUMNUMBEROFTHREADS)*/
# endif
        if ( ( self->length == a->dimension ) && ( a->dimension == b->dimension ) )
        {
            auto Integer i, j, k ;
            auto Real    sum ;
# ifdef USEOPENMP
            # pragma omp for schedule(dynamic)
# endif
            for ( i = 0 ; i < self->length ; i++ )
            {
                for ( k = 0 ; k < i ; k++ )
                {
                    sum = 0.0e+00 ;
                    for ( j = 0         ; j <= k           ; j++ ) sum += ( SymmetricMatrix_Item ( a, i, j ) * SymmetricMatrix_Item ( b, k, j ) - SymmetricMatrix_Item ( b, i, j ) * SymmetricMatrix_Item ( a, k, j ) ) ;
                    for ( j = ( k + 1 ) ; j <= i           ; j++ ) sum += ( SymmetricMatrix_Item ( a, i, j ) * SymmetricMatrix_Item ( b, j, k ) - SymmetricMatrix_Item ( b, i, j ) * SymmetricMatrix_Item ( a, j, k ) ) ;
                    for ( j = ( i + 1 ) ; j < self->length ; j++ ) sum += ( SymmetricMatrix_Item ( a, j, i ) * SymmetricMatrix_Item ( b, j, k ) - SymmetricMatrix_Item ( b, j, i ) * SymmetricMatrix_Item ( a, j, k ) ) ;
                    AntisymmetricMatrix_Item ( self, i, k ) = sum ;
                }
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of three symmetric matrices A * B * C - C * B * A.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorSSS ( AntisymmetricMatrix *self, const SymmetricMatrix *a, const SymmetricMatrix *b, const SymmetricMatrix *c, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) && ( c != NULL ) )
    {
        if ( ( self->length == a->dimension ) && ( a->dimension == b->dimension ) && ( b->dimension == c->dimension ) )
        {
            auto Integer i, ij, j, jk, k, kl, l ;
            auto Real    sum1, sum2, *tab = NULL, *tcb = NULL ;
            tab = Memory_Allocate_Array_Real ( self->length ) ;
            tcb = Memory_Allocate_Array_Real ( self->length ) ;
            if ( ( tab != NULL ) && ( tcb != NULL ) )
            {
                for ( i = 0 ; i < self->length ; i++ )
                {
                    for ( k = 0 ; k < self->length ; k++ )
                    {
                        sum1 = 0.0e+00 ;
                        sum2 = 0.0e+00 ;
                        for ( j = 0 ; j < self->length ; j++ )
                        {
                            if ( i >= j ) ij = ( i * ( i + 1 ) ) / 2 + j ;
                            else          ij = ( j * ( j + 1 ) ) / 2 + i ;
                            if ( j >= k ) jk = ( j * ( j + 1 ) ) / 2 + k ;
                            else          jk = ( k * ( k + 1 ) ) / 2 + j ;
                            sum1 += a->data[ij] * b->data[jk] ;
                            sum2 += c->data[ij] * b->data[jk] ;
                        }
                        tab[k] = sum1 ;
                        tcb[k] = sum2 ;
                    }
                    for ( l = 0 ; l < i ; l++ )
                    {
                        sum1 = 0.0e+00 ;
                        sum2 = 0.0e+00 ;
                        for ( k = 0 ; k < self->length ; k++ )
                        {
                            if ( k >= l ) kl = ( k * ( k + 1 ) ) / 2 + l ;
                            else          kl = ( l * ( l + 1 ) ) / 2 + k ;
                            sum1 += tab[k] * c->data[kl] ;
                            sum2 += tcb[k] * a->data[kl] ;
                        }
                        AntisymmetricMatrix_Item ( self, i, l ) = sum1 - sum2 ;
                    }
                }
            }
            Memory_Deallocate ( tab ) ;
            Memory_Deallocate ( tcb ) ;
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of three symmetric matrices and a transformation M^T * ( A * B * C - C * B * A ) * M.
! . The transformation matrix can be transposed, i.e. M * ( A * B * C - C * B * A ) * M^T.
! . Fast version that uses a lot of memory.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorTSSST (       AntisymmetricMatrix *self       ,
                                           const SymmetricMatrix     *a          ,
                                           const SymmetricMatrix     *b          ,
                                           const SymmetricMatrix     *c          ,
                                           const Real2DArray         *m          ,
                                           const Boolean              mTranspose ,
                                                 Real2DArray         *u          ,
                                                 Real2DArray         *v          ,
                                                 Real2DArray         *w          ,
                                                 Status              *status     )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) && ( c != NULL ) && ( m != NULL ) && ( u != NULL ) && ( v != NULL ) && ( w != NULL ) )
    {
        auto Boolean isOK ;
        isOK = ( a->dimension == b->dimension ) && ( a->dimension == c->dimension ) ;
        if ( mTranspose ) isOK = isOK && ( ( c->dimension == m->length1 ) && ( m->length0 == self->length ) ) ;
        else              isOK = isOK && ( ( c->dimension == m->length0 ) && ( m->length1 == self->length ) ) ;
        if ( isOK )
        {
            SymmetricMatrix_CopyToReal2DArray ( a, u, status ) ;
            SymmetricMatrix_CopyToReal2DArray ( b, v, status ) ;
            Real2DArray_MatrixMultiply ( False, False, 1.0e+00, u, v, 0.0e+00, w, status ) ;
            SymmetricMatrix_CopyToReal2DArray ( c, u, status ) ;
            Real2DArray_MatrixMultiply ( False, False, 1.0e+00, w, u, 0.0e+00, v, status ) ;
            if ( m->length0 == m->length1 )
            {
                Real2DArray_MatrixMultiply ( False,  mTranspose, 1.0e+00, v, m, 0.0e+00, u, status ) ;
                Real2DArray_MatrixMultiply ( !mTranspose, False, 1.0e+00, m, u, 0.0e+00, w, status ) ;
                AntisymmetricMatrix_CopyFromReal2DArray ( self, w, False, status ) ;
            }
            else
            {
                auto Integer n ;
                auto Real2DArray uSlice, wSlice ;
                if ( mTranspose ) n = m->length0 ;
                else              n = m->length1 ;
                Real2DArray_Slice ( u, 0, a->dimension, 1, 0, n, 1, &uSlice, status ) ;
                Real2DArray_Slice ( w, 0, n           , 1, 0, n, 1, &wSlice, status ) ;
                Real2DArray_MatrixMultiply ( False,  mTranspose, 1.0e+00, v, m      , 0.0e+00, &uSlice, status ) ;
                Real2DArray_MatrixMultiply ( !mTranspose, False, 1.0e+00, m, &uSlice, 0.0e+00, &wSlice, status ) ;
                AntisymmetricMatrix_CopyFromReal2DArray ( self, &wSlice, False, status ) ;
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy from a 2D array with antisymmetrization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CopyFromReal2DArray ( AntisymmetricMatrix *self, const Real2DArray *other, const Boolean scale, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer n ;
        n = self->length ;
        if ( ( n != other->length0 ) || ( n != other->length1 ) ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
        else
        {
            auto Integer i, j ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j < i ; j++ ) AntisymmetricMatrix_Item ( self, i, j ) = ( Real2DArray_Item ( other, i, j ) - Real2DArray_Item ( other, j, i ) ) ;
            }
            if ( scale ) AntisymmetricMatrix_Scale ( self, 0.5e+00 ) ;
/*
printf ( "\nASM in CopyFromReal2DArray:\n" ) ;
AntisymmetricMatrix_Print ( self ) ;
printf ( "\nR2D in CopyFromReal2DArray:\n" ) ;
Real2DArray_Print ( other ) ;
*/
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CopyTo ( const AntisymmetricMatrix *self, AntisymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) ) Real1DArray_CopyTo ( &(self->view), &(other->view), status ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy to a full 2D array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CopyToReal2DArray ( const AntisymmetricMatrix  *self, Real2DArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer n ;
        n = self->length ;
        if ( ( n != other->length0 ) || ( n != other->length1 ) ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
        else
        {
            auto Integer i, j ;
            auto Real    item ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    item = AntisymmetricMatrix_Item ( self, i, j ) ;
                    Real2DArray_Item ( other, i, j ) =  item ;
                    Real2DArray_Item ( other, j, i ) = -item ;
                }
                Real2DArray_Item ( other, i, i ) = 0.0e+00 ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Deallocate ( AntisymmetricMatrix **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        (*self)->view.data = NULL ;
        if ( (*self)->isOwner ) MEMORY_DEALLOCATE ( (*self)->data ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a column of the matrix.
! . To get a row scale by -1.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_GetColumn ( const AntisymmetricMatrix *self, const Integer n, Real1DArray *column, Status *status )
{
    if ( ( self != NULL ) && ( column != NULL ) )
    {
        if ( ( n < 0 ) || ( n >= self->length ) ) Status_Set ( status, Status_IndexOutOfRange ) ;
        else if ( self->length != column->length ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
        else
        {
            auto Integer i ;
            for ( i = 0     ; i < n            ; i++ ) Real1DArray_Item ( column, i ) = - AntisymmetricMatrix_Item ( self, n, i ) ;
            Real1DArray_Item ( column, n ) = 0.0e+00 ;
            for ( i = (n+1) ; i < self->length ; i++ ) Real1DArray_Item ( column, i ) =   AntisymmetricMatrix_Item ( self, i, n ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real AntisymmetricMatrix_GetItem ( const AntisymmetricMatrix *self, const Integer i, const Integer j, Status *status )
{
    Integer ij = -1 ;
    Real    sign = 0.0e+00, value = 0.0e+00 ;
    AntisymmetricMatrix_GetItemIndexAndSign ( self, i, j, &ij, &sign, status ) ;
    if ( ij >= 0 ) value = sign * self->data[ij] ;
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the index and sign of an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_GetItemIndexAndSign ( const AntisymmetricMatrix *self, const Integer i, const Integer j, Integer *index, Real *sign, Status *status )
{
    if ( self != NULL )
    {
        if ( ( i >= 0 ) && ( i < self->length ) && ( j >= 0 ) && ( j < self->length ) )
        {
            if ( i != j )
            {
                if ( i > j )
                {
                    (*index) = self->offset + ( ( i * ( i - 1 ) ) / 2 + j ) * self->stride ;
                    (*sign ) =  1.0e+00 ;
                }
                else
                {
                    (*index) = self->offset + ( ( j * ( j - 1 ) ) / 2 + i ) * self->stride ;
                    (*sign ) = -1.0e+00 ;
                }
            }
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Print ( const AntisymmetricMatrix *self )
{
    if ( self == NULL ) printf ( "Null antisymmetric matrix.\n" ) ;
    else
    {
        auto Integer i, j ;
        for ( i = 0 ; i < self->length ; i++ )
        {
            for ( j = 0 ; j < i ; j++ ) printf ( "%15.10f", AntisymmetricMatrix_Item ( self, i, j ) ) ;
            printf ( "   0.0000000000\n" ) ;
        }
        printf ( "\n" ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the items of a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Scale ( AntisymmetricMatrix *self, const Real value )
{
    if ( self != NULL ) Real1DArray_Scale ( &(self->view), value ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the items of a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Set ( AntisymmetricMatrix *self, const Real value )
{
    if ( self != NULL ) Real1DArray_Set ( &(self->view), value ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_SetItem ( AntisymmetricMatrix *self, const Integer i, const Integer j, const Real value, Status *status )
{
    Integer ij = -1 ;
    Real    sign = 0.0e+00 ;
    AntisymmetricMatrix_GetItemIndexAndSign ( self, i, j, &ij, &sign, status ) ;
    if ( ij >= 0 ) self->data[ij] = sign * value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform by a symmetric matrix (S * A * S).
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_SymmetricTransform  ( const AntisymmetricMatrix  *self, const SymmetricMatrix *matrix, AntisymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) )
    {
        if ( ( self->length == matrix->dimension ) && ( matrix->dimension == result->length ) )
        {
            auto Integer n ;
            auto Real1DArray *icolumn = NULL, *SAi = NULL, *tcolumn = NULL ;
            n = self->length ;
            icolumn = Real1DArray_Allocate ( n, status ) ;
            SAi     = Real1DArray_Allocate ( n, status ) ;
            tcolumn = Real1DArray_Allocate ( n, status ) ;
            if ( ( icolumn != NULL ) && ( SAi != NULL ) && ( tcolumn != NULL ) )
            {
                auto Integer i, k, l ;
                for ( i = 1 ; i < n ; i++ )
                {
                    SymmetricMatrix_GetColumn ( matrix, i, icolumn, NULL ) ;
                    for ( k = 0 ; k < n ; k++ )
                    {
                        AntisymmetricMatrix_GetColumn ( self, k, tcolumn, status ) ;
                        Real1DArray_Item ( SAi, k ) = Real1DArray_Dot ( icolumn, tcolumn, NULL ) ;
                    }
                    for ( l = 0 ; l < i ; l++ )
                    {
                        SymmetricMatrix_GetColumn ( matrix, l, tcolumn, NULL ) ;
                        AntisymmetricMatrix_Item ( result, i, l ) = Real1DArray_Dot ( SAi, tcolumn, NULL ) ;
                    }
                }
            }
            else Status_Set ( status, Status_MemoryAllocationFailure ) ;
            Real1DArray_Deallocate ( &icolumn ) ;
            Real1DArray_Deallocate ( &SAi     ) ;
            Real1DArray_Deallocate ( &tcolumn ) ;
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the trace of two matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real AntisymmetricMatrix_Trace2 ( const AntisymmetricMatrix *self, const AntisymmetricMatrix *other, Status *status )
{
    Real value = 0.0e+00 ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        value = - 2.0e+00 * Real1DArray_Dot ( &(self->view), &(other->view), status ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform by a matrix (B^T * A * B) or its transpose (B * A * B^T).
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Transform ( const AntisymmetricMatrix *self, const Real2DArray *matrix, const Boolean useTranspose, AntisymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) )
    {
        auto Real1DArray *BAik = NULL, *column = NULL ;
        BAik   = Real1DArray_Allocate ( self->length, status ) ;
        column = Real1DArray_Allocate ( self->length, status ) ;
        if ( ( BAik != NULL ) && ( column != NULL ) )
        {
            auto Integer     i, k, l ;
            auto Real1DArray iview, lview ;
            /* . (Bij * Ajk * (B^T)kl = Bij * Ajk * Blk). */
            if ( useTranspose )
            {
                if ( ( self->length == matrix->length1 ) && ( result->length == matrix->length0 ) )
                {
                    for ( i = 0 ; i < matrix->length0 ; i++ )
                    {
                        Real2DArray_RowSlice ( matrix, i, &iview, NULL ) ;
                        for ( k = 0 ; k < self->length ; k++ )
                        {
                            AntisymmetricMatrix_GetColumn ( self, k, column, status ) ;
                            Real1DArray_Item ( BAik, k ) = Real1DArray_Dot ( column, &iview, NULL ) ;
                        }
                        for ( l = 0 ; l < i ; l++ )
                        {
                            Real2DArray_RowSlice ( matrix, l, &lview, NULL ) ;
                            AntisymmetricMatrix_Item ( result, i, l ) = Real1DArray_Dot ( BAik, &lview, NULL ) ;
                        }
                    }
                }
                else Status_Set ( status, Status_DimensionError ) ;
            }
            /* . ((B^T)ij * Ajk * Bkl = Bji * Ajk * Bkl). */
            else
            {
                if ( ( self->length == matrix->length0 ) && ( result->length == matrix->length1 ) )
                {
                    for ( i = 0 ; i < matrix->length1 ; i++ )
                    {
                        Real2DArray_ColumnSlice ( matrix, i, &iview, NULL ) ;
                        for ( k = 0 ; k < self->length ; k++ )
                        {
                            AntisymmetricMatrix_GetColumn ( self, k, column, status ) ;
                            Real1DArray_Item ( BAik, k ) = Real1DArray_Dot ( column, &iview, NULL ) ;
                        }
                        for ( l = 0 ; l < i ; l++ )
                        {
                            Real2DArray_ColumnSlice ( matrix, l, &lview, NULL ) ;
                            AntisymmetricMatrix_Item ( result, i, l ) = Real1DArray_Dot ( BAik, &lview, NULL ) ;
                        }
                    }
                }
                else Status_Set ( status, Status_DimensionError ) ;
            }
        }
        else Status_Set ( status, Status_MemoryAllocationFailure ) ;
        Real1DArray_Deallocate ( &BAik   ) ;
        Real1DArray_Deallocate ( &column ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transpose the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Transpose ( AntisymmetricMatrix *self )
{
    if ( self != NULL ) Real1DArray_Scale ( &(self->view), -1.0e+00 ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a matrix view of a raw array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_ViewOfRaw ( AntisymmetricMatrix *self, const Integer length, const Integer stride, Real *data, const Integer rawstride, Status *status )
{
    if ( ( self != NULL ) && ( data != NULL ) )
    {
        auto Integer n ;
        auto Status  localstatus ;
        n = SliceIndices_GetLength ( 0, length, stride, length, &localstatus ) ;
        if ( ! Status_OK ( &localstatus ) ) Status_Set ( status, localstatus ) ;
        else
        {
            auto Integer fullstride ;
            fullstride = rawstride * stride ;
            /* . Basic data. */
            self->data    = data ;
            self->isOwner = False ;
            self->isView  = True  ;
            self->length  = n ;
            self->offset  = 0 ;
            self->size    = ( ( n * ( n - 1 ) ) * fullstride ) / 2 ;
            self->stride  = fullstride ;
            /*. View. */
            self->view.data    = data ;
            self->view.isOwner = False ;
            self->view.isView  = True  ;
            self->view.length  = ( n * ( n - 1 ) ) / 2 ;
            self->view.offset  = 0 ;
            self->view.size    = self->size ;
            self->view.stride  = fullstride ;
        }
    }
}
