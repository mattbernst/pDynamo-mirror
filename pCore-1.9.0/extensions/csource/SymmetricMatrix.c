/*------------------------------------------------------------------------------
! . File      : SymmetricMatrix.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Real symmetric matrices.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Notes:
!
!   In future it may be useful to have two strides, i.e. stride0 and stride1, which give the difference between row and column
!   elements, respectively. This would complicate the code but could be employed to have SM diagonal block slices of existing
!   SMs. Items in this scheme would be accessed as: i*s0 + (i*(i-1))*s1/2 + j*s1 (i >= j). 1D views would only then be possible
!   if s0 = s1.
!
!   It seems that most blas procedures do not allow strides for symmetric matrices.
!
!---------------------------------------------------------------------------------------------------------------------------------*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "cblas.h"
# include "f2clapack.h"

# include "Memory.h"
# include "SliceOld.h"
# include "SymmetricMatrix.h"

/*==================================================================================================================================
! . Standard procedures.
!=================================================================================================================================*/
# define __RANGE_CHECKING

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
! . The matrix is not initialized.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetricMatrix *SymmetricMatrix_Allocate ( const Integer n )
{
    SymmetricMatrix *self = NULL ;
    if ( n >= 0 )
    {
        self             = ( SymmetricMatrix * ) Memory_Allocate ( sizeof ( SymmetricMatrix ) ) ;
        self->dimensionP = n ;
        self->dimension  = self->dimensionP ;
        self->sizeP      = ( n * ( n + 1 ) ) / 2  ;
        self->size       = self->sizeP ;
        if ( n > 0 ) self->data = Memory_Allocate_Array_Real ( self->size ) ;
        else         self->data = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetricMatrix *SymmetricMatrix_AllocateN ( const Integer n, Status *status )
{
    SymmetricMatrix *self = NULL ;
    MEMORY_ALLOCATE ( self, SymmetricMatrix ) ;
    if ( self != NULL )
    {
        auto Integer length ;
        length = Maximum ( n, 0 ) ;
        self->dimension  = length ;
        self->dimensionP = length ;
        self->size       = ( length * ( length + 1 ) ) / 2 ;
        self->sizeP      = self->size ;
        self->data       = NULL ;
        if ( self->size > 0 )
        {
            MEMORY_ALLOCATEARRAY ( self->data, self->size, Real ) ;
            if ( self->data == NULL )
            {
                Status_Set ( status, Status_MemoryAllocationFailure ) ;
                SymmetricMatrix_Deallocate ( &self ) ;
            }
        }
        else if ( n < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Anticommutator of two symmetric matrices A * B + B * A.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_AnticommutatorSS ( SymmetricMatrix *self, const SymmetricMatrix *a, const SymmetricMatrix *b, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) )
    {
        if ( ( self->dimension == a->dimension ) && ( a->dimension == b->dimension ) )
        {
            auto Integer i, j, k ;
            auto Real    sum ;
            for ( i = 0 ; i < self->dimension ; i++ )
            {
                for ( k = 0 ; k <= i ; k++ )
                {
                    sum = 0.0e+00 ;
                    for ( j = 0         ; j <= k              ; j++ ) sum += ( SymmetricMatrix_Item ( a, i, j ) * SymmetricMatrix_Item ( b, k, j ) + SymmetricMatrix_Item ( b, i, j ) * SymmetricMatrix_Item ( a, k, j ) ) ;
                    for ( j = ( k + 1 ) ; j <= i              ; j++ ) sum += ( SymmetricMatrix_Item ( a, i, j ) * SymmetricMatrix_Item ( b, j, k ) + SymmetricMatrix_Item ( b, i, j ) * SymmetricMatrix_Item ( a, j, k ) ) ;
                    for ( j = ( i + 1 ) ; j < self->dimension ; j++ ) sum += ( SymmetricMatrix_Item ( a, j, i ) * SymmetricMatrix_Item ( b, j, k ) + SymmetricMatrix_Item ( b, j, i ) * SymmetricMatrix_Item ( a, j, k ) ) ;
                    SymmetricMatrix_Item ( self, i, k ) = sum ;
                }
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetricMatrix *SymmetricMatrix_Clone ( const SymmetricMatrix *self )
{
    SymmetricMatrix *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = SymmetricMatrix_Allocate ( self->dimension ) ;
        new->dimension = self->dimension ;
        new->size      = self->size      ;
        for ( i = 0 ; i < self->size ; i++ ) new->data[i] = self->data[i] ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy the elements of the first matrix by those of the second.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SymmetricMatrix_CopyTo ( const SymmetricMatrix *self, SymmetricMatrix *other )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( self->size == other->size )
        {
            auto Integer i ;
            for ( i = 0 ; i < self->size ; i++ ) other->data[i] = self->data[i] ;
        }
        else status = Status_DimensionError ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy from a full 2D array (and symmetrize at the same time).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_CopyFromReal2DArray ( SymmetricMatrix *self, const Real2DArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer n ;
        n = self->dimension ;
        if ( ( n != other->length0 ) || ( n != other->length1 ) ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
        else
        {
            auto Integer i, j ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j <= i ; j++ ) SymmetricMatrix_Item ( self, i, j ) = 0.5e+00 * ( Real2DArray_Item ( other, i, j ) + Real2DArray_Item ( other, j, i ) ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy to a full 2D array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_CopyToReal2DArray ( const SymmetricMatrix *self, Real2DArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer n ;
        n = self->dimension ;
        if ( ( n != other->length0 ) || ( n != other->length1 ) ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
        else
        {
            auto Integer i, j ;
            auto Real    item ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    item = SymmetricMatrix_Item ( self, i, j ) ;
                    Real2DArray_Item ( other, i, j ) = item ;
                    Real2DArray_Item ( other, j, i ) = item ;
                }
                Real2DArray_Item ( other, i, i ) = SymmetricMatrix_Item ( self, i, i ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Deallocate ( SymmetricMatrix **symmetricmatrix )
{
   if ( (*symmetricmatrix) != NULL )
   {
      Memory_Deallocate ( (*symmetricmatrix)->data ) ;
      Memory_Deallocate ( (*symmetricmatrix) ) ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dimension of the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
int SymmetricMatrix_Dimension ( const SymmetricMatrix *self )
{
   if ( self == NULL ) { return               0 ; }
   else                { return self->dimension ; }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a column (or a row) of the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_GetColumn ( const SymmetricMatrix *self, const Integer n, Real1DArray *column, Status *status )
{
    if ( ( self != NULL ) && ( column != NULL ) )
    {
        if ( ( n < 0 ) || ( n >= self->dimension ) ) Status_Set ( status, Status_IndexOutOfRange ) ;
        else if ( self->dimension != column->length ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
        else
        {
            auto Integer i ;
            for ( i = 0     ; i <= n              ; i++ ) Real1DArray_Item ( column, i ) = SymmetricMatrix_Item ( self, n, i ) ;
            for ( i = (n+1) ; i < self->dimension ; i++ ) Real1DArray_Item ( column, i ) = SymmetricMatrix_Item ( self, i, n ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_GetItem ( const SymmetricMatrix *self, const Integer i, const Integer j, Status *status )
{
    Real value = 0.0e+00 ;
    if ( self != NULL )
    {
        if ( ( i >= 0 ) && ( i < self->dimension ) && ( j >= 0 ) && ( j < self->dimension ) )
        {
            if ( i >= j ) value = SymmetricMatrix_Item ( self, i, j ) ;
            else          value = SymmetricMatrix_Item ( self, j, i ) ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Print ( SymmetricMatrix *self )
{
    if ( self == NULL ) printf ( "Null symmetric matrix.\n" ) ;
    else
    {
        auto Integer i, j ;
        for ( i = 0 ; i < self->dimension ; i++ )
        {
            for ( j = 0 ; j <= i ; j++ ) printf ( "%15.10f", SymmetricMatrix_Get_Component ( self, i, j ) ) ;
            printf ( "\n" ) ;
        }
        printf ( "\n" ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_SetItem ( SymmetricMatrix *self, const Integer i, const Integer j, const Real value, Status *status )
{
    if ( self != NULL )
    {
        if ( ( i >= 0 ) && ( i < self->dimension ) && ( j >= 0 ) && ( j < self->dimension ) )
        {
            if ( i >= j ) SymmetricMatrix_Item ( self, i, j ) = value ;
            else          SymmetricMatrix_Item ( self, j, i ) = value ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*==================================================================================================================================
! . Non-standard procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Add a scaled matrix to another.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SymmetricMatrix_AddScaledMatrix ( SymmetricMatrix *self, const Real alpha, const SymmetricMatrix *other )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( self->size == other->size )
        {
            cblas_daxpy ( self->size, alpha, other->data, 1, self->data, 1 ) ;
            status = Status_Success ;
        }
        else status = Status_DimensionError ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Decrement the elements of the first matrix by those of the second.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Decrement ( SymmetricMatrix *matrix1, SymmetricMatrix *matrix2 )
{
   if ( ( matrix1 != NULL ) && ( matrix2 != NULL ) )
   {
      auto Integer i ;
      for ( i = 0 ; i < matrix1->size ; i++ ) matrix1->data[i] -= matrix2->data[i] ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Diagonalization.
! . The input matrix is destroyed on output.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Diagonalize ( SymmetricMatrix *self, Real1DArray *eigenValues, Real2DArray *eigenVectors, Status *status )
{
    if ( ( self != NULL ) && ( eigenValues != NULL ) )
    {
        /* . Various checks - self and eigenValues need to be compact as does dimension 1 of eigenVectors. */
             if ( ! ( SymmetricMatrix_IsCompact ( self ) && Real1DArray_IsCompact ( eigenValues ) ) ) Status_Set ( status, Status_NonCompactArray     ) ;
        else if ( ( eigenVectors != NULL ) && ( ! Real2DArray_IsCompact ( eigenVectors, 1, NULL ) ) ) Status_Set ( status, Status_NonCompactDimension ) ;
        else
        {
            /* . LAPACK specific. */
            auto Real  *work = NULL ;
            auto integer  ifail = 0, *iwork = NULL, liwork, lwork, m, n ;
            n      = ( integer ) self->dimension ;
            liwork = 3 + 5 * n ;
            lwork  = 1 + 6 * n + n * n ;
            iwork  = ( integer * ) calloc ( liwork, sizeof ( integer ) ) ;
            work   = Memory_Allocate_Array_Real ( ( Integer )  lwork ) ;
            if ( ( iwork != NULL ) && ( work != NULL ) )
            {
                if ( eigenVectors == NULL )
                {
                    dspevd_ ( "N", "U", &n, self->data, Real1DArray_Data ( eigenValues ), NULL, NULL, work, &lwork, iwork, &liwork, &ifail ) ;
                }
                else
                {
                    m = ( integer ) eigenVectors->stride0 ;
                    dspevd_ ( "V", "U", &n, self->data, Real1DArray_Data ( eigenValues ), Real2DArray_Data ( eigenVectors ), &m, work, &lwork, iwork, &liwork, &ifail ) ;
                    Real2DArray_Transpose ( eigenVectors, status ) ;
                }
                if ( ifail != 0 ) { Status_Set ( status, Status_DiagonalizationFailure ) ; printf ( "\nDiagonalization Error = %d\n", ifail ) ; }
            }
            else Status_Set ( status, Status_MemoryAllocationFailure ) ;
            Memory_Deallocate ( iwork ) ;
            Memory_Deallocate (  work ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Partial diagonalization.
! . The input matrix is destroyed on output.
! . Note that the eigenVectors are stored rowwise in this procedure!
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_DiagonalizePartial ( SymmetricMatrix *self, Integer lower, Integer upper, Real1DArray *eigenValues, Real2DArray *eigenVectors, Status *status )
{
    if ( ( self != NULL ) && ( eigenValues != NULL ) )
    {
        auto Boolean doEigenVectors ;
        auto Integer n, numberOfEigenValues ;
        doEigenVectors = ( eigenVectors != NULL ) ;
        n              = self->dimension ;

        /* . Check that lower and upper are within range. */
        numberOfEigenValues = SliceIndices_GetLength ( lower, upper, 1, n, status ) ;
        if ( numberOfEigenValues > 0 )
        {
            /* . Various checks - self and eigenValues need to be compact as does dimension 1 of eigenVectors. */
                 if ( ! ( SymmetricMatrix_IsCompact ( self ) && Real1DArray_IsCompact ( eigenValues ) ) ) Status_Set ( status, Status_NonCompactArray     ) ;
            else if ( doEigenVectors  && ( ! Real2DArray_IsCompact ( eigenVectors, 1, NULL )          ) ) Status_Set ( status, Status_NonCompactDimension ) ;
            else if ( doEigenVectors  && ( ( eigenVectors->length0 < numberOfEigenValues ) || ( eigenVectors->length1 < n ) ) ) Status_Set ( status, Status_ArrayLengthMismatch ) ;
            else
            {
                auto Integer         dummy, il, info, iu ;
                auto Real            abstol ;
                auto Integer1DArray *ifail = NULL, *iwork = NULL ;
                auto Real1DArray    *rwork = NULL ;

                /* . Allocate space. */
                if ( doEigenVectors ) ifail = Integer1DArray_Allocate ( n, status ) ;
                iwork = Integer1DArray_Allocate ( 5 * n, status ) ;
                rwork = Real1DArray_Allocate    ( 8 * n, status ) ;
                if ( ( ( doEigenVectors && ( ifail != NULL ) ) || ( ! doEigenVectors ) ) && ( iwork != NULL ) && ( rwork != NULL ) )
                {
                    abstol = 2.0e+00 * dlamch_ ( "S" ) ; /* . The recommended value. */
                    info   = 0 ;
                    /* . Reset indices (Fortran convention). */
                    il = lower + 1 ;
                    iu = upper     ;
                    if ( doEigenVectors )
                    {
                        dspevx_ ( "V", "I", "U", &n, SymmetricMatrix_Data ( self ), NULL, NULL, &il, &iu, &abstol, &numberOfEigenValues, Real1DArray_Data ( eigenValues ), Real2DArray_Data ( eigenVectors ), &(eigenVectors->stride0),
                                                                                                                                    Real1DArray_Data ( rwork ), Integer1DArray_Data ( iwork ), Integer1DArray_Data ( ifail ), &info ) ;
                    }
                    else
                    {
                        dummy = 1 ;
                        dspevx_ ( "N", "I", "U", &n, SymmetricMatrix_Data ( self ), NULL, NULL, &il, &iu, &abstol, &numberOfEigenValues, Real1DArray_Data ( eigenValues ), NULL, &dummy, Real1DArray_Data ( rwork ), Integer1DArray_Data ( iwork ), NULL, &info ) ;
                    }
                    if ( info != 0 ) { Status_Set ( status, Status_DiagonalizationFailure ) ; printf ( "\nDiagonalization Error = %d\n", info ) ; }
                }
                else Status_Set ( status, Status_MemoryAllocationFailure ) ;
                Integer1DArray_Deallocate ( &ifail ) ;
                Integer1DArray_Deallocate ( &iwork ) ;
                Real1DArray_Deallocate    ( &rwork ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the value of an element.
! . Inefficient version.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_Get_Component ( SymmetricMatrix *symmetricmatrix, const Integer i, const Integer j )
{
   Real value = 0.0e+00 ;
   if ( symmetricmatrix != NULL )
   {
      auto Integer ij ;
      if ( i >= j ) ij = ( i * ( i + 1 ) ) / 2 + j ;
      else          ij = ( j * ( j + 1 ) ) / 2 + i ;
      if ( ij < symmetricmatrix->size ) value = symmetricmatrix->data[ij] ;
   }
   return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment the elements of the first matrix by those of the second.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Increment ( SymmetricMatrix *matrix1, SymmetricMatrix *matrix2 )
{
   if ( ( matrix1 != NULL ) && ( matrix2 != NULL ) )
   {
      auto Integer i ;
      for ( i = 0 ; i < matrix1->size ; i++ ) matrix1->data[i] += matrix2->data[i] ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment the value of a component.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_IncrementComponent ( SymmetricMatrix *self, const Integer i, const Integer j, const Real value )
{
    if ( self != NULL )
    {
        auto Integer ij ;
        if ( i >= j ) ij = ( i * ( i + 1 ) ) / 2 + j ;
        else          ij = ( j * ( j + 1 ) ) / 2 + i ;
        if ( ij < self->size ) self->data[ij] += value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment the diagonal elements of a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_IncrementDiagonal ( SymmetricMatrix *self, const Real value )
{
    if ( ( self != NULL ) && ( value != 0.0e+00 ) )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->dimension ; i++ ) self->data[(i*(i+3))/2] += value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment the diagonal elements given an array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_IncrementDiagonalFromArray ( const SymmetricMatrix *self, const Real1DArray *diagonal, Status *status )
{
    if ( ( self != NULL ) && ( diagonal != NULL ) )
    {
        auto Integer i, n ;
        n = Minimum ( self->dimension, diagonal->length ) ;
        for ( i = 0 ; i < n ; i++ ) self->data[(i*(i+3))/2] += Real1DArray_Item ( diagonal, i ) ;
        if ( self->dimension != diagonal->length ) Status_Set ( status, Status_InvalidDimension ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment the value of an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_IncrementItem ( SymmetricMatrix *self, const Integer i, const Integer j, const Real value, Status *status )
{
    if ( self != NULL )
    {
        auto Integer ij ;
        if ( i >= j ) ij = ( i * ( i + 1 ) ) / 2 + j ;
        else          ij = ( j * ( j + 1 ) ) / 2 + i ;
        if ( ( ij >= 0 ) && ( ij < self->size ) ) self->data[ij] += value ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment the elements of a diagonal block of a matrix.
! . The data in dblock is stored in upper triangular form.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Increment_DB ( SymmetricMatrix *symmetricmatrix, const Integer start, const Integer n, const Real *dblock )
{
   if ( ( symmetricmatrix != NULL ) && ( n > 0 ) )
   {
      auto Integer i, ij, j, k ;
      for ( i = start, k = 0 ; i <= ( start + n - 1 ) ; i++ )
      {
         ij = ( i * ( i + 1 ) ) / 2 + start ;
         for ( j = start ; j <= i ; ij++, j++, k++ ) symmetricmatrix->data[ij] += dblock[k] ;
# ifdef __RANGE_CHECKING
if ( ij > symmetricmatrix->size ) printf ( "SymmetricMatrix_Increment_DB: ij out of range: %d\n", ij ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment the elements of an on-diagonal block of a matrix.
! . The data in dblock is stored as a full square matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Increment_DBlockM ( SymmetricMatrix *symmetricmatrix, const Integer istart, const Integer ni, const Real *dblock )
{
   if ( ( symmetricmatrix != NULL ) && ( ni > 0 ) )
   {
      auto Integer i, ii, ij, j, n ;
      for ( i = 0 ; i < ni ; i++ )
      {
         ii = i + istart ;
         ij = ( ii * ( ii + 1 ) ) / 2 + istart ;
         n  = i * ni ;
         for ( j = 0 ; j <= i ; ij++, j++, n++ ) symmetricmatrix->data[ij] += dblock[n] ;
# ifdef __RANGE_CHECKING
if ( ij > symmetricmatrix->size ) printf ( "SymmetricMatrix_Increment_DBlockM: ij out of range: %d\n", ij ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment the elements of an off-diagonal block of a matrix.
! . The data in oblock is stored as a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Increment_OB ( SymmetricMatrix *symmetricmatrix, const Integer istart, const Integer ni,
                                                    const Integer jstart, const Integer nj, const Real *oblock )
{
   if ( ( symmetricmatrix != NULL ) && ( ni > 0 ) && ( nj > 0 ) )
   {
      auto Integer i, ij, j, n ;
      for ( i = istart, n = 0 ; i <= ( istart + ni - 1 ) ; i++ )
      {
         ij = ( i * ( i + 1 ) ) / 2 + jstart ;
         for ( j = jstart ; j <= ( jstart + nj - 1 ) ; ij++, j++, n++ ) symmetricmatrix->data[ij] += oblock[n] ;
# ifdef __RANGE_CHECKING
if ( ij > symmetricmatrix->size ) printf ( "SymmetricMatrix_Increment_OB: ij out of range: %d\n", ij ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment the elements of an off-diagonal block of a matrix.
! . The data in oblock is stored as a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Increment_OBlockM ( SymmetricMatrix *symmetricmatrix, const Integer istart, const Integer ni, const Integer jstart, const Integer nj, const Real *oblock )
{
   if ( ( symmetricmatrix != NULL ) && ( ni > 0 ) && ( nj > 0 ) )
   {
      auto Integer i, ij, j, n ;
      for ( i = istart, n = 0 ; i <= ( istart + ni - 1 ) ; i++ )
      {
         ij = ( i * ( i + 1 ) ) / 2 + jstart ;
         for ( j = jstart ; j <= ( jstart + nj - 1 ) ; ij++, j++, n++ ) symmetricmatrix->data[ij] += oblock[n] ;
# ifdef __RANGE_CHECKING
if ( ij > symmetricmatrix->size ) printf ( "SymmetricMatrix_Put_OBlockM: ij out of range: %d\n", ij ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy selected elements of the matrix to a square form.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_IndexedCopyToReal2DArray ( const SymmetricMatrix *self, const Integer1DArray *indices, Real2DArray *target, Status *status )
{
    if ( ( self != NULL ) && ( indices != NULL ) && ( target != NULL ) )
    {
        auto Integer i, j, m, n ;
        auto Real    v ;
        for ( i = 0 ; i < indices->length ; i++ )
        {
            m = Integer1DArray_Item ( indices, i ) ;
            for ( j = 0 ; j <= i ; j++ )
            {
                n = Integer1DArray_Item  ( indices, j ) ;
                v = SymmetricMatrix_Item ( self, m, n ) ;
                Real2DArray_Item ( target, i, j ) = v ;
                Real2DArray_Item ( target, j, i ) = v ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the inverse of the matrix from the eigenValues and eigenVectors.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define EIGENVALUE_TOLERANCE 1.0e-10
void SymmetricMatrix_Invert ( const SymmetricMatrix *self, const Real *tolerance, SymmetricMatrix *inverse, Status *status )
{
    if ( ( self != NULL ) && ( inverse != NULL ) )
    {
        if ( self->dimension == inverse->dimension )
        {
            auto Integer i, ij, j, n, u ;
            auto Real    fact, tol ;
            auto Real1DArray *eigenValues  ;
            auto Real2DArray *eigenVectors ;

            /* . Initialization. */
            n = self->dimension ;
            if ( tolerance == NULL ) tol = EIGENVALUE_TOLERANCE ;
            else                     tol = (*tolerance) ;

            /* . Allocate space. */
            eigenValues  = Real1DArray_Allocate ( n,    status ) ;
            eigenVectors = Real2DArray_Allocate ( n, n, status ) ;

            /* . Save self. */
            SymmetricMatrix_CopyTo ( self, inverse ) ;

            /* . Diagonalize the matrix. */
            SymmetricMatrix_Diagonalize ( inverse, eigenValues, eigenVectors, status ) ;

            /* . Construct the inverse matrix. */
            SymmetricMatrix_Set ( inverse, 0.0e+00 ) ;

            /* . Real loop over eigenVectors. */
            for ( u = 0 ; u < n ; u++ )
            {
                if ( fabs ( Real1DArray_Item ( eigenValues, u ) ) > tol )
                {
                    fact = 1.0e+00 / Real1DArray_Item ( eigenValues, u ) ;
                    for ( i = 0, ij = 0 ; i < n ; i++ )
                    {
                        for ( j = 0 ; j <= i ; ij++, j++ ) inverse->data[ij] += fact * Real2DArray_Item ( eigenVectors, i, u ) * Real2DArray_Item ( eigenVectors, j, u ) ;
                    }
                }
            }

            /* . Deallocate space. */
            Real1DArray_Deallocate ( &eigenValues  ) ;
            Real2DArray_Deallocate ( &eigenVectors ) ;
        }
        else Status_Set ( status, Status_InvalidDimension ) ;
    }
}
# undef EIGENVALUE_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Inverse square-root from eigenValues and eigenVectors.
! . The eigenValues are changed to their inverse square roots!
!---------------------------------------------------------------------------------------------------------------------------------*/
# define EIGENVALUE_TOLERANCE 1.0e-5
Status SymmetricMatrix_InverseSquareRoot ( SymmetricMatrix *self, Real1DArray *eigenValues, const Real2DArray *eigenVectors )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( eigenValues != NULL ) && ( eigenVectors != NULL ) )
    {
        if ( ( self->dimension == eigenVectors->length0 ) && ( eigenValues->length == eigenVectors->length1 ) )
        {
            auto Real f, *ivector ;
            auto Integer    i, n = self->dimension ;
            SymmetricMatrix_Set ( self, 0.0e+00 ) ;
            for ( i = 0 ; i < n ; i++ )
            {
                /* . Make the square root of the eigenValue and store. */
                f = sqrt ( Maximum ( Real1DArray_Item ( eigenValues, i ), 0.0e+00 ) ) ;
                if ( f > EIGENVALUE_TOLERANCE )
                {
                    f = 1.0e+00 / f ;
                    Real1DArray_Item ( eigenValues, i ) = f ;
                    /* . Update the overlap. */
                    ivector = Real2DArray_ItemPointer ( eigenVectors, 0, i ) ;
                    cblas_dspr ( CblasColMajor, CblasUpper, n, f, ivector, n, self->data ) ;
                }
            }
            status = Status_Success ;
        }
        else status = Status_DimensionError ;
    }
    return status ;
}
# undef EIGENVALUE_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for compactness.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SymmetricMatrix_IsCompact ( const SymmetricMatrix *self ) { return True ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for a diagonal matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SymmetricMatrix_IsDiagonal ( const SymmetricMatrix *self, const Real tolerance )
{
    Boolean isDiagonal = False ;
    if ( self != NULL )
    {
        auto Integer i, ij, j ;
        isDiagonal = True ;
        for ( i = 0, ij = 0 ; i < self->dimension ; i++, ij++ )
        {
            for ( j = 0 ; j < i ; ij++, j++ )
            {
                if ( fabs ( self->data[ij] ) > tolerance ) { isDiagonal = False ; break ; }
            }
            if ( ! isDiagonal ) break ;
        }
    }
    return isDiagonal ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for uniformness.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SymmetricMatrix_IsUniform ( const SymmetricMatrix *self ) { return True ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a symmetric matrix given a set of eigenvalues and eigenvectors.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Explicit case when both strides are 1? */
/*
Alternative/complementary strategies:
Use matrix/vector multiplication for j loop. However, conflicts possible with optimized blas if use threads.
Use matrix/matrix multiplication with two work spaces (although operation count is doubled):
Scale columns of matrix -> work1.
eigenvectors * work1 -> work2.
Symmetrize work2.
*/
void SymmetricMatrix_MakeFromEigensystem (       SymmetricMatrix *self            ,
                                           const Boolean          zeroMatrix      ,
                                           const Integer          numberOfVectors ,
                                           const Real1DArray     *eigenvalues     ,
                                           const Real2DArray     *eigenvectors    ,
                                                 Status          *status          )
{
    if ( zeroMatrix ) SymmetricMatrix_Set ( self, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( eigenvalues != NULL ) && ( eigenvectors != NULL ) && ( numberOfVectors != 0 ) )
    {
        auto Integer n, size ;
        size = eigenvectors->length0 ;
        if ( ( numberOfVectors < 0 ) || ( numberOfVectors > size ) ) n = size ;
        else n = numberOfVectors ;
        if ( ( n <= eigenvalues->length ) && ( self->dimension <= size ) && ( n <= eigenvectors->length1 ) )
        {
# ifdef USEOPENMP
            #pragma omp parallel shared ( status )
# endif
            {
                auto Integer i, j, o ;
                auto Real    sum, *work ;
                work = Real_Allocate ( n ) ;
                if ( work != NULL )
                {
# ifdef USEOPENMP
		    #pragma omp for schedule ( dynamic )
# endif
                    for ( i = 0 ; i < self->dimension ; i++ )
                    {
                        for ( o = 0 ; o < n ; o++ ) work[o] = Real1DArray_Item ( eigenvalues, o ) * Real2DArray_Item ( eigenvectors, i, o ) ;
                        for ( j = 0 ; j <= i ; j++ )
                        {
                            for ( o = 0, sum = 0.0e+00 ; o < n ; o++ ) sum += work[o] * Real2DArray_Item ( eigenvectors, j, o ) ;
                            SymmetricMatrix_Item ( self, i, j ) += sum ;
                        }
                    }
                    Real_Deallocate ( &work ) ;
                }
                else Status_Set ( status, Status_OutOfMemory ) ;
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the maximum absolute value.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_MaximumAbsoluteValue ( const SymmetricMatrix *self )
{
    Real maxval = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Real a ;
        auto Integer    i ;
        for ( i = 0 ; i < self->size ; i++ )
        {
            a = fabs ( self->data[i] ) ;
            if ( a > maxval ) maxval = a ;
        }
    }
    return maxval ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the product of two symmetric matrices (self * other).
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SymmetricMatrix_Multiply2 ( const SymmetricMatrix *self, const SymmetricMatrix *other, Real2DArray *result )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( other != NULL ) && ( result != NULL ) )
    {
        if ( ( self->dimension == other->dimension ) && ( result->length0 == self->dimension ) && ( result->length1 == self->dimension ) )
        {
            auto Real  sum ;
            auto Integer     i, ij, j, jk, k, n ;
            /* . Aij * Bjk = Cik. */
            n = self->dimension ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( k = 0 ; k < n ; k++ )
                {
                    for ( j = 0, sum = 0.0e+00 ; j < n ; j++ )
                    {
                        if ( i < j ) ij = ( j * ( j + 1 ) ) / 2 + i ;
                        else         ij = ( i * ( i + 1 ) ) / 2 + j ;
                        if ( k < j ) jk = ( j * ( j + 1 ) ) / 2 + k ;
                        else         jk = ( k * ( k + 1 ) ) / 2 + j ;
                        sum += self->data[ij] * other->data[jk] ;
                    }
                    Real2DArray_Item ( result, i, k ) = sum ;
                }
            }
            status = Status_Success ;
        }
        else status = Status_DimensionError ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the trace of the product of two symmetric matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_Multiply2_Trace ( const SymmetricMatrix *matrix1, const SymmetricMatrix *matrix2 )
{
   Real trace = 0.0e+00 ;
   if ( ( matrix1 != NULL ) && ( matrix2 != NULL ) )
   {
      auto Integer i, ij, j ;
      for ( i = 0, ij = 0 ; i < matrix1->dimension ; i++, ij++ )
      {
         for ( j = 0 ; j < i ; ij++, j++ ) trace += matrix1->data[ij] * matrix2->data[ij] ;
         trace += 0.5e+00 * matrix1->data[ij] * matrix2->data[ij] ;
      }
      trace *= 2.0e+00 ;
# ifdef __RANGE_CHECKING
if ( ( ij > matrix1->size ) || ( ij > matrix2->size ) ) printf ( "SymmetricMatrix_Multiply2_Trace: ij out of range: %d\n", ij ) ;
# endif
   }
   return trace ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the trace of the product of four symmetric matrices of the form A * S * B * S.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_MultiplyASBS_Trace ( const SymmetricMatrix *A, const SymmetricMatrix *B, const SymmetricMatrix *S, Status *status )
{
    Real trace = 0.0e+00 ;
    if ( ( A != NULL ) && ( B != NULL ) && ( S != NULL ) )
    {
        auto SymmetricMatrix *SBS ;
        SBS   = SymmetricMatrix_AllocateN ( A->dimension, status ) ;
        SymmetricMatrix_SymmetricTransform ( B, S, SBS, status ) ;
        trace = SymmetricMatrix_Multiply2_Trace ( A, SBS ) ;
        SymmetricMatrix_Deallocate ( &SBS ) ;
    }
    return trace ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Compute a canonical orthogonalizing transformation.
! . The matrix is destroyed (unless it is diagonal)!
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DIAGONAL_TOLERANCE   1.0e-10
# define EIGENVALUE_TOLERANCE 1.0e-10
Real2DArray *SymmetricMatrix_OrthogonalizingTransformation ( SymmetricMatrix *self, const Real *diagonaltolerance, const Real *eigenValuetolerance, Status *status )
{
    Real2DArray *matrix = NULL ;
    if ( self != NULL )
    {
        auto Integer      d, i, n ;
        auto Real         dtolerance, etolerance, evalue ;
        auto Real1DArray *eigenValues, icolumn, ncolumn ;

        /* . Initialization. */
        d = self->dimension ;
        if ( diagonaltolerance   == NULL ) dtolerance = DIAGONAL_TOLERANCE     ;
        else                               dtolerance = (*diagonaltolerance)   ;
        if ( eigenValuetolerance == NULL ) etolerance = EIGENVALUE_TOLERANCE   ;
        else                               etolerance = (*eigenValuetolerance) ;

        /* . Check for a diagonal matrix. */
        if ( SymmetricMatrix_IsDiagonal ( self, dtolerance ) )
        {
            /* . Find the number of diagonal elements within tolerance. */
            for ( i = 0, n = 0 ; i < d ; i++ )
            {
                evalue = fabs ( SymmetricMatrix_Get_Component ( self, i, i ) ) ;
                if ( evalue > etolerance ) n++ ;
            }

            /* . Allocate the matrix. */
            matrix = Real2DArray_Allocate ( d, n, NULL ) ;
            Real2DArray_Set ( matrix, 0.0e+00 ) ;

            /* . Set the diagonal elements. */
            for ( i = 0, n = 0 ; i < d ; i++ )
            {
                evalue = fabs ( SymmetricMatrix_Get_Component ( self, i, i ) ) ;
                if ( evalue > etolerance )
                {
                    evalue = 1.0e+00 / sqrt ( evalue ) ;
                    Real2DArray_Item ( matrix, i, n ) = evalue ;
                    n++ ;
                }
            }
        }
        /* . Compute the transformation. */
        else
        {
            /* . Diagonalize the matrix. */
            eigenValues = Real1DArray_Allocate ( d,    NULL ) ;
            matrix      = Real2DArray_Allocate ( d, d, NULL ) ;
            SymmetricMatrix_Diagonalize ( self, eigenValues, matrix, NULL ) ;

            /* . Loop over the eigenValues. */
            for ( i = 0, n = 0 ; i < d ; i++ )
            {
                evalue = fabs ( Real1DArray_Item ( eigenValues, i ) ) ;
                if ( evalue > etolerance )
                {
                    evalue = 1.0e+00 / sqrt ( evalue ) ;
                    /* . Scale the column. */
                    Real2DArray_ColumnSlice ( matrix, i, &icolumn, NULL ) ;
                    Real1DArray_Scale ( &icolumn, evalue ) ;
                    /* . Shift the column (down) if necessary. */
	            if ( n != i )
                    {
                        Real2DArray_ColumnSlice ( matrix, n, &ncolumn, NULL ) ;
                        Real1DArray_CopyTo ( &icolumn, &ncolumn, NULL ) ;
                    }
                    n++ ;
                }
            }

            /* . Change the dimensions of matrix if necessary. */
            if ( n < d )
            {
                matrix->length1 = n ;
                matrix->length  = matrix->length0 * n ;
            }

            /* . Finish up. */
            Real1DArray_Deallocate ( &eigenValues ) ;
        }
    }
    return matrix ;
}
# undef DIAGONAL_TOLERANCE
# undef EIGENVALUE_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Post-matrix multiply.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SymmetricMatrix_PostMatrixMultiply ( const SymmetricMatrix *self, const Real2DArray *matrix, const Boolean QTRANSPOSE, Real2DArray *result )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) )
    {
        auto Integer i ;
        if ( QTRANSPOSE )
        {
            if ( ( self->dimension == matrix->length1 ) && ( result->length0 == self->dimension ) && ( result->length1 == matrix->length0 ) )
            {
                for ( i = 0 ; i < matrix->length1 ; i++ )
                {
                    cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real2DArray_ItemPointer ( matrix, i, 0 ), matrix->stride1,
                                                                                         0.0e+00,  Real2DArray_ItemPointer ( result, 0, i ), result->stride0 ) ;
                }
                status = Status_Success ;
            }
            else status = Status_DimensionError ;
        }
        else
        {
            if ( ( self->dimension == matrix->length0 ) && ( result->length0 == self->dimension ) && ( result->length1 == matrix->length1 ) )
            {
                for ( i = 0 ; i < matrix->length1 ; i++ )
                {
                    cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real2DArray_ItemPointer ( matrix, 0, i ), matrix->stride0, 
                                                                                         0.0e+00,  Real2DArray_ItemPointer ( result, 0, i ), result->stride0 ) ;
                }
                status = Status_Success ;
            }
            else status = Status_DimensionError ;
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pre-matrix multiply.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_PreMatrixMultiply ( const SymmetricMatrix *self, const Real2DArray *matrix, const Boolean useTranspose, Real2DArray *result, Status *status )
{
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) )
    {
        auto Integer i ;
        if ( useTranspose )
        {
            if ( ( self->dimension == matrix->length0 ) && ( result->length0 == matrix->length1 ) && ( result->length1 == self->dimension ) )
            {
                for ( i = 0 ; i < matrix->length0 ; i++ )
                {
                    cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real2DArray_ItemPointer ( matrix, 0, i ), matrix->stride0,
                                                                                         0.0e+00,  Real2DArray_ItemPointer ( result, i, 0 ), result->stride1 ) ;
                }
                Status_Set ( status, Status_Success ) ;
            }
            else Status_Set ( status, Status_DimensionError ) ;
        }
        else
        {
            if ( ( self->dimension == matrix->length1 ) && ( result->length0 == matrix->length0 ) && ( result->length1 == self->dimension ) )
            {
                for ( i = 0 ; i < matrix->length0 ; i++ )
                {
                    cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real2DArray_ItemPointer ( matrix, i, 0 ), matrix->stride1,
                                                                                         0.0e+00,  Real2DArray_ItemPointer ( result, i, 0 ), result->stride1 ) ;
                }
                Status_Set ( status, Status_Success ) ;
            }
            else Status_Set ( status, Status_DimensionError ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Project out a set of vectors from a symmetric matrix. The projection is
! . done as follows: ( 1 - V V^T ) S ( 1 - V V^T ).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_ProjectOut ( SymmetricMatrix **self, const Real2DArray *vectors, Status *status )
{
    if ( ( (*self) != NULL ) && ( vectors != NULL ) )
    {
        if ( ( (*self)->dimension == vectors->length0 ) && ( vectors->length1 > 0 ) )
        {
            auto SymmetricMatrix *new, *projection ;
            projection = SymmetricMatrix_ProjectionMatrix ( vectors ) ;
            new        = SymmetricMatrix_AllocateN ( projection->dimension, status ) ;
            SymmetricMatrix_SymmetricTransform ( (*self), projection, new, NULL ) ;
            SymmetricMatrix_Deallocate ( self ) ;
            (*self) = new ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a projection matrix of the form ( 1 - V V^T ).
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetricMatrix *SymmetricMatrix_ProjectionMatrix ( const Real2DArray *vectors )
{
    SymmetricMatrix *self = NULL ;
    if ( vectors != NULL )
    {
        auto Integer i, ij, j, v ;
        auto Real    sum         ;
        self = SymmetricMatrix_Allocate ( vectors->length0 ) ;
        for ( i = 0, ij = 0 ; i < vectors->length0 ; i++ )
        {
            for ( j = 0 ; j <= i ; ij++, j++ )
            {
                for ( v = 0, sum = 0.0e+00 ; v < vectors->length1 ; v++ ) sum += Real2DArray_Item ( vectors, i, v ) * Real2DArray_Item ( vectors, j, v ) ;
                self->data[ij] = - sum ;
                if ( i == j ) self->data[ij] += 1.0e+00 ;
            }
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add alpha V V^T to a symmetric matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Raise ( SymmetricMatrix *self, const Real2DArray *vectors, const Real value, Status *status )
{
    if ( ( self != NULL ) && ( vectors != NULL ) && ( value != 0.0e+00 ) )
    {
        if ( ( self->dimension == vectors->length0 ) && ( vectors->length1 > 0 ) )
        {
            auto Integer i, ij, j, v ;
            auto Real    sum         ;
            for ( i = 0, ij = 0 ; i < vectors->length0 ; i++ )
            {
                for ( j = 0 ; j <= i ; ij++, j++ )
                {
                   for ( v = 0, sum = 0.0e+00 ; v < vectors->length1 ; v++ ) sum += Real2DArray_Item ( vectors, i, v ) * Real2DArray_Item ( vectors, j, v ) ;
                   self->data[ij] += value * sum ;
                }
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Rank-1 update.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . self needs to be compact. */
void SymmetricMatrix_Rank1Update ( SymmetricMatrix *self, const Real alpha, const Real1DArray *vector, Status *status )
{
    if ( ( self != NULL ) && ( vector != NULL ) ) cblas_dspr ( CblasColMajor, CblasUpper, self->dimension, alpha, Real1DArray_Data ( vector ), vector->stride, SymmetricMatrix_Data ( self ) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale all the elements of the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Scale ( SymmetricMatrix *symmetricmatrix, const Real value )
{
   if ( symmetricmatrix != NULL )
   {
      auto Integer i ;
      for ( i = 0 ; i < symmetricmatrix->size ; i++ ) symmetricmatrix->data[i] *= value ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the elements of a diagonal block of a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Scale_DB ( SymmetricMatrix *symmetricmatrix, const Integer start, const Integer n, const Real scale )
{
   if ( ( symmetricmatrix != NULL ) && ( n > 0 ) )
   {
      auto Integer i, ij, j ;
      for ( i = start ; i <= ( start + n - 1 ) ; i++ )
      {
         ij = ( i * ( i + 1 ) ) / 2 + start ;
         for ( j = start ; j <= i ; ij++, j++ ) symmetricmatrix->data[ij] *= scale ;
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale some diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Scale_Diagonal ( SymmetricMatrix *symmetricmatrix, const Integer start, const Integer n, const Real scale )
{
   if ( ( symmetricmatrix != NULL ) && ( n > 0 ) )
   {
      auto Integer i, ii ;
      for ( i = start ; i <= ( start + n - 1 ) ; i++ )
      {
         ii = ( i * ( i + 3 ) ) / 2 ;
         symmetricmatrix->data[ii] *= scale ;
# ifdef __RANGE_CHECKING
if ( ii >= symmetricmatrix->size ) printf ( "SymmetricMatrix_Scale_Diagonal: ii out of range: %d\n", ii ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the elements of an off-diagonal block of a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Scale_OB ( SymmetricMatrix *symmetricmatrix, const Integer istart, const Integer ni, const Integer jstart, const Integer nj, const Real scale )
{
   if ( ( symmetricmatrix != NULL ) && ( ni > 0 ) && ( nj > 0 ) )
   {
      auto Integer i, ij, j ;
      for ( i = istart ; i <= ( istart + ni - 1 ) ; i++ )
      {
         ij = ( i * ( i + 1 ) ) / 2 + jstart ;
         for ( j = jstart ; j <= ( jstart + nj - 1 ) ; ij++, j++ ) symmetricmatrix->data[ij] *= scale ;
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the off-diagonal elements of the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Scale_OD ( SymmetricMatrix *symmetricmatrix, const Real value )
{
   if ( symmetricmatrix != NULL )
   {
      auto Integer i, ij, j ;
      for ( i = 0, ij = 0 ; i < symmetricmatrix->dimension ; i++, ij++ )
      {
         for ( j = 0 ; j < i ; ij++, j++ ) symmetricmatrix->data[ij] *= value ;
      }
# ifdef __RANGE_CHECKING
if ( ij > symmetricmatrix->size ) printf ( "SymmetricMatrix_Scale_OD: ij out of range: %d\n", ij ) ;
# endif
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Assign a value to all the elemnts of the matrix.
! . Inefficient version.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set ( SymmetricMatrix *symmetricmatrix, const Real value )
{
   if ( symmetricmatrix != NULL )
   {
      auto Integer i ;
      for ( i = 0 ; i < symmetricmatrix->size ; i++ ) symmetricmatrix->data[i] = value ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the elements of a column of a matrix from  a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set_Column_From_Vector ( SymmetricMatrix *matrix, const Integer column, const Integer mstart,
                                                const Real1DArray *vector, const Integer vstart, const Integer n )
{
   if ( ( matrix != NULL ) && ( vector != NULL ) && ( n > 0 ) )
   {
      auto Integer i ;
      for ( i = 0 ; i < n ; i++ ) SymmetricMatrix_Set_Component ( matrix, column, mstart + i, Real1DArray_Item ( vector, vstart + i ) ) ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Assign a value to an element.
! . Inefficient version.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set_Component ( SymmetricMatrix *symmetricmatrix, const Integer i, const Integer j, const Real value )
{
   if ( symmetricmatrix != NULL )
   {
      auto Integer ij ;
      if ( i >= j ) ij = ( i * ( i + 1 ) ) / 2 + j ;
      else          ij = ( j * ( j + 1 ) ) / 2 + i ;
      if ( ij < symmetricmatrix->size ) symmetricmatrix->data[ij] = value ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the elements of an on-diagonal block of a matrix.
! . The data in dblock is stored as a full square matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set_DBlockM ( SymmetricMatrix *symmetricmatrix, const Integer istart, const Integer ni, const Real *dblock )
{
   if ( ( symmetricmatrix != NULL ) && ( ni > 0 ) )
   {
      auto Integer i, ii, ij, j, n ;
      for ( i = 0 ; i < ni ; i++ )
      {
         ii = i + istart ;
         ij = ( ii * ( ii + 1 ) ) / 2 + istart ;
         n  = i * ni ;
         for ( j = 0 ; j <= i ; ij++, j++, n++ ) symmetricmatrix->data[ij] = dblock[n] ;
# ifdef __RANGE_CHECKING
if ( ij > symmetricmatrix->size ) printf ( "SymmetricMatrix_Set_DBlockM: ij out of range: %d\n", ij ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the elements of a diagonal block of a matrix from a diagonal block of
! . another matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set_DB_From_SM_DB ( SymmetricMatrix *matrix1, const Integer start1,
                                         SymmetricMatrix *matrix2, const Integer start2,
                                                                        const Integer ncolumns )
{
   if ( ( matrix1 != NULL ) && ( matrix2 != NULL ) && ( ncolumns > 0 ) )
   {
      auto Integer i1, i2, ij1, ij2, m, n ;
      for ( i1 = start1, i2 = start2, n = 0 ; n < ncolumns ; i1++, i2++, n++ )
      {
         ij1 = ( i1 * ( i1 + 1 ) ) / 2 + start1 ;
         ij2 = ( i2 * ( i2 + 1 ) ) / 2 + start2 ;
         for ( m = 0 ; m <= n ; ij1++, ij2++, m++ ) matrix1->data[ij1] = matrix2->data[ij2] ;
# ifdef __RANGE_CHECKING
if ( ( ij1 > matrix1->size ) || ( ij2 > matrix2->size ) ) printf ( "SymmetricMatrix_Set_DB_From_SM_DB: ij1 or ij2 out of range: %d %d\n", ij1, ij2 ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Assign a value to some of the diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set_Diagonal ( SymmetricMatrix *symmetricmatrix, const Integer start, const Integer n, const Real value )
{
   if ( ( symmetricmatrix != NULL ) && ( n > 0 ) )
   {
      auto Integer i, ii ;
      for ( i = start ; i <= ( start + n - 1 ) ; i++ )
      {
         ii = ( i * ( i + 3 ) ) / 2 ;
         symmetricmatrix->data[ii] = value ;
# ifdef __RANGE_CHECKING
if ( ii >= symmetricmatrix->size ) printf ( "SymmetricMatrix_Set_Diagonal: ii out of range: %d\n", ii ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy a vector to the diagonal elements of a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set_Diagonal_V ( SymmetricMatrix *symmetricmatrix, const Integer start, const Integer n, const Real *values )
{
   if ( ( symmetricmatrix != NULL ) && ( n > 0 ) )
   {
      auto Integer i, ii, k ;
      for ( i = start, k = 0 ; i <= ( start + n - 1 ) ; i++, k++ )
      {
         ii = ( i * ( i + 3 ) ) / 2 ;
         symmetricmatrix->data[ii] = values[k] ;
# ifdef __RANGE_CHECKING
if ( ii >= symmetricmatrix->size ) printf ( "SymmetricMatrix_Set_Diagonal_V: ii out of range: %d\n", ii ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the elements of an off-diagonal block of a matrix.
! . The data in oblock is stored as a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set_OBlockM ( SymmetricMatrix *symmetricmatrix, const Integer istart, const Integer ni, const Integer jstart, const Integer nj, const Real *oblock )
{
   if ( ( symmetricmatrix != NULL ) && ( ni > 0 ) && ( nj > 0 ) )
   {
      auto Integer i, ij, j, n ;
      for ( i = istart, n = 0 ; i <= ( istart + ni - 1 ) ; i++ )
      {
         ij = ( i * ( i + 1 ) ) / 2 + jstart ;
         for ( j = jstart ; j <= ( jstart + nj - 1 ) ; ij++, j++, n++ ) symmetricmatrix->data[ij] = oblock[n] ;
# ifdef __RANGE_CHECKING
if ( ij > symmetricmatrix->size ) printf ( "SymmetricMatrix_Put_OBlockM: ij out of range: %d\n", ij ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the elements of an off-diagonal block of a matrix from an off-diagonal
! . block of another matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set_OB_From_SM_OB ( SymmetricMatrix *matrix1, const Integer istart1, const Integer jstart1,
                                         SymmetricMatrix *matrix2, const Integer istart2, const Integer jstart2,
                                                                                 const Integer ni, const Integer nj )
{
   if ( ( matrix1 != NULL ) && ( matrix2 != NULL ) && ( ni > 0 ) && ( nj > 0 ) )
   {
      auto Integer i1, i2, ij1, ij2, m, n ;
      for ( i1 = istart1, i2 = istart2, n = 0 ; n < ni ; i1++, i2++, n++ )
      {
         ij1 = ( i1 * ( i1 + 1 ) ) / 2 + jstart1 ;
         ij2 = ( i2 * ( i2 + 1 ) ) / 2 + jstart2 ;
         for ( m = 0 ; m < nj ; ij1++, ij2++, m++ ) matrix1->data[ij1] = matrix2->data[ij2] ;
# ifdef __RANGE_CHECKING
if ( ( ij1 > matrix1->size ) || ( ij2 > matrix2->size ) ) printf ( "SymmetricMatrix_Set_OB_From_SM_OB: ij1 or ij2 out of range: %d %d\n", ij1, ij2 ) ;
# endif
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Zero all the elements of the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set_Zero ( SymmetricMatrix *symmetricmatrix )
{
   if ( symmetricmatrix != NULL )
   {
      auto Integer i ;
      for ( i = 0 ; i < symmetricmatrix->size ; i++ ) symmetricmatrix->data[i] = 0.0e+00 ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Size of the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
int SymmetricMatrix_Size ( const SymmetricMatrix *self )
{
   if ( self == NULL ) { return          0 ; }
   else                { return self->size ; }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Solve a set of linear equations for a vector RHS.
! . On entry self and rhs contain the matrix and RHS whereas on exit self
! . is destroyed and rhs contains the solution.
! . Better memory management needed here.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _RCONDITIONTOLERANCE 1.0e-15
void SymmetricMatrix_SolveLinearEquations ( SymmetricMatrix *self, Real1DArray *rhs, Status *status )
{
    if ( ( self != NULL ) && ( rhs != NULL ) )
    {
        auto integer ifail = -1, n   ;
        n = self->dimension ;
        /* . Matrix must have same dimension as vector. */
        if ( n == rhs->length )
        {
            auto Real  *afp, *berr, *ferr, rcond, *t, *work, *x ;
            auto integer *ipiv, *iwork, nrhs ;
            /* . Allocate space. */
            nrhs  = 1 ;
            ipiv  = ( integer * ) calloc ( n, sizeof ( integer ) ) ;
            iwork = ( integer * ) calloc ( n, sizeof ( integer ) ) ;
            afp   = Memory_Allocate_Array_Real ( ( Integer ) ( ( n * ( n + 1 ) ) / 2 ) ) ;
            berr  = Memory_Allocate_Array_Real ( ( Integer )       n   ) ;
            ferr  = Memory_Allocate_Array_Real ( ( Integer )       n   ) ;
            work  = Memory_Allocate_Array_Real ( ( Integer ) ( 3 * n ) ) ;
            x     = Memory_Allocate_Array_Real ( ( Integer )       n   ) ;
	    /* . Transfer rhs to t (needed as no vector increment allowed). */
	    if ( rhs->stride == 1 )
	    {
	        t = rhs->data ;
            }
	    else
	    {
                t = Memory_Allocate_Array_Real ( ( Integer ) n ) ;
   	        cblas_dcopy ( n, rhs->data, rhs->stride, t, 1 ) ;
	    }
            /* . Solve. */
            dspsvx_ ( "N", "U", &n, &nrhs, self->data, afp, ipiv, t, &n, x, &n, &rcond, ferr, berr, work, iwork, &ifail ) ;

            /* . Copy solution to rhs. */
	    cblas_dcopy ( n, x, 1, rhs->data, rhs->stride ) ;

            /* . Check the condition number. */
            if ( fabs ( rcond ) < _RCONDITIONTOLERANCE ) ifail = n + 1 ;

            /* . Finish up. */
            Memory_Deallocate ( afp   ) ;
            Memory_Deallocate ( berr  ) ;
            Memory_Deallocate ( ferr  ) ;
            Memory_Deallocate ( ipiv  ) ;
            Memory_Deallocate ( iwork ) ;
            Memory_Deallocate ( work  ) ;
            Memory_Deallocate ( x     ) ;
	    if ( rhs->stride != 1 ) Memory_Deallocate ( t ) ;
        }
        /* . Status. */
             if ( ifail <  0 ) { Status_SafeSet ( status, Status_InvalidArgument ) ; }
        else if ( ifail == 0 ) { Status_SafeSet ( status, Status_Continue        ) ; }
        else                   { Status_SafeSet ( status, Status_SingularMatrix  ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sparsity (as a percentage).
! . Sparsity is defined as the number of items with a magnitude less than a given tolerance.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_Sparsity ( const SymmetricMatrix *self, const Real tolerance )
{
    Real sparsity = 1.0e+00 ;
    if ( ( self != NULL ) && ( self->dimension > 0 ) )
    {
        /* . Off-diagonal elements count twice. */
        auto Integer i, ij, j, n = 0 ;
        for ( i = 0, ij = 0 ; i < self->dimension ; i++, ij++ )
        {
            for ( j = 0 ; j < i ; ij++, j++ )
            {
                if ( fabs ( self->data[ij] ) <= tolerance ) n += 2 ;
            }
            if ( fabs ( self->data[ij] ) <= tolerance ) n += 1 ;
        }
        sparsity = ( ( Real ) n ) / ( ( Real ) ( self->dimension * self->dimension ) ) ;
    }
    return 100.0e+00 * sparsity ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Square-root from eigenValues and eigenVectors.
! . The eigenValues are changed to their square roots!
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SymmetricMatrix_SquareRoot ( SymmetricMatrix *self, Real1DArray *eigenValues, const Real2DArray *eigenVectors )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( eigenValues != NULL ) && ( eigenVectors != NULL ) )
    {
        if ( ( self->dimension == eigenVectors->length0 ) && ( eigenValues->length == eigenVectors->length1 ) )
        {
            auto Real f, *ivector ;
            auto Integer    i, n = self->dimension ;
            SymmetricMatrix_Set ( self, 0.0e+00 ) ;
            for ( i = 0 ; i < n ; i++ )
            {
                /* . Make the square root of the eigenValue and store. */
                f = sqrt ( Maximum ( Real1DArray_Item ( eigenValues, i ), 0.0e+00 ) ) ;
                Real1DArray_Item ( eigenValues, i ) = f ;
                /* . Update the overlap. */
                ivector = Real2DArray_ItemPointer ( eigenVectors, 0, i ) ;
                cblas_dspr ( CblasColMajor, CblasUpper, n, f, ivector, n, self->data ) ;
            }
            status = Status_Success ;
        }
        else status = Status_DimensionError ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the product of two symmetric matrices as another symmetric matrix.
! . Equivalent to ( A * B + B * A ) / 2.
!--------------------------------------------------------------------- -------*/
Status SymmetricMatrix_SymmetricMultiply2 ( const SymmetricMatrix *self, const SymmetricMatrix *other, SymmetricMatrix *result )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( other != NULL ) && ( result != NULL ) )
    {
        if ( ( self->dimension == other->dimension ) && ( result->dimension == self->dimension ) )
        {
            auto Real  sum ;
            auto Integer     i, ij, ik, j, jk, k, n ;
            /* . ( Aij * Bjk + Bij * Ajk ) / 2 = Cik. */
            n = self->dimension ;
            for ( i = ik = 0 ; i < n ; i++ )
            {
                for ( k = 0 ; k <= i ; ik++, k++ )
                {
                    for ( j = 0, sum = 0.0e+00 ; j < n ; j++ )
                    {
                        if ( i < j ) ij = ( j * ( j + 1 ) ) / 2 + i ;
                        else         ij = ( i * ( i + 1 ) ) / 2 + j ;
                        if ( k < j ) jk = ( j * ( j + 1 ) ) / 2 + k ;
                        else         jk = ( k * ( k + 1 ) ) / 2 + j ;
                        sum += ( self->data[ij] * other->data[jk] + other->data[ij] * self->data[jk] ) ;
                    }
                    result->data[ik] = 0.5e+00 * sum ;
                }
            }
            status = Status_Success ;
        }
        else status = Status_DimensionError ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform a matrix by another symmetric matrix (i.e. calculate B * A * B).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_SymmetricTransform ( const SymmetricMatrix *A, const SymmetricMatrix *B, SymmetricMatrix *BAB, Status *status )
{
    if ( ( A != NULL ) && ( B != NULL ) && ( BAB != NULL ) )
    {
        auto Integer      n ;
        auto Real1DArray *t ;
        n = A->dimension ;
        t = Real1DArray_Allocate ( n, status ) ;
        if ( t != NULL )
        {
            auto Integer i, ij, il, j, jk, k, kl, l ;
            auto Real    sum ;
            for ( i = il = 0 ; i < n ; i++ )
            {
                for ( k = 0 ; k < n ; k++ )
                {
                    sum = 0.0e+00 ;
                    for ( j = 0 ; j < n ; j++ )
                    {
                        if ( i >= j ) ij = ( i * ( i + 1 ) ) / 2 + j ;
                        else          ij = ( j * ( j + 1 ) ) / 2 + i ;
                        if ( j >= k ) jk = ( j * ( j + 1 ) ) / 2 + k ;
                        else          jk = ( k * ( k + 1 ) ) / 2 + j ;
                        sum += B->data[ij] * A->data[jk] ;
                    }
                    Real1DArray_Item ( t, k ) = sum ;
                }
                for ( l = 0 ; l <= i ; il++, l++ )
                {
                    sum = 0.0e+00 ;
                    for ( k = 0 ; k < n ; k++ )
                    {
                        if ( k >= l ) kl = ( k * ( k + 1 ) ) / 2 + l ;
                        else          kl = ( l * ( l + 1 ) ) / 2 + k ;
                        sum += Real1DArray_Item ( t, k ) * B->data[kl] ;
                    }
                    BAB->data[il] = sum ;
                }
            }
            Real1DArray_Deallocate ( &t ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the trace.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_Trace ( const SymmetricMatrix *self )
{
    Real trace = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer i, ii ;
        for ( i = 0 ; i < self->dimension ; i++ )
        {
            ii = ( i * ( i + 3 ) ) / 2 ;
            trace += self->data[ii] ;
        }
    }
    return trace ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the trace of the product of two matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_Trace2 ( const SymmetricMatrix *self, const SymmetricMatrix *other, Status *status )
{
    Real value = 0.0e+00 ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( self->dimension == other->dimension )
        {
            auto Integer i, ij, j ;
            for ( i = 0, ij = 0 ; i < self->dimension ; i++, ij++ )
            {
               for ( j = 0 ; j < i ; ij++, j++ ) value += self->data[ij] * other->data[ij] ;
               value += 0.5e+00 * self->data[ij] * other->data[ij] ;
            }
            value *= 2.0e+00 ;
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transfer the elements of the second matrix to those of the first.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Transfer ( SymmetricMatrix *matrix1, SymmetricMatrix *matrix2 )
{
   if ( ( matrix1 != NULL ) && ( matrix2 != NULL ) )
   {
      auto Integer i ;
      for ( i = 0 ; i < matrix1->size ; i++ ) matrix1->data[i] = matrix2->data[i] ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform by a matrix (B^T * A * B) or its transpose (B * A * B^T).
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SymmetricMatrix_Transform ( const SymmetricMatrix *self, const Real2DArray *matrix, const Boolean QTRANSPOSE, SymmetricMatrix *result )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) )
    {
        auto Real *tmpvec = NULL ;
        /* . Allocate space. */
        tmpvec = Memory_Allocate_Array_Real_Initialize ( self->dimension, 0.0e+00 ) ;
        if ( tmpvec == NULL ) status = Status_OutOfMemory ;
        else
        {
            auto Integer i, il, l ;
            if ( QTRANSPOSE )
            {
                if ( ( self->dimension == matrix->length1 ) && ( result->dimension == matrix->length0 ) )
                {
                    for ( i = 0, il = 0 ; i < matrix->length0 ; i++ )
                    {
                        cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real2DArray_ItemPointer ( matrix, i, 0 ), matrix->stride1, 0.0e+00, tmpvec, 1 ) ;
                        for ( l = 0 ; l <= i ; il++, l++ ) result->data[il] = cblas_ddot ( self->dimension, tmpvec, 1, Real2DArray_ItemPointer ( matrix, l, 0 ), matrix->stride1 ) ;
                    }
                    status = Status_Success ;
                }
                else status = Status_DimensionError ;
            }
            else
            {
                if ( ( self->dimension == matrix->length0 ) && ( result->dimension == matrix->length1 ) )
                {
                    for ( i = 0, il = 0 ; i < matrix->length1 ; i++ )
                    {
                        cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real2DArray_ItemPointer ( matrix, 0, i ), matrix->stride0, 0.0e+00, tmpvec, 1 ) ;
                        for ( l = 0 ; l <= i ; il++, l++ ) result->data[il] = cblas_ddot ( self->dimension, tmpvec, 1, Real2DArray_ItemPointer ( matrix, 0, l ), matrix->stride0 ) ;
                    }
                    status = Status_Success ;
                }
                else status = Status_DimensionError ;
            }
            /* . Finish up. */
            Memory_Deallocate ( tmpvec ) ;
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform by a matrix (B^T * A * B) or its transpose (B * A * B^T).
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . self needs to be compact. */
Status SymmetricMatrix_TransformN ( const SymmetricMatrix *self, const Real2DArray *matrix, const Boolean useTranspose, SymmetricMatrix *result )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) )
    {
        auto Real1DArray *tmpvec = NULL ;
        /* . Allocate space. */
        tmpvec = Real1DArray_Allocate ( self->dimension, &status ) ;
        Real1DArray_Set ( tmpvec, 0.0e+00 ) ;
        if ( tmpvec == NULL ) Status_Set ( &status, Status_MemoryAllocationFailure ) ;
        else
        {
            auto Integer i, il, l ;
            if ( useTranspose )
            {
                if ( ( self->dimension == matrix->length1 ) && ( result->dimension == matrix->length0 ) )
                {
                    for ( i = 0, il = 0 ; i < matrix->length0 ; i++ )
                    {
                        cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real2DArray_ItemPointer ( matrix, i, 0 ), matrix->stride1, 0.0e+00, Real1DArray_Data ( tmpvec ), tmpvec->stride ) ;

                        for ( l = 0 ; l <= i ; il++, l++ ) result->data[il] = cblas_ddot ( self->dimension, Real1DArray_Data ( tmpvec ), tmpvec->stride, Real2DArray_ItemPointer ( matrix, l, 0 ), matrix->stride1 ) ;
                    }
                    Status_Set ( &status, Status_Continue ) ;
                }
                else Status_Set ( &status, Status_DimensionError ) ;
            }
            else
            {
                if ( ( self->dimension == matrix->length0 ) && ( result->dimension == matrix->length1 ) )
                {
                    for ( i = 0, il = 0 ; i < matrix->length1 ; i++ )
                    {
                        cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real2DArray_ItemPointer ( matrix, 0, i ), matrix->stride0, 0.0e+00, Real1DArray_Data ( tmpvec ), tmpvec->stride ) ;
                        for ( l = 0 ; l <= i ; il++, l++ ) result->data[il] = cblas_ddot ( self->dimension, Real1DArray_Data ( tmpvec ), tmpvec->stride, Real2DArray_ItemPointer ( matrix, 0, l ), matrix->stride0 ) ;
                    }
                    Status_Set ( &status, Status_Continue ) ;
                }
                else Status_Set ( &status, Status_DimensionError ) ;
            }
            /* . Finish up. */
            Real1DArray_Deallocate ( &tmpvec ) ;
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform a symmetric matrix, S, by the columns of a  matrix, M, to give Mi^T * S * Ml.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_TransformByMatrixColumns ( SymmetricMatrix *self, Real2DArray *matrix, const Integer i, const Integer l )
{
    Real value = 0.0 ;
    if ( ( self != NULL ) && ( matrix != NULL ) )
    {
        auto Integer j, jk, k ;
        auto Real  mki, sum, *tmpvec ;

        /* . Allocate space. */
        tmpvec = Memory_Allocate_Array_Real_Initialize ( matrix->length0, 0.0e+00 ) ;

        /* . Loop over the ith row of M^T. */
        /* . Form the (ik)th elements of (M^T)ij * Sjk. */
        for ( k = 0, jk = 0 ; k < matrix->length0 ; k++ )
        {
            mki = Real2DArray_Item ( matrix, k, i ) ;
            sum = 0.0e+00 ;
            for ( j = 0 ; j < k ; j++, jk++ )
            {
                tmpvec[j] += self->data[jk] * mki ;
                sum       += self->data[jk] * Real2DArray_Item ( matrix, j, i ) ;
            }
            tmpvec[k] = sum + self->data[jk] * mki ;
            jk++ ;
        }

        /* . Calculate the (il)th elements of the transformed matrix. */
        value = cblas_ddot ( matrix->length0, tmpvec, 1, Real2DArray_ItemPointer ( matrix, 0, l ), matrix->stride0 ) ;

        /* . Deallocate space. */
        Memory_Deallocate ( tmpvec ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform a symmetric matrix, S, by a matrix, M, to give M^T * S * M.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Transform_In_Place ( SymmetricMatrix *self, Real2DArray *matrix )
{
   if ( ( self != NULL ) && ( matrix != NULL ) )
   {
      auto Integer i, il, l ;
      auto Real  *result, *tmpvec ;

      /* . Allocate space. */
      result = Memory_Allocate_Array_Real ( self->sizeP ) ;
      tmpvec = Memory_Allocate_Array_Real_Initialize ( matrix->length0, 0.0e+00 ) ;

      /* . Loop over the rows of M^T. */
      for ( i = 0, il = 0 ; i < matrix->length1 ; i++ )
      {
         /* . Form the (ik)th elements of (M^T)ij * Sjk. */
         cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real2DArray_ItemPointer ( matrix, 0, i ), matrix->stride0, 0.0e+00, tmpvec, 1 ) ;
/*
         for ( k = 0, jk = 0 ; k < matrix->length0 ; k++ )
         {
            mki = Real2DArray_Item ( matrix, k, i ) ;
            sum = 0.0e+00 ;
            for ( j = 0 ; j < k ; j++, jk++ )
            {
               tmpvec[j] += self->data[jk] * mki ;
               sum       += self->data[jk] * Real2DArray_Item ( matrix, j, i ) ;
            }
            tmpvec[k] = sum + self->data[jk] * mki ;
            jk++ ;
         }

void SymmetricMatrix_Multiply_By_Vector ( const SymmetricMatrix *matrix, const Vector *other, Vector *result )

      cblas_dspmv ( CblasColMajor, CblasUpper, matrix->dimension, 1.0e+00, matrix->data, Real1DArray_Item ( other, 0 ), 1, 0.0e+00,  Real1DArray_Item ( result, 0 ), 1 ) ;
*/

         /* . Fill the (il)th elements of the transformed matrix. */
         for ( l = 0 ; l <= i ; il++, l++ ) result[il] = cblas_ddot ( matrix->length0, tmpvec, 1, Real2DArray_ItemPointer ( matrix, 0, l ), matrix->stride0 ) ;
# ifdef __RANGE_CHECKING
if ( il > self->size ) printf ( "SymmetricMatrix_Transform_In_Place: il out of range: %d\n", il ) ;
# endif
      }

      /* . Deallocate space. */
      Memory_Deallocate ( self->data ) ;
      Memory_Deallocate ( tmpvec ) ;

      /* . Juggle the variables in S. */
      self->data      = result ;
      self->dimension = matrix->length1 ;
      self->size      = ( matrix->length1 * ( matrix->length1 + 1 ) ) / 2 ;
   }
}
# undef __RANGE_CHECKING

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unweighting.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Unweight ( SymmetricMatrix *self )
{
    if ( self != NULL )
    {
        /* . Off-diagonal elements only. */
        auto Integer i, ij, j ;
        for ( i = 0, ij = 0 ; i < self->dimension ; i++, ij++ )
        {
            for ( j = 0 ; j < i ; ij++, j++ ) self->data[ij] *= 0.5e+00 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Update a symmetric matrix using an appropriate updating formula.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define UPDATE_TOLERANCE 1.0e-08
Status SymmetricMatrix_Update ( SymmetricMatrix *self, const Real1DArray *dx, const Real1DArray *dg, const SYMMETRICMATRIXUPDATING_OPTION option, const Real *tolerance )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( dx != NULL ) && ( dg != NULL ) )
    {
        auto Real       aa, aaxx, ax, fac, gx, mix, tol, xx ;
	auto Real1DArray *a = NULL ;
        /* . Get the tolerance. */
        if ( tolerance == NULL ) tol = UPDATE_TOLERANCE ;
        else                     tol = (*tolerance)     ;
        /* . Allocate space. */
	a = Real1DArray_Allocate ( self->dimension, NULL ) ;
        if ( a != NULL )
        {
            /* . Get A * dx. */
            cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real1DArray_Data ( dx ), dx->stride, 0.0e+00, Real1DArray_Data ( a ), a->stride ) ;

            /* . For non-BFGS options get g - A * dx. */
            if ( option != SYMMETRICMATRIXUPDATING_BFGS )
            {
                Real1DArray_Scale ( a, -1.0e+00 ) ;
                Real1DArray_Add   ( a, dg, NULL ) ;
            }

            /* . Calculate some dot products. */
            ax = Real1DArray_Dot ( a, dx, NULL ) ;

            /* . BFGS option. */
            if ( option == SYMMETRICMATRIXUPDATING_BFGS )
            {
                gx = Real1DArray_Dot ( dg, dx, NULL ) ;
	        if ( fabs ( ax ) > tol ) cblas_dspr ( CblasColMajor, CblasUpper, self->dimension, ( - 1.0e+00 / ax ), Real1DArray_Data ( a  ), a->stride , self->data ) ;
	        if ( fabs ( gx ) > tol ) cblas_dspr ( CblasColMajor, CblasUpper, self->dimension, (   1.0e+00 / gx ), Real1DArray_Data ( dg ), dg->stride, self->data ) ;
            }
            /* . Other options. */
            else
            {

                /* . Calculate some dot products. */
                aa = Real1DArray_Dot (  a,  a, NULL ) ;
                xx = Real1DArray_Dot ( dx, dx, NULL ) ;
                aaxx = aa * xx ;

                /* . Bofill or MS. */
                if ( ( option == SYMMETRICMATRIXUPDATING_BOFILL ) || ( option == SYMMETRICMATRIXUPDATING_MS ) )
                {
                    if ( option == SYMMETRICMATRIXUPDATING_BOFILL )
                    {
                        if ( fabs ( aaxx ) > tol ) fac = ax / aaxx ;
                        else                       fac = 0.0e+00   ;
                    }
                    else
                    {
                        if ( fabs ( ax  ) > tol ) fac = 1.0e+00 / ax ;
                        else                      fac = 0.0e+00      ;
                    }
                    cblas_dspr ( CblasColMajor, CblasUpper, self->dimension, fac, Real1DArray_Data ( a ), a->stride, self->data ) ;
                }

                /* . Bofill or Powell. */
                if ( ( option == SYMMETRICMATRIXUPDATING_BOFILL ) || ( option == SYMMETRICMATRIXUPDATING_POWELL ) )
                {
                    /* . Get the mixing factor. */
                    if ( ( option == SYMMETRICMATRIXUPDATING_BOFILL ) && ( fabs ( aaxx ) > tol ) ) mix = 1.0e+00 - ax * ax / aaxx ;
                    else                                                                           mix = 1.0e+00                  ;
                    if ( fabs ( xx * xx ) > tol ) cblas_dspr  ( CblasColMajor, CblasUpper, self->dimension, ( - mix * ax / ( xx * xx ) ),
                                                                                                   Real1DArray_Data ( dx ), dx->stride, self->data ) ;
                    if ( fabs ( xx      ) > tol ) cblas_dspr2 ( CblasColMajor, CblasUpper, self->dimension, (   mix     /    xx        ),
                                                                Real1DArray_Data ( a ), a->stride, Real1DArray_Data ( dx ), dx->stride, self->data ) ;
                }
            }

            /* . Finish up. */
	    Real1DArray_Deallocate ( &a ) ;
            status = Status_Success ;
        }
        else status = Status_OutOfMemory ;
    }
    return status ;
}
# undef UPDATE_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Update the hessian using the BFGS formula.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define UPDATE_TOLERANCE 1.0e-08
void SymmetricMatrix_BFGS_Update ( SymmetricMatrix *hessian, const Real1DArray *dx, const Real1DArray *dg )
{
    if ( ( hessian != NULL ) && ( dx != NULL ) && ( dg != NULL ) )
    {
	auto Real         dxdg, dxhdx ;
	auto Real1DArray *hdx ;
	hdx = Real1DArray_Allocate ( hessian->dimension, NULL ) ;
	cblas_dspmv ( CblasColMajor, CblasUpper, hessian->dimension, 1.0e+00, hessian->data, Real1DArray_Data ( dx ), dx->stride, 0.0e+00,  Real1DArray_Data ( hdx ), hdx->stride ) ;
	dxhdx = Real1DArray_Dot  ( dx, hdx, NULL ) ;
	dxdg  = Real1DArray_Dot  ( dx, dg,  NULL ) ;
	if ( fabs ( dxdg )  > UPDATE_TOLERANCE ) cblas_dspr ( CblasColMajor, CblasUpper, hessian->dimension, (   1.0e+00 / dxdg  ), Real1DArray_Data ( dg  ), dg->stride , hessian->data ) ;
	if ( fabs ( dxhdx ) > UPDATE_TOLERANCE ) cblas_dspr ( CblasColMajor, CblasUpper, hessian->dimension, ( - 1.0e+00 / dxhdx ), Real1DArray_Data ( hdx ), hdx->stride, hessian->data ) ;
	Real1DArray_Deallocate ( &hdx ) ;
    }
}
# undef UPDATE_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Multiply by a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . self needs to be compact. */
void SymmetricMatrix_VectorMultiply ( const SymmetricMatrix *self, const Real1DArray *other, Real1DArray *result, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && ( result != NULL ) )
    {
        if ( ( self->dimension == other->length ) && ( other->length == result->length ) )
        {
            cblas_dspmv ( CblasColMajor, CblasUpper, self->dimension, 1.0e+00, self->data, Real1DArray_Data ( other ), other->stride, 0.0e+00, Real1DArray_Data ( result ), result->stride ) ;
        }
        else Status_Set ( status, Status_InvalidDimension ) ;
    }
}
/*
void SymmetricMatrix_Multiply_By_Vector ( const SymmetricMatrix *matrix, const Vector *other, Vector *result )
{
   if ( ( matrix != NULL ) && ( other != NULL ) && ( result != NULL ) )
   {
      cblas_dspmv ( CblasColMajor, CblasUpper, matrix->dimension, 1.0e+00, matrix->data, Real1DArray_Item ( other, 0 ), 1, 0.0e+00,  Real1DArray_Item ( result, 0 ), 1 ) ;
   }
}
*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a matrix view of a raw array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_ViewOfRaw ( SymmetricMatrix *self, const Integer length, const Integer stride, Real *data, const Integer rawstride, Status *status )
{
    if ( ( self != NULL ) && ( data != NULL ) )
    {
        auto Integer n ;
        auto Status  localStatus ;
        Status_Set ( &localStatus, Status_Continue ) ;
        n = SliceIndices_GetLength ( 0, length, stride, length, &localStatus ) ;
        if ( ! Status_OK ( &localStatus ) ) Status_Set ( status, localStatus ) ;
        else
        {
/*
            auto Integer fullstride ;
            fullstride = rawstride * stride ;
*/
            self->data       = data ;
            self->dimension  = n    ;
            self->size       = ( n * ( n + 1 ) ) / 2 ;
            self->dimensionP = self->dimension ;
            self->sizeP      = self->size      ;
            /* . Basic data. */
/*
            self->data    = data ;
            self->isOwner = False ;
            self->isView  = True  ;
            self->length  = n ;
            self->offset  = 0 ;
            self->size    = ( n * ( n + 1 ) ) * fullstride / 2 ;
            self->stride  = fullstride ;
*/
            /*. View. */
/*
            self->view.data    = data ;
            self->view.isOwner = False ;
            self->view.isView  = True  ;
            self->view.length  = ( n * ( n + 1 ) ) / 2 ;
            self->view.offset  = 0 ;
            self->view.size    = self->size ;
            self->view.stride  = fullstride ;
*/
        }
    }
}
