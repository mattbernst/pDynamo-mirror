/*------------------------------------------------------------------------------
! . File      : Real2DArray.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . 2-D arrays.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "cblas.h"
# include "f2clapack.h"

# include "Boolean1DArray.h"
# include "Macros.h"
# include "Memory.h"
# include "Real2DArray.h"
# include "SliceOld.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*# define CHECKMULTIPLICATIONS*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a 1D slice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_1DSlice ( const Real2DArray *self, const Integer start0, const Integer stop0, const Integer stride0, const Integer start1, const Integer stop1, const Integer stride1, Real1DArray *slice, Status *status )
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
            slice->offset  = Real2DArray_ItemIndex ( self, start0, start1 ) ;
            slice->size    = self->size ;
            if ( length0 == 1 ) slice->stride = self->stride1 * stride1 ;
            else                slice->stride = self->stride0 * stride0 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Absolute maximum.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real2DArray_AbsoluteMaximum ( const Real2DArray *self )
{
    Real value = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        auto Real1DArray sview ;
        for ( i = 0 ; i < self->length0 ; i++ )
        {
            Real2DArray_RowSlice ( self, i, &sview, NULL ) ;
            value = Maximum ( value, Real1DArray_AbsoluteMaximum ( &sview ) ) ;
        }
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add a scaled array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_AddScaledArray ( Real2DArray *self, const Real value, const Real2DArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && ( value != 0.0e+00 ) )
    {
        if ( ( self->length0 == other->length0 ) && ( self->length1 == other->length1 ) )
        {
            auto Integer i ;
            auto Real1DArray oview, sview ;
            for ( i = 0 ; i < self->length0 ; i++ )
            {
                Real2DArray_RowSlice ( other, i, &oview, NULL ) ;
                Real2DArray_RowSlice ( self , i, &sview, NULL ) ;
                Real1DArray_AddScaledArray ( &sview, value, &oview, NULL ) ;
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real2DArray *Real2DArray_Allocate ( const Integer length0, const Integer length1, Status *status )
{
    Real2DArray *self = NULL ;
    MEMORY_ALLOCATE ( self, Real2DArray ) ;
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
            MEMORY_ALLOCATEARRAY ( self->data, self->length, Real ) ;
            if ( self->data == NULL ) Real2DArray_Deallocate ( &self ) ;
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
Real2DArray *Real2DArray_Clone ( const Real2DArray *self, Status *status )
{
    Real2DArray *clone = NULL ;
    if ( self != NULL )
    {
        clone = Real2DArray_Allocate ( self->length0, self->length1, status ) ;
        Real2DArray_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form a 1-D array from the dot products of the columns of 2 2-D arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_ColumnDotProducts ( const Boolean      initialize ,
                                     const Real2DArray *a          ,
                                     const Real2DArray *b          ,
                                           Real1DArray *c          )
{
    Integer p ;
    Real1DArray sliceA, sliceB ;
    if ( initialize ) Real1DArray_Set ( c, 0.0e+00 ) ;
    for ( p = 0 ; p < c->length ; p++ )
    {
        Real2DArray_ColumnSlice ( a, p, &sliceA, NULL ) ;
        Real2DArray_ColumnSlice ( b, p, &sliceB, NULL ) ;
        Real1DArray_Item ( c, p ) += Real1DArray_Dot ( &sliceA, &sliceB, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a slice of a column.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_ColumnSlice ( const Real2DArray *self, const Integer column, Real1DArray *slice, Status *status )
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
! . Fill a matrix by row from a vector.
! . The vector will be the same size or bigger than the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Real2DArray_CopyFromArrayByRow ( Real2DArray *self, const Real1DArray *vector, const Selection *selection )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( vector != NULL ) )
    {
	auto Integer d0, d1, l1 ;
	d0 = self->length0 ;
	d1 = self->length1 ;
        l1 = self->stride1 ;
        /* . Null selection implies all rows to be copied. */
	if ( selection == NULL )
	{
	    if ( self->length <= vector->length )
            {
        	if ( self->length > 0 )
        	{
		    auto Integer i, n ;
		    /* . Loop over rows. */
		    for ( i = 0, n = 0 ; i < d0 ; i++, n += d1 ) cblas_dcopy ( d1, Real1DArray_ItemPointer ( vector, n ), vector->stride, Real2DArray_RowPointer ( self, i ), l1 ) ;
        	}
        	status = Status_Success ;
            }
            else status = Status_DimensionError ;
	}
	/* . Selection present. */
	else
	{
	    /* . Checks are made on space during the loop and not beforehand. */
	    status = Status_Success ;
	    if ( self->length > 0 )
	    {
	        auto Integer i, n, s ;
	        /* . Loop over rows. */
		for ( s = 0, n = 0 ; s < selection->nindices ; s++ )
		{
                    i = selection->indices[s] ;
		    if ( ( s >= 0 ) && ( s < d0 ) )
		    {
		        /* . Copy row. */
			if ( ( n + d1 ) > vector->length )
			{
			    status = Status_DimensionError ;
			    break ;
			}
			else
			{
			    cblas_dcopy ( d1, Real1DArray_ItemPointer ( vector, n ), vector->stride, Real2DArray_RowPointer ( self, i ), l1 ) ;
			    n += d1 ;
			}
		    }
		}
	    }
	}
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
! . As much data as possible is copied.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_CopyTo ( const Real2DArray *self, Real2DArray *other, Status *status )
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
            auto Integer i ;
            for ( i = 0 ; i < n0 ; i++ ) cblas_dcopy ( n1, Real2DArray_ItemPointer ( self, i, 0 ), self->stride1, Real2DArray_ItemPointer ( other, i, 0 ), other->stride1 ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Empty a matrix by row into a vector.
! . The vector needs to be big enough to hold all components.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Real2DArray_CopyToArrayByRow ( const Real2DArray *self, Real1DArray *vector, const Selection *selection )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( vector != NULL ) )
    {
	auto Integer d0, d1, l1 ;
	d0 = self->length0 ;
	d1 = self->length1 ;
	l1 = self->stride1 ;
        /* . Null selection implies all rows to be copied. */
	if ( selection == NULL )
	{
	    if ( self->length <= vector->length )
            {
        	if ( self->length > 0 )
        	{
		    auto Integer i, n ;
		    /* . Loop over rows. */
		    for ( i = 0, n = 0 ; i < d0 ; i++, n += d1 ) cblas_dcopy ( d1, Real2DArray_RowPointer ( self, i ), l1, Real1DArray_ItemPointer ( vector, n ), vector->stride ) ;
        	}
        	status = Status_Success ;
            }
            else status = Status_DimensionError ;
	}
	/* . Selection present. */
	else
	{
	    /* . Checks are made on space during the loop and not beforehand. */
	    status = Status_Success ;
	    if ( self->length > 0 )
	    {
	        auto Integer i, n, s ;
	        /* . Loop over rows. */
		for ( s = 0, n = 0 ; s < selection->nindices ; s++ )
		{
                    i = selection->indices[s] ;
		    if ( ( s >= 0 ) && ( s < d0 ) )
		    {
		        /* . Copy row. */
			if ( ( n + d1 ) > vector->length )
			{
			    status = Status_DimensionError ;
			    break ;
			}
			else
			{
			    cblas_dcopy ( d1, Real2DArray_RowPointer ( self, i ), l1, Real1DArray_ItemPointer ( vector, n ), vector->stride ) ;
			    n += d1 ;
			}
		    }
		}
	    }
	}
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_Deallocate ( Real2DArray **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        if ( (*self)->isOwner ) MEMORY_DEALLOCATE ( (*self)->data ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the determinant of a square matrix using LUP factorization.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . self is destroyed as it is overwritten with the LUP factorization.
*/
Real Real2DArray_Determinant ( Real2DArray *self, Status *status )
{
    Real determinant = 0.0e+00 ;
    if ( self != NULL )
    {
        /* . The matrix must be square. */
        if ( self->length0 == self->length1 )
        {
            auto Integer         current, cycleLength, i, n, p, start, t ;
            auto Real            parity ;
            auto Boolean1DArray *isVisited ;
            auto Integer1DArray *permutation, *pivots ;

            /* . Allocate space. */
            n = self->length0 ;
            isVisited   = Boolean1DArray_Allocate ( n, status ) ;
            permutation = Integer1DArray_Allocate ( n, status ) ;
            pivots      = Integer1DArray_Allocate ( n, status ) ;

            /* . Do the factorization. */
            Real2DArray_LUPFactorizationInPlace ( self, pivots, status ) ;

            /* . Find the determinant of U. */
            determinant = 1.0e+00 ;
            for ( i = 0 ; i < n ; i++ ) determinant *= Real2DArray_Item ( self, i, i ) ;

            /* . Convert pivots to a permutation. */
            for ( i = 0 ; i < n   ; i++ ) Integer1DArray_Item ( permutation, i ) = i ;
            for ( i = 0 ; i < n-1 ; i++ )
            {
                p = Integer1DArray_Item ( pivots,      i ) - 1 ;
                t = Integer1DArray_Item ( permutation, p ) ;
                Integer1DArray_Item ( permutation, p ) = Integer1DArray_Item ( permutation, i ) ;
                Integer1DArray_Item ( permutation, i ) = t ;
            }

            /* . Find the parity of the pivots. */
            parity = 1.0e+00 ;
            Boolean1DArray_Set ( isVisited, False ) ;
            for ( i = 0 ; i < n ; i++ )
            {
                if ( ! Boolean1DArray_Item ( isVisited, i ) )
                {
                    current     = Integer1DArray_Item ( permutation, i ) ;
                    cycleLength = 1 ;
                    start       = i ;
                    Boolean1DArray_Item ( isVisited, i ) = True ;
                    while ( current != start )
                    {
                        Boolean1DArray_Item ( isVisited, current ) = True ;
                        current      = Integer1DArray_Item ( permutation, current ) ;
                        cycleLength += 1 ;
                    }
                    if ( cycleLength % 2 == 0 ) parity *= -1.0e+00 ;
                }
            }
            determinant *= parity ;

            /* . Finish up. */
            Boolean1DArray_Deallocate ( &isVisited   ) ;
            Integer1DArray_Deallocate ( &permutation ) ;
            Integer1DArray_Deallocate ( &pivots      ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return determinant ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Element-wise exponentiation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_Exp ( Real2DArray *self )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i, j ;
        auto Real    v ;
        for ( i = 0 ; i < self->length0 ; i++ )
        {
            for ( j = 0 ; j < self->length1 ; j++ )
            {
                v = Real2DArray_Item ( self, i, j ) ;
                if ( v > Real_MaximumExponent ) v = Real_Huge ;
                else                            v = exp ( v ) ;
                Real2DArray_Item ( self, i, j ) = v ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Frobenius norm.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real2DArray_FrobeniusNorm ( const Real2DArray *self )
{
    Real norm = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i ;
        auto Real1DArray sview ;
        for ( i = 0 ; i < self->length0 ; i++ )
        {
            Real2DArray_RowSlice ( self, i, &sview, NULL ) ;
            norm += Real1DArray_Dot ( &sview, &sview, NULL ) ;
        }
        norm = sqrt ( norm ) ;
    }
    return norm ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real2DArray_GetItem ( const Real2DArray *self, const Integer i, const Integer j, Status *status )
{
    Real value = 0.0e+00 ;
    if ( self != NULL )
    {
        if ( self->length > 0 )
        {
            if ( ( i >= 0 ) && ( i < self->length0 ) && ( j >= 0 ) && ( j < self->length1 ) ) value = Real2DArray_Item ( self, i, j ) ;
            else Status_Set ( status, Status_IndexOutOfRange ) ;
        }
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gram-Schmidt orthogonalize a set of vectors stored columnwise in-place.
! . A modified (as opposed to classical) iterative algorithm is employed.
! . The number of orthogonal vectors is returned (<= old number).
! . There is an option to treat the first numberConstant vectors as already orthogonalized.
! . The tolerance is the size of norm2 per element.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DEBUG_GSO
# define DEFAULT_TOLERANCE 1.0e-10
Integer Real2DArray_GramSchmidtOrthogonalize ( const Real2DArray *self, const Integer *maximumIterations, const Integer *numberConstant, const Real *tolerance, Status *status )
{
    Integer numberOrthogonalized = 0 ;
    if ( self != NULL )
    {
        auto Integer nstart, nvectors ;

        /* . Get various indices. */
        nvectors = self->length1 ;
        if ( numberConstant != NULL ) nstart = Maximum ( 0, (*numberConstant) ) ;
        else                          nstart = 0 ;
        if ( nstart < nvectors )
        {
            auto Integer     i, iteration, j, ncurrent, niterations ;
            auto Real        delta, factor ;
            auto Real1DArray ivector, jvector ;

            /* . Get the number of iterations. */
            if ( maximumIterations != NULL ) niterations = Maximum ( 1, (*maximumIterations) ) ;
            else                             niterations = 1 ;

            /* . Get the tolerance for normalization. */
            if ( tolerance == NULL ) delta = DEFAULT_TOLERANCE     ;
            else                     delta = fabs ( (*tolerance) ) ;
            delta *= sqrt ( ( Real ) self->length0 ) ;

            /* . Loop over vectors to be orthogonalized. */
            for ( i = ncurrent = nstart ; i < nvectors ; i++ )
            {
                Real2DArray_ColumnSlice ( self, i, &ivector, status ) ;

                /* . Loop over iterations. */
                for ( iteration = 0 ; iteration < niterations ; iteration++ )
                {
	            /* . Loop over vectors to be orthogonalized against. */
                    for ( j = 0 ; j < ncurrent ; j++ )
                    {
                        Real2DArray_ColumnSlice    ( self, j,  &jvector, status ) ;
                        factor = Real1DArray_Dot   ( &ivector, &jvector, status ) ;
                        Real1DArray_AddScaledArray ( &ivector, -factor, &jvector, status ) ;
                    }
                }

                /* . Normalization. If OK, scale and move the vector if necessary. */
                factor = Real1DArray_Norm2 ( &ivector ) ;
/*printf ( "\nVector, iteration, norm2 = %d %d %20.15f\n", i, iteration, factor ) ;*/
                if ( factor > delta )
                {
                    Real1DArray_Scale ( &ivector, 1.0e+00 / factor ) ;
                    if ( ncurrent != i )
                    {
                        Real2DArray_ColumnSlice ( self, ncurrent,  &jvector, status ) ;
                        Real1DArray_CopyTo      ( &ivector, &jvector, status ) ;
                    }
                    ncurrent++ ;
                    numberOrthogonalized++ ;
                }
            }

# ifdef DEBUG_GSO
{
    auto Real deviation, tolerance ;
    auto Real2DArray vectors ;
    tolerance = DEFAULT_TOLERANCE ;
    Real2DArray_Slice ( self, 0, self->length0, 1, 0, ncurrent, 1, &vectors, status ) ;
    if ( ! Real2DArray_IsOrthogonal ( &vectors, &tolerance, &deviation, status ) )
    {
/*printf ( "\nStart, Stop, Found = %d %d %d\n", nstart, nvectors, numberOrthogonalized ) ;*/
        printf ( "\nGram-Schmidt orthogonalized vectors not orthogonal. Deviation = %20.15f\n", deviation ) ;
    }
}
# endif
        }
    }
    return numberOrthogonalized ;
}
# undef DEBUG_GSO
# undef DEFAULT_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for compactness.
! . -1 indicates the full array.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean Real2DArray_IsCompact ( const Real2DArray *self, const Integer dimension, Status *status )
{
    Boolean isCompact = True ;
    if ( self != NULL )
    {
             if ( dimension == 0 ) isCompact = ( self->stride0 == self->length1 ) ;
        else if ( dimension == 1 ) isCompact = ( self->stride1 == 1             ) ;
        else
        {
            if ( dimension != -1 ) Status_Set ( status, Status_InvalidDimension ) ;
            isCompact = ( self->stride0 == self->length1 ) && ( self->stride1 == 1 ) ;
        }
    }
    return isCompact ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for orthogonality (check procedure).
!---------------------------------------------------------------------------------------------------------------------------------*/
# define ORTHOGONALITY_TOLERANCE 1.0e-10
Boolean Real2DArray_IsOrthogonal ( const Real2DArray *self, const Real *tolerance, Real *deviation, Status *status )
{
    Boolean isOrthogonal = False ;
    if ( deviation != NULL ) (*deviation) = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer      i ;
        auto Real         absmax, tol ;
        auto Real2DArray *r ;

        /* . Allocate space. */
        r = Real2DArray_Allocate ( self->length1, self->length1, status ) ;
        if ( r != NULL )
        {
            /* . Get the tolerance. */
            if ( tolerance == NULL ) tol = ORTHOGONALITY_TOLERANCE ;
            else                     tol = (*tolerance) ;

            /* . Get self^T * self - I. */
            Real2DArray_MatrixMultiply ( True, False, 1.0e+00, self, self, 0.0e+00, r, status ) ;
            for ( i = 0 ; i < r->length0 ; i++ ) Real2DArray_Item ( r, i, i ) -= 1.0e+00 ;

            /* . Find maximum deviation. */
            absmax       = Real2DArray_AbsoluteMaximum ( r ) ;
            isOrthogonal = ( absmax <= tol ) ;
            if ( deviation != NULL ) (*deviation) = absmax ;
        }

        /* . Finish up. */
        Real2DArray_Deallocate ( &r ) ;
    }
    return isOrthogonal ;
}
# undef ORTHOGONALITY_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for a symmetric matrix (check procedure).
!---------------------------------------------------------------------------------------------------------------------------------*/
# define SYMMETRIC_TOLERANCE 1.0e-10
Boolean Real2DArray_IsSymmetric ( const Real2DArray *self, const Real *tolerance, Real *deviation )
{
    Boolean isSymmetric = False ;
    if ( deviation != NULL ) (*deviation) = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->length0 == self->length1 ) )
    {
        auto Integer i, j ;
        auto Real    difference, tol ;

        /* . Get the tolerance. */
        if ( tolerance == NULL ) tol = SYMMETRIC_TOLERANCE ;
        else                     tol = (*tolerance) ;

        /* . Loop over indices. */
        difference = 0.0e+00 ;
        for ( i = 0 ; i < self->length0 ; i++ )
        {
            for ( j = 0 ; j < i ; j++ ) difference = Maximum ( difference, fabs ( Real2DArray_Item ( self, i, j ) - Real2DArray_Item ( self, j, i ) ) ) ;
        }

        /* . Find maximum difference. */
        isSymmetric = ( difference <= tol ) ;
        if ( deviation != NULL ) (*deviation) = difference ;
    }
    return isSymmetric ;
}
# undef SYMMETRIC_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for uniformness.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean Real2DArray_IsUniform ( const Real2DArray *self )
{
    Boolean isUniform = True ;
    if ( self != NULL ) isUniform = ( self->stride0 == self->length1 * self->stride1 ) ;
    return isUniform ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Real2DArray_Length ( const Real2DArray *self, const Integer dimension )
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
! . LUP factorization of a matrix in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . self is overwritten with L in the lower triangle and U in the upper triangle with diagonal.
! . The diagonal elements of L are 1.
*/
void Real2DArray_LUPFactorizationInPlace ( Real2DArray *self, Integer1DArray *pivots, Status *status )
{
    if ( ( self != NULL ) && ( pivots != NULL ) )
    {
        auto Integer info, m, n ;

        /* . Array dimensions. */
        info = 0 ;
        m    = self->length0 ;
        n    = self->length1 ;
        if ( ( Integer1DArray_Length ( pivots ) >= Minimum ( m, n ) ) && ( self->stride1 == 1 ) )
        {
            Integer1DArray_Set ( pivots, -1 ) ;
            /* . Reverse m and n due to C <-> Fortran ordering. */
            dgetrf_( &n, &m, Real2DArray_Data ( self ), &(self->stride0), Integer1DArray_Data ( pivots ), &info ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Matrix-matrix multiply. A straight interface to cblas_dgemm. C = alpha A B + beta C (A and B can be transposed).
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_MatrixMultiply ( const Boolean aTranspose, const Boolean bTranspose, const Real alpha, const Real2DArray *a, const Real2DArray *b, const Real beta, Real2DArray *c, Status *status )
{
    if ( ( a != NULL ) && ( b != NULL ) && ( c != NULL ) )
    {
        /* . All matrices need to be compact in dimension 1. */
        if ( ! ( Real2DArray_IsCompact ( a, 1, NULL ) && Real2DArray_IsCompact ( b, 1, NULL ) && Real2DArray_IsCompact ( c, 1, NULL ) ) ) Status_Set ( status, Status_NonCompactDimension ) ;
        else
        {
            auto Integer  k, m, n ;
            m = c->length0 ;
            n = c->length1 ;
	    if ( aTranspose ) k = a->length0 ;
	    else              k = a->length1 ;
            if ( ( ( ! aTranspose ) && ( ! bTranspose ) && ( m == a->length0 ) && ( n == b->length1 ) && ( a->length1 == b->length0 ) ) ||
                 ( ( ! aTranspose ) && (   bTranspose ) && ( m == a->length0 ) && ( n == b->length0 ) && ( a->length1 == b->length1 ) ) ||
                 ( (   aTranspose ) && ( ! bTranspose ) && ( m == a->length1 ) && ( n == b->length1 ) && ( a->length0 == b->length0 ) ) ||
                 ( (   aTranspose ) && (   bTranspose ) && ( m == a->length1 ) && ( n == b->length0 ) && ( a->length0 == b->length1 ) ) )
	    {
                auto enum CBLAS_TRANSPOSE aT, bT ;
# ifdef CHECKMULTIPLICATIONS
auto Integer i, j, l ;
auto Real    deviation, sum ;
auto Real2DArray *store ;
deviation = 0.0e+00 ;
store     = Real2DArray_Allocate ( c->length0, c->length1, status ) ;
Real2DArray_CopyTo ( c, store, status ) ;
# endif
                if ( aTranspose ) aT = CblasTrans   ;
                else              aT = CblasNoTrans ;
                if ( bTranspose ) bT = CblasTrans   ;
                else              bT = CblasNoTrans ;
                cblas_dgemm ( CblasRowMajor, aT, bT, m, n, k, alpha, Real2DArray_Data ( a ), a->stride0, Real2DArray_Data ( b ), b->stride0, beta, Real2DArray_Data ( c ), c->stride0 ) ;
# ifdef CHECKMULTIPLICATIONS
if ( aTranspose )
{
    if ( bTranspose )
    {
        for ( i = 0 ; i < c->length0 ; i++ )
        {
            for ( j = 0 ; j < c->length1 ; j++ )
            {
                sum = beta * Real2DArray_Item ( store, i, j ) - Real2DArray_Item ( c, i, j ) ;
                for ( l = 0 ; l < k ; l++ ) sum += alpha * Real2DArray_Item ( a, l, i ) * Real2DArray_Item ( b, j, l ) ;
                deviation = Maximum ( deviation, fabs ( sum ) ) ;
            }
        }
    }
    else
    {
        for ( i = 0 ; i < c->length0 ; i++ )
        {
            for ( j = 0 ; j < c->length1 ; j++ )
            {
                sum = beta * Real2DArray_Item ( store, i, j ) - Real2DArray_Item ( c, i, j ) ;
                for ( l = 0 ; l < k ; l++ ) sum += alpha * Real2DArray_Item ( a, l, i ) * Real2DArray_Item ( b, l, j ) ;
                deviation = Maximum ( deviation, fabs ( sum ) ) ;
            }
        }
    }
}
else
{
    if ( bTranspose )
    {
        for ( i = 0 ; i < c->length0 ; i++ )
        {
            for ( j = 0 ; j < c->length1 ; j++ )
            {
                sum = beta * Real2DArray_Item ( store, i, j ) - Real2DArray_Item ( c, i, j ) ;
                for ( l = 0 ; l < k ; l++ ) sum += alpha * Real2DArray_Item ( a, i, l ) * Real2DArray_Item ( b, j, l ) ;
                deviation = Maximum ( deviation, fabs ( sum ) ) ;
            }
        }
    }
    else
    {
        for ( i = 0 ; i < c->length0 ; i++ )
        {
            for ( j = 0 ; j < c->length1 ; j++ )
            {
                sum = beta * Real2DArray_Item ( store, i, j ) - Real2DArray_Item ( c, i, j ) ;
                for ( l = 0 ; l < k ; l++ ) sum += alpha * Real2DArray_Item ( a, i, l ) * Real2DArray_Item ( b, l, j ) ;
                deviation = Maximum ( deviation, fabs ( sum ) ) ;
            }
        }
    }
}
if ( deviation > 1.0e-10 ) printf ( "\nReal2DArray Matrix-Matrix Deviation (Transposes %u %u) = %.10f\n", aTranspose, bTranspose, deviation ) ;
Real2DArray_Deallocate ( &store ) ;
# endif
            }
            else Status_Set ( status, Status_InvalidDimension ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Multiply two arrays item by item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_Multiply ( Real2DArray *self, const Real2DArray *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( ( self->length0 == other->length0 ) && ( self->length1 == other->length1 ) )
        {
            auto Integer i, j ;
            for ( i = 0 ; i < self->length0 ; i++ )
            {
                for ( j = 0 ; j < self->length1 ; j++ ) Real2DArray_Item ( self , i, j ) *= Real2DArray_Item ( other, i, j ) ;
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_Print ( const Real2DArray *self )
{
    if ( self == NULL ) printf ( "Null real 2-D array.\n" ) ;
    else
    {
        auto Integer i, j ;
        for ( i = 0 ; i < self->length0 ; i++ )
        {
            for ( j = 0 ; j < self->length1 ; j++ ) printf ( "%15.10f", Real2DArray_Item ( self, i, j ) ) ;
            printf ( "\n" ) ;
        }
        printf ( "\n" ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Project a matrix from a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_ProjectOutOf1DArray ( const Real2DArray *self, Real1DArray *vector, Status *status )
{
    if ( ( self != NULL ) && ( vector != NULL ) )
    {
        auto Real1DArray *Pv ;
        Pv = Real1DArray_Allocate ( self->length1, status ) ;
        Real2DArray_VectorMultiply ( True ,  1.0e+00, self,  vector, 0.0e+00,     Pv, status ) ;
        Real2DArray_VectorMultiply ( False, -1.0e+00, self,  Pv    , 1.0e+00, vector, status ) ;
        Real1DArray_Deallocate ( &Pv ) ;

# ifdef CHECKMULTIPLICATIONS
{
/* . Check orthogonality of self. */
auto Integer i, j ;
auto Real deviation, sum ;
auto Real2DArray *r ;
r = Real2DArray_Allocate ( self->length1, self->length1, status ) ;
/* . Get U^T * U - I. */
Real2DArray_MatrixMultiply ( True, False, 1.0e+00, self, self, 0.0e+00, r, NULL ) ;
for ( i = 0 ; i < r->length0 ; i++ ) Real2DArray_Item ( r, i, i ) -= 1.0e+00 ;
deviation = Real2DArray_AbsoluteMaximum ( r ) ;
if ( deviation > 1.0e-10 ) printf ( "\nReal2DArray_ProjectOutOf1DArray> Maximum orthogonality deviation = %20.10f", deviation ) ;
Real2DArray_Deallocate ( &r ) ;
/* . Check projection. */
deviation = 0.0e+00 ;
for ( i = 0 ; i < self->length1 ; i++ )
{
    sum = 0.0e+00 ;
    for ( j = 0 ; j < self->length0 ; j++ ) sum += Real2DArray_Item ( self, j, i ) * Real1DArray_Item ( vector, j ) ;
    deviation = Maximum ( deviation, fabs ( sum ) ) ;
}
if ( deviation > 1.0e-10 ) printf ( "\nReal2DArray_ProjectOutOf1DArray> Maximum projection deviation = %20.10f", deviation ) ;

}
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pruning by row.
! . Selection is prescreened for the indices that are in range.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real2DArray *Real2DArray_PruneByRow ( const Real2DArray *self, const Selection *selection )
{
    Real2DArray *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Integer i, m, n ;
        /* . Find the number of indices in range. */
        for ( i = 0, n = 0 ; i < selection->nindices ; i++ )
        {
            m = selection->indices[i] ;
            if ( ( m >= 0 ) && ( m < self->length0 ) ) n += 1 ;
        }
        /* . Create the matrix. */
        new = Real2DArray_Allocate ( n, self->length1, NULL ) ;
	if ( ( new != NULL ) && ( new->data != NULL ) )
	{
            for ( i = 0, n = 0 ; i < selection->nindices ; i++ )
            {
        	m = selection->indices[i] ;
        	if ( ( m >= 0 ) && ( m < self->length0 ) )
        	{
                    cblas_dcopy ( self->length1, Real2DArray_RowPointer ( self, m ), self->stride1, Real2DArray_RowPointer ( new, n ), new->stride1 ) ;
                    n += 1 ;
        	}
            }
	}
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
! . This is only done for the slowest changing dimension.
! . This can only be done if the array owns the data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_Resize ( Real2DArray *self, const Integer length0, const Real *initializer, Status *status )
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
                    auto Real   *new = NULL ;
                    /* . The minimal number of items needed is, in fact, offset + length0 * stride0 - ( stride0 - 1 ).
                    ! . However, the longer, stride-consistent, value is left here for the moment just in case. */
                    n = self->offset + length0 * self->stride0 ;
                    MEMORY_REALLOCATEARRAY ( new, self->data, n, Real ) ;
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
                                for ( j = 0 ; j < self->length1 ; j++ ) Real2DArray_Item ( self, i, j ) = (*initializer) ;
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
! . Get the RMS value of the components.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real2DArray_RootMeanSquare ( const Real2DArray *self )
{
    Real alpha = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->data != NULL ) )
    {
        auto Integer i ;
        /* . Loop over rows. */
        for ( i = 0 ; i < self->length0 ; i++ ) alpha += cblas_ddot ( self->length1, Real2DArray_RowPointer ( self, i ), self->stride1,
                                                                                     Real2DArray_RowPointer ( self, i ), self->stride1 ) ;
        alpha = sqrt ( alpha / ( ( Real ) self->length ) ) ;
    }
    return alpha ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy a row from one array to another.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_RowCopyTo ( const Real2DArray  *self, const Integer m, Real2DArray *other, const Integer n, Status *status )
{
    if ( ( self != NULL )  && ( other != NULL ) )
    {
        auto Real1DArray sliceM, sliceN ;
        Real2DArray_RowSlice ( self , m, &sliceM, status ) ;
        Real2DArray_RowSlice ( other, n, &sliceN, status ) ;
        Real1DArray_CopyTo ( &sliceM, &sliceN, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a slice of a row.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_RowSlice ( const Real2DArray  *self, const Integer row, Real1DArray *slice, Status *status )
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
! . Scale all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_Scale ( Real2DArray *self, Real value )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i, j ;
        for ( i = 0 ; i < self->length0 ; i++ )
        {
            for ( j = 0 ; j < self->length1 ; j++ ) Real2DArray_Item ( self, i, j ) *= value ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_Set ( Real2DArray *self, Real value )
{
    if ( ( self != NULL ) && ( self->length > 0 ) )
    {
        auto Integer i, j ;
        for ( i = 0 ; i < self->length0 ; i++ )
        {
            for ( j = 0 ; j < self->length1 ; j++ ) Real2DArray_Item ( self, i, j ) = value ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the components of a matrix by row.
! . Selection indices out of range are ignored.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_SetByRow ( Real2DArray *self, const Selection *selection, const Real alpha )
{
    if ( ( self != NULL ) && ( selection != NULL ) && ( self->data != NULL ) )
    {
        auto Integer i, j, s ;
        for ( i = 0 ; i < selection->nindices ; i++ )
        {
            s = selection->indices[i] ;
            if ( ( s >= 0 ) && ( s < self->length0 ) )
            {
                for ( j = 0 ; j < self->length1 ; j++ ) Real2DArray_Item ( self, s, j ) = alpha ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_SetItem ( Real2DArray *self, const Integer i, const Integer j, const Real value, Status *status )
{
    if ( self != NULL )
    {
        if ( self->length > 0 )
        {
            if ( ( i >= 0 ) && ( i < self->length0 ) && ( j >= 0 ) && ( j < self->length1 ) ) Real2DArray_Item ( self, i, j ) = value ;
            else Status_Set ( status, Status_IndexOutOfRange ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a 2D slice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_Slice ( const Real2DArray *self, const Integer start0, const Integer stop0, const Integer stride0, const Integer start1, const Integer stop1, const Integer stride1, Real2DArray *slice, Status *status )
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
            slice->offset  = Real2DArray_ItemIndex ( self, start0, start1 ) ;
            slice->size    = self->size ;
            slice->stride0 = self->stride0 * stride0 ;
            slice->stride1 = self->stride1 * stride1 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Trace.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Real2DArray_Trace ( const Real2DArray *self, const Real2DArray *other, Status *status )
{
    Real trace = 0.0e+00 ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( ( self->length0 == other->length0 ) && ( self->length1 == other->length1 ) )
        {
            auto Integer i ;
            for ( i = 0 ; i < self->length0 ; i++ ) trace += cblas_ddot ( self->length1, Real2DArray_ItemPointer ( self , i, 0 ), self->stride1  ,
                                                                                         Real2DArray_ItemPointer ( other, i, 0 ), other->stride1 ) ;
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
    return trace ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transpose a matrix in-place.
! . This is always possible for square arrays and for non-square arrays that are uniform.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The function pi. */
# define IndexFunction( pIn, r, c, pOut ) { t = pIn / r ; pOut = c * ( pIn - r * t ) + t ; }

void Real2DArray_Transpose ( Real2DArray *self, Status *status )
{
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        /* . Square arrays. */
        if ( self->length0 == self->length1 )
        {
            auto Integer i, ij, j, ji ;
            auto Real    t ;
            /* . Loop over rows and columns. */
            for ( i = 0 ; i < self->length0 ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    ij = Real2DArray_ItemIndex ( self, i, j ) ;
                    ji = Real2DArray_ItemIndex ( self, j, i ) ;
                    t  = Real2DArray_ItemByIndex ( self, ij ) ;
                    Real2DArray_ItemByIndex ( self, ij ) = Real2DArray_ItemByIndex ( self, ji ) ;
                    Real2DArray_ItemByIndex ( self, ji ) = t ;
                }
            }
        }
        /* . Non-square uniform arrays. */
        else if ( Real2DArray_IsUniform ( self ) )
        {
            /*
            !
            ! . From E. G. Cate & D. W. Twigg, "Algorithm 513", ACM Transactions in Math. Software, 3, 104-110, 1977.
            !
            ! . Indexing functions (without offset and stride):
            !
            !   index[i,j] = ( i * c + j ) -> transposed[j,i] = ( j * r + i )
            !
            !   pi    ( index     [i,j] ) = transposed[j,i] => pi    ( a ) = ( a * r ) mod ( c * r )
            !   pi^-1 ( transposed[j,i] ) = index     [i,j] => pi^-1 ( a ) = ( a * c )
            !
            */
            Boolean1DArray *isMoved ;
            Integer         columns, cycles, iwork, last, maximumP0, moved, offset, p0, p1, p2, q0, q1, q2, rows, size, stride, t ;
            Real            B, C ;

            /* . Initialization. */
            columns = Real2DArray_Columns ( self ) ;
            rows    = Real2DArray_Rows    ( self ) ;
            offset  = self->offset  ;
            stride  = self->stride1 ;

            /* . Try and allocate space - not necessary but increases efficiency. */
            iwork   = ( columns + rows + 1 ) / 2 ;
            isMoved = Boolean1DArray_Allocate ( iwork, NULL ) ;
            if ( isMoved == NULL ) iwork = 0 ;
            else                   Boolean1DArray_Set ( isMoved, False ) ;

            /* . Find the number of items (fixed points) which do not need to be moved. */
            moved = 2 ;
            if ( ( rows > 2 ) && ( columns > 2 ) ) moved += ( Integer_GCD ( rows - 1, columns - 1 ) - 1 ) ;

            /* . Initial values for the search. */
            size = rows * columns ;
            last = size - 1 ;
            p0   = 1 ;

            /* . Loop until all items moved. */
            cycles = 0 ;
            while ( moved < size )
            {
                /* . Find the start of the next cycle. */
                if ( cycles > 0 )
                {
                    while ( True )
                    {
                        maximumP0 = last - p0 ;
                        p0       += 1 ;
                        IndexFunction ( p0, rows, columns, p2 ) ;
                        if ( p0 > maximumP0 ) { Status_Set ( status, Status_LogicError ) ; break ; }
                        else if ( p0 != p2 )
                        {
                            if ( p0 >= iwork )
                            {
                                while ( ( p2 > p0 ) && ( p2 < maximumP0 ) )
                                {
                                    p1 = p2 ;
                                    IndexFunction ( p1, rows, columns, p2 ) ;
                                }
                                if ( p2 == p0 ) break ;
                            }
                            else if ( ! Boolean1DArray_Item ( isMoved, p0 ) ) break ;
                        }
                    }
                }

                /* . Rearrange the items of a loop and its companion loop. */
                p1 = p0 ;
                q0 = last - p0 ;
                q1 = q0 ;
                B  = Real2DArray_ItemByIndex ( self, offset + p1 * stride ) ;
                C  = Real2DArray_ItemByIndex ( self, offset + q1 * stride ) ;
                while ( True )
                {
                    IndexFunction ( p1, rows, columns, p2 ) ;
                    q2 = last - p2 ;
                    if ( p1 < iwork ) Boolean1DArray_Item ( isMoved, p1 ) = True ;
                    if ( q1 < iwork ) Boolean1DArray_Item ( isMoved, q1 ) = True ;
                    moved += 2 ;
                    if ( p2 == p0 )
                    {
                        Real2DArray_ItemByIndex ( self, offset + p1 * stride ) = B ;
                        Real2DArray_ItemByIndex ( self, offset + q1 * stride ) = C ;
                        break ;
                    }
                    if ( p2 == q0 )
                    {
                        Real2DArray_ItemByIndex ( self, offset + p1 * stride ) = C ;
                        Real2DArray_ItemByIndex ( self, offset + q1 * stride ) = B ;
                        break ;
                    }
                    Real2DArray_ItemByIndex ( self, offset + p1 * stride ) = Real2DArray_ItemByIndex ( self, offset + p2 * stride ) ;
                    Real2DArray_ItemByIndex ( self, offset + q1 * stride ) = Real2DArray_ItemByIndex ( self, offset + q2 * stride ) ;
                    p1 = p2 ;
                    q1 = q2 ;
                }

                /* . Finish up. */
                cycles += 1 ;
                /*printf ( "Loop information ( cycle, tag, start, moved, size ): %5d %5d %5d %5d", cycles, p0, moved, size ) ; */
            }

            /* . Finish up. */
            Boolean1DArray_Deallocate ( &isMoved ) ;

            /* . Reset the view variables. */
            t = self->length0 ;
            self->length0 = self->length1 ;
            self->length1 = t ;
            self->stride0 = t * stride ;
            self->stride1 = stride ;
        }
        else Status_Set ( status, Status_NonUniformArray ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Matrix-vector multiply. A straight interface to cblas_dgemv. y = alpha A x + beta y (A can be transposed).
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_VectorMultiply ( const Boolean aTranspose, const Real alpha, const Real2DArray *a, const Real1DArray *x, const Real beta, Real1DArray *y, Status *status )
{
    if ( ( a != NULL ) && ( x != NULL ) && ( y != NULL ) )
    {
        /* . The matrix needs to be compact in dimension 1. */
        if ( ! Real2DArray_IsCompact ( a, 1, NULL ) ) Status_Set ( status, Status_NonCompactDimension ) ;
        else
        {
            auto Integer m, n ;
            m = a->length0 ;
            n = a->length1 ;
            if ( ( ( ! aTranspose ) && ( n == x->length ) && ( m == y->length ) ) ||
                 ( (   aTranspose ) && ( m == x->length ) && ( n == y->length ) ) )
            {
                auto enum CBLAS_TRANSPOSE aT ;
# ifdef CHECKMULTIPLICATIONS
auto Integer i, j ;
auto Real    deviation, sum ;
auto Real1DArray *store ;
deviation = 0.0e+00 ;
store     = Real1DArray_Allocate ( y->length, NULL ) ;
Real1DArray_CopyTo ( y, store, NULL ) ;
# endif
                if ( aTranspose ) aT = CblasTrans   ;
                else              aT = CblasNoTrans ;
                cblas_dgemv ( CblasRowMajor, aT, m, n, alpha, Real2DArray_Data ( a ), a->stride0, Real1DArray_Data ( x ), x->stride, beta, Real1DArray_Data ( y ), y->stride ) ;
# ifdef CHECKMULTIPLICATIONS
if ( aTranspose )
{
for ( i = 0 ; i < y->length ; i++ )
{
    sum = beta * Real1DArray_Item ( store, i ) - Real1DArray_Item ( y, i ) ;
    for ( j = 0 ; j < x->length ; j++ ) sum += alpha * Real2DArray_Item ( a, j, i ) * Real1DArray_Item ( x, j ) ;
    deviation = Maximum ( deviation, fabs ( sum ) ) ;
}
}
else
{
for ( i = 0 ; i < y->length ; i++ )
{
    sum = beta * Real1DArray_Item ( store, i ) - Real1DArray_Item ( y, i ) ;
    for ( j = 0 ; j < x->length ; j++ ) sum += alpha * Real2DArray_Item ( a, i, j ) * Real1DArray_Item ( x, j ) ;
    deviation = Maximum ( deviation, fabs ( sum ) ) ;
}
}
if ( deviation > 1.0e-10 ) printf ( "\nReal2DArray Matrix-Vector Deviation (Transpose %u) = %.10f\n", aTranspose, deviation ) ;
Real1DArray_Deallocate ( &store ) ;
# endif
            }
            else Status_Set ( status, Status_InvalidDimension ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create an array view of a raw array.
! . The input offset, length and stride values are final.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Real2DArray_ViewOfRaw ( Real2DArray *self, const Integer offset, const Integer length0, const Integer stride0, const Integer length1, const Integer stride1, Real *data, const Integer size, Status *status )
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
