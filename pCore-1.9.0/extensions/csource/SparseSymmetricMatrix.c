/*------------------------------------------------------------------------------
! . File      : SparseSymmetricMatrix.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Real sparse symmetric matrices.
!=================================================================================================================================*/

/*
! . There are many alternative ways of storing items in the matrix. Currently all diagonal items
! . are stored followed by all off-diagonal items. This is done because:
!
! * The separate storage of diagonal and off-diagonal items facilitates matrix-vector multiplication.
! * The storage of all diagonal items (as opposed to only non-zero ones) is required for some
!   operations (e.g. Cholesky decomposition) and also facilitates some others (e.g. row iteration).
!
! . Equally, at the moment, both diagonal and off-diagonal items employ the same item type although
! . this involves redundant extra storage for diagonal items.
!
! . A matrix is canonical when its off-diagonal items are put in ascending order of (i,j) pairs
! . with i > j.
!
! . These choices may change if matrices with many zero diagonal elements are to be routinely treated.
*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "Memory.h"
# include "SparseSymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer Value_Compare ( const void *vItem1, const void *vItem2 ) ;

static void SparseSymmetricMatrix_IndexItems              ( SparseSymmetricMatrix *self ) ;
static void SparseSymmetricMatrix_InitializeDiagonalItems ( SparseSymmetricMatrix *self ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Minimum size allocation is length. */
SparseSymmetricMatrix *SparseSymmetricMatrix_Allocate ( const Integer length, const Integer size, Status *status )
{
    SparseSymmetricMatrix *self = NULL ;
    if ( length < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
    else
    {
        MEMORY_ALLOCATE ( self, SparseSymmetricMatrix ) ;
        if ( self != NULL )
        {
            /* . Scalars. */
            self->isCanonical            = False  ;
            self->length                 = length ;
            self->maximumNonZeroRowItems = 0      ;
            self->numberOfItems          = length ;
            self->size                   = Maximum ( size, length ) ;
            /* . Arrays. */
            self->rowIndex               = NULL ;
            self->items                  = NULL ;
            /* . Allocation. */
            MEMORY_ALLOCATEARRAY ( self->items, self->size, SparseSymmetricMatrixItem ) ;
            self->rowIndex = Integer1DArray_Allocate ( length+1, status ) ;
            if ( ( self->items == NULL ) || ( self->rowIndex == NULL ) ) SparseSymmetricMatrix_Deallocate ( &self ) ;
            else
            {
                Integer1DArray_Set ( self->rowIndex, -1 ) ;
                SparseSymmetricMatrix_InitializeDiagonalItems ( self ) ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    return self ;
}


/*----------------------------------------------------------------------------------------------------------------------------------
! . Append an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Use only for off-diagonal items? */

void SparseSymmetricMatrix_AppendItem ( SparseSymmetricMatrix *self, const Integer i, const Integer j, const Real value, Status *status )
{
    if ( self != NULL )
    {
	auto Integer n ;
        n = self->numberOfItems ;
        /* . General checking. */
        if ( ( i < 0 ) || ( i >= self->length ) || ( j < 0 ) || ( j >= self->length ) ) Status_Set ( status, Status_IndexOutOfRange ) ;
        /* . Diagonal items. */
	else if ( i == j )
	{
	    if ( self->items[i].value == 0.0e+00 ) self->items[i].value = value ;
	    else Status_Set ( status, Status_DuplicateSparseArrayItems ) ;
	}
        /* . Overflow. */
        else if ( n >= self->size ) Status_Set ( status, Status_ArrayOverFlow ) ;
	/* . Off-diagonal items. */
        else
        {
	    self->isCanonical    = False ;
            self->items[n].i     =     i ;
            self->items[n].j     =     j ;
            self->items[n].next  =    -1 ;
            self->items[n].value = value ;
            self->numberOfItems ++ ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply an incomplete Cholesky decomposition to a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . The equation M * x = b is solved where M is an incomplete Cholesky decomposition
! . computed by ComputeIncompleteCholeskyDecomposition.
*/

void SparseSymmetricMatrix_ApplyIncompleteCholeskyDecomposition ( const SparseSymmetricMatrix *self, const Real1DArray *b, Real1DArray *x )
{
    if ( self != NULL )
    {
        auto Integer c, i, j, n ;
	auto Real    sum ;

        /* . Length. */
        n = self->length ;

        /* . Initialization. */
        Real1DArray_CopyTo ( b, x, NULL ) ;

        /* . Compute x = P^-1 * b = P^T * b. */
/*
        if ( self->permutation == NULL ) Real1DArray_CopyTo ( b, x, NULL ) ;
        else
        {
            for ( i = 0 ; i < n ; i++ ) Real1DArray_Item ( x, i ) = Real1DArray_Item ( b, Integer1DArray_Item ( self->permutation, i ) ) ;
        }
*/

        /* . Compute x = L^-1 * x. */
        for ( i = 1 ; i < n ; i++ )
        {
	    sum = 0.0e+00 ;
            for ( j = Integer1DArray_Item ( self->rowIndex, i ) ; j < Integer1DArray_Item ( self->rowIndex, i+1 ) ; j++ )
            {
                c = self->items[j].j ;
                sum += ( self->items[j].value * Real1DArray_Item ( x, c ) ) ;
            }
	    Real1DArray_Item ( x, i ) -= sum ;
        }

        /* . x = D^-1 * x. */
        for ( i = 0 ; i < n ; i++ ) Real1DArray_Item ( x, i ) *= self->items[i].value ;

        /* . Compute x = (L^T)^-1 * x. */
        for ( i = n-1 ; i > 0 ; i-- )
        {
            for ( j = Integer1DArray_Item ( self->rowIndex, i ) ; j < Integer1DArray_Item ( self->rowIndex, i+1 ) ; j++ )
            {
                c = self->items[j].j ;
                Real1DArray_Item ( x, c ) -= ( self->items[j].value * Real1DArray_Item ( x, i ) ) ;
            }
        }

        /* . Compute x = (P^T)^-1 * b = P * x. */
/*
        Real1DArray_Permute ( x, self->permutation, NULL ) ;
*/
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Canonicalize the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_Canonicalize ( SparseSymmetricMatrix *self, Status *status )
{
    if ( ( self != NULL ) && ( ! self->isCanonical ) )
    {
        /* . There are off-diagonal items. */
	if ( self->numberOfItems > self->length )
	{
            auto Integer i, j, l, n ;

            /* . Ensure i > j for each off-diagonal item. */
            for ( l = self->length ; l < self->numberOfItems ; l++ )
            {
        	i = self->items[l].i ;
        	j = self->items[l].j ;
        	if ( i < j ) { self->items[l].i = j ; self->items[l].j = i ; }
            }

            /* . Put all off-diagonal items in order. */
	    n = self->numberOfItems - self->length ;
            qsort ( ( void * ) SparseSymmetricMatrix_ItemPointer ( self, self->length ), ( size_t ) n, sizeof ( SparseSymmetricMatrixItem ), ( void * ) Value_Compare ) ;

            /* . Remove items with duplicate indices. */
            for ( l = n = self->length+1 ; l < self->numberOfItems ; l++ )
            {
        	i = self->items[l].i ;
        	j = self->items[l].j ;
        	if ( ( i != self->items[n-1].i ) || ( j != self->items[n-1].j ) )
        	{
                    if ( n != l )
                    {
                	self->items[n].i     = self->items[l].i ;
                	self->items[n].j     = self->items[l].j ;
                	self->items[n].value = self->items[l].value ;
                    }
                    n++ ;
        	}
            }
            if ( n < self->numberOfItems ) Status_Set ( status, Status_DuplicateSparseArrayItems ) ;
            self->numberOfItems = n ;
        }

        /* . Index the items. */
	SparseSymmetricMatrix_IndexItems ( self ) ;

        /* . Finish up. */
        self->isCanonical = True ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Clear all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_Clear ( SparseSymmetricMatrix *self )
{
    if ( self != NULL )
    {
        self->isCanonical            = False  ;
        self->maximumNonZeroRowItems =     0  ;
        self->numberOfItems          = self->length ;
        SparseSymmetricMatrix_InitializeDiagonalItems ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
! . A minimum size array is allocated.
!---------------------------------------------------------------------------------------------------------------------------------*/
SparseSymmetricMatrix *SparseSymmetricMatrix_Clone ( const SparseSymmetricMatrix *self, Status *status )
{
    SparseSymmetricMatrix *clone = NULL ;
    if ( self != NULL )
    {
        clone = SparseSymmetricMatrix_Allocate ( self->length, self->numberOfItems, status ) ;
        SparseSymmetricMatrix_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Compute a modified or unmodified incomplete Cholesky decomposition in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . The decomposition that is computed is A = M + R where M = P*L*D*L^T*P^T.
! . Here L is lower triangular with unit diagonal, D is diagonal and P a permutation.
! . It is stored as (L-I) + D. D can also be inverted D^-1.
! . A non-zero alpha permits the modified ICD to be computed (0 <= alpha <= 1).
!
! . At the moment no fill-in is attempted so P is the identity.
!
! . For fill-in, Markowitz strategy chooses current row which has least number of non-zero elements.
! . When there are ties - choose row of lowest index.
*/

# define _SmallFactor 1.0e-02

void SparseSymmetricMatrix_ComputeIncompleteCholeskyDecomposition ( SparseSymmetricMatrix *self, const Real alpha, Integer *numberOfModifiedPivots, Status *status )
{
    if ( self != NULL )
    {
        auto Boolean                              doFillIn ;
        auto Boolean1DArray                      *visitedIndices ;
	auto Integer                              i, j, jActive, k, kActive, kIndex, l, lActive, lower, n, next, numberActive, pivotRow ;
        auto Integer1DArray                      *columnIndices, *itemIndices ;
        auto Real                                 f, fillIn, pivot, reciprocalPivot, sum ;
        auto SparseSymmetricMatrixItem           *item ;
	auto SparseSymmetricMatrixRowItemIterator iterator ;

        /* . Initialization. */
	doFillIn = ( alpha != 0.0e+00 ) ;
        n        = self->length ;
	(*numberOfModifiedPivots) = 0 ;
        SparseSymmetricMatrix_Canonicalize ( self, status ) ;

        /* . Allocate space. */
        columnIndices  = Integer1DArray_Allocate ( self->maximumNonZeroRowItems, status ) ;
        itemIndices    = Integer1DArray_Allocate ( self->maximumNonZeroRowItems, status ) ;
	visitedIndices = Boolean1DArray_Allocate ( self->maximumNonZeroRowItems, status ) ;
/*
        iPivot        = Integer1DArray_Allocate ( n, status ) ;
        Integer1DArray_Set ( iPivot, -1 ) ;
*/

        /* . Perform the decomposition. */
        if ( Status_OK ( status ) )
        {
	    /* . Loop over all rows. */
            for ( i = 0 ; i < n ; i++ )
            {

                /* . Select pivot row. */
	        /* . The Markowitz strategy selects the active row with the smallest
	        !  . number of non-zero elements. In the event of a tie, the row of
	        !  . lowest index wins.
	        */

                /* . Use the next row. */
                pivotRow = i ;

                /* . Save iPivot - the row is now inactive. */
/*
	        Integer1DArray_Item ( iPivot, pivotRow ) = i ;
*/

                /* . Find the indices of the non-zero active elements in the row. */
	        /* . The column indices will always be in ascending order. */
	        numberActive = 0 ;
	        sum          = 0.0e+00 ;
                SparseSymmetricMatrixRowItemIterator_Initialize ( &iterator, self, pivotRow, status ) ;
	        while ( ( next = SparseSymmetricMatrixRowItemIterator_Next ( &iterator ) ) >= 0 )
	        {
	            item = &(self->items[next]) ;
		    if ( item->i == pivotRow ) j = item->j ;
		    else                       j = item->i ;
/*
		    if ( Integer1DArray_Item ( iPivot, j ) < 0 )
*/
		    if ( j > i )
		    {
		        Integer1DArray_Item ( columnIndices , numberActive ) = j     ;
		        Integer1DArray_Item ( itemIndices   , numberActive ) = next  ;
		        Boolean1DArray_Item ( visitedIndices, numberActive ) = False ;
   	                numberActive += 1 ;
		        sum          += fabs ( item->value ) ;
		    }
	        }

                /* . Find the reciprocal of the diagonal element. */
	        /* . A check is made to ensure positive-definiteness. */
                pivot = self->items[pivotRow].value ;
                if ( pivot <= _SmallFactor * sum )
	        {
                    self->items[pivotRow].value = sum ;
                    if ( sum == 0.0e+00 ) self->items[pivotRow].value = 1.0e+00 ;
                    (*numberOfModifiedPivots) += 1 ;
	        }
                reciprocalPivot = 1.0e+00 / self->items[pivotRow].value ;
/*
printf ( "\nRow Number Active: %d %d %25.15f %25.15f\n", i, numberActive, sum, reciprocalPivot ) ;
Integer1DArray_Print ( columnIndices ) ;
Integer1DArray_Print ( itemIndices   ) ;
*/
                /* . Loop over non-zero active rows. */
                for ( jActive = 0 ; jActive < numberActive ; jActive++ )
	        {

                    /* . Initialization - Aij/Aii. */
                    j = Integer1DArray_Item ( columnIndices, jActive ) ;
                    f = reciprocalPivot * self->items[Integer1DArray_Item ( itemIndices, jActive )].value ;
/*
printf ( "\nJ Active>: %d %d %25.15f\n", jActive, j, f ) ;
*/
                    /* . Fill and Markowitz manipulations here. */

                    /* . Loop over non-zero active items of row j - including the diagonal. */
		    /* . k always increases. */
                    lower = 0 ;
                    SparseSymmetricMatrixRowItemIterator_Initialize ( &iterator, self, j, status ) ;
	            while ( ( next = SparseSymmetricMatrixRowItemIterator_Next ( &iterator ) ) >= 0 )
		    {
	                item = &(self->items[next]) ;
		        if ( item->i == j ) k = item->j ;
		        else                k = item->i ;
/*
printf ( "\nk = %d\n", k ) ;
*/
		        /* . Search to see if k also occurs in the active rows of i. */
		        kActive = -1 ;
		        for ( lActive = lower ; lActive < numberActive ; lActive++ )
		        {
			    l = Integer1DArray_Item ( columnIndices, lActive ) ;
		                 if ( l >  k ) break ;
			    else if ( l == k ) { kActive = lActive ; lower = lActive + 1 ; break ; }
			    else lower = lActive ; /* . l < k. */
		        }
		        /*
		        ! . If the column is active and also occurs in the active columns of i,
		        ! . modify Ajk by - Aji*Aki/Aii and flag the item as having been used.
		        */
/*
		        if ( ( Integer1DArray_Item ( iPivot, k ) < 0 ) && ( kActive >= 0 ) )
*/
		        if ( ( k > i ) && ( kActive >= 0 ) )
		        {
		            kIndex      = Integer1DArray_Item ( itemIndices, kActive ) ;
/*printf ( "\nMatch k = %d %d %d %25.15f\n", k, kActive, kIndex, - ( f * self->items[kIndex].value ) ) ;*/
			    item->value -= ( f * self->items[kIndex].value ) ;
                            if ( doFillIn ) Boolean1DArray_Item ( visitedIndices, kActive ) = True ;
		        }
		    }

                    /* . Check for fill-in or a modified ICD. */
                    if ( doFillIn )
		    {
                        /* . Reloop over active items. */
	                for ( kActive = 0 ; kActive < numberActive ; kActive++ )
		        {
                	    /* . Reactivate already visited items. */
			    if ( Boolean1DArray_Item ( visitedIndices, kActive ) )
			    {
			        Boolean1DArray_Item ( visitedIndices, kActive ) = False ;
			    }
			    /* . Item not visited. */
			    /*
			    ! . This means that the pivot row has a non-zero entry but row j does not.
			    ! . This will create fill-in or a correction to the corresponding diagonal
			    ! . items.
			    */
        		    else
			    {
		                /* . Calculate fill-in if appropriate. */
                	        k = Integer1DArray_Item ( columnIndices, kActive ) ;
                                if ( j >= k )
			        {
                        	    fillIn = -f * self->items[Integer1DArray_Item ( itemIndices, kActive )].value ; /* . - Aji*Aki/Aii. */

                        	    /* . Fill-in manipulations - add into Ajk. */

                        	    /* . Add in the correction to the diagonals (modified ICD). */
                        	    fillIn *= alpha ;
                        	    self->items[j].value += fillIn ;
                        	    self->items[k].value += fillIn ;
                                }
        		    }
                        }
		    }
	        }
            }

            /* . Replace diagonal elements by their inverses. */
	    for ( i = 0 ; i < n ; i++ )
	    {
	        self->items[i].value = 1.0e+00 / self->items[i].value ;
	    }

            /* . Scale all lower triangular elements in a column by the column diagonal. */
	    for ( i = n ; i < self->numberOfItems ; i++ )
	    {
	        item = &(self->items[i]) ;
	        item->value *= self->items[item->j].value ;
	    }

            /* . Invert iPivot. */
/*
	    Integer1DArray_CopyTo ( iPivot, iWork, NULL ) ;
	    for ( i = 0 ; i < n ; i++ ) Integer1DArray_Item ( iPivot, Integer1DArray_Item ( iWork, i ) ) = i ;
*/
        }

        /* . Finish up - recanonicalization unnecesary in principle. */
        Integer1DArray_Deallocate ( &columnIndices  ) ;
        Integer1DArray_Deallocate ( &itemIndices    ) ;
        Boolean1DArray_Deallocate ( &visitedIndices ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
! . As much off-diagonal data as possible is copied.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_CopyTo ( const SparseSymmetricMatrix *self, SparseSymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        /* . Matrix dimension. */
        if ( self->length == other->length )
        {
            auto Integer l, n ;
            if ( self->numberOfItems <= other->size ) n = self->numberOfItems ;
            else
            {
                n = other->size ;
                Status_Set ( status, Status_ArrayOverFlow ) ;
            }
            for ( l = 0 ; l < n ; l++ )
            {
                other->items[l].i     = self->items[l].i     ;
                other->items[l].j     = self->items[l].j     ;
                other->items[l].value = self->items[l].value ;
            }
            other->numberOfItems = n ;
            if ( self->isCanonical ) SparseSymmetricMatrix_Canonicalize ( other, status ) ;
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_Deallocate ( SparseSymmetricMatrix **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Integer1DArray_Deallocate ( &((*self)->rowIndex) ) ;
        MEMORY_DEALLOCATE         (   (*self)->items     ) ;
        MEMORY_DEALLOCATE         (   (*self)            ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Retrieve diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_GetDiagonal ( const SparseSymmetricMatrix *self, Real1DArray *diagonal, Status *status )
{
    Real1DArray_Set ( diagonal, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( diagonal != NULL ) )
    {
        if ( self->length == diagonal->length )
        {
            auto Integer l ;
            for ( l = 0 ; l < self->length ; l++ ) Real1DArray_Item ( diagonal, l ) = self->items[l].value ;
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Index the items.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void SparseSymmetricMatrix_IndexItems ( SparseSymmetricMatrix *self )
{
    if ( self != NULL )
    {
        auto Integer l, last, n ;

        /* . Find the maximum number of non-zero items in a row. */
        /* . Initialization - use rowIndex as workspace. */
        Integer1DArray_Set ( self->rowIndex, 1 ) ;
        for ( l = self->length ; l < self->numberOfItems ; l++ )
        {
            Integer1DArray_Item ( self->rowIndex, self->items[l].i ) += 1 ;
            Integer1DArray_Item ( self->rowIndex, self->items[l].j ) += 1 ;
        }
        self->maximumNonZeroRowItems = Integer1DArray_Maximum ( self->rowIndex ) ;

        /* . Set up the row index - do in reverse order. */
        /* . Initialization. */
        Integer1DArray_Set ( self->rowIndex, -1 ) ;

        /* . Upper-triangular items. */
        for ( l = self->numberOfItems - 1 ; l >= self->length ; l-- )
        {
            last = Integer1DArray_Item ( self->rowIndex, self->items[l].j ) ;
            self->items[l].next = last ;
            Integer1DArray_Item ( self->rowIndex, self->items[l].j ) = l ;
        }

        /* . Diagonal items. */
        for ( l = 0 ; l < self->length ; l++ )
        {
            last = Integer1DArray_Item ( self->rowIndex, self->items[l].i ) ;
            self->items[l].next = last ;
            Integer1DArray_Item ( self->rowIndex, self->items[l].i ) = l ;
        }

        /* . Row index. */
        Integer1DArray_Set ( self->rowIndex, 0 ) ;
        for ( l = self->length ; l < self->numberOfItems ; l++ )
        {
            Integer1DArray_Item ( self->rowIndex, self->items[l].i ) += 1 ;
        }
	last = self->length ;
        for ( l = 0 ; l < self->length ; l++ )
        {
	    n = Integer1DArray_Item ( self->rowIndex, l ) ;
            Integer1DArray_Item ( self->rowIndex, l ) = last ;
	    last += n ;
        }
        Integer1DArray_Item ( self->rowIndex, self->length ) = last ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialize the diagonal items.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void SparseSymmetricMatrix_InitializeDiagonalItems ( SparseSymmetricMatrix *self )
{
    if ( ( self != NULL ) && ( self->items != NULL ) )
    {
        auto Integer l ;
        for ( l = 0 ; l < self->length ; l++ )
        {
            self->items[l].i     =  l ;
            self->items[l].j     =  l ;
            self->items[l].next  = -1 ;
            self->items[l].value = 0.0e+00 ;
        }
    }
}


/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a diagonal preconditioner from the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Small 1.0e-06
void SparseSymmetricMatrix_MakeDiagonalPreconditioner ( const SparseSymmetricMatrix *self, Real1DArray *preconditioner, const Real *tolerance, Status *status )
{
    if ( ( self != NULL ) && ( preconditioner != NULL ) )
    {
        auto Integer i ;
        auto Real    t, v ;
        if ( tolerance == NULL ) t = _Small ;
        else                     t = (*tolerance ) ;
        SparseSymmetricMatrix_GetDiagonal ( self, preconditioner, status ) ;
        for ( i = 0 ; i < preconditioner->length ; i++ )
        {
            v = Maximum ( fabs ( Real1DArray_Item ( preconditioner, i ) ), t ) ;
            Real1DArray_Item ( preconditioner, i ) = 1.0e+00 / v ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_Print ( const SparseSymmetricMatrix  *self )
{
    if ( self == NULL ) printf ( "\nNull sparse symmetric matrix.\n" ) ;
    else
    {
        auto Integer l, start ;

        /* . Items. */
        printf ( "\nItems (index, i, j, next-in-column, value):\n" ) ;
        for ( l = 0 ; l < self->numberOfItems ; l++ ) printf ( "%10d %10d %10d %10d %15.10f\n", l, self->items[l].i, self->items[l].j, self->items[l].next, self->items[l].value ) ;

        /* . Row index. */
        if ( self->isCanonical )
        {
            printf ( "\nRow Index (row, start, number of items):\n" ) ;
            for ( l = 0 ; l < self->length ; l++ )
            {
                start = Integer1DArray_Item ( self->rowIndex, l ) ;
                printf ( "%10d %10d %10d\n", l, start, Integer1DArray_Item ( self->rowIndex, l+1 ) - start ) ;
            }
        }

        /* . Other data. */
        printf ( "\nOther Data (length, size, numberOfItems, maximumNonZeroRowItems): %10d %10d %10d %10d\n", self->length, self->size, self->numberOfItems, self->maximumNonZeroRowItems ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Matrix-vector multiplication.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_VectorMultiply ( const SparseSymmetricMatrix *self, const Real1DArray *x, Real1DArray *y, Status *status )
{
    Real1DArray_Set ( y, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( x != NULL ) && ( y != NULL ) )
    {
        if ( ( self->length == x->length ) && ( self->length == y->length ) )
        {
            auto Integer i, j, l ;
            auto Real    value ;
            for ( l = 0 ; l < self->length ; l++ ) Real1DArray_Item ( y, l ) += self->items[l].value * Real1DArray_Item ( x, l ) ;
            for ( l = self->length ; l < self->numberOfItems ; l++ )
            {
                i     = self->items[l].i ;
                j     = self->items[l].j ;
                value = self->items[l].value ;
                Real1DArray_Item ( y, i ) += value * Real1DArray_Item ( x, j ) ;
                Real1DArray_Item ( y, j ) += value * Real1DArray_Item ( x, i ) ;
            }
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*==================================================================================================================================
! . Row iterator procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrixRowItemIterator_Initialize ( SparseSymmetricMatrixRowItemIterator *self, SparseSymmetricMatrix *target, const Integer row, Status *status )
{
    if ( ( self != NULL ) && ( target != NULL ) && ( target->isCanonical ) )
    {
        if ( ( row >= 0 ) && ( row < target->length ) )
        {
            auto Integer n ;
            /* . Initialization. */
            self->inLT    = True   ;
	    self->current = -1     ;
            self->ltLast  = Integer1DArray_Item ( target->rowIndex, row + 1 ) ;
            self->row     = row    ;
            self->target  = target ;
	    /* . Find the first item in row, either in lower triangle or the diagonal. */
	    n = self->ltLast - Integer1DArray_Item ( target->rowIndex, row ) ;
	    if ( n > 0 ) self->current = Integer1DArray_Item ( target->rowIndex, row ) ;
	    else       { self->current = row ; self->inLT = False ; }
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Next item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer SparseSymmetricMatrixRowItemIterator_Next ( SparseSymmetricMatrixRowItemIterator *self )
{
    Integer current = -1 ;
    if ( ( self != NULL ) && ( self->current != -1 ) )
    {
        /* . Set current. */
        current = self->current ;
	/* . Find next item. */
        if ( self->inLT )
        {
            /* . Move to diagonal. */
	    if ( current == ( self->ltLast - 1 ) ) { self->current = self->row ; self->inLT = False ; }
            /* . Continue in lower triangle. */
	    else self->current += 1 ;
        }
        /* . Upper triangle. */
        else self->current = self->target->items[current].next ;
/*
if ( ( current >= self->target->numberOfItems ) || ( self->current >= self->target->numberOfItems ) )
{
printf ( "\nInformation: %d %d %d %d %d %d %d\n", current, self->current, self->target->length, self->inLT, self->ltLast, self->row, self->target->numberOfItems ) ;
}
*/
    }
    return current ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/* . This is applied to off-diagonal items only. */
/* . For this to work, i >= j. */
static Integer Value_Compare ( const void *vItem1, const void *vItem2 )
{
    Integer i ;
    SparseSymmetricMatrixItem *item1, *item2 ;
    item1 = ( SparseSymmetricMatrixItem * ) vItem1 ;
    item2 = ( SparseSymmetricMatrixItem * ) vItem2 ;
         if ( item1->i < item2->i ) i = -1 ;
    else if ( item1->i > item2->i ) i =  1 ;
    else if ( item1->j < item2->j ) i = -1 ;
    else if ( item1->j > item2->j ) i =  1 ;
    else i = 0 ;
    return i ;
}
