/*------------------------------------------------------------------------------
! . File      : SparseSymmetricMatrix.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SPARSESYMMETRICMATRIX
# define _SPARSESYMMETRICMATRIX

# include "Boolean.h"
# include "Boolean1DArray.h"
# include "Integer.h"
# include "Integer1DArray.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The item type. */
typedef struct {
    Integer i     ;
    Integer j     ;
    Integer next  ; /* . Next in column. */
    Real    value ;
} SparseSymmetricMatrixItem ;

/* . The matrix type. */
typedef struct {
    Boolean                    isCanonical            ;
    Integer                    length                 ;
    Integer                    maximumNonZeroRowItems ;
    Integer                    numberOfItems          ;
    Integer                    size                   ;
    Integer1DArray            *rowIndex               ;
    SparseSymmetricMatrixItem *items                  ;
} SparseSymmetricMatrix ;

/* . The row iterator type. */
typedef struct {
    Boolean inLT    ;
    Integer current ;
    Integer ltLast  ;
    Integer row     ;
    SparseSymmetricMatrix *target ;
} SparseSymmetricMatrixRowItemIterator ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to a particular item. */
# define SparseSymmetricMatrix_ItemPointer( self, i ) ( &((self)->items[i]) )

/* . A pointer to the start of the data. */
# define SparseSymmetricMatrix_Items( self ) ( &((self)->items[0]) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern SparseSymmetricMatrix *SparseSymmetricMatrix_Allocate                               ( const Integer length, const Integer size, Status *status ) ;
extern void                   SparseSymmetricMatrix_AppendItem                             (       SparseSymmetricMatrix  *self, const Integer i, const Integer j, const Real value, Status *status ) ;
extern void                   SparseSymmetricMatrix_ApplyIncompleteCholeskyDecomposition   ( const SparseSymmetricMatrix  *self, const Real1DArray *b, Real1DArray *x ) ;
extern void                   SparseSymmetricMatrix_Canonicalize                           (       SparseSymmetricMatrix  *self, Status *status ) ;
extern void                   SparseSymmetricMatrix_Clear                                  (       SparseSymmetricMatrix  *self ) ;
extern SparseSymmetricMatrix *SparseSymmetricMatrix_Clone                                  ( const SparseSymmetricMatrix  *self, Status *status ) ;
extern void                   SparseSymmetricMatrix_ComputeIncompleteCholeskyDecomposition (       SparseSymmetricMatrix  *self, const Real alpha, Integer *numberOfModifiedPivots, Status *status ) ;
extern void                   SparseSymmetricMatrix_CopyTo                                 ( const SparseSymmetricMatrix  *self, SparseSymmetricMatrix *other, Status *status ) ;
extern void                   SparseSymmetricMatrix_Deallocate                             (       SparseSymmetricMatrix **self ) ;
extern void                   SparseSymmetricMatrix_GetDiagonal                            ( const SparseSymmetricMatrix  *self, Real1DArray *diagonal, Status *status ) ;
extern void                   SparseSymmetricMatrix_MakeDiagonalPreconditioner             ( const SparseSymmetricMatrix  *self, Real1DArray *preconditioner, const Real *tolerance, Status *status ) ;
extern void                   SparseSymmetricMatrix_Print                                  ( const SparseSymmetricMatrix  *self ) ;
extern void                   SparseSymmetricMatrix_VectorMultiply                         ( const SparseSymmetricMatrix  *self, const Real1DArray *x, Real1DArray *y, Status *status ) ;

/* . Iterators. */
extern void    SparseSymmetricMatrixRowItemIterator_Initialize ( SparseSymmetricMatrixRowItemIterator *self, SparseSymmetricMatrix *target, const Integer row, Status *status ) ;
extern Integer SparseSymmetricMatrixRowItemIterator_Next       ( SparseSymmetricMatrixRowItemIterator *self ) ;

/*
! . Options for canonicalize - keep, remove, etc.?
! . For remove - set small values to zero.
extern Real                   SparseSymmetricMatrix_AbsoluteMaximum  ( const SparseSymmetricMatrix  *self ) ;
extern Real                   SparseSymmetricMatrix_GetItem          (       SparseSymmetricMatrix  *self, const Integer i, const Integer j, Status *status ) ;
extern void                   SparseSymmetricMatrix_RemoveSmallItems (       SparseSymmetricMatrix  *self, const Real tolerance ) ;
extern void                   SparseSymmetricMatrix_Resize           (       SparseSymmetricMatrix  *self, Status *status ) ;
extern void                   SparseSymmetricMatrix_Scale            (	    SparseSymmetricMatrix  *self, const Real value ) ;
extern void                   SparseSymmetricMatrix_Set              (	    SparseSymmetricMatrix  *self, const Real value ) ;
extern void                   SparseSymmetricMatrix_SetItem          (       SparseSymmetricMatrix  *self, const Integer i, const Integer j, const Real value, Status *status ) ;
extern Real                   SparseSymmetricMatrix_Sparsity         (       SparseSymmetricMatrix  *self, const Real tolerance ) ;
*/

# endif
