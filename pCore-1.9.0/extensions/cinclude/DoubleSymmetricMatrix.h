/*------------------------------------------------------------------------------
! . File      : DoubleSymmetricMatrix.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _DOUBLESYMMETRICMATRIX
# define _DOUBLESYMMETRICMATRIX

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The array type. */
typedef struct {
    Boolean  isOwner  ;
    Integer  length   ; /* . Number of items. */
    Integer  length0  ; /* . Length of each dimension. */
    Integer  length01 ; /* . Length of first two and last two dimensions. */
    Integer  size     ; /* . Overall data size. */
    Real    *data     ;
} DoubleSymmetricMatrix ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Access an item (i >= j, k >= l, ij >= kl). */
/*
# define IJIndex(i,j) ( ( (i) * ( i + 1 ) ) / 2 + j )
# define DoubleSymmetricMatrix_Item(i,j,k,l) IJIndex ( IJIndex ( i, j ), IJIndex ( k, l ) )
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern DoubleSymmetricMatrix *DoubleSymmetricMatrix_Allocate      ( const Integer length, Status *status ) ;
extern DoubleSymmetricMatrix *DoubleSymmetricMatrix_Clone         ( const DoubleSymmetricMatrix  *self, Status *status ) ;
extern void                   DoubleSymmetricMatrix_CopyTo        ( const DoubleSymmetricMatrix  *self, DoubleSymmetricMatrix *other, Status *status ) ;
extern void                   DoubleSymmetricMatrix_Deallocate    (       DoubleSymmetricMatrix **self ) ;
extern Real                   DoubleSymmetricMatrix_GetItem       ( const DoubleSymmetricMatrix  *self, const Integer i, const Integer j, const Integer k, const Integer l, Status *status ) ;
extern void                   DoubleSymmetricMatrix_IncrementItem ( const DoubleSymmetricMatrix  *self, const Integer i, const Integer j, const Integer k, const Integer l, const Real value, Status *status ) ;
extern Integer                DoubleSymmetricMatrix_Index         ( const Integer i, const Integer j, const Integer k, const Integer l ) ;
extern Integer                DoubleSymmetricMatrix_Length        ( const DoubleSymmetricMatrix  *self ) ;
extern void                   DoubleSymmetricMatrix_Print         ( const DoubleSymmetricMatrix  *self ) ;
extern void                   DoubleSymmetricMatrix_Set           (       DoubleSymmetricMatrix  *self, Real value ) ;
extern void                   DoubleSymmetricMatrix_SetItem       ( const DoubleSymmetricMatrix  *self, const Integer i, const Integer j, const Integer k, const Integer l, const Real value, Status *status ) ;
extern void                   DoubleSymmetricMatrix_Unweight      (       DoubleSymmetricMatrix  *self ) ;

# endif
