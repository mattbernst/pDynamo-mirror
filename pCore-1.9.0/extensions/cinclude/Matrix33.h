/*------------------------------------------------------------------------------
! . File      : Matrix33.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _MATRIX33
# define _MATRIX33

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Real2DArray.h"
# include "Status.h"
# include "Vector3.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . Type. */
# define Matrix33 Real2DArray

/* . Procedures. */
# define Matrix33_AddScaledMatrix33        Real2DArray_AddScaledArray
# define Matrix33_Clone                    Real2DArray_Clone
# define Matrix33_CopyFromArray            Real2DArray_CopyFromArrayByRow
# define Matrix33_CopyTo                   Real2DArray_CopyTo
# define Matrix33_CopyToArray              Real2DArray_CopyToArrayByRow
# define Matrix33_Deallocate               Real2DArray_Deallocate
# define Matrix33_GetItem                  Real2DArray_GetItem
# define Matrix33_GramSchmidtOrthogonalize Real2DArray_GramSchmidtOrthogonalize
# define Matrix33_Length                   Real2DArray_Length
# define Matrix33_ProjectOutOfVector       Real2DArray_ProjectOutOfVector
# define Matrix33_PruneByRow               Real2DArray_PruneByRow
# define Matrix33_RootMeanSquare           Real2DArray_RootMeanSquare
# define Matrix33_Scale                    Real2DArray_Scale
# define Matrix33_Set                      Real2DArray_Set
# define Matrix33_SetItem                  Real2DArray_SetItem
# define Matrix33_Transpose                Real2DArray_Transpose

/* . Macros. */
# define Matrix33_Data                     Real2DArray_Data        
# define Matrix33_Item                     Real2DArray_Item        
# define Matrix33_ItemIndex                Real2DArray_ItemIndex   
# define Matrix33_ItemPointer              Real2DArray_ItemPointer 
# define Matrix33_RowPointer               Real2DArray_RowPointer 

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Matrix33_GetRow( self, i, x, y, z ) \
         { \
	   x = Matrix33_Item ( self, i, 0 ) ; \
	   y = Matrix33_Item ( self, i, 1 ) ; \
	   z = Matrix33_Item ( self, i, 2 ) ; \
         }

# define Matrix33_IncrementRow( self, i, a, b, c ) \
         { \
	    Matrix33_Item ( self, i, 0 ) += a ; \
	    Matrix33_Item ( self, i, 1 ) += b ; \
	    Matrix33_Item ( self, i, 2 ) += c ; \
         }

# define Matrix33_SetRow( self, i, x, y, z ) \
         { \
	    Matrix33_Item ( self, i, 0 ) = x ; \
	    Matrix33_Item ( self, i, 1 ) = y ; \
	    Matrix33_Item ( self, i, 2 ) = z ; \
         }

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern Matrix33 *Matrix33_Allocate               ( void ) ;
extern void      Matrix33_ApplyToVector3         ( const Matrix33  *self, Vector3 *vector3 ) ;
extern void      Matrix33_AngleAxisFromRotation  ( const Matrix33  *self, const Real *tolerance, Real *angle, Vector3 *axis, Status *status ) ;
extern Real      Matrix33_Determinant            ( const Matrix33  *self ) ;
extern void      Matrix33_InverseDerivative      ( const Matrix33  *self, const Integer i, const Integer j, Matrix33 *other ) ;
extern void      Matrix33_Invert                 (       Matrix33  *self, const Matrix33 *other ) ;
extern Boolean   Matrix33_IsEqual                ( const Matrix33  *self, const Matrix33 *other ) ;
extern Boolean   Matrix33_IsIdentity             ( const Matrix33  *self ) ;
extern Boolean   Matrix33_IsImproperRotation     ( const Matrix33  *self ) ;
extern Boolean   Matrix33_IsOrthogonal           ( const Matrix33  *self ) ;
extern Boolean   Matrix33_IsProperRotation       ( const Matrix33  *self ) ;
extern void      Matrix33_PostMultiplyBy         (       Matrix33  *self, const Matrix33 *other ) ;
extern void      Matrix33_PreMultiplyBy          (       Matrix33  *self, const Matrix33 *other ) ;
extern Status    Matrix33_Reflection             (       Matrix33 **self, const Vector3 *normal ) ;
extern Status    Matrix33_RotationAboutAxis      (       Matrix33 **self, const Real angle, const Real x, const Real y, const Real z ) ;
extern Status    Matrix33_RotationFromQuaternion (       Matrix33 **self, const Real q0, const Real q1, const Real q2, const Real q3 ) ;

# endif
