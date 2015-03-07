/*------------------------------------------------------------------------------
! . File      : Vector3.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _VECTOR3
# define _VECTOR3

# include "Boolean.h"
# include "Definitions.h"
# include "Real1DArray.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Type. */
# define Vector3 Real1DArray

/* . Procedures. */
# define Vector3_AbsoluteMaximum      Real1DArray_AbsoluteMaximum
# define Vector3_AbsoluteMaximumIndex Real1DArray_AbsoluteMaximumIndex
# define Vector3_Add                  Real1DArray_Add
# define Vector3_AddScalar            Real1DArray_AddScalar
# define Vector3_AddScaledVector      Real1DArray_AddScaledArray
# define Vector3_Clone                Real1DArray_Clone
# define Vector3_CopyTo               Real1DArray_CopyTo
# define Vector3_Deallocate           Real1DArray_Deallocate
# define Vector3_Divide               Real1DArray_Divide
# define Vector3_Dot                  Real1DArray_Dot
# define Vector3_GetItem              Real1DArray_GetItem
# define Vector3_Length               Real1DArray_Length
# define Vector3_Maximum              Real1DArray_Maximum
# define Vector3_Minimum              Real1DArray_Minimum
# define Vector3_Multiply             Real1DArray_Multiply
# define Vector3_Norm2                Real1DArray_Norm2
# define Vector3_Normalize            Real1DArray_Normalize
# define Vector3_RootMeanSquare       Real1DArray_RootMeanSquare
# define Vector3_Scale                Real1DArray_Scale
# define Vector3_Set                  Real1DArray_Set
# define Vector3_SetItem              Real1DArray_SetItem
# define Vector3_Sum                  Real1DArray_Sum
# define Vector3_ViewOfRaw            Real1DArray_ViewOfRaw

/* . Macros. */
# define Vector3_Data                 Real1DArray_Data
# define Vector3_Item                 Real1DArray_Item
# define Vector3_ItemPointer          Real1DArray_ItemPointer

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Vector3 *Vector3_Allocate     ( void ) ;
extern void     Vector3_CrossProduct (       Vector3 *self, const Vector3 *other ) ;
extern Boolean  Vector3_IsEqual      ( const Vector3 *self, const Vector3 *other ) ;
extern Boolean  Vector3_IsNull       ( const Vector3 *self ) ;

# endif
