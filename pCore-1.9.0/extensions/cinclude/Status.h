/*------------------------------------------------------------------------------
! . File      : Status.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . The status type is used to indicate program state.
! . The messages differ in severity - some may be recoverable and others not.
!=================================================================================================================================*/
# ifndef _STATUS
# define _STATUS

# include "Definitions.h"

/* . Check and rationalize. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef enum {
    Status_Normal                    =  0 , /* . Everything OK. */
    Status_ArrayInvalidDimension     =  1 ,
    Status_ArrayInvalidExtent        =  2 ,
    Status_ArrayInvalidIndex         =  3 ,
    Status_ArrayInvalidRank          =  4 ,
    Status_ArrayInvalidSlice         =  5 ,
    Status_ArrayNonConformableShapes =  6 ,
    Status_ArrayNonConformableSizes  =  7 ,
    Status_ArrayNonReshapeable       =  8 ,
    Status_ArrayReadOnly             =  9 ,
    Status_LogicError                = 10 , /* . Too vague? */
    Status_MemoryAllocationFailure   = 11 ,
    Status_ArrayDimensionMismatch    = 12 , /* . From here not rationalized. */
    Status_ArrayLengthMismatch       = 13 ,
    Status_ArrayOutOfBounds          = 14 ,
    Status_ArrayOverFlow             = 15 ,
    Status_Continue                  = 16 ,
    Status_DiagonalizationFailure    = 17 ,
    Status_DimensionError            = 18 ,
    Status_DivideByZero              = 19 ,
    Status_DuplicateSparseArrayItems = 20 ,
    Status_IndexOutOfRange           = 21 ,
    Status_InvalidArgument           = 22 ,
    Status_InvalidArrayExtent        = 23 ,
    Status_InvalidArrayOperation     = 24 ,
    Status_InvalidArraySize          = 25 ,
    Status_InvalidDimension          = 26 ,
    Status_InvalidRank               = 27 ,
    Status_InvalidPermutation        = 28 ,
    Status_InvalidSlice              = 29 ,
    Status_InvalidSliceIndices       = 30 ,
    Status_LinearEquationFailure     = 31 ,
    Status_NegativeArrayLength       = 32 ,
    Status_NonCompactArray           = 33 ,
    Status_NonCompactDimension       = 34 ,
    Status_NonConformableArrays      = 35 ,
    Status_NonPositiveDefiniteMatrix = 36 ,
    Status_NonReshapeableArray       = 37 ,
    Status_NonUniformArray           = 38 ,
    Status_Null                      = 39 ,
    Status_NullVector                = 40 ,
    Status_OutOfMemory               = 41 ,
    Status_OverFlow                  = 42 ,
    Status_ReadOnlyArray             = 43 ,
    Status_SingularMatrix            = 44 ,
    Status_Success                   = 45 ,
    Status_ThreadError               = 46 ,
    Status_ValueError                = 47
} Status ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Check to see if everything is normal. */
# define Status_IsNormal( self ) ( (*self) == Status_IsNormal )

/* . Old. */
# define Status_IsOK( self ) ( ( self) == Status_Success  ) /* . Redefine with Status_Normal. */
# define Status_OK( self )   ( (*self) == Status_Continue )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern CSize Status_CSize   ( void ) ;
extern void  Status_SafeSet ( Status *self, Status status ) ;
extern void  Status_Set     ( Status *self, Status status ) ;

# endif
