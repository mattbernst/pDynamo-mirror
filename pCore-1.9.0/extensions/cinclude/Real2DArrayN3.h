/*------------------------------------------------------------------------------
! . File      : Real2DArrayN3.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _REAL2DARRAYN3
# define _REAL2DARRAYN3

# include "RealNDArray.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Real2DArray_DecrementRowN3( self, i, xij, yij, zij ) \
        { \
	    Real2DArray_Item ( self, i, 0 ) -= xij ; \
	    Real2DArray_Item ( self, i, 1 ) -= yij ; \
	    Real2DArray_Item ( self, i, 2 ) -= zij ; \
        }

# define Real2DArray_DifferenceRowN3( self, i, j, xij, yij, zij ) \
        { \
	    xij = Real2DArray_Item ( self, i, 0 ) - Real2DArray_Item ( self, j, 0 ) ; \
	    yij = Real2DArray_Item ( self, i, 1 ) - Real2DArray_Item ( self, j, 1 ) ; \
	    zij = Real2DArray_Item ( self, i, 2 ) - Real2DArray_Item ( self, j, 2 ) ; \
        }

# define Real2DArray_GetRowN3( self, i, x, y, z ) \
        { \
	   x = Real2DArray_Item ( self, i, 0 ) ; \
	   y = Real2DArray_Item ( self, i, 1 ) ; \
	   z = Real2DArray_Item ( self, i, 2 ) ; \
        }

# define Real2DArray_IncrementRowN3( self, i, xij, yij, zij ) \
        { \
	    Real2DArray_Item ( self, i, 0 ) += xij ; \
	    Real2DArray_Item ( self, i, 1 ) += yij ; \
	    Real2DArray_Item ( self, i, 2 ) += zij ; \
        }

# define Real2DArray_ScaleRowN3( self, i, value ) \
        { \
	    Real2DArray_Item ( self, i, 0 ) *= value ; \
	    Real2DArray_Item ( self, i, 1 ) *= value ; \
	    Real2DArray_Item ( self, i, 2 ) *= value ; \
        }

# define Real2DArray_SetRowN3( self, i, x, y, z ) \
        { \
	    Real2DArray_Item ( self, i, 0 ) = x ; \
	    Real2DArray_Item ( self, i, 1 ) = y ; \
	    Real2DArray_Item ( self, i, 2 ) = z ; \
        }

# endif
