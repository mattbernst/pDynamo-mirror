/*------------------------------------------------------------------------------
! . File      : Integer2DArrayN3.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _INTEGER2DARRAYN3
# define _INTEGER2DARRAYN3

# include "IntegerNDArray.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Integer2DArray__DecrementRowN3( self, i, xij, yij, zij ) \
        { \
	    Integer2DArray_Item ( self, i, 0 ) -= xij ; \
	    Integer2DArray_Item ( self, i, 1 ) -= yij ; \
	    Integer2DArray_Item ( self, i, 2 ) -= zij ; \
        }

# define Integer2DArray_DifferenceRowN3( self, i, j, xij, yij, zij ) \
        { \
	    xij = Integer2DArray_Item ( self, i, 0 ) - Integer2DArray_Item ( self, j, 0 ) ; \
	    yij = Integer2DArray_Item ( self, i, 1 ) - Integer2DArray_Item ( self, j, 1 ) ; \
	    zij = Integer2DArray_Item ( self, i, 2 ) - Integer2DArray_Item ( self, j, 2 ) ; \
        }

# define Integer2DArray_GetRowN3( self, i, x, y, z ) \
        { \
	   x = Integer2DArray_Item ( self, i, 0 ) ; \
	   y = Integer2DArray_Item ( self, i, 1 ) ; \
	   z = Integer2DArray_Item ( self, i, 2 ) ; \
        }

# define Integer2DArray_IncrementRowN3( self, i, xij, yij, zij ) \
        { \
	    Integer2DArray_Item ( self, i, 0 ) += xij ; \
	    Integer2DArray_Item ( self, i, 1 ) += yij ; \
	    Integer2DArray_Item ( self, i, 2 ) += zij ; \
        }

# define Integer2DArray_ScaleRowN3( self, i, value ) \
        { \
	    Integer2DArray_Item ( self, i, 0 ) *= value ; \
	    Integer2DArray_Item ( self, i, 1 ) *= value ; \
	    Integer2DArray_Item ( self, i, 2 ) *= value ; \
        }

# define Integer2DArray_SetRowN3( self, i, x, y, z ) \
        { \
	    Integer2DArray_Item ( self, i, 0 ) = x ; \
	    Integer2DArray_Item ( self, i, 1 ) = y ; \
	    Integer2DArray_Item ( self, i, 2 ) = z ; \
        }

# endif
