/*------------------------------------------------------------------------------
! . File      : Macros.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . General macros.
!=================================================================================================================================*/
# ifndef _MACROS
# define _MACROS

# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Concatenation. */
# define _MakeToken0( a, b ) a ## b
# define _MakeToken(  a, b ) _MakeToken0 ( a, b )

/* . Even and odd. */
# define IsEven( n ) ( ( n % 2 ) == 0 )
# define IsOdd( n )  ( ( n % 2 ) != 0 )

/* . Maximum and minimum functions. */
# define Maximum( a, b ) ( ( a ) > ( b ) ? ( a ) : ( b ) )
# define Minimum( a, b ) ( ( a ) < ( b ) ? ( a ) : ( b ) )

/* . Modulo function returning positive number always. */
# define Modulo( a, b ) ( ( ( a < 0 ) ? ( ( a % b ) + b ) : a ) % b )

/* . Rounding function to return the nearest integer from a real. */
# define Round( a ) ( ( ( a ) >= 0 ) ? ( Integer ) ( ( a ) + 0.5 ) : ( Integer ) ( ( a ) - 0.5 ) )

/* . Return the sign of a number. */
# define Sign( a ) ( ( a ) >= 0.0 ? 1 : -1 )

/* . Square function. */
# define Square( a ) ( ( a ) * ( a ) )

# endif
