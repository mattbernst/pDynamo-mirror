/*------------------------------------------------------------------------------
! . File      : Real.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . The default real type (at least 64 bit floating point).
!=================================================================================================================================*/
# ifndef _REAL
# define _REAL

# include <float.h>

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Type. */
typedef double Real ;

/* . Maximum and minimum reals. */
# define Real_Maximum  DBL_MAX
# define Real_Minimum -DBL_MAX

/* . These need redoing. */
/* . The largest real - used only for setting invalid calculations. */
# define Real_Huge 1.0e+300

/* . The maximum exponent for exp (). */
# define Real_MaximumExponent 700.0e+00

/* . Some approximation to the smallest normal number. */
# define Real_SmallestNormalNumber 2.5e-308

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real Real_SafeMinimum  ( void ) ;
extern Real Real_UnitRoundOff ( void ) ;

# endif
