/*------------------------------------------------------------------------------
! . File      : Definitions.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This file contains general definitions used throughout the cDynamo program.
!=================================================================================================================================*/
# ifndef _DEFINITIONS
# define _DEFINITIONS

# include <limits.h>
# include <stddef.h>
# include <stdint.h>

/* . The standard cardinal and integer need to be at least 32 bits. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Primitive value types.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Default types. */
# include "Boolean.h"
# include "Cardinal.h"
# include "Integer.h"
# include "Real.h"

/* . Specific types. */
typedef uint16_t Cardinal16 ;
typedef uint32_t Cardinal32 ;
typedef uint64_t Cardinal64 ;
typedef  int16_t Integer16  ;
typedef  int32_t Integer32  ;
typedef  int64_t Integer64  ;
typedef float    Real32     ;
typedef double   Real64     ;

typedef char     Character  ;
typedef size_t   CSize      ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
# include "Macros.h"

# endif
