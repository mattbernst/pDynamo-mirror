/*------------------------------------------------------------------------------
! . File      : RysQuadrature.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _RYSQUADRATURE
# define _RYSQUADRATURE

# include "Definitions.h"

/*------------------------------------------------------------------------------
! . Parameters.
!-----------------------------------------------------------------------------*/
/* . The maximum number of rys roots. */
# define MAXRYS 9

/*------------------------------------------------------------------------------
! . Rys Quadrature Data.
!-----------------------------------------------------------------------------*/
/* . The Rys quadrature type. */
typedef struct
{
   Real roots  [MAXRYS] ;
   Real weights[MAXRYS] ;
} RysQuadrature ;

/*------------------------------------------------------------------------------
! . Procedures.
!-----------------------------------------------------------------------------*/
extern void RysQuadrature_Roots ( RysQuadrature *roots, Integer nroots, Real x ) ;

# endif
