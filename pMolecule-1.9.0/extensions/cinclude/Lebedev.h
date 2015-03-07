/*------------------------------------------------------------------------------
! . File      : Lebedev.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _LEBEDEVLAIKOV
# define _LEBEDEVLAIKOV

# include "Definitions.h"

/*------------------------------------------------------------------------------
! . Procedures.
!-----------------------------------------------------------------------------*/
extern Integer LebedevLaikov_Angular_Momentum_Value ( const Integer npts   ) ;
extern Integer LebedevLaikov_Number_Of_Points       ( const Integer lvalue ) ;
extern Integer LebedevLaikov_Points                 ( const Integer N, Real *X, Real *Y, Real *Z, Real *W ) ;

# endif
