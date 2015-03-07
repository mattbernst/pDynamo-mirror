/*------------------------------------------------------------------------------
! . File      : RandomNumberDistribution.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _RANDOMNUMBERDISTRIBUTIONS
# define _RANDOMNUMBERDISTRIBUTIONS

# include "RandomNumberGenerator.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Gaussian/Normal distribution. */
extern Real RandomNumberDistribution_GaussianBoxMueller ( RandomNumberGenerator *rng, const Real mu, const Real sigma ) ;

# endif
