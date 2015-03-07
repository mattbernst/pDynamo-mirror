/*------------------------------------------------------------------------------
! . File      : Correlation.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _CORRELATION
# define _CORRELATION

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real1DArray *Correlation_MakeDotProduct ( Real2DArray *x, Real2DArray *y, Real1DArray *c, const Boolean useFFT, const Boolean normalize, const Boolean removeMean, const Integer tCorrelation, const Real *tolerance, Status *status ) ;
extern Real1DArray *Correlation_MakeSimple     ( Real1DArray *x, Real1DArray *y, Real1DArray *c, const Boolean useFFT, const Boolean normalize, const Boolean removeMean, const Integer tCorrelation, const Real *tolerance, Status *status ) ;

# endif
