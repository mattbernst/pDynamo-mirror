/*------------------------------------------------------------------------------
! . File      : GaussianBasisGrid.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _GAUSSIANBASISGRID
# define _GAUSSIANBASISGRID

# include "Coordinates3.h"
# include "GaussianBasis.h"
# include "GridFunctionDataBlock.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_GridPointValues ( const GaussianBasis         *iBasis ,
                                     const Real                  *rI     ,
                                     const Coordinates3          *rG     ,
                                           GridFunctionDataBlock *data   ) ;
# endif
