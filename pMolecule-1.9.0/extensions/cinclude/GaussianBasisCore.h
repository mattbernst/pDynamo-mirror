/*------------------------------------------------------------------------------
! . File      : GaussianBasisCore.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _GAUSSIANBASISCORE
# define _GAUSSIANBASISCORE

# include <stdio.h>

# include "Definitions.h"
# include "GaussianBasis.h"
# include "Real2DArray.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void         GaussianBasis_Normalize             (       GaussianBasis *self      ) ;
extern void         GaussianBasis_ScaleExponents        (       GaussianBasis *self, ... ) ;
extern void         GaussianBasis_UnnormalizePrimitives (       GaussianBasis *self      ) ;

extern Real2DArray *GaussianBasis_2Coulomb              ( const GaussianBasis *ibasis, const Real *ri, const GaussianBasis *jbasis, const Real *rj ) ;
extern Real2DArray *GaussianBasis_2Overlap              ( const GaussianBasis *ibasis, const Real *ri, const GaussianBasis *jbasis, const Real *rj ) ;
extern void         GaussianBasis_2OverlapD             ( const GaussianBasis *ibasis, const Real *ri, const GaussianBasis *jbasis, const Real *rj,
                                                                                        Real2DArray **sijx, Real2DArray **sijy, Real2DArray **sijz ) ;
extern void         GaussianBasis_2OverlapFD            ( const GaussianBasis  *ibasis ,
                                                          const Real           *ri     ,
                                                          const GaussianBasis  *jbasis ,
                                                          const Real           *rj     ,
                                                                Real2DArray   **sij    ,
                                                                Real2DArray   **sijx   ,
                                                                Real2DArray   **sijy   ,
                                                                Real2DArray   **sijz   ) ;

# endif
