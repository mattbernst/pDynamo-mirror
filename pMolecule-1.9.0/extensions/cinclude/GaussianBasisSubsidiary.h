/*------------------------------------------------------------------------------
! . File      : GaussianBasisSubsidiary.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _GAUSSIANBASISSUBSIDIARY
# define _GAUSSIANBASISSUBSIDIARY

# include "Definitions.h"

/*------------------------------------------------------------------------------
! . Procedures.
!-----------------------------------------------------------------------------*/
extern void Subsidiary_Integral_Derivative2 ( const Real *x, const Real *y, const Real *z, const Real a, const Integer ni, const Integer nj, const Integer jdim, Real *xd, Real *yd, Real *zd ) ;
extern void Subsidiary_Integral_Derivative3 ( const Real *x, const Real *y, const Real *z, Real *xg, Real *yg, Real *zg,
                                                             Real *xh, Real *yh, Real *zh, const Real ag, const Real ah,
                                                                   const Integer ni, const Integer nj, const Integer nf,
                                                                                 const Integer dim1, const Integer dim2,
                                                                            const Integer ddim1, const Integer ddim2 ) ;
extern void Subsidiary_Integral_Dipole      ( Real *x, Real *y, Real *z, const Real aa, const Real *r0, const Real *ri, const Real *rj, const Real *center, const Integer ni, const Integer nj ) ;
extern void Subsidiary_Integral_Kinetic     ( const Real *x, const Real *y, const Real *z, Real *xt, Real *yt, Real *zt, const Real aj, const Integer ni, const Integer nj, const Integer jdimo, const Integer jdimt ) ;
extern void Subsidiary_Integral_Nuclear2C   ( const Integer iangmom, const Integer jangmom, const Real b00, const Real b10, const Real bp01, const Real f00, const Real xc00,
                                                                                      const Real xcp00, const Real yc00, const Real ycp00, const Real zc00, const Real zcp00,
                                                                                                                   const Integer jdim, Real *xint, Real *yint, Real *zint ) ;
extern void Subsidiary_Integral_Nuclear3C   ( const Integer ni, const Integer nj, const Integer nk, const Boolean qij0, const Boolean qij1, const Boolean qn0, const Boolean qn1,
                                              const Real b00, const Real b10, const Real bp01, const Real dxij, const Real dyij, const Real dzij,
                                              const Real f00, const Real xc00, const Real xcp00, const Real  yc00, const Real  ycp00, const Real zc00,
                                              const Real zcp00, const Integer jdim1, const Integer jdim2, Real *xint, Real *yint, Real *zint ) ;
extern void Subsidiary_Integral_Overlap2    ( Real *x, Real *y, Real *z, const Real aa, const Real *r0, const Real *ri, const Real *rj, const Integer ni, const Integer nj ) ;
extern void Subsidiary_Integral_Overlap3    ( Real *x, Real *y, Real *z, const Real aa, const Real *r0, const Real *ri, const Real *rj, const Real *rk, const Integer ni, const Integer nj, const Integer nk ) ;

# endif
