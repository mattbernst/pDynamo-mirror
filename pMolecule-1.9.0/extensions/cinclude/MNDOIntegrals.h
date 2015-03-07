/*------------------------------------------------------------------------------
! . File      : MNDOIntegrals.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _MNDOINTEGRALS
# define _MNDOINTEGRALS

# include "BlockStorage.h"
# include "Definitions.h"
# include "MNDOParameters.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDOIntegrals_AddInOneCenterTEIs         ( const MNDOParameters  *self, const Integer istart, BlockStorage *twoelectronintegrals ) ;
extern void MNDOIntegrals_MakeFockG                  ( BlockStorage *twoelectronintegrals, SymmetricMatrix *densitya, SymmetricMatrix *focka, SymmetricMatrix *densityb, SymmetricMatrix *fockb ) ;
extern void MNDOIntegrals_MolecularFrame2CIntegrals  ( const MNDOParameters *idata, const Integer istart, const Real *xi, const MNDOParameters *jdata, const Integer jstart, const Real *xj, Real *enuc,
                                                                                                                              Real1DArray *e1b, Real1DArray *e2a, BlockStorage *twoelectronintegrals ) ;
extern void MNDOIntegrals_MolecularFrame2CIntegralsD ( const MNDOParameters *idata, const Integer ifirst, const Real *xi, const MNDOParameters *jdata, const Integer jfirst, const Real *xj,
                                                                                                const Real1DArray *dOneI, const Real1DArray *dOneJ, const Real2DArray *dTwoIJ, Real *eng ) ;

# endif
