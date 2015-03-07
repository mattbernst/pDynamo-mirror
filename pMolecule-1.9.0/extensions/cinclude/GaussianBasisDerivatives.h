/*------------------------------------------------------------------------------
! . File      : GaussianBasisDerivatives.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _GAUSSIANBASISDERIVATIVES
# define _GAUSSIANBASISDERIVATIVES

# include "Coordinates3.h"
# include "Definitions.h"
# include "GaussianBasis.h"
# include "QCAtomContainer.h"
# include "Real1DArray.h"
# include "SymmetricMatrix.h"

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern void GaussianBasis_Electron_FitD     ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                                      SymmetricMatrix *odensityt, const Real1DArray *fpotential, const Real1DArray *wvector,
                                                                                                              Coordinates3 *qcGradients3 ) ;
extern void GaussianBasis_Electron_NuclearD ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                                                                  SymmetricMatrix *odensityt, Coordinates3 *qcGradients3 ) ;
extern void GaussianBasis_Fit_FitD          ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                                                                  const Real1DArray *fpotential, const Real1DArray *wvector,
                                                                                                              Coordinates3 *qcGradients3 ) ;
extern void GaussianBasis_Kinetic_2OverlapD ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                                                         const SymmetricMatrix *odensityt, const SymmetricMatrix *odensityw,
                                                                                                              Coordinates3 *qcGradients3 ) ;
extern void GaussianBasis_Nuclear_NuclearD  ( const QCAtomContainer *qcAtoms, Coordinates3 *qccoordinates3, Coordinates3 *qcGradients3 ) ;

# endif
