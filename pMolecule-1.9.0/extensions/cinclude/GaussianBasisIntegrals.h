/*------------------------------------------------------------------------------
! . File      : GaussianBasisIntegrals.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _GAUSSIANBASISINTEGRALS
# define _GAUSSIANBASISINTEGRALS

# include "BlockStorage.h"
# include "Coordinates3.h"
# include "Definitions.h"
# include "GaussianBasis.h"
# include "Macros.h"
# include "QCAtomContainer.h"
# include "QCOnePDM.h"
# include "Real1DArray.h"
# include "SymmetricMatrix.h"

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_Dipole  ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                             const Real *center, SymmetricMatrix **dipx, SymmetricMatrix **dipy , SymmetricMatrix **dipz ) ;
extern void GaussianBasis_Electron_Fit     ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                                                                                                        BlockStorage **fitintegrals ) ;
extern void GaussianBasis_Electron_Nuclear ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                                                                                                  SymmetricMatrix *onelectronmatrix ) ;
extern void GaussianBasis_Fit_Fit          ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                                                            const Real1DArray *fselfoverlap, SymmetricMatrix **inversefitmatrix ) ;
extern Real GaussianBasis_Fock             ( const Integer nfbasis, QCOnePDM *density, BlockStorage *fitintegrals, Real1DArray *fpotential,
                                                                                                              SymmetricMatrix *inversefitmatrix ) ;
extern void GaussianBasis_Kinetic_2Overlap ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Coordinates3 *qccoordinates3,
                                                                                          SymmetricMatrix **kinetic, SymmetricMatrix **overlap ) ;
extern Real GaussianBasis_Nuclear_Nuclear  ( const QCAtomContainer *qcAtoms, Coordinates3 *qccoordinates3 ) ;
extern void GaussianBasis_Point_Electron   ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, const Coordinates3 *qccoordinates3, const Coordinates3 *points,
                                                                             const SymmetricMatrix *odensityt, Real1DArray *potentials ) ;
extern void GaussianBasis_Point_Nuclear    ( const QCAtomContainer *qcAtoms, const Real1DArray *znuclear, const Coordinates3 *qccoordinates3, const Coordinates3 *points, Real1DArray *potentials ) ;
extern void GaussianBasis_SelfOverlap      ( const QCAtomContainer *qcAtoms, const QCParameter *qcParameters, Real1DArray **fselfoverlap ) ;

# endif
