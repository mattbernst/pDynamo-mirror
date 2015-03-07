/*------------------------------------------------------------------------------
! . File      : QCModelDFT.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCMODELDFT
# define _QCMODELDFT

/*# define TESTGRADIENTS*/

# include "BlockStorage.h"
# include "Coordinates3.h"
# include "Definitions.h"
# include "DFTFunctionalModel.h"
# include "DFTGrid.h"
# include "QCAtomContainer.h"
# include "QCChargeModelOptions.h"
# include "QCOnePDM.h"
# include "QCModelDFTState.h"
# include "QCOnePDM.h"
# include "QCParameters.h"
# include "Real1DArray.h"
# include "SymmetricMatrix.h"
# include "Vector3.h"

# ifdef TESTGRADIENTS
# include "Memory.h"
# include "RysQuadrature.h"
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The DFT QC model type. */
typedef struct {
    Boolean               inCore                    ;
    Boolean               isSpinRestricted          ;
    Boolean               keepOrbitalData           ;
    Boolean               linkAtomRatio             ;
    Integer               numberFractionalAlphaHOOs ;
    Integer               numberFractionalAlphaLUOs ;
    Integer               numberFractionalBetaHOOs  ;
    Integer               numberFractionalBetaLUOs  ;
    Integer               numberFractionalHOOs      ;
    Integer               numberFractionalLUOs      ;
    Real                  fermiBroadening           ;
    DFTFunctionalModel   *functionalModel           ;
    DFTGridAccuracy       accuracy                  ;
    QCChargeModel         qcChargeModel             ;
    QCOnePDMOccupancyType occupancyType             ;
} QCModelDFT ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern QCModelDFT  *QCModelDFT_Allocate                          ( void ) ;
extern void         QCModelDFT_AtomicCharges                     ( const QCModelDFT  *self, QCModelDFTState *qcState, const QCChargeModel *qcchargemodel, const Boolean QSPIN, Real1DArray *qcCharges ) ;
extern QCModelDFT  *QCModelDFT_Clone                             ( const QCModelDFT  *self ) ;
extern void         QCModelDFT_Deallocate                        (       QCModelDFT **self ) ;
extern Vector3     *QCModelDFT_DipoleMoment                      ( const QCModelDFT  *self, QCModelDFTState *qcState, const Coordinates3 *coordinates3, const Vector3 *center ) ;
extern void         QCModelDFT_Fock                              ( const QCModelDFT  *self, QCModelDFTState *qcState, Real *eelectronic ) ;
extern void         QCModelDFT_Gradients                         ( const QCModelDFT  *self, QCModelDFTState *qcState ) ;
extern void         QCModelDFT_GridPointDensities                ( const QCModelDFT  *self, QCModelDFTState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data, const Boolean QSPIN ) ;
extern void         QCModelDFT_GridPointOrbitals                 ( const QCModelDFT  *self, QCModelDFTState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data, const Integer norbitals, Integer *orbitals, const Boolean QALPHA ) ;
extern void         QCModelDFT_GridPointPotentials               ( const QCModelDFT  *self, QCModelDFTState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data ) ;
extern void         QCModelDFT_Integrals                         ( const QCModelDFT  *self, QCModelDFTState *qcState ) ;
extern void         QCModelDFT_MakeOrthogonalizingTransformation ( const QCModelDFT  *self, QCModelDFTState *qcState ) ;
extern void         QCModelDFT_MayerBondOrders                   ( const QCModelDFT  *self, QCModelDFTState *qcState, const QCChargeModel *qcchargemodel, SymmetricMatrix *bondorders, Real1DArray *charges, Real1DArray *freevalence, Real1DArray *totalvalence ) ;
extern Real1DArray *QCModelDFT_OrbitalEnergies                   ( const QCModelDFT  *self, QCModelDFTState *qcState, const Boolean QALPHAORBITALS, Integer *homo, Integer *lumo ) ;
extern Real2DArray *QCModelDFT_Orbitals                          ( const QCModelDFT  *self, QCModelDFTState *qcState, const Boolean QALPHA ) ;
extern void         QCModelDFT_QCMMFock                          ( const QCModelDFT  *self, QCModelDFTState *qcState, Real *eelectronic ) ;
extern void         QCModelDFT_QCMMMakeWeightedDensity           ( const QCModelDFT  *self, QCModelDFTState *qcState ) ;

# ifdef TESTGRADIENTS
extern void QCModelDFT_GradientsTest ( QCAtomContainer    *qcAtoms         ,
                                       QCParameter        *qcParameters    ,
                                       Coordinates3       *qccoordinates   ,
                                       DFTFunctionalModel *functionalModel ,
                                       DFTGrid            *dftgrid         ,
                                       QCOnePDM           *densityp        ,
                                       QCOnePDM           *densityq        ,
                                       Real1DArray        *fpotential      ,
                                       SymmetricMatrix    *wdensity        ,
                                       Coordinates3       *gradients       ) ;
# endif

# endif
