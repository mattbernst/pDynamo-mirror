/*------------------------------------------------------------------------------
! . File      : QCModelMNDO.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCMODELMNDO
# define _QCMODELMNDO

# include "ChargeConstraintContainer.h"
# include "Coordinates3.h"
# include "DefineStatements.h"
# include "Definitions.h"
# include "Integer1DArray.h"
# include "Matrix33.h"
# ifdef MNDOCI
# include "MNDOCIModel.h"
# endif /*MNDOCI*/
# include "QCChargeModelOptions.h"
# include "QCModelMNDOState.h"
# include "QCOnePDM.h"
# include "Real1DArray.h"
# include "Selection.h"
# include "SymmetricMatrix.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The MNDO QC model type. */
typedef struct {
    Boolean               isSpinRestricted          ;
    Boolean               keepOrbitalData           ;
    Boolean               linkAtomRatio             ;
    Real                  fermiBroadening           ;
    Integer               numberFractionalAlphaHOOs ;
    Integer               numberFractionalAlphaLUOs ;
    Integer               numberFractionalBetaHOOs  ;
    Integer               numberFractionalBetaLUOs  ;
    Integer               numberFractionalHOOs      ;
    Integer               numberFractionalLUOs      ;
# ifdef MNDOCI
    MNDOCIModel          *ciModel                   ;
# endif /*MNDOCI*/
    QCChargeModel         qcChargeModel             ;
    QCOnePDMOccupancyType occupancyType             ;
} QCModelMNDO ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern QCModelMNDO *QCModelMNDO_Allocate            ( void ) ;
extern void         QCModelMNDO_CCFockLowdin        ( const QCModelMNDO  *self, QCModelMNDOState *qcState, const ChargeConstraintContainer *chargeConstraints, const Real1DArray *lambdas, Status *status ) ;
extern Real         QCModelMNDO_CCFunction          ( const QCModelMNDO               *self              ,
                                                            QCModelMNDOState          *qcState           ,
                                                      const ChargeConstraintContainer *chargeConstraints ,
                                                      const Real1DArray               *lambdas           ,
                                                      const Boolean                    buildFock         ,
                                                            Real1DArray               *gradients         ,
                                                            SymmetricMatrix           *hessian           ,
                                                            Status                    *status            ) ;
extern void         QCModelMNDO_CCHessian           ( const QCModelMNDO  *self, QCModelMNDOState *qcState, const ChargeConstraintContainer *chargeConstraints, SymmetricMatrix *hessian, Status *status ) ;
# ifdef MNDOCI
extern void         QCModelMNDO_CIEnergy            ( const QCModelMNDO  *self, QCModelMNDOState *state ) ;
extern Real1DArray *QCModelMNDO_CIStateCharacters   ( const QCModelMNDO  *self, QCModelMNDOState *state, const Matrix33 *rotation, const Integer1DArray *mapping,
                                                                             const Selection *stateIndices, const Boolean includeCoreOrbitals, Status *status ) ;
# endif /*MNDOCI*/
extern QCModelMNDO *QCModelMNDO_Clone               ( const QCModelMNDO  *self ) ;
extern void         QCModelMNDO_Deallocate          (       QCModelMNDO **self ) ;
extern Vector3     *QCModelMNDO_DipoleMoment        ( const QCModelMNDO  *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3, const Vector3 *center ) ;
extern void         QCModelMNDO_Fock                ( const QCModelMNDO  *self, QCModelMNDOState *qcState, Real *eelectronic ) ;
extern void         QCModelMNDO_Gradients           ( const QCModelMNDO  *self, QCModelMNDOState *qcState ) ;
extern void         QCModelMNDO_GridPointDensities  ( const QCModelMNDO  *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data, const Boolean QSPIN ) ;
extern void         QCModelMNDO_GridPointOrbitals   ( const QCModelMNDO  *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data, const Integer norbitals, Integer *orbitals, const Boolean QALPHA ) ;
extern void         QCModelMNDO_GridPointPotentials ( const QCModelMNDO  *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data ) ;
extern void         QCModelMNDO_Integrals           ( const QCModelMNDO  *self, QCModelMNDOState *qcState ) ;
extern void         QCModelMNDO_LocalizeOrbitals    ( const QCModelMNDO  *self, QCModelMNDOState      *state, const Boolean doQ, const Integer startOrbital, const Integer stopOrbital, Status *status ) ;
extern void         QCModelMNDO_LowdinCharges       ( const QCModelMNDO  *self, QCModelMNDOState *qcState, const Boolean QSPIN, Real1DArray *qcCharges ) ;
extern void         QCModelMNDO_MayerBondOrders     ( const QCModelMNDO  *self, QCModelMNDOState *qcState, SymmetricMatrix *bondorders, Real1DArray *charges, Real1DArray *freevalence, Real1DArray *totalvalence ) ;
extern Real1DArray *QCModelMNDO_OrbitalCharacters   ( const QCModelMNDO *self, QCModelMNDOState *state, const Matrix33 *rotation, const Integer1DArray *mapping,
                                                                                  const Boolean useDensityP, const Selection *orbitalIndices, Status *status ) ;
extern Real1DArray *QCModelMNDO_OrbitalEnergies     ( const QCModelMNDO  *self, QCModelMNDOState *qcState, const Boolean QALPHAORBITALS, Integer *homo, Integer *lumo ) ;
extern Real2DArray *QCModelMNDO_Orbitals            ( const QCModelMNDO  *self, QCModelMNDOState *qcState, const Boolean QALPHA ) ;
extern void         QCModelMNDO_RotateOrbital       ( const QCModelMNDO *self, const QCModelMNDOState *state, const Matrix33 *rotation, const Integer1DArray *mapping,
                                                                                                             const Real1DArray *inOrbital, Real1DArray *outOrbital ) ;
extern void         QCModelMNDO_QCMMFockLowdin      ( const QCModelMNDO  *self, QCModelMNDOState *qcState, Real *eelectronic ) ;
extern void         QCModelMNDO_TransformIntegrals  ( Real2DArray **integrals, const Real2DArray *ic2o, const Real2DArray *jc2o ) ;

# endif
