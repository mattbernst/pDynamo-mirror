/*------------------------------------------------------------------------------
! . File      : QCModelMNDOState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCMODELMNDOSTATE
# define _QCMODELMNDOSTATE

# include "BlockStorage.h"
# include "Coordinates3.h"
# include "DefineStatements.h"
# include "Definitions.h"
# include "Integer.h"
# ifdef MNDOCI
# include "MNDOCIState.h"
# endif /*MNDOCI*/
# include "MNDOParameters.h"
# include "QCAtomContainer.h"
# include "QCOnePDM.h"
# include "QCMMInteractionState.h"
# include "QCParameters.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "SymmetricMatrix.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The MNDO QC model state type. */
typedef struct {
    /* . Scalar data. */
    Boolean               isConverged          ;
    Integer               ncycles              ;
    Real                  eelectronic          ;
    Real                  enuclear             ;
    Real                  eocc                 ;
    Real                  eoei                 ;
    Real                  eqcmm                ;
    Real                  eqcqc                ;
    Real                  etei                 ;
    /* . Aliases - persistent. */
    QCAtomContainer      *qcAtoms              ;
    QCParameter          *qcParameters         ;
    /* . Aliases - temporary. */
    Coordinates3         *coordinates3         ;
    Coordinates3         *gradients3           ;
    QCMMInteractionState *qcmmstate            ;
    /* . Allocated data - persistent. */
    QCOnePDM             *densityp             ;
    QCOnePDM             *densityq             ;
    Real2DArray          *c2o                  ;
    /* . Allocated data - temporary. */
    BlockStorage         *twoelectronintegrals ;
    Coordinates3         *qccoordinates3       ;
    Real2DArray          *qptransformation     ;
    SymmetricMatrix      *oneelectronmatrix    ;
    /* . Backup data. */
    SymmetricMatrix      *fockA                ;
    SymmetricMatrix      *fockB                ;
# ifdef MNDOCI
    /* . CI state. */
    MNDOCIState          *cistate              ;
# endif /*MNDOCI*/
} QCModelMNDOState ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern QCModelMNDOState *QCModelMNDOState_Allocate                   ( void ) ;
extern void              QCModelMNDOState_BackupFock                 ( QCModelMNDOState  *self, Status *status ) ;
extern void              QCModelMNDOState_Deallocate                 ( QCModelMNDOState **self ) ;
extern Real              QCModelMNDOState_EnergyFromFock             ( const QCModelMNDOState *self ) ;
extern void              QCModelMNDOState_Finalize                   ( QCModelMNDOState  *self, const Boolean QKEEPDATA ) ;
extern void              QCModelMNDOState_GetGradientDensityTerms    ( QCModelMNDOState  *self, const MNDOParameters *idata, const Integer ifirst, const MNDOParameters *jdata, const Integer jfirst,
                                                                                                                                         const SymmetricMatrix *dtotal, const SymmetricMatrix *dspin,
                                                                                                                                   Real1DArray **dOneI, Real1DArray **dOneJ, Real2DArray **dTwoIJ ) ;
extern void              QCModelMNDOState_Initialize                 ( QCModelMNDOState  *self, Coordinates3 *coordinates3, QCMMInteractionState *qcmmstate, Coordinates3 *gradients3 ) ;
extern Real              QCModelMNDOState_MakeDensities              ( QCModelMNDOState  *self ) ;
extern void              QCModelMNDOState_MakeOrbitalTransformations ( QCModelMNDOState  *self ) ;
extern void              QCModelMNDOState_RestoreFock                ( QCModelMNDOState  *self ) ;
extern QCModelMNDOState *QCModelMNDOState_Setup                      ( QCAtomContainer *qcAtoms, QCParameter *qcParameters, const Real alphaCharge, const Real betaCharge,
                                                                                                const QCOnePDMOccupancyType occupancyType, const Boolean isSpinRestricted,
                                                                                         const Integer numberFractionalHOOs     , const Integer numberFractionalLUOs     ,
                                                                                         const Integer numberFractionalAlphaHOOs, const Integer numberFractionalAlphaLUOs,
                                                                                         const Integer numberFractionalBetaHOOs , const Integer numberFractionalBetaLUOs ,
                                                                                                                                            const Real fermiBroadening ) ;

# endif
