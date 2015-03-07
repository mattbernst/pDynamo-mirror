/*------------------------------------------------------------------------------
! . File      : QCModelDFTState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCMODELDFTSTATE
# define _QCMODELDFTSTATE

# include "BlockStorage.h"
# include "Coordinates3.h"
# include "Definitions.h"
# include "DFTGrid.h"
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
/* . The DFT QC model state type. */
typedef struct {
    /* . Scalar data. */
    Boolean               isConverged          ;
    Integer               ncycles              ;
# ifdef USEOPENMP
    Integer               numberOfThreads      ;
# endif
    Real                  eelectronic          ;
    Real                  enuclear             ;
    Real                  eocc                 ;
    Real                  eoei                 ;
    Real                  eqcmm                ;
    Real                  eqcqc                ;
    Real                  equad                ;
    Real                  etei                 ;
    Real                  rhoquad              ;
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
    Real2DArray          *o2c                  ;
    /* . Allocated data - temporary. */
    BlockStorage         *fitintegrals         ;
    Coordinates3         *qccoordinates3       ;
# ifdef USEOPENMP
    Coordinates3        **qcGradients3Array    ;
# endif
    DFTGrid              *dftgrid              ;
    Real1DArray          *fpotential           ;
    Real1DArray          *fselfoverlap         ;
    Real1DArray          *overlapeigenvalues   ;
    Real1DArray          *wvector              ;
    Real2DArray          *lowdintransformation ;
    Real2DArray          *orthogonalizer       ;
    Real2DArray          *overlapeigenvectors  ;
# ifdef USEOPENMP
    SymmetricMatrix     **fockAArray           ;
    SymmetricMatrix     **fockBArray           ;
# endif
    SymmetricMatrix      *inversefitmatrix     ;
    SymmetricMatrix      *oneelectronmatrix    ;
    SymmetricMatrix      *overlap              ;
    SymmetricMatrix      *wdensity             ;
} QCModelDFTState ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifdef USEOPENMP
extern void             QCModelDFTState_AccumulateFock       ( QCModelDFTState  *self, SymmetricMatrix **fockArray ) ;
extern void             QCModelDFTState_AccumulateGradients3 ( QCModelDFTState  *self ) ;
# endif
extern QCModelDFTState *QCModelDFTState_Allocate             ( void ) ;
extern void             QCModelDFTState_Deallocate           ( QCModelDFTState **self ) ;
extern void             QCModelDFTState_Finalize             ( QCModelDFTState  *self, const Boolean QKEEPDATA ) ;
extern void             QCModelDFTState_Initialize           ( QCModelDFTState  *self, Coordinates3 *coordinates3, QCMMInteractionState *qcmmstate, Coordinates3 *gradients3 ) ;
# ifdef USEOPENMP
extern void             QCModelDFTState_InitializeFockArray  ( QCModelDFTState  *self, QCOnePDM *qcOnePDM, SymmetricMatrix **fockArray ) ;
# endif
extern QCModelDFTState *QCModelDFTState_Setup                (       QCAtomContainer      *qcAtoms                   ,
                                                                     QCParameter          *qcParameters              ,
                                                               const Real                  alphaCharge               ,
                                                               const Real                  betaCharge                ,
                                                               const QCOnePDMOccupancyType occupancyType             ,
                                                               const Boolean               isSpinRestricted          ,
                                                               const Integer               numberFractionalHOOs      ,
                                                               const Integer               numberFractionalLUOs      ,
                                                               const Integer               numberFractionalAlphaHOOs ,
                                                               const Integer               numberFractionalAlphaLUOs ,
                                                               const Integer               numberFractionalBetaHOOs  ,
                                                               const Integer               numberFractionalBetaLUOs  ,
                                                               const Real                  fermiBroadening           ) ;
# endif
