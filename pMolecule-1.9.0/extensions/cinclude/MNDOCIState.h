/*------------------------------------------------------------------------------
! . File      : MNDOCIState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _MNDOCISTATE
# define _MNDOCISTATE

# include "BlockStorage.h"
# include "Boolean.h"
# include "DoubleSymmetricMatrix.h"
# include "Integer.h"
# include "Integer1DArray.h"
# include "Integer2DArray.h"
# include "JDEigenvalueSolver.h"
# include "QCOnePDM.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "RealNDArray.h"
# include "SparseSymmetricMatrix.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The MNDO CI configuration type. */
typedef struct {
    Integer         nalphas ;
    Integer         nspqr   ;
    Real            spin    ;
    Integer1DArray *alphas  ;
    Integer1DArray *betas   ;
    Integer1DArray *spqr    ;
} MNDOCIConfiguration ;

/* . The MNDO CI state type. */
typedef struct {
    Boolean                  doFull                       ;
    Boolean                  doSparse                     ;
    Boolean                  fractionallyOccupiedInactive ;
    Boolean                  orbitalDegeneracies          ;
    Boolean                  rootNotFound                 ;
    Boolean                  usePreconditioning           ;
    Integer                  ciMatrixNonZero              ;
    Integer                  localizeStart                ;
    Integer                  localizeStop                 ;
    Integer                  nactive                      ;
    Integer                  nconfigurations              ;
    Integer                  ncore                        ;
    Integer                  nelectrons                   ;
    Integer                  norbitals                    ;
    Integer                  numberOfStates               ;
    Integer                  root                         ;
    Real                     baseline                     ;
    Real                     ciEnergy                     ;
    Real                     ciMatrixSparsity             ;
    Real                     ecore                        ;
    Real                     requiredSpin                 ;
    Real                     rootenergy                   ;
    /* . Configurations. */
    MNDOCIConfiguration     *configurations               ;
    /* . Aliases. */
    Real1DArray              *energies                    ;
    Real1DArray              *occupancies                 ;
    Real2DArray              *orbitals                    ;
    BlockStorage             *twoelectronintegrals        ;
    SymmetricMatrix          *oneelectronmatrix           ;
    /* . Allocated data. */
    DoubleSymmetricMatrix    *moteis                      ;
    DoubleSymmetricMatrix    *twopdm                      ;
    Real1DArray              *checkEnergies               ;
    Real1DArray              *ciEnergies                  ;
    Real1DArray              *ciMatrixPreconditioner      ;
    Real1DArray              *ciVector                    ;
    Real1DArray              *ocore                       ;
    Real1DArray              *spins                       ;
    Real2DArray              *checkVectors                ;
    Real2DArray              *ciVectors                   ;
    Real2DArray              *kpa                         ;
    Real2DArray              *kpaMO                       ;
    Real2DArray              *motei34                     ;
    RealNDArray              *motei234                    ;
    SparseSymmetricMatrix    *ciMatrixSparse              ;
    SymmetricMatrix          *ciMatrixFull                ;
    SymmetricMatrix          *fcore                       ;
    SymmetricMatrix          *fcoreMO                     ;
    SymmetricMatrix          *onepdma                     ;
    SymmetricMatrix          *onepdmb                     ;
    SymmetricMatrix          *onepdmMOa                   ;
    SymmetricMatrix          *onepdmMOb                   ;
    SymmetricMatrix          *pcore                       ;
    SymmetricMatrix          *spinDensity                 ;
    /* . Slice. */
    Real2DArray               activemos                   ;
    /* . Solver. */
    JDEigenvalueSolverReport  eigenvalueSolverReport      ;
    JDEigenvalueSolverState  *eigenvalueSolverState       ;
    JDEigenvalueSolverTarget  eigenvalueSolverTarget      ;
    /* . CI gradient items. */
    Boolean                   ciGradientError             ;
    Boolean                   doGradients                 ;
    Integer                   cphfErrorFlag               ;
    Integer                   cphfIterations              ;
    Integer                   numberDegenerateRedundant   ;
    Integer                   numberNonRedundant          ;
    Integer                   numberRedundant             ;
    Integer2DArray           *indicesNR                   ;
    Integer2DArray           *indicesR                    ;
    Real1DArray              *aDiagonal                   ;
    Real1DArray              *qNR                         ;
    Real1DArray              *qR                          ;
    Real1DArray              *preconditioner              ;
    Real1DArray              *tpdm1                       ;
    Real2DArray              *tpdm2                       ;
    RealNDArray              *tpdm3                       ;
    SymmetricMatrix          *onepdm                      ;
    SymmetricMatrix          *onepdmHF                    ;
    SymmetricMatrix          *onepdmMO                    ;
    SymmetricMatrix          *work1                       ;
    SymmetricMatrix          *work2                       ;
    SymmetricMatrix          *zMatrix                     ;
} MNDOCIState ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern MNDOCIState *MNDOCIState_Allocate                ( const Integer nactive, const Integer nconfigurations, Status *status ) ;
extern void         MNDOCIState_CalculateKPA            (       MNDOCIState  *self ) ;
extern void         MNDOCIState_CalculateZMatrix        (       MNDOCIState  *self ) ;
extern void         MNDOCIState_Deallocate              (       MNDOCIState **self ) ;
extern void         MNDOCIState_DiagonalizeCIMatrix     (       MNDOCIState  *self, const JDEigenvalueSolver *eigenvalueSolver, Status *status ) ;
extern void         MNDOCIState_Finalize                (       MNDOCIState  *self, const Boolean keepWavefunction ) ;
extern void         MNDOCIState_FourIndexTransformation (       MNDOCIState  *self ) ;
extern void         MNDOCIState_GetCIMatrixSparsity     (       MNDOCIState  *self ) ;
extern Boolean      MNDOCIState_Initialize              (       MNDOCIState  *self, const Boolean doGradients, QCOnePDM *densityp ) ;
extern void         MNDOCIState_MakeCIMatrix            (       MNDOCIState  *self ) ;
extern void         MNDOCIState_MakeDensities           (       MNDOCIState  *self, QCOnePDM *densityp, QCOnePDM *densityq ) ;
extern MNDOCIState *MNDOCIState_MakeFull                ( const Integer nactive, const Integer nup, const Integer ndown, Status *status ) ;
extern MNDOCIState *MNDOCIState_MakeSinglesDoubles      ( const Boolean doSingles, const Boolean doDoubles, const Boolean doAllStates, const Integer nactive, const Integer nclosed, const Integer nopen, Status *status ) ;
extern MNDOCIState *MNDOCIState_MakeUserSpecified       ( const Integer2DArray *microstates, const Integer activeOrbitals, const Integer activeElectrons, Status *status ) ;

# endif
