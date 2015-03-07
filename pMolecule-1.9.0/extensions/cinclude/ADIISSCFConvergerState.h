/*------------------------------------------------------------------------------
! . File      : ADIISSCFConvergerState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _ADIISSCFCONVERGERSTATE
# define _ADIISSCFCONVERGERSTATE

# include "AntisymmetricMatrix.h"
# include "Definitions.h"
# include "QCOnePDM.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The ADIIS SCF converger iteration type. */
typedef enum {
    ADIISSCFConvergerIterationType_ADIIS = 0,
    ADIISSCFConvergerIterationType_DIIS  = 1,
    ADIISSCFConvergerIterationType_Null  = 2,
    ADIISSCFConvergerIterationType_ODA   = 3
} ADIISSCFConvergerIterationType ;

/* . The storage type for past iterations. */
typedef struct {
    AntisymmetricMatrix *error   ;
    SymmetricMatrix     *density ;
    SymmetricMatrix     *fock    ;
} ADIISSCFConvergerOnePDM ;

/* . The ADIIS SCF converger state type. */
typedef struct {
    ADIISSCFConvergerIterationType iterationType  ;
    Boolean                        isConverged    ;
    Integer                        diisActive     ;
    Integer                        history        ;
    Integer                        iteration      ;
    Integer                        maximumHistory ;
    Real                           diisError      ;
    Real                           energyChange   ;
    Real                           odaOldEnergy   ;
    Real                           oldEnergy      ;
    Real                           odaMu          ;
    Real                           rmsDifference  ;
    /* . Aliases. */
    QCOnePDM                      *densityP       ;
    QCOnePDM                      *densityQ       ;
    Real2DArray                   *orthogonalizer ;
    SymmetricMatrix               *overlap        ;
    /* . Allocated arrays. */
    ADIISSCFConvergerOnePDM       *densitiesP     ;
    ADIISSCFConvergerOnePDM       *densitiesQ     ;
    ADIISSCFConvergerOnePDM       *odaDensityP    ;
    ADIISSCFConvergerOnePDM       *odaDensityQ    ;
    AntisymmetricMatrix           *asmWork        ;
    Integer1DArray                *historyIndices ;
    Real1DArray                   *adiisAlphas    ;
    Real1DArray                   *adiisGradients ;
    Real1DArray                   *energies       ;
    Real2DArray                   *adiisTracesP   ;
    Real2DArray                   *adiisTracesQ   ;
    Real2DArray                   *diisTraces     ;
    SymmetricMatrix               *adiisHessian   ;
} ADIISSCFConvergerState ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void                    ADIISSCFConvergerState_ADIISSetUp   (       ADIISSCFConvergerState  *self, const Boolean useEDIIS ) ;
extern void                    ADIISSCFConvergerState_ADIISUpdate  (       ADIISSCFConvergerState  *self ) ;
extern ADIISSCFConvergerState *ADIISSCFConvergerState_Allocate     ( void ) ;
extern Boolean                 ADIISSCFConvergerState_Converged    ( const ADIISSCFConvergerState  *self ) ;
extern void                    ADIISSCFConvergerState_Deallocate   (       ADIISSCFConvergerState **self ) ;
extern void                    ADIISSCFConvergerState_DIISIterate  (       ADIISSCFConvergerState  *self ) ;
extern Real                    ADIISSCFConvergerState_DIISSave     (       ADIISSCFConvergerState  *self, const Real energy ) ;
extern void                    ADIISSCFConvergerState_ODAIterate   (       ADIISSCFConvergerState  *self, const Real energy, const Real minimumMu ) ;
extern void                    ADIISSCFConvergerState_ODASave      (       ADIISSCFConvergerState  *self ) ;
extern ADIISSCFConvergerState *ADIISSCFConvergerState_SetUp        ( const Boolean QUSERCA, const Integer maximumHistory, QCOnePDM *densityp, QCOnePDM *densityq,
                                                                                        SymmetricMatrix *overlap, Real2DArray *orthogonalizer, Status *status ) ;

# endif
