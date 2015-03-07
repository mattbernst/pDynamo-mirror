/*------------------------------------------------------------------------------
! . File      : QCOnePDM.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCONEPDM
# define _QCONEPDM

# include "Boolean.h"
# include "Definitions.h"
# include "GridFunctionDataBlock.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Density type. */
typedef enum {
    QCOnePDMDensityType_Alpha = 0,
    QCOnePDMDensityType_Beta  = 1,
    QCOnePDMDensityType_Spin  = 2,
    QCOnePDMDensityType_Total = 3
} QCOnePDMDensityType ;

/* . Occupancy type. */
typedef enum {
    QCOnePDMOccupancyType_Cardinal           = 0,
    QCOnePDMOccupancyType_FractionalFixed    = 1,
    QCOnePDMOccupancyType_FractionalVariable = 2
} QCOnePDMOccupancyType ;

/* . The QC one particle density matrix type. */
typedef struct {
   Boolean                isValid         ;
   Integer                numberOccupied  ;
   Integer                numberOrbitals  ;
   Real                   fermiBroadening ;
   Real                   fermiEnergy     ;
   Real                   occupancyEnergy ;
   Real                   occupancyFactor ; /* . 1 or 2 depending on density type. */
   Real                   totalCharge     ;
   QCOnePDMDensityType    densityType     ;
   QCOnePDMOccupancyType  occupancyType   ;
   Real1DArray           *energies        ;
   Real1DArray           *occupancies     ; /* . Full occupancy. */
   Real2DArray           *orbitals        ;
   SymmetricMatrix       *density         ; /* . Full density. */
   SymmetricMatrix       *fock            ;
} QCOnePDM ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern QCOnePDM        *QCOnePDM_Allocate                 ( const Integer length, Status *status ) ;
extern void             QCOnePDM_AllocateOrbitalData      (       QCOnePDM  *self, Real2DArray *orthogonalizer, Status *status ) ;
extern QCOnePDM        *QCOnePDM_Clone                    ( const QCOnePDM  *self, Status *status ) ;
extern void             QCOnePDM_AlphaBetaFromTotalSpin   (       QCOnePDM  *self, QCOnePDM *other ) ;
extern void             QCOnePDM_AlphaBetaToTotalSpin     (       QCOnePDM  *self, QCOnePDM *other ) ;
extern void             QCOnePDM_Deallocate               (       QCOnePDM **self ) ;
extern void             QCOnePDM_DeallocateOrbitalData    (       QCOnePDM  *self ) ;
extern void             QCOnePDM_DensityGridValues        ( const QCOnePDM  *self, const GridFunctionDataBlock *basisData, Real1DArray *rho, Status *status ) ;
extern QCOnePDM        *QCOnePDM_FromDiagonalGuess        ( const QCOnePDMDensityType densityType, const QCOnePDMOccupancyType occupancyType, const Integer length, const Real totalCharge,
                                                                                                   const Integer numberFractionalHOOs, const Integer numberFractionalLUOs, Status *status ) ;
extern Integer          QCOnePDM_HOMOIndex                ( const QCOnePDM  *self ) ;
extern Integer          QCOnePDM_LUMOIndex                ( const QCOnePDM  *self ) ;
extern Real             QCOnePDM_Make                     (       QCOnePDM  *self, SymmetricMatrix *scratch, Status *status ) ;
extern void             QCOnePDM_MakeElementary           ( SymmetricMatrix *self, const Integer noccupied, const Real1DArray *occupancies, const Real2DArray *orbitals, Status *status ) ;
extern Real             QCOnePDM_MakeFromFock             (       QCOnePDM  *self, Real2DArray *orthogonalizer, Status *status ) ;
extern SymmetricMatrix *QCOnePDM_MakeWeighted             ( const QCOnePDM  *self, const QCOnePDM *other, Status *status ) ;
extern void             QCOnePDM_Reorthogonalize          (       QCOnePDM  *self, Real2DArray *orthogonalizer, Status *status ) ;
extern Boolean          QCOnePDM_ResetDensityFromDensity  (       QCOnePDM  *self, SymmetricMatrix *density, SymmetricMatrix *overlap, Status *status ) ;
extern Boolean          QCOnePDM_ResetDensityFromOrbitals (       QCOnePDM  *self, Real2DArray *orbitals, Status *status ) ;
extern void             QCOnePDM_SpinExpectationValues    (       QCOnePDM  *self, QCOnePDM *other, SymmetricMatrix *overlap, Real *Sz, Real *S2 ) ;

# endif
