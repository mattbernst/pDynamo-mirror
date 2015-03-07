/*------------------------------------------------------------------------------
! . File      : MNDOCIModel.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _MNDOCIMODEL
# define _MNDOCIMODEL

# include "Boolean.h"
# include "Integer.h"
# include "Integer2DArray.h"
# include "JDEigenvalueSolver.h"
# include "MNDOCIState.h"
# include "QCOnePDM.h"
# include "Real.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The MNDO CI algorithm type. */
typedef enum {
    MNDOCIAlgorithm_Direct = 0, /* . Not yet implemented. */
    MNDOCIAlgorithm_Full   = 1,
    MNDOCIAlgorithm_Sparse = 2
} MNDOCIAlgorithm ;

/* . The MNDO CI method type. */
typedef enum {
    MNDOCIMethod_Doubles        = 0,
    MNDOCIMethod_Full           = 1,
    MNDOCIMethod_Singles        = 2,
    MNDOCIMethod_SinglesDoubles = 3,
    MNDOCIMethod_UserSpecified  = 4
} MNDOCIMethod ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The MNDO CI model type. */
typedef struct {
    Boolean            checkAlgorithm               ;
    Boolean            doAllStates                  ; /* . Unused for the moment. */
    Boolean            identifyRootSpin             ;
    Boolean            localizeOrbitals             ;
    Integer            activeElectrons              ;
    Integer            activeOrbitals               ;
    Integer            localizeStart                ;
    Integer            localizeStop                 ;
    Integer            minimalMultiplicity          ;
    Integer            numberOfStates               ;
    Integer            requiredRoot                 ;
    Integer            rootMultiplicity             ;
    Real               degeneracyTolerance          ;
    Real               fractionalOccupancyTolerance ;
    Real               spinTolerance                ;
    MNDOCIAlgorithm    algorithm                    ;
    MNDOCIMethod       method                       ;
    Integer2DArray    *microstates                  ;
    JDEigenvalueSolver eigenvalueSolver             ;
} MNDOCIModel ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern MNDOCIModel *MNDOCIModel_Allocate   ( Status *status ) ;
extern MNDOCIModel *MNDOCIModel_Clone      ( const MNDOCIModel  *self, Status *status ) ;
extern void         MNDOCIModel_Deallocate (       MNDOCIModel **self ) ;
extern void         MNDOCIModel_Energy     ( const MNDOCIModel  *self, MNDOCIState *state, const Real electronic, const Real enuclear, QCOnePDM *densityp, QCOnePDM *densityq, Status *status ) ;
extern MNDOCIState *MNDOCIModel_MakeState  ( const MNDOCIModel  *self, const Integer multiplicity, const Integer nalpha, const Integer nbeta, const Integer norbitals, Status *status ) ;

# endif
