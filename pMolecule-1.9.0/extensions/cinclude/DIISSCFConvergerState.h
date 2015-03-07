/*------------------------------------------------------------------------------
! . File      : DIISSCFConvergerState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _DIISSCFCONVERGERSTATE
# define _DIISSCFCONVERGERSTATE

# include "AntisymmetricMatrix.h"
# include "Definitions.h"
# include "QCOnePDM.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "SymmetricMatrix.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The DIIS SCF converger type. */
typedef struct {
    Boolean               isConverged    ;
    Boolean               QDAMP          ;
    Boolean               QDIIS          ;
    Boolean               QRCA           ;
    Boolean               QRCAPREDICTED  ;
    Integer               iteration      ;
    Integer               matnum         ;
    Integer               ndiis          ;
    Real                  damp           ;
    Real                  deavg          ;
    Real                  deltaeold      ;
    Real                  diiserror      ;
    Real                  dmpsav         ;
    Real                  eold           ;
    Real                  rcamu          ;
    Real                  rmsdifference  ;
    Integer              *matind         ;
/*    AntisymmetricMatrix  *asmwork        ;*/
    AntisymmetricMatrix **errorp         ;
    AntisymmetricMatrix **errorq         ;
    QCOnePDM             *densityp       ;
    QCOnePDM             *densityq       ;
    Real2DArray          *bcoeff         ;
    Real2DArray          *orthogonalizer ;
    Real2DArray         **work           ;
    SymmetricMatrix      *dfockp1        ;
    SymmetricMatrix      *dfockp2        ;
    SymmetricMatrix      *dfockq1        ;
    SymmetricMatrix      *dfockq2        ;
    SymmetricMatrix      *overlap        ;
    SymmetricMatrix      *rdensityp      ;
    SymmetricMatrix      *rfockp         ;
    SymmetricMatrix      *rdensityq      ;
    SymmetricMatrix      *rfockq         ;
    SymmetricMatrix     **xdensityp      ;
    SymmetricMatrix     **xdensityq      ;
    SymmetricMatrix     **xfockp         ;
    SymmetricMatrix     **xfockq         ;
} DIISSCFConvergerState ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern DIISSCFConvergerState *DIISSCFConvergerState_Allocate     ( void ) ;
extern Boolean                DIISSCFConvergerState_Converged    ( const DIISSCFConvergerState  *self ) ;
extern void                   DIISSCFConvergerState_Deallocate   (       DIISSCFConvergerState **self ) ;
extern void                   DIISSCFConvergerState_GetTableData ( const DIISSCFConvergerState  *self, Integer *iteration, Real *deltaeold, Real *rmsdifference ) ;
extern DIISSCFConvergerState *DIISSCFConvergerState_SetUp        ( const Boolean QUSERCA, const Integer ndiis, QCOnePDM *densityp, QCOnePDM *densityq, SymmetricMatrix *overlap, Real2DArray *orthogonalizer ) ;

# endif
