/*------------------------------------------------------------------------------
! . File      : DFTIntegratorDataBlock.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _DFTINTEGRATORDATABLOCK
# define _DFTINTEGRATORDATABLOCK

# include "Boolean.h"
# include "Integer.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The DFT integration block data view type. */
typedef struct {
    Real1DArray dRhoX         ;
    Real1DArray dRhoY         ;
    Real1DArray dRhoZ         ;
    Real1DArray laplacianRho  ;
    Real1DArray rho           ; 
    Real1DArray sigma         ;
    Real1DArray tau           ;
    Real1DArray vLaplacianRho ;
    Real1DArray vRho          ;
    Real1DArray vSigma        ;
    Real1DArray vTau          ;
} DFTIntegratorDataBlockView ;

/* . The DFT integration block data type. */
typedef struct {
    Boolean      hasLocalData        ;
    Integer      numberOfPoints      ;
    Real1DArray *eXC                 ;
    Real1DArray *localEXC            ;
    Real2DArray *dRhoX               ;
    Real2DArray *dRhoY               ;
    Real2DArray *dRhoZ               ;
    Real2DArray *localVLaplacianRho  ;
    Real2DArray *localVRho           ;
    Real2DArray *localVSigma         ;
    Real2DArray *localVTau           ;
    Real2DArray *laplacianRho        ;
    Real2DArray *rho                 ;
    Real2DArray *sigma               ;
    Real2DArray *tau                 ;
    Real2DArray *vLaplacianRho       ;
    Real2DArray *vRho                ;
    Real2DArray *vSigma              ;
    Real2DArray *vTau                ;
    DFTIntegratorDataBlockView viewP ;
    DFTIntegratorDataBlockView viewQ ;
    Real1DArray  sigmaPQ             ;
    Real1DArray  vSigmaPQ            ;
} DFTIntegratorDataBlock ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void                    DFTIntegratorDataBlock_Accumulate     ( DFTIntegratorDataBlock  *self ) ;
extern DFTIntegratorDataBlock *DFTIntegratorDataBlock_Allocate       ( const Integer numberOfFunctionals ,
                                                                       const Integer numberOfPoints      ,
                                                                       const Boolean hasSigma            ,
                                                                       const Boolean hasLaplacian        ,
                                                                       const Boolean hasTau              ,
                                                                       const Boolean isSpinRestricted    ,
                                                                             Status  *status             ) ;
extern void                    DFTIntegratorDataBlock_Deallocate     ( DFTIntegratorDataBlock **self ) ;
extern void                    DFTIntegratorDataBlock_Initialize     ( DFTIntegratorDataBlock  *self ) ;
extern void                    DFTIntegratorDataBlock_InitializeView ( DFTIntegratorDataBlock  *self, const Integer c, DFTIntegratorDataBlockView *view ) ;

# endif
