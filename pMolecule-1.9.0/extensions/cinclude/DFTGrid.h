/*------------------------------------------------------------------------------
! . File      : DFTGrid.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _DFTGRID
# define _DFTGRID

# include "Coordinates3.h"
# include "Definitions.h"
# include "DFTGridWeights.h"
# include "GridFunctionDataBlock.h"
# include "List.h"
# include "QCAtomContainer.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The number of accuracy classes. */
# define NDFTGRID_ACCURACY 5

/* . DFT grid accuracy definitions. */
typedef enum {
    DFTGridAccuracy_VeryLow  = 0,
    DFTGridAccuracy_Low      = 1,
    DFTGridAccuracy_Medium   = 2,
    DFTGridAccuracy_High     = 3,
    DFTGridAccuracy_VeryHigh = 4
} DFTGridAccuracy ;

/* . The type for storing grid points. */
typedef struct {
    Integer       atom           ;
    Integer       numberOfPoints ;
    Coordinates3 *coordinates3   ;
    Real1DArray  *weights        ;
    GridFunctionDataBlock *functionData ;
} DFTGridPointBlock ;

/* . The grid type. */
typedef struct {
    DFTGridAccuracy     accuracy        ;
    Integer             blockSize       ;
    Integer             numberOfPoints  ;
    Integer             numberOfRecords ;
    Real                bfTolerance     ;
    Real                rhoTolerance    ;
    DFTGridPointBlock **records         ;
    DFTGridWeights     *weights         ;
    List               *points          ;
} DFTGrid ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern DFTGrid           *DFTGrid_Allocate               ( const DFTGridAccuracy accuracy ) ;
extern DFTGrid           *DFTGrid_Construct              ( const DFTGridAccuracy accuracy, const QCAtomContainer *qcAtoms, const Coordinates3 *qcCoordinates3 ) ;
extern void               DFTGrid_Deallocate             (       DFTGrid **self ) ;
extern void               DFTGrid_DeallocateFunctionData (       DFTGrid  *self ) ;
extern Real               DFTGrid_FunctionDataSize       (       DFTGrid  *self ) ;
extern Boolean            DFTGrid_HasFunctionData        (       DFTGrid  *self ) ; 
extern DFTGridPointBlock *DFTGrid_Iterate                (       DFTGrid  *self ) ; 
extern void               DFTGrid_MakeRecords            (       DFTGrid  *self ) ;
extern Integer            DFTGrid_NumberOfFunctionValues (       DFTGrid  *self ) ;
extern Integer            DFTGrid_NumberOfRecords        (       DFTGrid  *self ) ; 

# endif
