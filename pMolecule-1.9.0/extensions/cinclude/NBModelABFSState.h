/*------------------------------------------------------------------------------
! . File      : NBModelABFSState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _NBMODELABFSSTATE
# define _NBMODELABFSSTATE

# include "Coordinates3.h"
# include "Definitions.h"
# include "ImageList.h"
# include "LJParameterContainer.h"
# include "MMAtomContainer.h"
# include "PairList.h"
# include "PairListGenerator.h"
# include "QCAtomContainer.h"
# include "QCMMLinkAtomCouplingOptions.h"
# include "Real1DArray.h"
# include "RegularGrid.h"
# include "RegularGridOccupancy.h"
# include "Selection.h"
# include "SelectionContainer.h"
# include "SymmetricMatrix.h"
# include "SymmetryParameters.h"
# include "SymmetryParameterGradients.h"
# include "Transformation3Container.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The NB model state type. */
typedef struct {
    /*. Options. */
    Boolean                     isNew                       ;
    Boolean                     useCentering                ;
    Boolean                     useGridSearch               ;
    Integer                     numberOfCalls               ;
# ifdef USEOPENMP
    Integer                     numberOfThreads             ;
# endif
    Integer                     numberOfUpdates             ;
    QCMMLinkAtomCoupling        qcmmCoupling                ;
    /* . Cutoff values. */
    Real                        listCutoff                  ;
    Real                        outerCutoff                 ;
    /* . Energies. */
    Real                        eimmmel                     ;
    Real                        eimmmlj                     ;
    Real                        eimqcmmlj                   ;
    Real                        eimqcqclj                   ;
    Real                        emmel                       ;
    Real                        emmel14                     ;
    Real                        emmlj                       ;
    Real                        emmlj14                     ;
    Real                        eqcmmlj                     ;
    Real                        eqcmmlj14                   ;
    /* . Statistics. */
    Real                        numberOfMMMMPairs           ;
    Real                        numberOfQCMMElPairs         ;
    Real                        numberOfQCMMLJPairs         ;
    Real                        numberOfMMMMImageImages     ;
    Real                        numberOfMMMMImagePairs      ;
    Real                        numberOfQCMMElImageImages   ;
    Real                        numberOfQCMMElImagePairs    ;
    Real                        numberOfQCMMLJImageImages   ;
    Real                        numberOfQCMMLJImagePairs    ;
    Real                        numberOfQCQCElImageImages   ;
    Real                        numberOfQCQCElImagePairs    ;
    Real                        numberOfQCQCLJImageImages   ;
    Real                        numberOfQCQCLJImagePairs    ;
    /* . Allocated arrays. */
    Coordinates3               *mmCoordinates3              ;
    Coordinates3               *qcCoordinates3              ;
    Coordinates3               *isolateTranslations3        ;
    Coordinates3               *referenceCoordinates3       ;
# ifdef USEOPENMP
    Coordinates3              **gradientsArray              ;
    Coordinates3              **threadGradients             ;
# endif
    PairList                   *qcmmExclusions              ;
    Real1DArray                *mmCharges                   ;
    Selection                  *freeSelection               ;
    Selection                  *mmSelection                 ;
    Selection                  *qcbSelection                ;
    Selection                  *qcpSelection                ;
    SelectionContainer         *isolates                    ;
    SymmetryParameters         *referenceSymmetryParameters ;
    /* . Lists. */
    ImageList                  *inbmmmm                     ;
    ImageList                  *inbqcmmel                   ;
    ImageList                  *inbqcmmlj                   ;
    ImageList                  *inbqcqcel                   ;
    ImageList                  *inbqcqclj                   ;
    PairList                   *nbmmmm                      ;
    PairList                   *nbmmmm14                    ;
    PairList                   *nbqcmmel                    ;
    PairList                   *nbqcmmel14                  ;
    PairList                   *nbqcmmlj                    ;
    PairList                   *nbqcmmlj14                  ;
    /* . Aliases. */
    Coordinates3               *coordinates3                ;
    PairList                   *exclusions                  ;
    Coordinates3               *gradients3                  ;
    Coordinates3               *inputCoordinates3           ;
    Coordinates3               *mmGradients3                ;
    PairList                   *interactions14              ;
    LJParameterContainer       *ljParameters                ;
    LJParameterContainer       *ljParameters14              ;
    MMAtomContainer            *mmAtoms                     ;
    QCAtomContainer            *qcAtoms                     ;
    Real1DArray                *qcCharges                   ;
    Real1DArray                *qcmmPotentials              ;
    Selection                  *fixedAtoms                  ;
    SymmetricMatrix            *qcqcPotentials              ;
    SymmetryParameters         *symmetryParameters          ;
    SymmetryParameterGradients *symmetryParameterGradients  ;
    Transformation3Container   *transformations             ;
    /* . Grid search. */
    RegularGrid                *mmGrid                      ;
    RegularGrid                *qcGrid                      ;
    RegularGridOccupancy       *mmOccupancy                 ;
    RegularGridOccupancy       *qcOccupancy                 ;
} NBModelABFSState ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifdef USEOPENMP
extern void              NBModelABFSState_AccumulateGradients    (       NBModelABFSState           *self                       ) ;                             
# endif
extern NBModelABFSState *NBModelABFSState_Allocate               ( const Integer                     n                          ) ;                                     
extern void              NBModelABFSState_Deallocate             (       NBModelABFSState          **self                       ) ;                    
extern void              NBModelABFSState_GridFinalize           (       NBModelABFSState           *self                       ) ;                    
extern void              NBModelABFSState_GridInitialize         (       NBModelABFSState           *self                       ,
                                                                   const PairListGenerator          *generator                  ,
                                                                         Status                     *status                     ) ; 
extern void              NBModelABFSState_Initialize             (       NBModelABFSState           *self                       ,                
                                                                         Coordinates3               *coordinates3               ,                
                                                                         SymmetryParameters         *symmetryParameters         ,                
                                                                         Coordinates3               *gradients3                 ,                
                                                                         SymmetryParameterGradients *symmetryParameterGradients ) ;              
extern void              NBModelABFSState_InitializeCoordinates3 (       NBModelABFSState           *self                       ,
                                                                   const Boolean                     doUpdate                   ) ;
extern NBModelABFSState *NBModelABFSState_SetUp                  (       MMAtomContainer            *mmAtoms                    ,                  
                                                                         QCAtomContainer            *qcAtoms                    ,                
                                                                         Selection                  *fixedAtoms                 ,                
                                                                         PairList                   *exclusions                 ,                
                                                                         PairList                   *interactions14             ,                
                                                                         LJParameterContainer       *ljParameters               ,                
                                                                         LJParameterContainer       *ljParameters14             ,                
                                                                         Real1DArray                *qcCharges                  ,                
                                                                         Real1DArray                *qcmmPotentials             ,                
                                                                         SymmetricMatrix            *qcqcPotentials             ,                
                                                                         Transformation3Container   *transformations            ,                
                                                                   const QCMMLinkAtomCoupling        qcmmCoupling               ) ;              
extern void              NBModelABFSState_SetUpCentering         (       NBModelABFSState           *self                       ,
                                                                   const Boolean                     useCentering               ,
                                                                         Status                     *status                     ) ;
extern void              NBModelABFSState_StatisticsAccumulate   (       NBModelABFSState           *self                       ) ;
extern void              NBModelABFSState_StatisticsInitialize   (       NBModelABFSState           *self                       ) ;

# endif
