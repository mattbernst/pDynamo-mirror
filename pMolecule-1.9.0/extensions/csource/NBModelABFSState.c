/*------------------------------------------------------------------------------
! . File      : NBModelABFSState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module implements ABFS non-bonding interaction state procedures.
!=================================================================================================================================*/

# include "ExecutionEnvironment.h"
# include "Memory.h"
# include "NBModelABFSState.h"

# ifdef USEOPENMP
/*----------------------------------------------------------------------------------------------------------------------------------
! . Accumulation - could certainly be better here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFSState_AccumulateGradients ( NBModelABFSState *self )
{
    if ( ( self->numberOfThreads > 1 ) && ( self->gradients3 != NULL ) )
    {
        #pragma omp parallel
        {
	    auto Integer r, t ;
            auto Real    gX, gY, gZ ;
            #pragma omp for schedule(dynamic)
            for ( r = 0 ; r < Coordinates3_Rows ( self->gradients3 ) ; r++ )
            {
                for ( t = 1 ; t < self->numberOfThreads ; t++ )
                {
                    Coordinates3_GetRow       ( self->threadGradients[t] , r, gX, gY, gZ ) ;
                    Coordinates3_IncrementRow ( self->gradients3         , r, gX, gY, gZ ) ;
                }
            }
        }
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
NBModelABFSState *NBModelABFSState_Allocate ( const Integer n )
{
    NBModelABFSState *self = NULL ;
    if ( n > 0 )
    {
        self = ( NBModelABFSState * ) Memory_Allocate ( sizeof ( NBModelABFSState ) ) ;
        if ( self != NULL )
        {
            /* . Options. */
            self->isNew                       = True  ; /* . The state is new. */
            self->useCentering                = False ;
            self->useGridSearch               = False ;
            self->qcmmCoupling                = QCMMLinkAtomCoupling_MM ;
            /* . Cutoff values. */
            self->listCutoff                  = 0.0e+00 ;
            self->outerCutoff                 = 0.0e+00 ;
            /* . Energies. */
            self->eimmmel                     = 0.0e+00 ;
            self->eimmmlj                     = 0.0e+00 ;
            self->eimqcmmlj                   = 0.0e+00 ;
            self->eimqcqclj                   = 0.0e+00 ;
            self->emmel                       = 0.0e+00 ;
            self->emmel14                     = 0.0e+00 ;
            self->emmlj                       = 0.0e+00 ;
            self->emmlj14                     = 0.0e+00 ;
            self->eqcmmlj                     = 0.0e+00 ;
            self->eqcmmlj14                   = 0.0e+00 ;
            /* . Statistics. */
            NBModelABFSState_StatisticsInitialize ( self ) ;
            /* . Lists. */
            self->inbmmmm                     = NULL ; 
            self->inbqcmmel                   = NULL ; 
            self->inbqcmmlj                   = NULL ; 
            self->inbqcqcel                   = NULL ; 
            self->inbqcqclj                   = NULL ; 
            self->nbmmmm                      = NULL ; 
            self->nbmmmm14                    = NULL ; 
            self->nbqcmmel                    = NULL ; 
            self->nbqcmmel14                  = NULL ; 
            self->nbqcmmlj                    = NULL ; 
            self->nbqcmmlj14                  = NULL ; 
            /* . Pointers. */
            self->coordinates3                = NULL ;
            self->exclusions                  = NULL ;
            self->fixedAtoms                  = NULL ;
            self->freeSelection               = NULL ;
            self->gradients3                  = NULL ;
            self->inputCoordinates3           = NULL ;
            self->interactions14              = NULL ;
            self->isolates                    = NULL ;
            self->isolateTranslations3        = NULL ;
            self->ljParameters                = NULL ;
            self->ljParameters14              = NULL ;
            self->mmAtoms                     = NULL ;
            self->mmSelection                 = NULL ;
            self->qcbSelection                = NULL ;
            self->qcAtoms                     = NULL ;
            self->qcCharges                   = NULL ;
            self->qcCoordinates3              = NULL ;
            self->qcmmExclusions              = NULL ;
            self->qcmmPotentials              = NULL ;
            self->qcpSelection                = NULL ;
            self->qcqcPotentials              = NULL ;
            self->symmetryParameters          = NULL ;
            self->symmetryParameterGradients  = NULL ;
            self->transformations             = NULL ;
            /* . QC/MM arrays. */
            self->mmCharges                   = NULL ;
            self->mmCoordinates3              = NULL ;
            self->mmGradients3                = NULL ;
            /* . Arrays to be allocated. */
            self->referenceCoordinates3       = Coordinates3_Allocate ( n )      ;
            self->referenceSymmetryParameters = SymmetryParameters_Allocate (  ) ;
            /* . Grid search. */
            self->mmGrid                      = NULL ;
            self->qcGrid                      = NULL ;
            self->mmOccupancy                 = NULL ;
            self->qcOccupancy                 = NULL ;
# ifdef USEOPENMP
            self->numberOfThreads             =   -1 ;
            self->gradientsArray              = NULL ;
            self->threadGradients             = NULL ;
# endif
            /* . Deallocate if there is not enough memory. */
            if ( ( self->referenceCoordinates3 == NULL ) || ( self->referenceSymmetryParameters == NULL ) ) NBModelABFSState_Deallocate ( &self ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFSState_Deallocate ( NBModelABFSState **self )
{
    if ( (*self) != NULL )
    {
        /* . Objects. */
        Coordinates3_Deallocate       ( &((*self)->qcCoordinates3             ) ) ;
        Coordinates3_Deallocate       ( &((*self)->referenceCoordinates3      ) ) ;
        PairList_Deallocate           ( &((*self)->qcmmExclusions             ) ) ;
        Selection_Deallocate          ( &((*self)->freeSelection              ) ) ;
        Selection_Deallocate          ( &((*self)->mmSelection                ) ) ;
        Selection_Deallocate          ( &((*self)->qcbSelection               ) ) ;
        Selection_Deallocate          ( &((*self)->qcpSelection               ) ) ;
        SymmetryParameters_Deallocate ( &((*self)->referenceSymmetryParameters) ) ;
        /* . Centering. */
        if ( (*self)->useCentering )
        {
            Coordinates3_Deallocate       ( &((*self)->isolateTranslations3) ) ;
            Coordinates3_Deallocate       ( &((*self)->coordinates3        ) ) ;
            SelectionContainer_Deallocate ( &((*self)->isolates            ) ) ;
        }
        /* . QC/MM data. */
        Real1DArray_Deallocate        ( &((*self)->mmCharges                  ) ) ;
        if ( (*self)->qcmmCoupling != QCMMLinkAtomCoupling_MM )
        {
            Coordinates3_Deallocate   ( &((*self)->mmCoordinates3             ) ) ;
            Coordinates3_Deallocate   ( &((*self)->mmGradients3               ) ) ;
        }
        /* . Lists. */
        ImageList_Deallocate ( &((*self)->inbmmmm   ) ) ;
        ImageList_Deallocate ( &((*self)->inbqcmmel ) ) ;
        ImageList_Deallocate ( &((*self)->inbqcmmlj ) ) ;
        ImageList_Deallocate ( &((*self)->inbqcqcel ) ) ;
        ImageList_Deallocate ( &((*self)->inbqcqclj ) ) ;
        PairList_Deallocate  ( &((*self)->nbmmmm    ) ) ;
        PairList_Deallocate  ( &((*self)->nbmmmm14  ) ) ;
        PairList_Deallocate  ( &((*self)->nbqcmmel  ) ) ;
        PairList_Deallocate  ( &((*self)->nbqcmmel14) ) ;
        PairList_Deallocate  ( &((*self)->nbqcmmlj  ) ) ;
        PairList_Deallocate  ( &((*self)->nbqcmmlj14) ) ;
# ifdef USEOPENMP
        /* . Thread data. */
        if ( ( (*self)->numberOfThreads > 0 ) && ( (*self)->threadGradients != NULL ) )
        {
            auto Integer t ;
            for ( t = 1 ; t < (*self)->numberOfThreads ; t++ ) Coordinates3_Deallocate ( &((*self)->threadGradients[t]) ) ;
            MEMORY_DEALLOCATE ( (*self)->threadGradients ) ;
        }
# endif
        /* . Object itself. */
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Finalization procedure for grid updating.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFSState_GridFinalize ( NBModelABFSState *self )
{
    if ( ( self != NULL ) && self->useGridSearch )
    {
        RegularGrid_Deallocate          ( &(self->mmGrid     ) );
        RegularGrid_Deallocate          ( &(self->qcGrid     ) );
        RegularGridOccupancy_Deallocate ( &(self->mmOccupancy) );
        RegularGridOccupancy_Deallocate ( &(self->qcOccupancy) );
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization procedure for grid updating.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFSState_GridInitialize ( NBModelABFSState *self, const PairListGenerator *generator, Status *status )
{
    if ( self != NULL )
    {
        /* . Determine whether to use grid search updating if the state is new. */
        if ( self->isNew ) self->useGridSearch = PairListGenerator_DetermineMethod ( generator, self->coordinates3, NULL ) ;
        /* . Calculate the grids if necessary. */
        if ( self->useGridSearch )
        {
            if ( self->mmGrid == NULL )
            {
                Coordinates3_MakeGridAndOccupancy ( self->coordinates3, NULL, generator->cellSize, &(self->mmGrid), &(self->mmOccupancy), status ) ;
            }
            if ( ( self->qcAtoms != NULL ) && ( QCAtomContainer_Size ( self->qcAtoms ) > 1 ) && ( self->qcGrid == NULL ) )
            {
                Coordinates3_MakeGridAndOccupancy ( self->qcCoordinates3, NULL, generator->cellSize, &(self->qcGrid), &(self->qcOccupancy), status ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization for an energy/gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFSState_Initialize ( NBModelABFSState           *self                       ,
                                   Coordinates3               *coordinates3               ,
                                   SymmetryParameters         *symmetryParameters         ,
                                   Coordinates3               *gradients3                 ,
                                   SymmetryParameterGradients *symmetryParameterGradients )
{
    if ( self != NULL )
    {
        /* . Statistics. */
        self->numberOfCalls += 1 ;
        /* . Energies. */
        self->eimmmel   = 0.0e+00 ;
        self->eimmmlj   = 0.0e+00 ;
        self->eimqcmmlj = 0.0e+00 ;
        self->eimqcqclj = 0.0e+00 ;
        self->emmel     = 0.0e+00 ;
        self->emmel14   = 0.0e+00 ;
        self->emmlj     = 0.0e+00 ;
        self->emmlj14   = 0.0e+00 ;
        self->eqcmmlj   = 0.0e+00 ;
        self->eqcmmlj14 = 0.0e+00 ;
        /* . Coordinates, gradients and symmetry parameters. */
        self->gradients3                 = gradients3                 ;
        self->inputCoordinates3          = coordinates3               ;
        self->symmetryParameters         = symmetryParameters         ;
        self->symmetryParameterGradients = symmetryParameterGradients ;
        /* . QC/MM quantities. */
        Real1DArray_Set     ( self->qcCharges,      0.0e+00 ) ;
        Real1DArray_Set     ( self->qcmmPotentials, 0.0e+00 ) ;
        SymmetricMatrix_Set ( self->qcqcPotentials, 0.0e+00 ) ;
# ifdef USEOPENMP
        /* . Thread data. */
        if ( gradients3 == NULL ) self->gradientsArray = NULL ;
        else
        {
            auto Integer t ;
            self->threadGradients[0] = gradients3 ;
            for ( t = 1 ; t < self->numberOfThreads ; t++ ) Coordinates3_Set ( self->threadGradients[t], 0.0e+00 ) ;
            self->gradientsArray = self->threadGradients ;
        }
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization of coordinate-dependent quantities.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFSState_InitializeCoordinates3 ( NBModelABFSState *self, const Boolean doUpdate )
{
    if ( self != NULL )
    {
        /* . Set the coordinates. */
        /* . Use centering and update the isolate translations if necessary. */
        if ( self->useCentering )
        {
            Coordinates3_CopyTo ( self->inputCoordinates3, self->coordinates3, NULL ) ;
            if ( doUpdate )
            {
                SymmetryParameters_CenterCoordinates3ByIsolate ( self->symmetryParameters, self->isolates, NULL, self->coordinates3 ) ;
                Coordinates3_CopyTo ( self->coordinates3, self->isolateTranslations3, NULL ) ;
                Coordinates3_AddScaledArray ( self->isolateTranslations3, -1.0e+00, self->inputCoordinates3, NULL ) ;
            }
            else { Coordinates3_AddScaledArray ( self->coordinates3, 1.0e+00, self->isolateTranslations3, NULL ) ; }
        }
        /* . Do not use centering. */
        else { self->coordinates3 = self->inputCoordinates3 ; }
        /* . QC/MM quantities. */
        if ( self->qcmmCoupling == QCMMLinkAtomCoupling_MM )
        {
            self->mmCoordinates3 = self->coordinates3 ;
            self->mmGradients3   = self->gradients3   ;
        }
        else
        {
            QCAtomContainer_GetMMCoordinates3 ( self->qcAtoms, self->coordinates3, self->mmCoordinates3 ) ;
            Coordinates3_Set ( self->mmGradients3, 0.0e+00 ) ;
        }
        /* . QC coordinates. */
        QCAtomContainer_GetCoordinates3 ( self->qcAtoms, self->coordinates3, False, &(self->qcCoordinates3) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setup the state.
!---------------------------------------------------------------------------------------------------------------------------------*/
NBModelABFSState *NBModelABFSState_SetUp (       MMAtomContainer          *mmAtoms         ,
                                                 QCAtomContainer          *qcAtoms         ,
                                                 Selection                *fixedAtoms      ,
                                                 PairList                 *exclusions      ,
                                                 PairList                 *interactions14  ,
                                                 LJParameterContainer     *ljParameters    ,
                                                 LJParameterContainer     *ljParameters14  ,
                                                 Real1DArray              *qcCharges       ,
                                                 Real1DArray              *qcmmPotentials  ,
                                                 SymmetricMatrix          *qcqcPotentials  ,
                                                 Transformation3Container *transformations ,
                                           const QCMMLinkAtomCoupling      qcmmCoupling    )
{
    NBModelABFSState *self = NULL ;
    if ( mmAtoms != NULL )
    {
        auto Boolean QOK    ;
        auto Integer n      ;
        auto Status  status ;
        n    = mmAtoms->natoms ;
        self = NBModelABFSState_Allocate ( n ) ;
        if ( self != NULL )
        {
            /* . Make sure that the pairlists are in the correct format. */
            status = SelfPairList_MakeConnections ( exclusions,     n ) ; QOK =        ( status != Status_OutOfMemory ) ;
            status = SelfPairList_MakeConnections ( interactions14, n ) ; QOK = QOK && ( status != Status_OutOfMemory ) ;

            /* . Selections. */
            /* . Fixed/free atoms. */
            if ( fixedAtoms != NULL ) self->freeSelection = Selection_Complement ( fixedAtoms, n ) ;

            /* . MM and QC atom selections. */
            /* . The MM selection can be NULL in which case all atoms are MM. */
            if ( qcAtoms != NULL )
            {
                /* . Selections. */
                self->qcbSelection = QCAtomContainer_MakeFullSelection ( qcAtoms ) ;
                self->qcpSelection = QCAtomContainer_MakePureSelection ( qcAtoms ) ;
                self->mmSelection  = Selection_Complement ( self->qcpSelection, n ) ;
                /* . Other data. */
                self->qcCoordinates3 = Coordinates3_Allocate ( qcAtoms->natoms ) ;
                self->qcmmExclusions = SelfPairList_ToCrossPairList ( exclusions, self->qcbSelection, self->mmSelection, NULL, True, True, True, False ) ;
                status = CrossPairList_MakeConnections ( self->qcmmExclusions, qcAtoms->natoms ) ;
                /* . Check memory. */
                QOK = QOK && ( self->mmSelection    != NULL ) &&
                             ( self->qcbSelection   != NULL ) &&
                             ( self->qcCoordinates3 != NULL ) &&
                             ( self->qcpSelection   != NULL ) &&
                             ( ( self->qcmmExclusions != NULL ) || ( exclusions == NULL ) ) && ( status != Status_OutOfMemory ) ;
                /* . MM coupling. */
                self->qcmmCoupling = qcmmCoupling ;
                if ( qcmmCoupling == QCMMLinkAtomCoupling_MM )
                {
                    self->mmCharges = MMAtomContainer_GetCharges ( mmAtoms, True, NULL ) ;
                    QOK = QOK && ( self->mmCharges != NULL ) ;
                }
                /* . RC/RD coupling. */
                else
                {
                    self->mmCharges = QCAtomContainer_GetMMCharges ( qcAtoms, mmAtoms, ( qcmmCoupling == QCMMLinkAtomCoupling_RC ) ) ;
                    if ( self->mmCharges != NULL )
                    {
                        self->mmCoordinates3 = Coordinates3_Allocate ( self->mmCharges->length ) ;
                        self->mmGradients3   = Coordinates3_Allocate ( self->mmCharges->length ) ;
                    }
                    QOK = QOK && ( self->mmCharges != NULL ) && ( self->mmCoordinates3 != NULL ) && ( self->mmGradients3 != NULL ) ;
                }
            }
            /* . Assign aliases. */
            self->exclusions      = exclusions      ;
            self->fixedAtoms      = fixedAtoms      ;
            self->interactions14  = interactions14  ;
            self->ljParameters    = ljParameters    ;
            self->ljParameters14  = ljParameters14  ;
            self->mmAtoms         = mmAtoms         ;
            self->qcAtoms         = qcAtoms         ;
            self->qcCharges       = qcCharges       ;
            self->qcmmPotentials  = qcmmPotentials  ;
            self->qcqcPotentials  = qcqcPotentials  ;
            self->transformations = transformations ;
# ifdef USEOPENMP
            /* . Thread data. */
            {
                auto Integer t ;
                self->gradientsArray  = NULL ;
                self->numberOfThreads = omp_get_max_threads ( ) ;
                MEMORY_ALLOCATEARRAY ( self->threadGradients, self->numberOfThreads, Coordinates3 * ) ;
                if ( self->threadGradients == NULL ) QOK = False ;
                else
                {
                    self->threadGradients[0] = NULL ;
                    for ( t = 1 ; t < self->numberOfThreads ; t++ )
                    {
                        self->threadGradients[t] = Coordinates3_Allocate ( n ) ;
                        QOK = QOK && ( self->threadGradients[t] != NULL ) ;
                    }
                }
            }
# endif
            /* . Check that everything is OK. */
            if ( ! QOK ) NBModelABFSState_Deallocate ( &self ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setup centering for systems with symmetry.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFSState_SetUpCentering ( NBModelABFSState *self, const Boolean useCentering, Status *status )
{
    if ( ( self != NULL ) && ( self->exclusions != NULL ) && ( self->transformations != NULL ) && useCentering )
    {
        Integer n = self->mmAtoms->natoms ;
        SelectionContainer *nbIsolates = NULL ;
        /* . Create isolates. */
        nbIsolates = SelfPairList_ToIsolateSelectionContainer ( self->exclusions, n ) ; /* . MM isolates. */
        SelectionContainer_MergeIsolates  ( nbIsolates, self->qcbSelection ) ;          /* . All QC atoms must stay together. */
        SelectionContainer_RemoveIsolates ( nbIsolates, self->fixedAtoms   ) ;          /* . Fixed atoms must not move. */
        /* . Setup centering. */
        if ( nbIsolates == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ; 
        else if ( nbIsolates->nitems > 1 )
        {
            self->coordinates3         = Coordinates3_Allocate ( n ) ;
            self->isolateTranslations3 = Coordinates3_Allocate ( n ) ;
            if ( ( self->coordinates3 == NULL ) || ( self->isolateTranslations3 == NULL ) ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
            else
            {
                self->isolates     = nbIsolates ;
                self->useCentering = True       ;
            }
        }
        else { SelectionContainer_Deallocate ( &nbIsolates ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Accumulate statistics counters after an update.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFSState_StatisticsAccumulate ( NBModelABFSState *self )
{
    if ( self != NULL )
    {
        self->numberOfUpdates += 1 ;
        if ( self->nbmmmm   != NULL ) self->numberOfMMMMPairs   += self->nbmmmm->npairs   ;
        if ( self->nbqcmmel != NULL ) self->numberOfQCMMElPairs += self->nbqcmmel->npairs ;
        if ( self->nbqcmmlj != NULL ) self->numberOfQCMMLJPairs += self->nbqcmmlj->npairs ;
        if ( self->inbmmmm  != NULL )
        {
            self->numberOfMMMMImageImages   += ImageList_NumberOfImages ( self->inbmmmm   ) ;
            self->numberOfMMMMImagePairs    += ImageList_NumberOfPairs  ( self->inbmmmm   ) ;
        }
        if ( self->inbqcmmel != NULL )
        {
            self->numberOfQCMMElImageImages += ImageList_NumberOfImages ( self->inbqcmmel ) ;
            self->numberOfQCMMElImagePairs  += ImageList_NumberOfPairs  ( self->inbqcmmel ) ;
        }
        if ( self->inbqcmmlj != NULL )
        {
            self->numberOfQCMMLJImageImages += ImageList_NumberOfImages ( self->inbqcmmlj ) ;
            self->numberOfQCMMLJImagePairs  += ImageList_NumberOfPairs  ( self->inbqcmmlj ) ;
        }
        if ( self->inbqcqcel != NULL )
        {
            self->numberOfQCQCElImageImages += ImageList_NumberOfImages ( self->inbqcqcel ) ;
            self->numberOfQCQCElImagePairs  += ImageList_NumberOfPairs  ( self->inbqcqcel ) ;
        }
        if ( self->inbqcqclj != NULL )
        {
            self->numberOfQCQCLJImageImages += ImageList_NumberOfImages ( self->inbqcqclj ) ;
            self->numberOfQCQCLJImagePairs  += ImageList_NumberOfPairs  ( self->inbqcqclj ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialize statistics counters.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFSState_StatisticsInitialize ( NBModelABFSState *self )
{
    if ( self != NULL )
    {
        self->numberOfCalls              = 0   ;
        self->numberOfUpdates            = 0   ;
        self->numberOfMMMMPairs          = 0.0 ;
        self->numberOfQCMMElPairs        = 0.0 ;
        self->numberOfQCMMLJPairs        = 0.0 ;
        self->numberOfMMMMImageImages    = 0.0 ;
        self->numberOfMMMMImagePairs     = 0.0 ;
        self->numberOfQCMMElImageImages  = 0.0 ;
        self->numberOfQCMMElImagePairs   = 0.0 ;
        self->numberOfQCMMLJImageImages  = 0.0 ;
        self->numberOfQCMMLJImagePairs   = 0.0 ;
        self->numberOfQCQCElImageImages  = 0.0 ;
        self->numberOfQCQCElImagePairs   = 0.0 ;
        self->numberOfQCQCLJImageImages  = 0.0 ;
        self->numberOfQCQCLJImagePairs   = 0.0 ;
    }
}
