/*------------------------------------------------------------------------------
! . File      : DFTIntegrator.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module handles DFT integration.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "DFTIntegrator.h"
# include "DFTIntegratorDataBlock.h"
# include "ExecutionEnvironment.h"
# include "GridFunctionDataBlock.h"
# include "Memory.h"

# define DFTGRIDWEIGHTDERIVATIVES

/* . Arrays are defined as number of basis functions (columns) * number of grid points (rows). */

/*
! . Care is needed not to mix OpenMP and parallel Atlas - in the general case will never want to do this anyway.
! . Can add compilation flag to prevent this - i.e. don't use OpenMP where parallel BLAS is being used.
*/

/*
!
! . Formulae for rho-dependent terms:
!    Rp     = Sum_mn Bmp Bnp Pmn
!    Del Rp = Sum_mn ( Xmp Bnp + Bmp Xnp ) Pmn, etc.
!    Sp     = Del Rp . Del Rp
!    Lp     = Sum_mn ( XXmp Bnp + 2 Xmp Xnp + Bmp XXnp ) Pmn + ...
!    Tp     = 1/2 Sum_mn ( Xmp Xnp + Ymp Ynp + Zmp Znp ) Pmn
!
! . Formulae for derivatives (all v-terms weighted by Wp):
!    Rp      -     Sum_p Bmp Bnp vRp
!    Sp      - 2 * Sum_p ( Xmp Bnp + Bmp Xnp ) * dRhoXp         * vSp + ...
!    Sp (ab) -     Sum_p ( Xmp Bnp + Bmp Xnp ) * dRhoXp (other) * vSp + ...
!    Lp      -     Sum_p ( XXmp Bnp + 2 Xmp Xnp + Bmp XXnp ) vLp + ...
!    Tp      - 1/2 Sum_p ( Xmp Xnp + Ymp Ynp + Zmp Znp ) vTp
!
! . Key:
!    B                      - basis function values
!    P                      - density matrix
!    X, Y, Z                - derivative values of basis functions
!    XX, XY, XZ, YY, YZ, ZZ - second derivative values of basis functions
!
!    m, n refer to basis functions, p to points.
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_Fock ( const Boolean                hasSigma     ,
                                 const Boolean                hasLaplacian ,
                                 const Boolean                hasTau       ,
                                 const GridFunctionDataBlock *basisData    ,
                                 const Real1DArray           *dRhoX        ,
                                 const Real1DArray           *dRhoY        ,
                                 const Real1DArray           *dRhoZ        ,
                                 const Real1DArray           *vRho         ,
                                 const Real1DArray           *vSigma       ,
                                 const Real1DArray           *vLaplacianRho,
                                 const Real1DArray           *vTau         ,
                                       SymmetricMatrix       *fock         ,
                                       Real2DArray           *work2D       ) ;

static void DFTIntegrator_FockSigma ( const Boolean                hasSigma  ,
                                      const Boolean                scale     ,
                                      const GridFunctionDataBlock *basisData ,
                                      const Real1DArray           *dRhoX     ,
                                      const Real1DArray           *dRhoY     ,
                                      const Real1DArray           *dRhoZ     ,
                                      const Real1DArray           *vSigma    ,
                                            SymmetricMatrix       *fock      ,
                                            Real2DArray           *work2D    ) ;

static void DFTIntegrator_FormReducedDensity ( const Integer1DArray   *indices        ,
                                               const SymmetricMatrix  *density        ,
                                                     Real2DArray     **reducedDensity ,
                                                     Status           *status         ) ;

static void DFTIntegrator_Gradients ( const Boolean                hasSigma      ,
                                      const Boolean                hasLaplacian  ,
                                      const Boolean                hasTau        ,
                                      const Integer1DArray        *atomIndices   ,
                                      const GridFunctionDataBlock *basisData     ,
                                      const Real1DArray           *dRhoX         ,
                                      const Real1DArray           *dRhoY         ,
                                      const Real1DArray           *dRhoZ         ,
                                      const Real1DArray           *vRho          ,
                                      const Real1DArray           *vSigma        ,
                                      const Real1DArray           *vLaplacianRho ,
                                      const Real1DArray           *vTau          ,
                                      const Real2DArray           *density       ,
                                      const Integer                gridAtom      ,
                                            Coordinates3          *gradients3    ,
                                            Real2DArray           *temp2D        ,
                                            Real2DArray           *work2D        ) ;

static void DFTIntegrator_GradientsSigma ( const Boolean                hasSigma    ,
                                           const Boolean                isSelf      ,
                                           const Integer1DArray        *atomIndices ,
                                           const GridFunctionDataBlock *basisData   ,
                                           const Real1DArray           *dRhoX       ,
                                           const Real1DArray           *dRhoY       ,
                                           const Real1DArray           *dRhoZ       ,
                                           const Real1DArray           *vSigma      ,
                                           const Real2DArray           *density     ,
                                           const Integer                gridAtom    ,
                                                 Coordinates3          *gradients3  ,
                                                 Real2DArray           *temp2D      ,
                                                 Real2DArray           *work2D      ) ;

static void DFTIntegrator_GridPointRho ( const Boolean                hasSigma     ,
                                         const Boolean                hasLaplacian ,
                                         const Boolean                hasTau       ,
                                         const GridFunctionDataBlock *basisData    ,
                                         const Real2DArray           *density      ,
                                               Real1DArray           *rho          ,
                                               Real1DArray           *dRhoX        ,
                                               Real1DArray           *dRhoY        ,
                                               Real1DArray           *dRhoZ        ,
                                               Real1DArray           *sigma        ,
                                               Real1DArray           *laplacianRho ,
                                               Real1DArray           *tau          ,
                                               Real1DArray           *work1D       ,
                                               Real2DArray           *work2D       ) ;

static void DFTIntegrator_GridPointSigma ( const Boolean      hasSigma ,
                                           const Real1DArray *dRhoXa   ,
                                           const Real1DArray *dRhoYa   ,
                                           const Real1DArray *dRhoZa   ,
                                           const Real1DArray *dRhoXb   ,
                                           const Real1DArray *dRhoYb   ,
                                           const Real1DArray *dRhoZb   ,
                                                 Real1DArray *sigma    ) ;

static void UGradientContributions ( const Integer1DArray *atomIndices ,
                                     const Integer1DArray *indices     ,
                                     const Real2DArray    *a           ,
                                     const Real2DArray    *x           ,
                                     const Real2DArray    *y           ,
                                     const Real2DArray    *z           ,
                                     const Integer         gridAtom    ,
                                           Coordinates3   *gradients3  ) ;

static void UReal2DArray_ColumnAddScaledArray ( const Real2DArray *a       ,
                                                const Real1DArray *weights ,
                                                      Real2DArray *b       ) ;

static void UReal2DArray_ColumnScale ( const Real1DArray *a ,
                                             Real2DArray *b ) ;

static void USymmetricMatrix_DotProductIncrement ( const Integer1DArray  *indices ,
                                                   const Real2DArray     *a       ,
                                                   const Real2DArray     *b       ,
                                                         SymmetricMatrix *c       ) ;

/*==================================================================================================================================
! . Integrate over a grid.
!=================================================================================================================================*/
extern void DFTIntegrator_Integrate ( const DFTFunctionalModel *functionalModel    ,
                                            DFTGrid            *grid               ,
                                      const QCAtomContainer    *qcAtoms            ,
                                      const QCParameter        *qcParameters       ,
                                            Coordinates3       *qcCoordinates3     ,
                                      const QCOnePDM           *densityP           ,
                                      const QCOnePDM           *densityQ           ,
                                      const Boolean             inCore             ,
                                      const Boolean             isSpinUnrestricted ,
                                            Real               *eQuad              ,
                                            Real               *rhoQuad            ,
# ifdef USEOPENMP
                                            SymmetricMatrix   **fockA              ,
                                            SymmetricMatrix   **fockB              ,
                                            Coordinates3      **gradients3         ,
# else
                                            SymmetricMatrix    *fockA              ,
                                            SymmetricMatrix    *fockB              ,
                                            Coordinates3       *gradients3         ,
# endif
                                            Status             *status             )
{
    if ( eQuad   != NULL ) (*eQuad  ) = 0.0e+00 ;
    if ( rhoQuad != NULL ) (*rhoQuad) = 0.0e+00 ;
    if ( ( functionalModel != NULL ) && ( grid != NULL ) && ( qcAtoms != NULL ) && ( qcParameters != NULL ) && ( qcCoordinates3 != NULL ) && ( densityP != NULL ) )
    {
        auto Boolean         determineFunctionData, doFock, doGradients, storeFunctionData ;
        auto Integer         order ;
        auto Real            eXCTotal = 0.0e+00, rhoTotal = 0.0e+00 ;
        auto Integer1DArray *atomIndices     = NULL ;

        /* . Initialization. */
        DFTGrid_MakeRecords ( grid ) ;
        doFock      = ( fockA      != NULL ) ;
        doGradients = ( gradients3 != NULL ) ;
        order       = functionalModel->order ;
        if ( doGradients )
        {
            order += 1 ;
            atomIndices = QCAtomContainer_OrbitalBasisAtomIndices ( qcAtoms, status ) ;
            DFTGrid_DeallocateFunctionData ( grid ) ;
        }
        determineFunctionData = ( ! inCore      ) || ( inCore && ( grid->records[0]->functionData == NULL ) ) ;
        storeFunctionData     = ( ! doGradients ) && ( inCore && ( grid->records[0]->functionData == NULL ) ) ;

# ifdef USEOPENMP
        # pragma omp parallel shared(status)
# endif
        {
            auto Integer gridAtom, r ;
            auto Status  localStatus = Status_Continue ;
            auto Coordinates3                  *coordinates3    = NULL, *localGradients3  = NULL ;
            auto DFTGridPointBlock             *block           = NULL ;
            auto DFTGridWeightsDerivativesWork *weightsWork     = NULL ;
            auto DFTIntegratorDataBlock        *rhoData         = NULL ;
            auto DFTIntegratorDataBlockView    *rhoDataP        = NULL , *rhoDataQ        = NULL ;
            auto GridFunctionDataBlock         *basisData       = NULL ;
            auto Real1DArray                   *weights         = NULL , *work1D          = NULL ;
            auto Real2DArray                   *reducedDensityP = NULL , *reducedDensityQ = NULL, *temp2D = NULL, *work2D = NULL ;
            auto SymmetricMatrix               *localFockA      = NULL , *localFockB      = NULL ;
# ifdef USEOPENMP
            auto Integer thisThread = omp_get_thread_num ( ) ;
# endif
            /* . Assign local Fock and gradients variables. */
            if ( doFock )
            {
# ifdef USEOPENMP
                localFockA = fockA[thisThread] ;
                if ( isSpinUnrestricted ) localFockB = fockB[thisThread] ;
# else
                localFockA = fockA ;
                if ( isSpinUnrestricted ) localFockB = fockB ;
# endif
            }
            if ( doGradients )
            {
# ifdef USEOPENMP
                localGradients3 = gradients3[thisThread] ;
# else
                localGradients3 = gradients3 ;
# endif
            }

            /* . Loop over the grid point blocks. */
# ifdef USEOPENMP
            #pragma omp for reduction(+:eXCTotal,rhoTotal) schedule(dynamic)
# endif
            for ( r = 0 ; r < grid->numberOfRecords ; r++ )
            {
                /* . Get the block and associated data. */
                block        = grid->records[r]    ;
                gridAtom     = block->atom         ;
                coordinates3 = block->coordinates3 ;
                weights      = block->weights      ;
                /* . Determine basis function values and their derivatives at the grid points. */
                if ( determineFunctionData )
                {
                    /* . Allocate a data block. */
                    if ( ( basisData == NULL ) || ( basisData->numberOfPoints != block->numberOfPoints ) )
                    {
                        GridFunctionDataBlock_Deallocate ( &basisData ) ;
                        basisData = GridFunctionDataBlock_Allocate ( qcAtoms->nobasisw, block->numberOfPoints, order, &localStatus ) ;
                    }
                    else GridFunctionDataBlock_Resize ( basisData, qcAtoms->nobasisw, &localStatus ) ;
                    QCAtomContainer_OrbitalBasisGridPointValues ( qcAtoms, qcParameters, qcCoordinates3, coordinates3, True, &(grid->bfTolerance), basisData, &localStatus ) ;
                }
                /* . Retrieve function data. */
                else basisData = block->functionData ;
                if ( ( basisData == NULL ) || ( localStatus != Status_Continue ) || ( basisData->numberOfFunctions <= 0 ) ) goto EndOfLoop ;

                /* . Ensure that there is an integration data block of the correct size. */
                if ( ( rhoData == NULL ) || ( rhoData->numberOfPoints != block->numberOfPoints ) )
                {
                    DFTIntegratorDataBlock_Deallocate ( &rhoData ) ;
                    rhoData = DFTIntegratorDataBlock_Allocate ( functionalModel->numberOfFunctionals ,
                                                                block->numberOfPoints                ,
                                                                functionalModel->hasSigma            ,
                                                                functionalModel->hasLaplacian        ,
                                                                functionalModel->hasTau              ,
                                                                functionalModel->isSpinRestricted    ,
                                                                &localStatus                         ) ;
                    if ( rhoData == NULL ) goto EndOfLoop ;
                    rhoDataP = &(rhoData->viewP) ;
                    rhoDataQ = &(rhoData->viewQ) ;
                }
                /* . Allocate scratch space. */
                if ( ( work2D == NULL ) || ( work2D->length0 != basisData->numberOfFunctions ) || ( work2D->length1 != rhoData->numberOfPoints ) )
                {
                    Real2DArray_Deallocate ( &work2D ) ;
                    work2D = Real2DArray_Allocate ( basisData->numberOfFunctions, rhoData->numberOfPoints, &localStatus ) ;
                    if ( work2D == NULL ) goto EndOfLoop ;
                }
                if ( ( functionalModel->hasLaplacian || functionalModel->hasTau ) && ( ( work1D == NULL ) || ( work1D->length != rhoData->numberOfPoints ) ) )
                {
                    Real1DArray_Deallocate ( &work1D ) ;
                    work1D = Real1DArray_Allocate ( rhoData->numberOfPoints, &localStatus ) ;
                    if ( work1D == NULL ) goto EndOfLoop ;
                }
                if ( doGradients && functionalModel->hasSigma && ( ( temp2D == NULL ) || ( temp2D->length0 != basisData->numberOfFunctions ) || ( temp2D->length1 != rhoData->numberOfPoints ) ) )
                {
                    Real2DArray_Deallocate ( &temp2D ) ;
                    temp2D = Real2DArray_Allocate ( basisData->numberOfFunctions, rhoData->numberOfPoints, &localStatus ) ;
                    if ( temp2D == NULL ) goto EndOfLoop ;
                }

                /* . Evaluate the densities and associated quantities at the grid points. */
                DFTIntegrator_FormReducedDensity ( basisData->indices, densityP->density, &reducedDensityP, &localStatus ) ;
                if ( reducedDensityP == NULL ) goto EndOfLoop ;
                DFTIntegrator_GridPointRho   ( functionalModel->hasSigma     ,
                                               functionalModel->hasLaplacian ,
                                               functionalModel->hasTau       ,
                                               basisData                     ,
                                               reducedDensityP               ,  
                                               &(rhoDataP->rho         )     ,  
                                               &(rhoDataP->dRhoX       )     ,
                                               &(rhoDataP->dRhoY       )     ,
                                               &(rhoDataP->dRhoZ       )     ,
                                               &(rhoDataP->sigma       )     ,
                                               &(rhoDataP->laplacianRho)     ,
                                               &(rhoDataP->tau         )     ,
                                               work1D                        ,  
                                               work2D                        ) ;
                if ( isSpinUnrestricted )
                {
                    DFTIntegrator_FormReducedDensity ( basisData->indices, densityQ->density, &reducedDensityQ, &localStatus ) ;
                    if ( reducedDensityQ == NULL ) goto EndOfLoop ;
                    DFTIntegrator_GridPointRho   ( functionalModel->hasSigma     ,
                                                   functionalModel->hasLaplacian ,
                                                   functionalModel->hasTau       ,
                                                   basisData                     ,
                                                   reducedDensityQ               ,  
                                                   &(rhoDataQ->rho         )     ,  
                                                   &(rhoDataQ->dRhoX       )     ,
                                                   &(rhoDataQ->dRhoY       )     ,
                                                   &(rhoDataQ->dRhoZ       )     ,
                                                   &(rhoDataQ->sigma       )     ,
                                                   &(rhoDataQ->laplacianRho)     ,
                                                   &(rhoDataQ->tau         )     ,
                                                   work1D                        ,  
                                                   work2D                        ) ;
                    DFTIntegrator_GridPointSigma ( functionalModel->hasSigma ,
                                                   &(rhoDataP->dRhoX )       ,
                                                   &(rhoDataP->dRhoY )       ,
                                                   &(rhoDataP->dRhoZ )       ,
                                                   &(rhoDataQ->dRhoX )       ,
                                                   &(rhoDataQ->dRhoY )       ,
                                                   &(rhoDataQ->dRhoZ )       ,
                                                   &(rhoData->sigmaPQ)       ) ;
                }

                /* . Skip the block if all densities are insignificant. */
                if ( Real2DArray_AbsoluteMaximum ( rhoData->rho ) <= grid->rhoTolerance )  goto EndOfLoop ;

                /* . Evaluate the functional terms. */
                DFTFunctionalModel_Evaluate ( functionalModel, rhoData ) ;

                /* . Accumulation and weighting of the integration data. */
                Real1DArray_Multiply ( &(rhoDataP->vRho), weights, NULL ) ;
                if ( functionalModel->hasLaplacian ) Real1DArray_Multiply ( &(rhoDataP->vLaplacianRho), weights, NULL ) ;
                if ( functionalModel->hasSigma     ) Real1DArray_Multiply ( &(rhoDataP->vSigma       ), weights, NULL ) ;
                if ( functionalModel->hasTau       ) Real1DArray_Multiply ( &(rhoDataP->vTau         ), weights, NULL ) ;
                if ( isSpinUnrestricted )
                {
                    Real1DArray_AddScaledArray ( &(rhoDataP->rho ), 1.0e+00, &(rhoDataQ->rho), NULL ) ;
                    Real1DArray_Multiply       ( &(rhoDataQ->vRho),          weights         , NULL ) ;
                    if ( functionalModel->hasLaplacian ) Real1DArray_Multiply ( &(rhoDataQ->vLaplacianRho), weights, NULL ) ;
                    if ( functionalModel->hasSigma     ) Real1DArray_Multiply ( &(rhoDataQ->vSigma       ), weights, NULL ) ;
                    if ( functionalModel->hasSigma     ) Real1DArray_Multiply ( &(rhoData->vSigmaPQ      ), weights, NULL ) ;
                    if ( functionalModel->hasTau       ) Real1DArray_Multiply ( &(rhoDataQ->vTau         ), weights, NULL ) ;
                }
                /* . Total energy and density - eXC is multiplied by rhoP which is now the total density. */
                Real1DArray_Multiply ( rhoData->eXC, &(rhoDataP->rho), NULL ) ;
                eXCTotal += Real1DArray_Dot ( rhoData->eXC    , weights, NULL ) ;
                rhoTotal += Real1DArray_Dot ( &(rhoDataP->rho), weights, NULL ) ;

                /* . Fock terms. */
                if ( doFock )
                {
                    DFTIntegrator_Fock ( functionalModel->hasSigma     ,
                                         functionalModel->hasLaplacian ,
                                         functionalModel->hasTau       ,
                                         basisData                     ,
                                         &(rhoDataP->dRhoX        )    ,
                                         &(rhoDataP->dRhoY        )    ,
                                         &(rhoDataP->dRhoZ        )    ,
                                         &(rhoDataP->vRho         )    ,
                                         &(rhoDataP->vSigma       )    ,
                                         &(rhoDataP->vLaplacianRho)    ,
                                         &(rhoDataP->vTau         )    ,
                                         localFockA                    ,
                                         work2D                        ) ;
                    if ( isSpinUnrestricted )
                    {
                        DFTIntegrator_Fock ( functionalModel->hasSigma     ,
                                             functionalModel->hasLaplacian ,
                                             functionalModel->hasTau       ,
                                             basisData                     ,
                                             &(rhoDataQ->dRhoX        )    ,
                                             &(rhoDataQ->dRhoY        )    ,
                                             &(rhoDataQ->dRhoZ        )    ,
                                             &(rhoDataQ->vRho         )    ,
                                             &(rhoDataQ->vSigma       )    ,
                                             &(rhoDataQ->vLaplacianRho)    ,
                                             &(rhoDataQ->vTau         )    ,
                                             localFockB                    ,
                                             work2D                        ) ;
                        DFTIntegrator_FockSigma ( functionalModel->hasSigma ,
                                                  False                     ,
                                                  basisData                 ,
                                                  &(rhoDataQ->dRhoX  )      ,
                                                  &(rhoDataQ->dRhoY  )      ,
                                                  &(rhoDataQ->dRhoZ  )      ,
                                                  &(rhoData->vSigmaPQ)      ,
                                                  localFockA                ,
                                                  work2D                    ) ;
                        DFTIntegrator_FockSigma ( functionalModel->hasSigma ,
                                                  False                     ,
                                                  basisData                 ,
                                                  &(rhoDataP->dRhoX  )      ,
                                                  &(rhoDataP->dRhoY  )      ,
                                                  &(rhoDataP->dRhoZ  )      ,
                                                  &(rhoData->vSigmaPQ)      ,
                                                  localFockB                ,
                                                  work2D                    ) ;
                    }
                }

                /* . Gradient terms. */
                if ( doGradients )
                {
                    /* . Direct terms. */
                    DFTIntegrator_Gradients ( functionalModel->hasSigma     ,
                                              functionalModel->hasLaplacian ,
                                              functionalModel->hasTau       ,
                                              atomIndices                   ,
                                              basisData                     ,
                                              &(rhoDataP->dRhoX        )    ,
                                              &(rhoDataP->dRhoY        )    ,
                                              &(rhoDataP->dRhoZ        )    ,
                                              &(rhoDataP->vRho         )    ,
                                              &(rhoDataP->vSigma       )    ,
                                              &(rhoDataP->vLaplacianRho)    ,
                                              &(rhoDataP->vTau         )    ,
                                              reducedDensityP               ,
                                              gridAtom                      ,
                                              localGradients3               ,
                                              temp2D                        ,
                                              work2D                        ) ;
                    if ( isSpinUnrestricted )
                    {   
                        DFTIntegrator_Gradients ( functionalModel->hasSigma     ,
                                                  functionalModel->hasLaplacian ,
                                                  functionalModel->hasTau       ,
                                                  atomIndices                   ,
                                                  basisData                     ,
                                                  &(rhoDataQ->dRhoX        )    ,
                                                  &(rhoDataQ->dRhoY        )    ,
                                                  &(rhoDataQ->dRhoZ        )    ,
                                                  &(rhoDataQ->vRho         )    ,
                                                  &(rhoDataQ->vSigma       )    ,
                                                  &(rhoDataQ->vLaplacianRho)    ,
                                                  &(rhoDataQ->vTau         )    ,
                                                  reducedDensityQ               ,
                                                  gridAtom                      ,
                                                  localGradients3               ,
                                                  temp2D                        ,
                                                  work2D                        ) ;
                        DFTIntegrator_GradientsSigma ( functionalModel->hasSigma ,
                                                       False                ,
                                                       atomIndices          ,
                                                       basisData            ,
                                                       &(rhoDataQ->dRhoX )  ,
                                                       &(rhoDataQ->dRhoY )  ,
                                                       &(rhoDataQ->dRhoZ )  ,
                                                       &(rhoData->vSigmaPQ) ,
                                                       reducedDensityP      ,
                                                       gridAtom             ,
                                                       localGradients3      ,
                                                       temp2D               ,
                                                       work2D               ) ;
                        DFTIntegrator_GradientsSigma ( functionalModel->hasSigma ,
                                                       False                ,
                                                       atomIndices          ,
                                                       basisData            ,
                                                       &(rhoDataP->dRhoX )  ,
                                                       &(rhoDataP->dRhoY )  ,
                                                       &(rhoDataP->dRhoZ )  ,
                                                       &(rhoData->vSigmaPQ) ,
                                                       reducedDensityQ      ,
                                                       gridAtom             ,
                                                       localGradients3      ,
                                                       temp2D               ,
                                                       work2D               ) ;
                    }

# ifdef DFTGRIDWEIGHTDERIVATIVES
                    /* . Weight terms. */
                    if ( weightsWork == NULL )
                    {
                        weightsWork = DFTGridWeightsDerivativesWork_Allocate ( grid->weights, &localStatus ) ;
                        if ( weightsWork == NULL ) goto EndOfLoop ;
                    }
                    DFTGridWeights_Derivatives ( grid->weights         ,
                                                 gridAtom              ,
                                                 block->numberOfPoints ,
                                                 coordinates3          ,
                                                 weights               ,
                                                 rhoData->eXC          ,
                                                 localGradients3       ,
                                                 weightsWork           ) ;
# endif
                }

                /* . End of loop. */
                EndOfLoop:

                /* . Store function data. */
                if ( storeFunctionData ) { block->functionData = basisData ; basisData = NULL ; }
            }

            /* . Deallocate space. */
            DFTGridWeightsDerivativesWork_Deallocate ( &weightsWork ) ;
            DFTIntegratorDataBlock_Deallocate        ( &rhoData     ) ;
            Real1DArray_Deallocate    ( &work1D          ) ;
            Real2DArray_Deallocate    ( &reducedDensityP ) ;
            Real2DArray_Deallocate    ( &reducedDensityQ ) ;
            Real2DArray_Deallocate    ( &temp2D          ) ;
            Real2DArray_Deallocate    ( &work2D          ) ;
            if ( doGradients || ( ! inCore ) ) GridFunctionDataBlock_Deallocate ( &basisData ) ;

            /* . Finish up. */
            if ( localStatus != Status_Continue )
            {
# ifdef USEOPENMP
		#pragma omp critical
# endif
                Status_Set ( status, localStatus ) ;
            }
        }

        /* . Deallocate space. */
        Integer1DArray_Deallocate ( &atomIndices ) ;

        /* . Finish up. */
        if ( eQuad   != NULL ) (*eQuad  ) = eXCTotal ;
        if ( rhoQuad != NULL ) (*rhoQuad) = rhoTotal ;
    }
}

/*==================================================================================================================================
! . Local procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Contributions to a Fock matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_Fock ( const Boolean                hasSigma     ,
                                 const Boolean                hasLaplacian ,
                                 const Boolean                hasTau       ,
                                 const GridFunctionDataBlock *basisData    ,
                                 const Real1DArray           *dRhoX        ,
                                 const Real1DArray           *dRhoY        ,
                                 const Real1DArray           *dRhoZ        ,
                                 const Real1DArray           *vRho         ,
                                 const Real1DArray           *vSigma       ,
                                 const Real1DArray           *vLaplacianRho,
                                 const Real1DArray           *vTau         ,
                                       SymmetricMatrix       *fock         ,
                                       Real2DArray           *work2D       )
{
    /* . Aliases. */
    Integer1DArray *indices = basisData->indices ;
    Real2DArray *b   = basisData->f   ,
                *bX  = basisData->fX  ,
                *bY  = basisData->fY  ,
                *bZ  = basisData->fZ  ,
                *bXX = basisData->fXX ,
                *bYY = basisData->fYY ,
                *bZZ = basisData->fZZ ;

   /* . Rho. */
    Real2DArray_CopyTo       ( b, work2D, NULL ) ;
    UReal2DArray_ColumnScale ( vRho, work2D ) ;
    USymmetricMatrix_DotProductIncrement ( indices, b, work2D, fock ) ;

    /* . Sigma. */
    DFTIntegrator_FockSigma ( hasSigma, True, basisData, dRhoX, dRhoY, dRhoZ, vSigma, fock, work2D ) ;

    /* . Laplacian. */
    if ( hasLaplacian )
    {
        /* . First derivative contribution. */
        Real2DArray_CopyTo       ( bX, work2D, NULL ) ;
        UReal2DArray_ColumnScale ( vLaplacianRho, work2D ) ;
        Real2DArray_Scale        ( work2D, 2.0e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bX, work2D, fock ) ;
        Real2DArray_CopyTo       ( bY, work2D, NULL ) ;
        UReal2DArray_ColumnScale ( vLaplacianRho, work2D ) ;
        Real2DArray_Scale        ( work2D, 2.0e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bY, work2D, fock ) ;
        Real2DArray_CopyTo       ( bZ, work2D, NULL ) ;
        UReal2DArray_ColumnScale ( vLaplacianRho, work2D ) ;
        Real2DArray_Scale        ( work2D, 2.0e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bZ, work2D, fock ) ;
        /* . Second derivative contribution. */
        Real2DArray_CopyTo         ( bXX, work2D, NULL ) ;
        Real2DArray_AddScaledArray ( work2D, 1.0e+00, bYY, NULL ) ;
        Real2DArray_AddScaledArray ( work2D, 1.0e+00, bZZ, NULL ) ;
        UReal2DArray_ColumnScale   ( vLaplacianRho, work2D ) ;
        USymmetricMatrix_DotProductIncrement ( indices, b, work2D, fock ) ;
        USymmetricMatrix_DotProductIncrement ( indices, work2D, b, fock ) ;
    }

    /* . Tau. */
    if ( hasTau )
    {
        Real2DArray_CopyTo       ( bX, work2D, NULL ) ;
        UReal2DArray_ColumnScale ( vTau, work2D ) ;
        Real2DArray_Scale        ( work2D, 0.5e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bX, work2D, fock ) ;
        Real2DArray_CopyTo       ( bY, work2D, NULL ) ;
        UReal2DArray_ColumnScale ( vTau, work2D ) ;
        Real2DArray_Scale        ( work2D, 0.5e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bY, work2D, fock ) ;
        Real2DArray_CopyTo       ( bZ, work2D, NULL ) ;
        UReal2DArray_ColumnScale ( vTau, work2D ) ;
        Real2DArray_Scale        ( work2D, 0.5e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bZ, work2D, fock ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sigma contributions to a Fock matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_FockSigma ( const Boolean                hasSigma  ,
                                      const Boolean                scale     ,
                                      const GridFunctionDataBlock *basisData ,
                                      const Real1DArray           *dRhoX     ,
                                      const Real1DArray           *dRhoY     ,
                                      const Real1DArray           *dRhoZ     ,
                                      const Real1DArray           *vSigma    ,
                                            SymmetricMatrix       *fock      ,
                                            Real2DArray           *work2D    )
{
    if ( hasSigma )
    {
        /* . Aliases. */
        auto Integer1DArray *indices = basisData->indices ;
        auto Real2DArray    *b  = basisData->f  ,
                            *bX = basisData->fX ,
                            *bY = basisData->fY ,
                            *bZ = basisData->fZ ;
        Real2DArray_CopyTo                ( bX, work2D, NULL ) ;
        UReal2DArray_ColumnScale          ( dRhoX , work2D ) ;
        UReal2DArray_ColumnAddScaledArray ( bY, dRhoY, work2D ) ;
        UReal2DArray_ColumnAddScaledArray ( bZ, dRhoZ, work2D ) ;
        UReal2DArray_ColumnScale          ( vSigma, work2D ) ;
        if ( scale ) Real2DArray_Scale    ( work2D, 2.0e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, b, work2D, fock ) ;
        USymmetricMatrix_DotProductIncrement ( indices, work2D, b, fock ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form a square reduced density matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_FormReducedDensity ( const Integer1DArray   *indices        ,
                                               const SymmetricMatrix  *density        ,
                                                     Real2DArray     **reducedDensity ,
                                                     Status           *status         )
{
    Integer      n = indices->length ;
    Real2DArray *new = NULL, *old = (*reducedDensity) ;
    if ( ( old != NULL ) && ( old->length0 == n ) && ( old->length1 == n ) ) new = old ;
    else
    {
        Real2DArray_Deallocate ( reducedDensity ) ;
        new = Real2DArray_Allocate ( n, n, status ) ;
    }
    if ( new != NULL ) SymmetricMatrix_IndexedCopyToReal2DArray ( density, indices, new, NULL ) ;
    (*reducedDensity) = new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Contributions to the gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_Gradients ( const Boolean                hasSigma      ,
                                      const Boolean                hasLaplacian  ,
                                      const Boolean                hasTau        ,
                                      const Integer1DArray        *atomIndices   ,
                                      const GridFunctionDataBlock *basisData     ,
                                      const Real1DArray           *dRhoX         ,
                                      const Real1DArray           *dRhoY         ,
                                      const Real1DArray           *dRhoZ         ,
                                      const Real1DArray           *vRho          ,
                                      const Real1DArray           *vSigma        ,
                                      const Real1DArray           *vLaplacianRho ,
                                      const Real1DArray           *vTau          ,
                                      const Real2DArray           *density       ,
                                      const Integer                gridAtom      ,
                                            Coordinates3          *gradients3    ,
                                            Real2DArray           *temp2D        ,
                                            Real2DArray           *work2D        )
{
    /* . Aliases. */
    auto Integer1DArray *indices = basisData->indices ;
    auto Real2DArray *b   = basisData->f   ,
                     *bX  = basisData->fX  ,
                     *bY  = basisData->fY  ,
                     *bZ  = basisData->fZ  ,
                     *bXX = basisData->fXX ,
                     *bXY = basisData->fXY ,
                     *bXZ = basisData->fXZ ,
                     *bYY = basisData->fYY ,
                     *bYZ = basisData->fYZ ,
                     *bZZ = basisData->fZZ ;

    /* . Rho. */
    Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, b, 0.0e+00, work2D, NULL ) ;
    UReal2DArray_ColumnScale   ( vRho, work2D ) ;
    Real2DArray_Scale          ( work2D, 2.0e+00 ) ;
    UGradientContributions     ( atomIndices, indices, work2D, bX, bY, bZ, gridAtom, gradients3 ) ;

    /* . Sigma. */
    if ( hasSigma )
    {
        DFTIntegrator_GradientsSigma ( hasSigma, True, atomIndices, basisData, dRhoX, dRhoY, dRhoZ, vSigma, density, gridAtom, gradients3, temp2D, work2D ) ;
    }

    /* . Laplacian. */
    if ( hasLaplacian )
    {
        /* . More aliases. */
        auto Real2DArray *bXXX = basisData->fXXX ,
                         *bXXY = basisData->fXXY ,
                         *bXXZ = basisData->fXXZ ,
                         *bXYY = basisData->fXYY ,
                         *bXYZ = basisData->fXYZ ,
                         *bXZZ = basisData->fXZZ ,
                         *bYYY = basisData->fYYY ,
                         *bYYZ = basisData->fYYZ ,
                         *bYZZ = basisData->fYZZ ,
                         *bZZZ = basisData->fZZZ ,
                         *sum, *sumX, *sumY, *sumZ ;
        /* . Aliases - some of the integrals are destroyed. */
        sum = bXYZ ; sumX = bXXX ; sumY = bYYY ; sumZ = bZZZ ;
        /* . First derivative contribution. */
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, bX, 0.0e+00, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vLaplacianRho, work2D ) ;
        Real2DArray_Scale          ( work2D, 4.0e+00 ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXX, bXY, bXZ, gridAtom, gradients3 ) ;
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, bY, 0.0e+00, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vLaplacianRho, work2D ) ;
        Real2DArray_Scale          ( work2D, 4.0e+00 ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXY, bYY, bYZ, gridAtom, gradients3 ) ;
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, bZ, 0.0e+00, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vLaplacianRho, work2D ) ;
        Real2DArray_Scale          ( work2D, 4.0e+00 ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXZ, bYZ, bZZ, gridAtom, gradients3 ) ;
        /* . Second derivative contribution - 1. */
        Real2DArray_CopyTo         ( bXX, sum, NULL ) ;
        Real2DArray_AddScaledArray ( sum, 1.0e+00, bYY, NULL ) ;
        Real2DArray_AddScaledArray ( sum, 1.0e+00, bZZ, NULL ) ;
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, sum, 0.0e+00, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vLaplacianRho, work2D ) ;
        Real2DArray_Scale          ( work2D, 2.0e+00 ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bX, bY, bZ, gridAtom, gradients3 ) ;
        /* . Second derivative contribution - 2. */
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, b, 0.0e+00, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vLaplacianRho, work2D ) ;
        Real2DArray_Scale          ( work2D, 2.0e+00 ) ;
        Real2DArray_AddScaledArray ( sumX, 1.0e+00, bXYY, NULL ) ;
        Real2DArray_AddScaledArray ( sumX, 1.0e+00, bXZZ, NULL ) ;
        Real2DArray_AddScaledArray ( sumY, 1.0e+00, bXXY, NULL ) ;
        Real2DArray_AddScaledArray ( sumY, 1.0e+00, bYZZ, NULL ) ;
        Real2DArray_AddScaledArray ( sumZ, 1.0e+00, bXXZ, NULL ) ;
        Real2DArray_AddScaledArray ( sumZ, 1.0e+00, bYYZ, NULL ) ;
        UGradientContributions     ( atomIndices, indices, work2D, sumX, sumY, sumZ, gridAtom, gradients3 ) ;
    }

    /* . Tau. */
    if ( hasTau )
    {
        /* . No scaling as factors of 2 and 1/2. */
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, bX, 0.0e+00, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vTau, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXX, bXY, bXZ, gridAtom, gradients3 ) ;
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, bY, 0.0e+00, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vTau, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXY, bYY, bYZ, gridAtom, gradients3 ) ;
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, bZ, 0.0e+00, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vTau, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXZ, bYZ, bZZ, gridAtom, gradients3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sigma contributions to the gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_GradientsSigma ( const Boolean                hasSigma    ,
                                           const Boolean                isSelf      ,
                                           const Integer1DArray        *atomIndices ,
                                           const GridFunctionDataBlock *basisData   ,
                                           const Real1DArray           *dRhoX       ,
                                           const Real1DArray           *dRhoY       ,
                                           const Real1DArray           *dRhoZ       ,
                                           const Real1DArray           *vSigma      ,
                                           const Real2DArray           *density     ,
                                           const Integer                gridAtom    ,
                                                 Coordinates3          *gradients3  ,
                                                 Real2DArray           *temp2D      ,
                                                 Real2DArray           *work2D      )
{
    if ( hasSigma )
    {
        auto Real factor ;
        /* . Aliases. */
        auto Integer1DArray *indices = basisData->indices ;
        auto Real2DArray *b   = basisData->f   ,
                         *bX  = basisData->fX  ,
                         *bY  = basisData->fY  ,
                         *bZ  = basisData->fZ  ,
                         *bXX = basisData->fXX ,
                         *bXY = basisData->fXY ,
                         *bXZ = basisData->fXZ ,
                         *bYY = basisData->fYY ,
                         *bYZ = basisData->fYZ ,
                         *bZZ = basisData->fZZ ;
        if ( isSelf ) factor = 4.0e+00 ;
        else          factor = 2.0e+00 ;
        /* . Contribution 1. */
        Real2DArray_CopyTo         ( bX, temp2D, NULL ) ;
        UReal2DArray_ColumnScale   ( dRhoX, temp2D ) ;
        UReal2DArray_ColumnAddScaledArray ( bY, dRhoY, temp2D ) ;
        UReal2DArray_ColumnAddScaledArray ( bZ, dRhoZ, temp2D ) ;
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, temp2D, 0.0e+00, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vSigma, work2D ) ;
        Real2DArray_Scale          ( work2D, factor ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bX, bY, bZ, gridAtom, gradients3 ) ;
        /* . Set up for contributions 2, 3 and 4. */
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, density, b, 0.0e+00, temp2D, NULL ) ;
        UReal2DArray_ColumnScale   ( vSigma, temp2D ) ;
        Real2DArray_Scale          ( temp2D, factor ) ;
        /* . Contribution 2. */
        Real2DArray_CopyTo         ( temp2D, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( dRhoX, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXX, bXY, bXZ, gridAtom, gradients3 ) ;
        /* . Contribution 3. */
        Real2DArray_CopyTo         ( temp2D, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( dRhoY, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXY, bYY, bYZ, gridAtom, gradients3 ) ;
        /* . Contribution 4. */
        Real2DArray_CopyTo         ( temp2D, work2D, NULL ) ;
        UReal2DArray_ColumnScale   ( dRhoZ, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXZ, bYZ, bZZ, gridAtom, gradients3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of a single density and associated quantities at the grid points.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_GridPointRho ( const Boolean                hasSigma     ,
                                         const Boolean                hasLaplacian ,
                                         const Boolean                hasTau       ,
                                         const GridFunctionDataBlock *basisData    ,
                                         const Real2DArray           *density      ,
                                               Real1DArray           *rho          ,
                                               Real1DArray           *dRhoX        ,
                                               Real1DArray           *dRhoY        ,
                                               Real1DArray           *dRhoZ        ,
                                               Real1DArray           *sigma        ,
                                               Real1DArray           *laplacianRho ,
                                               Real1DArray           *tau          ,
                                               Real1DArray           *work1D       ,
                                               Real2DArray           *work2D       )
{
    /* . Aliases. */
    Real2DArray *b   = basisData->f   ,
                *bX  = basisData->fX  ,
                *bY  = basisData->fY  ,
                *bZ  = basisData->fZ  ,
                *bXX = basisData->fXX ,
                *bYY = basisData->fYY ,
                *bZZ = basisData->fZZ ;

    /* . Rho. */
    Real2DArray_MatrixMultiply    ( False, False, 1.0e+00, density, b, 0.0e+00, work2D, NULL ) ;
    Real2DArray_ColumnDotProducts ( True, b, work2D, rho ) ;

    /* . Sigma. */
    if ( hasSigma )
    {
        Real2DArray_ColumnDotProducts ( True, bX, work2D, dRhoX ) ;
        Real2DArray_ColumnDotProducts ( True, bY, work2D, dRhoY ) ;
        Real2DArray_ColumnDotProducts ( True, bZ, work2D, dRhoZ ) ;
        DFTIntegrator_GridPointSigma ( True, dRhoX, dRhoY, dRhoZ, dRhoX, dRhoY, dRhoZ, sigma ) ;
        Real1DArray_Scale ( dRhoX, 2.0e+00 ) ;
        Real1DArray_Scale ( dRhoY, 2.0e+00 ) ;
        Real1DArray_Scale ( dRhoZ, 2.0e+00 ) ;
        Real1DArray_Scale ( sigma, 4.0e+00 ) ;
    }

    /* . Laplacian or tau. */
    if ( hasLaplacian || hasTau )
    {
        Real2DArray_MatrixMultiply    ( False, False, 1.0e+00, density, bX, 0.0e+00, work2D, NULL ) ;
        Real2DArray_ColumnDotProducts ( True , bX, work2D, work1D ) ;
        Real2DArray_MatrixMultiply    ( False, False, 1.0e+00, density, bY, 0.0e+00, work2D, NULL ) ;
        Real2DArray_ColumnDotProducts ( False, bY, work2D, work1D ) ;
        Real2DArray_MatrixMultiply    ( False, False, 1.0e+00, density, bZ, 0.0e+00, work2D, NULL ) ;
        Real2DArray_ColumnDotProducts ( False, bZ, work2D, work1D ) ;
    }

    /* . Laplacian. */
    if ( hasLaplacian )
    {
        Real1DArray_CopyTo ( work1D, laplacianRho, NULL ) ;
        Real2DArray_ColumnDotProducts ( False, bXX, work2D, laplacianRho ) ;
        Real2DArray_ColumnDotProducts ( False, bYY, work2D, laplacianRho ) ;
        Real2DArray_ColumnDotProducts ( False, bZZ, work2D, laplacianRho ) ;
        Real1DArray_Scale  ( laplacianRho, 2.0e+00 ) ;
    }

    /* . Tau. */
    if ( hasTau )
    {
        Real1DArray_CopyTo ( work1D, tau, NULL ) ;
        Real1DArray_Scale ( tau, 0.5e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the cross-sigma values for two densities at the grid points.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_GridPointSigma ( const Boolean      hasSigma ,
                                           const Real1DArray *dRhoXa   ,
                                           const Real1DArray *dRhoYa   ,
                                           const Real1DArray *dRhoZa   ,
                                           const Real1DArray *dRhoXb   ,
                                           const Real1DArray *dRhoYb   ,
                                           const Real1DArray *dRhoZb   ,
                                                 Real1DArray *sigma    )
{
    if ( hasSigma )
    {
        auto Integer p ;
        for ( p = 0 ; p < sigma->length ; p++ )
        {
            Real1DArray_Item ( sigma, p ) = Real1DArray_Item ( dRhoXa, p ) * Real1DArray_Item ( dRhoXb, p ) +
                                            Real1DArray_Item ( dRhoYa, p ) * Real1DArray_Item ( dRhoYb, p ) +
                                            Real1DArray_Item ( dRhoZa, p ) * Real1DArray_Item ( dRhoZb, p ) ;
        }
    }
}

/*==================================================================================================================================
! . Local utilities - may be generalized and moved.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Utility procedure for determining gradient contributions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void UGradientContributions ( const Integer1DArray *atomIndices ,
                                     const Integer1DArray *indices     ,
                                     const Real2DArray    *a           ,
                                     const Real2DArray    *x           ,
                                     const Real2DArray    *y           ,
                                     const Real2DArray    *z           ,
                                     const Integer         gridAtom    ,
                                           Coordinates3   *gradients3  )
{
    Integer     c, i, m ;
    Real        gX, gY, gZ ;
    Real1DArray sliceA, sliceB ;
    for ( i = 0 ; i < indices->length ; i++ )
    {
        m = Integer1DArray_Item ( indices, i ) ;
        c = Integer1DArray_Item ( atomIndices, m ) ;
        /* . Calculation the contributions. */
        Real2DArray_RowSlice ( a     , i, &sliceA, NULL ) ;
        Real2DArray_RowSlice ( x     , i, &sliceB, NULL ) ;
        gX = Real1DArray_Dot ( &sliceA  , &sliceB, NULL ) ;
        Real2DArray_RowSlice ( y     , i, &sliceB, NULL ) ;
        gY = Real1DArray_Dot ( &sliceA  , &sliceB, NULL ) ;
        Real2DArray_RowSlice ( z     , i, &sliceB, NULL ) ;
        gZ = Real1DArray_Dot ( &sliceA  , &sliceB, NULL ) ;
        /* . Add in the contributions. */
        Coordinates3_DecrementRow ( gradients3, c, gX, gY, gZ ) ;
# ifdef DFTGRIDWEIGHTDERIVATIVES
        /* . This is the dE/drg term (i.e. the derivative with respect to the grid point which belongs to gridAtom). */
        Coordinates3_IncrementRow ( gradients3, gridAtom, gX, gY, gZ ) ;
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add scaled array by column for a 2-D array.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void UReal2DArray_ColumnAddScaledArray ( const Real2DArray *a       ,
                                                const Real1DArray *weights ,
                                                      Real2DArray *b       )
{
    Integer     i ;
    Real        w ;
    Real1DArray sliceA, sliceB ;
    for ( i = 0 ; i < weights->length ; i++ )
    {
        w = Real1DArray_Item ( weights, i ) ;
        Real2DArray_ColumnSlice ( a, i, &sliceA, NULL ) ;
        Real2DArray_ColumnSlice ( b, i, &sliceB, NULL ) ;
        Real1DArray_AddScaledArray ( &sliceB, w, &sliceA, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the columns of a 2-D array.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void UReal2DArray_ColumnScale ( const Real1DArray *a ,
                                             Real2DArray *b )
{
    Integer     i ;
    Real        w ;
    Real1DArray slice ;
    for ( i = 0 ; i < a->length ; i++ )
    {
        w = Real1DArray_Item ( a, i ) ;
        Real2DArray_ColumnSlice ( b, i, &slice, NULL ) ;
        Real1DArray_Scale ( &slice, w ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment selected entries of a symmetric matrix using the dot products of the rows of two 2-D arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void USymmetricMatrix_DotProductIncrement ( const Integer1DArray  *indices ,
                                                   const Real2DArray     *a       ,
                                                   const Real2DArray     *b       ,
                                                         SymmetricMatrix *c       )
{
    Integer     i, j, m, n ;
    Real1DArray sliceA, sliceB ;
    for ( i = 0 ; i < indices->length ; i++ )
    {
        m = Integer1DArray_Item ( indices, i ) ;
        Real2DArray_RowSlice ( a, i, &sliceA, NULL ) ;
        for ( j = 0 ; j <= i ; j++ )
        {
            n = Integer1DArray_Item ( indices, j ) ;
            Real2DArray_RowSlice ( b, j, &sliceB, NULL ) ;
            SymmetricMatrix_Item ( c, m, n ) += Real1DArray_Dot ( &sliceA, &sliceB, NULL ) ;
        }
    }
}
