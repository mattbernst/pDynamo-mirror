/*------------------------------------------------------------------------------
! . File      : DFTIntegratorDataBlock.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module defines a block data structure that is needed for DFT integration.
!=================================================================================================================================*/

# include "DFTIntegratorDataBlock.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Accumulation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTIntegratorDataBlock_Accumulate ( DFTIntegratorDataBlock *self )
{
    if ( ( self != NULL ) && ( self->hasLocalData ) )
    {
        Real1DArray_AddScaledArray ( self->eXC          , 1.0e+00, self->localEXC          , NULL ) ;
        Real2DArray_AddScaledArray ( self->vLaplacianRho, 1.0e+00, self->localVLaplacianRho, NULL ) ;
        Real2DArray_AddScaledArray ( self->vRho         , 1.0e+00, self->localVRho         , NULL ) ;
        Real2DArray_AddScaledArray ( self->vSigma       , 1.0e+00, self->localVSigma       , NULL ) ;
        Real2DArray_AddScaledArray ( self->vTau         , 1.0e+00, self->localVTau         , NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTIntegratorDataBlock *DFTIntegratorDataBlock_Allocate ( const Integer numberOfFunctionals ,
                                                          const Integer numberOfPoints      ,
                                                          const Boolean hasSigma            ,
                                                          const Boolean hasLaplacian        ,
                                                          const Boolean hasTau              ,
                                                          const Boolean isSpinRestricted    ,
                                                                Status  *status             )
{
    DFTIntegratorDataBlock *self = NULL ;
    MEMORY_ALLOCATE ( self, DFTIntegratorDataBlock ) ;
    if ( self != NULL )
    {
        auto Integer n ;
        n = Maximum ( numberOfPoints    , 0 ) ;
        self->hasLocalData   = ( numberOfFunctionals > 1 ) ;
        self->numberOfPoints = n ;
        /* . Basic initialization. */
        self->eXC                = NULL ;
        self->localEXC           = NULL ;
        self->dRhoX              = NULL ;
        self->dRhoY              = NULL ;
        self->dRhoZ              = NULL ;
        self->localVLaplacianRho = NULL ;
        self->localVRho          = NULL ;
        self->localVSigma        = NULL ;
        self->localVTau          = NULL ;
        self->laplacianRho       = NULL ;
        self->rho                = NULL ;
        self->sigma              = NULL ;
        self->tau                = NULL ;
        self->vLaplacianRho      = NULL ;
        self->vRho               = NULL ;
        self->vSigma             = NULL ;
        self->vTau               = NULL ;
        /* . Allocation. */ 
        if ( n > 0 )
        {
            auto Boolean isOK ;
            auto Integer c, d ;
            /* . Determine array sizes. */
            if ( isSpinRestricted ) { c = 1 ; d = 1 ; }
            else                    { c = 2 ; d = 3 ; }
            /* . Allocation. */
            self->eXC    = Real1DArray_Allocate ( n,    status ) ;
            self->rho    = Real2DArray_Allocate ( n, c, status ) ;
            self->vRho   = Real2DArray_Allocate ( n, c, status ) ;
            if ( self->hasLocalData )
            {
                self->localEXC  = Real1DArray_Allocate ( n,    status ) ;
                self->localVRho = Real2DArray_Allocate ( n, c, status ) ;
            }
            else { self->localEXC = self->eXC ; self->localVRho = self->vRho ; }
            isOK = ( ( self->eXC       != NULL ) &&
                     ( self->localEXC  != NULL ) &&
                     ( self->localVRho != NULL ) &&
                     ( self->rho       != NULL ) &&
                     ( self->vRho      != NULL ) ) ;
            if ( hasSigma )
            {
                self->dRhoX       = Real2DArray_Allocate ( n, c, status ) ;
                self->dRhoY       = Real2DArray_Allocate ( n, c, status ) ;
                self->dRhoZ       = Real2DArray_Allocate ( n, c, status ) ;
                self->sigma       = Real2DArray_Allocate ( n, d, status ) ;
                self->vSigma      = Real2DArray_Allocate ( n, d, status ) ;
                if ( self->hasLocalData ) self->localVSigma = Real2DArray_Allocate ( n, d, status ) ;
                else                      self->localVSigma = self->vSigma ;
                isOK = isOK && ( ( self->dRhoX       != NULL ) &&
                                 ( self->dRhoY       != NULL ) &&
                                 ( self->dRhoZ       != NULL ) &&
                                 ( self->localVSigma != NULL ) &&
                                 ( self->sigma       != NULL ) &&
                                 ( self->vSigma      != NULL ) ) ;
            }
            if ( hasLaplacian )
            {
                self->laplacianRho  = Real2DArray_Allocate ( n, c, status ) ;
                self->vLaplacianRho = Real2DArray_Allocate ( n, c, status ) ;
                if ( self->hasLocalData ) self->localVLaplacianRho = Real2DArray_Allocate ( n, c, status ) ;
                else                      self->localVLaplacianRho = self->vLaplacianRho ;
                isOK = isOK && ( ( self->localVLaplacianRho != NULL ) &&
                                 ( self->laplacianRho       != NULL ) &&
                                 ( self->vLaplacianRho      != NULL ) ) ;
            }
            if ( hasTau )
            {
                self->tau  = Real2DArray_Allocate ( n, c, status ) ;
                self->vTau = Real2DArray_Allocate ( n, c, status ) ;
                if ( self->hasLocalData ) self->localVTau = Real2DArray_Allocate ( n, c, status ) ;
                else                      self->localVTau = self->vTau ;
                isOK = isOK && ( ( self->localVTau != NULL ) &&
                                 ( self->tau       != NULL ) &&
                                 ( self->vTau      != NULL ) ) ;
            }
            if ( isOK )
            {   
                DFTIntegratorDataBlock_InitializeView ( self, 0, &(self->viewP) ) ;
                if ( ! isSpinRestricted )
                {
                    DFTIntegratorDataBlock_InitializeView ( self, 1, &(self->viewQ) ) ;
                    Real2DArray_ColumnSlice ( self->sigma , 1, &(self->sigmaPQ ), NULL ) ;
                    Real2DArray_ColumnSlice ( self->vSigma, 1, &(self->vSigmaPQ), NULL ) ;
                }
            }
            else DFTIntegratorDataBlock_Deallocate ( &self ) ;
        }
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTIntegratorDataBlock_Deallocate ( DFTIntegratorDataBlock **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Real1DArray_Deallocate ( &((*self)->eXC          ) ) ;
        Real2DArray_Deallocate ( &((*self)->dRhoX        ) ) ;
        Real2DArray_Deallocate ( &((*self)->dRhoY        ) ) ;
        Real2DArray_Deallocate ( &((*self)->dRhoZ        ) ) ;
        Real2DArray_Deallocate ( &((*self)->laplacianRho ) ) ;
        Real2DArray_Deallocate ( &((*self)->rho          ) ) ;
        Real2DArray_Deallocate ( &((*self)->sigma        ) ) ;
        Real2DArray_Deallocate ( &((*self)->tau          ) ) ;
        Real2DArray_Deallocate ( &((*self)->vLaplacianRho) ) ;
        Real2DArray_Deallocate ( &((*self)->vRho         ) ) ;
        Real2DArray_Deallocate ( &((*self)->vSigma       ) ) ;
        Real2DArray_Deallocate ( &((*self)->vTau         ) ) ;
        if ( (*self)->hasLocalData )
        {
            Real1DArray_Deallocate ( &((*self)->localEXC          ) ) ;
            Real2DArray_Deallocate ( &((*self)->localVLaplacianRho) ) ;
            Real2DArray_Deallocate ( &((*self)->localVRho         ) ) ;
            Real2DArray_Deallocate ( &((*self)->localVSigma       ) ) ;
            Real2DArray_Deallocate ( &((*self)->localVTau         ) ) ;
        }
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization., DFTIntegratorDataBlockView *view
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTIntegratorDataBlock_Initialize ( DFTIntegratorDataBlock *self )
{
    if ( ( self != NULL ) && ( self->hasLocalData ) )
    {
        Real1DArray_Set ( self->eXC          , 0.0e+00 ) ;
        Real2DArray_Set ( self->vLaplacianRho, 0.0e+00 ) ;
        Real2DArray_Set ( self->vRho         , 0.0e+00 ) ;
        Real2DArray_Set ( self->vSigma       , 0.0e+00 ) ;
        Real2DArray_Set ( self->vTau         , 0.0e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTIntegratorDataBlock_InitializeView ( DFTIntegratorDataBlock *self, const Integer c, DFTIntegratorDataBlockView *view )
{
    if ( ( self != NULL ) && ( view != NULL ) )
    {
        auto Integer d ;
        d = ( c == 0 ? 0 : 2 ) ;
        Real2DArray_ColumnSlice ( self->dRhoX        , c, &(view->dRhoX        ), NULL ) ;
        Real2DArray_ColumnSlice ( self->dRhoY        , c, &(view->dRhoY        ), NULL ) ;
        Real2DArray_ColumnSlice ( self->dRhoZ        , c, &(view->dRhoZ        ), NULL ) ;
        Real2DArray_ColumnSlice ( self->laplacianRho , c, &(view->laplacianRho ), NULL ) ;
        Real2DArray_ColumnSlice ( self->rho          , c, &(view->rho          ), NULL ) ;
        Real2DArray_ColumnSlice ( self->sigma        , d, &(view->sigma        ), NULL ) ;
        Real2DArray_ColumnSlice ( self->tau          , c, &(view->tau          ), NULL ) ;
        Real2DArray_ColumnSlice ( self->vLaplacianRho, c, &(view->vLaplacianRho), NULL ) ;
        Real2DArray_ColumnSlice ( self->vRho         , c, &(view->vRho         ), NULL ) ;
        Real2DArray_ColumnSlice ( self->vSigma       , d, &(view->vSigma       ), NULL ) ;
        Real2DArray_ColumnSlice ( self->vTau         , c, &(view->vTau         ), NULL ) ;
    }
}
