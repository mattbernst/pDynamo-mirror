/*------------------------------------------------------------------------------
! . File      : DFTFunctionalModel.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module defines the DFT functional model. It is an interface to the libxc library.
!=================================================================================================================================*/

# include "DFTFunctionalModel.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTFunctionalModel *DFTFunctionalModel_Allocate ( const Integer numberOfFunctionals, Status *status )
{
    DFTFunctionalModel *self = NULL ;
    MEMORY_ALLOCATE ( self, DFTFunctionalModel ) ;
    if ( self != NULL )
    {
        auto Integer n ;
        n = Maximum ( numberOfFunctionals, 0 ) ;
        /* . Basic data. */
        self->functionals         = NULL  ;
        self->hasLaplacian        = False ;
        self->hasSigma            = False ;
        self->hasTau              = False ;
        self->isSpinRestricted    = True  ;
        self->numberOfFunctionals = n     ;
        self->order               = -1    ;
        /* . Functionals array. */
        if ( n > 0 )
        {
            MEMORY_ALLOCATEARRAY ( self->functionals, n, xc_func_type ) ;
            if ( self->functionals == NULL ) DFTFunctionalModel_Deallocate ( &self ) ;
        }
        else if ( n < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTFunctionalModel *DFTFunctionalModel_Clone ( const DFTFunctionalModel *self, Status *status )
{
    DFTFunctionalModel *clone = NULL ;
    if ( self != NULL )
    {
        auto Integer f ;
        auto Integer1DArray *ids ;
        ids   = Integer1DArray_Allocate ( self->numberOfFunctionals, status ) ;
        for ( f = 0 ; f < self->numberOfFunctionals ; f++ ) Integer1DArray_Item ( ids, f ) = self->functionals[f].info->number ;
        clone = DFTFunctionalModel_MakeFromIDs ( ids, self->isSpinRestricted, status ) ;
        Integer1DArray_Deallocate ( &ids ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTFunctionalModel_Deallocate ( DFTFunctionalModel **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        auto Integer f ;
        for ( f = 0 ; f < (*self)->numberOfFunctionals ; f++ ) xc_func_end ( &((*self)->functionals[f]) ) ;
        MEMORY_DEALLOCATE ( (*self)->functionals ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluation of the energy density and its first derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTFunctionalModel_Evaluate ( const DFTFunctionalModel *self, DFTIntegratorDataBlock *data )
{
    if ( ( self != NULL ) && ( data != NULL ) )
    {
        auto Integer f ;
        auto xc_func_type *functional ;
        DFTIntegratorDataBlock_Initialize ( data ) ;
        for ( f = 0 ; f < self->numberOfFunctionals ; f++ )
        {
            functional = &(self->functionals[f]) ;
            switch ( functional->info->family )
            {
                case XC_FAMILY_LDA:
                    xc_lda_exc_vxc ( functional                           ,
                                     data->numberOfPoints                 ,
                                     Real2DArray_Data ( data->rho       ) ,
                                     Real1DArray_Data ( data->localEXC  ) ,
                                     Real2DArray_Data ( data->localVRho ) ) ;
                    break ;
                case XC_FAMILY_GGA:
                case XC_FAMILY_HYB_GGA:
                    xc_gga_exc_vxc ( functional                             ,
                                     data->numberOfPoints                   ,
                                     Real2DArray_Data ( data->rho         ) ,
                                     Real2DArray_Data ( data->sigma       ) ,
                                     Real1DArray_Data ( data->localEXC    ) ,
                                     Real2DArray_Data ( data->localVRho   ) ,
                                     Real2DArray_Data ( data->localVSigma ) ) ;
                    break ;
                case XC_FAMILY_MGGA:
                case XC_FAMILY_HYB_MGGA:
                    xc_mgga_exc_vxc ( functional                                    ,
                                      data->numberOfPoints                          ,
                                      Real2DArray_Data ( data->rho                ) ,
                                      Real2DArray_Data ( data->sigma              ) ,
                                      Real2DArray_Data ( data->laplacianRho       ) ,
                                      Real2DArray_Data ( data->tau                ) ,
                                      Real1DArray_Data ( data->localEXC           ) ,
                                      Real2DArray_Data ( data->localVRho          ) ,
                                      Real2DArray_Data ( data->localVSigma        ) ,
                                      Real2DArray_Data ( data->localVLaplacianRho ) ,
                                      Real2DArray_Data ( data->localVTau          ) ) ;
            }
            DFTIntegratorDataBlock_Accumulate ( data ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor given an array of functional IDs.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTFunctionalModel *DFTFunctionalModel_MakeFromIDs ( const Integer1DArray *ids, const Boolean isSpinRestricted, Status *status )
{
    DFTFunctionalModel *self = NULL ;
    if ( ( ids != NULL ) && ( Integer1DArray_Length ( ids ) > 0 ) )
    {
        self = DFTFunctionalModel_Allocate ( Integer1DArray_Length ( ids ), status ) ;
        if ( self != NULL )
        {
            auto Boolean isOK ;
            auto Integer f, failures, id, spin ;
            self->isSpinRestricted = isSpinRestricted ;
            if ( isSpinRestricted ) spin = XC_UNPOLARIZED ;
            else                    spin = XC_POLARIZED   ;
            for ( f = failures = 0 ; f < self->numberOfFunctionals ; f++ )
            {
                id   = Integer1DArray_Item ( ids, f ) ;
                isOK = xc_func_init ( &(self->functionals[f]), id, spin ) >= 0 ;
                if ( isOK )
                {
                    switch ( self->functionals[f].info->family )
                    {
                        case XC_FAMILY_LDA:
                            self->order        = Maximum ( self->order, 0 ) ;
                            break ;
                        case XC_FAMILY_GGA:
                        case XC_FAMILY_HYB_GGA:
                            self->hasSigma     = True ;
                            self->order        = Maximum ( self->order, 1 ) ;
                            break ;
                        case XC_FAMILY_MGGA:
                        case XC_FAMILY_HYB_MGGA:
                            self->hasLaplacian = True ;
                            self->hasSigma     = True ;
                            self->hasTau       = True ;
                            self->order        = Maximum ( self->order, 2 ) ;
                    }
                }
                else failures += 1 ;
            }
            if ( failures > 0 )
            {
                DFTFunctionalModel_Deallocate ( &self ) ;
                Status_Set ( status, Status_InvalidArgument ) ;
            }
        }
        else Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    return self ;
}
