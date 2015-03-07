/*------------------------------------------------------------------------------
! . File      : QCMMInteractionState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module implements a container for quantities needed for the calculation of QC/MM interactions.
!=================================================================================================================================*/

# include "Memory.h"
# include "QCMMInteractionState.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCMMInteractionState *QCMMInteractionState_Allocate ( const Integer n, const Boolean includeQCQC, Status *status )
{
    QCMMInteractionState *self = NULL ;
    if ( n > 0 )
    {
        self = ( QCMMInteractionState * ) Memory_Allocate ( sizeof ( QCMMInteractionState ) ) ;
        if ( self != NULL )
        {
            self->qcCharges      = Real1DArray_Allocate ( n, status ) ;
            self->qcmmPotentials = Real1DArray_Allocate ( n, status ) ;
            if ( includeQCQC ) self->qcqcPotentials = SymmetricMatrix_Allocate ( n ) ;
            else               self->qcqcPotentials = NULL ;
            /* . Deallocate if there is not enough memory. */
            if ( ( self->qcCharges == NULL ) || ( self->qcmmPotentials == NULL ) || ( includeQCQC && ( self->qcqcPotentials == NULL ) ) )
            {
                QCMMInteractionState_Deallocate ( &self ) ;
            }
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate the QCQC potentials array for an already existing state.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Nothing is done if the array already exists. */
void QCMMInteractionState_AllocateQCQCPotentials ( QCMMInteractionState *self, Status *status )
{
    if ( ( self != NULL ) && ( self->qcqcPotentials == NULL ) )
    {
        auto Integer n ;
        n = self->qcCharges->length ;
        self->qcqcPotentials = SymmetricMatrix_Allocate ( n ) ;
        if ( self->qcqcPotentials == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCMMInteractionState_Deallocate ( QCMMInteractionState **self )
{
    if ( (*self) != NULL )
    {
        Real1DArray_Deallocate     ( &((*self)->qcCharges)      ) ;
        Real1DArray_Deallocate     ( &((*self)->qcmmPotentials) ) ;
        SymmetricMatrix_Deallocate ( &((*self)->qcqcPotentials) ) ;
        Memory_Deallocate          ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization for an energy/gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCMMInteractionState_Initialize ( QCMMInteractionState *self )
{
    if ( self != NULL )
    {
        Real1DArray_Set     ( self->qcCharges,      0.0e+00 ) ;
        Real1DArray_Set     ( self->qcmmPotentials, 0.0e+00 ) ;
        SymmetricMatrix_Set ( self->qcqcPotentials, 0.0e+00 ) ;
    }
}
