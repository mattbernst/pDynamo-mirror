/*------------------------------------------------------------------------------
! . File      : NBModelFullState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . This module implements full non-bonding interaction state procedures.
!=============================================================================*/

# include "Memory.h"
# include "NBModelFullState.h"

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
NBModelFullState *NBModelFullState_Allocate ( const Integer n )
{
    NBModelFullState *self = NULL ;
    if ( n > 0 )
    {
        self = ( NBModelFullState * ) Memory_Allocate ( sizeof ( NBModelFullState ) ) ;
        if ( self != NULL )
        {
            /* . Scalars. */
            self->emmel          = 0.0e+00 ;
            self->emmel14        = 0.0e+00 ;
            self->emmlj          = 0.0e+00 ;
            self->emmlj14        = 0.0e+00 ;
            self->eqcmmlj        = 0.0e+00 ;
            self->eqcmmlj14      = 0.0e+00 ;
            /* . Pointers. */
            self->coordinates3   = NULL ;
            self->exclusions     = NULL ;
            self->gradients3     = NULL ;
            self->interactions14 = NULL ;
            self->ljParameters   = NULL ;
            self->ljParameters14 = NULL ;
            self->mmAtoms        = NULL ;
            self->qcAtoms        = NULL ;
            self->qcCharges      = NULL ;
            self->qcmmPotentials = NULL ;
            /* . Arrays to be allocated. */
            self->QE14    = Memory_Allocate_Array_Boolean_Initialize ( n, False ) ;
            self->QFREE   = Memory_Allocate_Array_Boolean_Initialize ( n, True  ) ;
            self->QINCL   = Memory_Allocate_Array_Boolean_Initialize ( n, True  ) ;
            self->QMM     = Memory_Allocate_Array_Boolean_Initialize ( n, True  ) ;
            self->baindex = Memory_Allocate_Array_Integer_Initialize ( n, -1    ) ;
            /* . Deallocate if there is not enough memory. */
            if ( ( self->QE14 == NULL ) || ( self->QFREE == NULL ) || ( self->QINCL == NULL ) || ( self->QMM == NULL ) || ( self->baindex == NULL ) ) NBModelFullState_Deallocate ( &self ) ;
        }
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void NBModelFullState_Deallocate ( NBModelFullState **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate_Boolean ( &((*self)->QE14)    ) ;
        Memory_Deallocate_Boolean ( &((*self)->QFREE)   ) ;
        Memory_Deallocate_Boolean ( &((*self)->QINCL)   ) ;
        Memory_Deallocate_Boolean ( &((*self)->QMM)     ) ;
        Memory_Deallocate_Integer ( &((*self)->baindex) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Initialization for an energy/gradient calculation.
!-----------------------------------------------------------------------------*/
void NBModelFullState_Initialize ( NBModelFullState *self, Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    if ( self != NULL )
    {
        /* . Energies. */
        self->emmel     = 0.0e+00 ;
        self->emmel14   = 0.0e+00 ;
        self->emmlj     = 0.0e+00 ;
        self->emmlj14   = 0.0e+00 ;
        self->eqcmmlj   = 0.0e+00 ;
        self->eqcmmlj14 = 0.0e+00 ;
        /* . Coordinates and gradients. */
        self->coordinates3 = coordinates3 ;
        self->gradients3   = gradients3   ;
        /* . QC/MM quantities. */
        Real1DArray_Set ( self->qcCharges,      0.0e+00 ) ;
        Real1DArray_Set ( self->qcmmPotentials, 0.0e+00 ) ;
    }
}

/*------------------------------------------------------------------------------
! . Setup the state.
!-----------------------------------------------------------------------------*/
NBModelFullState *NBModelFullState_Setup ( MMAtomContainer *mmAtoms, QCAtomContainer *qcAtoms, Selection *fixedatoms, PairList *exclusions, PairList *interactions14,
                                                                                            LJParameterContainer *ljParameters, LJParameterContainer *ljParameters14,
                                                                                                                           Real1DArray *qcCharges, Real1DArray *qcmmPotentials )
{
    NBModelFullState *self = NULL ;
    if ( mmAtoms != NULL )
    {
        auto Integer i, n ;
        n = mmAtoms->natoms ;
        self = NBModelFullState_Allocate ( n ) ;
        if ( self != NULL )
        {
            /* . Make sure that the pairlists are in the correct format. */
            SelfPairList_MakeConnections ( exclusions,     n ) ;
            SelfPairList_MakeConnections ( interactions14, n ) ;

            /* . Modify the boolean arrays. */
            /* . Fixed atoms. */
            if ( fixedatoms != NULL )
            {
                for ( n = 0 ; n < fixedatoms->nindices ; n++ ) self->QFREE[fixedatoms->indices[n]] = False ;
            }
            /* . The MM atoms. */
            if ( qcAtoms != NULL )
            {
                for ( n = 0 ; n < qcAtoms->natoms ; n++ )
                {
                    i = qcAtoms->data[n].index ;
                    if ( qcAtoms->data[n].QBOUNDARY ) self->baindex[i] = n     ; /* . Save the boundary atom QC index.    */
                    else                              self->QMM[i]     = False ; /* . Only non-boundary atoms are not MM. */
                }
            }

            /* . Assign aliases. */
            self->exclusions     = exclusions     ;
            self->interactions14 = interactions14 ;
            self->ljParameters   = ljParameters   ;
            self->ljParameters14 = ljParameters14 ;
            self->mmAtoms        = mmAtoms        ;
            self->qcAtoms        = qcAtoms        ;
            self->qcCharges      = qcCharges      ;
            self->qcmmPotentials = qcmmPotentials ;
        }
    }
    return self ;
}
