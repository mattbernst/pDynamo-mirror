/*------------------------------------------------------------------------------
! . File      : NBModelORCAState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . This module implements full non-bonding interaction state procedures which
! . are compatible with the ORCA program.
!=============================================================================*/

# include "Memory.h"
# include "NBModelORCAState.h"

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
NBModelOrcaState *NBModelORCAState_Allocate ( const Integer n )
{
    NBModelOrcaState *self = NULL ;
    if ( n > 0 )
    {
        self = ( NBModelOrcaState * ) Memory_Allocate ( sizeof ( NBModelOrcaState ) ) ;
        if ( self != NULL )
        {
            /* . Options. */
            self->qcmmcoupling     = QCMMLinkAtomCoupling_RC ;
            /* . Scalars. */
            self->emmel            = 0.0e+00 ;
            self->emmel14          = 0.0e+00 ;
            self->emmlj            = 0.0e+00 ;
            self->emmlj14          = 0.0e+00 ;
            self->eqcmmlj          = 0.0e+00 ;
            self->eqcmmlj14        = 0.0e+00 ;
            /* . Pointers. */
            self->coordinates3     = NULL ;
            self->exclusions       = NULL ;
            self->gradients3       = NULL ;
            self->interactions14   = NULL ;
            self->ljParameters     = NULL ;
            self->ljParameters14   = NULL ;
            self->mmAtoms          = NULL ;
            self->qcAtoms          = NULL ;
            /* . Arrays to be allocated. */
            self->QE14             = Memory_Allocate_Array_Boolean_Initialize ( n, False ) ;
            self->QFREE            = Memory_Allocate_Array_Boolean_Initialize ( n, True  ) ;
            self->QINCL            = Memory_Allocate_Array_Boolean_Initialize ( n, True  ) ;
            self->QMM              = Memory_Allocate_Array_Boolean_Initialize ( n, True  ) ;
            /* . MM atom arrays. */
            self->mmcoordinates3   = NULL ;
            self->mmCharges        = NULL ;
            self->mmexclusions     = NULL ;
            self->mmgradients3     = NULL ;
            self->mminteractions14 = NULL ;
            /* . Deallocate if there is not enough memory. */
            if ( ( self->QE14 == NULL ) || ( self->QFREE == NULL ) || ( self->QINCL == NULL ) || ( self->QMM == NULL ) ) NBModelORCAState_Deallocate ( &self ) ;
        }
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void NBModelORCAState_Deallocate ( NBModelOrcaState **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate_Boolean ( &((*self)->QE14)  ) ;
        Memory_Deallocate_Boolean ( &((*self)->QFREE) ) ;
        Memory_Deallocate_Boolean ( &((*self)->QINCL) ) ;
        Memory_Deallocate_Boolean ( &((*self)->QMM)   ) ;
        Coordinates3_Deallocate ( &((*self)->mmcoordinates3)   ) ;
        Coordinates3_Deallocate   ( &((*self)->mmgradients3)     ) ;
        PairList_Deallocate     ( &((*self)->mmexclusions)     ) ;
        PairList_Deallocate     ( &((*self)->mminteractions14) ) ;
        Real1DArray_Deallocate  ( &((*self)->mmCharges)        ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Finalization for an energy/gradient calculation.
!-----------------------------------------------------------------------------*/
void NBModelORCAState_Finalize ( NBModelOrcaState *self )
{
    if ( self != NULL )
    {
        /* . MM atoms. */
        QCAtomContainer_SetMMGradients3 ( self->qcAtoms, self->mmgradients3, self->gradients3 ) ;
    }
}

/*------------------------------------------------------------------------------
! . Initialization for an energy/gradient calculation.
!-----------------------------------------------------------------------------*/
void NBModelORCAState_Initialize ( NBModelOrcaState *self, Coordinates3 *coordinates3, Coordinates3 *gradients3 )
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
        /* . MM atoms. */
        QCAtomContainer_GetMMCoordinates3 ( self->qcAtoms, coordinates3, self->mmcoordinates3 ) ;
        Coordinates3_Set ( self->mmgradients3, 0.0e+00 ) ;
    }
}

/*------------------------------------------------------------------------------
! . Setup the state.
!-----------------------------------------------------------------------------*/
NBModelOrcaState *NBModelORCAState_Setup ( MMAtomContainer *mmAtoms, QCAtomContainer *qcAtoms, Selection *fixedatoms, PairList *exclusions, PairList *interactions14,
                                                                                                                LJParameterContainer *ljParameters, LJParameterContainer *ljParameters14,
                                                                                                                                                          const QCMMLinkAtomCoupling qcmmcoupling )
{
    NBModelOrcaState *self = NULL ;
    if ( mmAtoms != NULL )
    {
        auto Boolean QOK = True ;
        auto Integer i, n ;
        n = mmAtoms->natoms ;
        self = NBModelORCAState_Allocate ( n ) ;
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

            /* . There are QC atoms. */
            if ( qcAtoms != NULL )
            {
                /* . Set MM flags and modify boundary atom data. */
                for ( n = 0 ; n < qcAtoms->natoms ; n++ )
                {
                    i = qcAtoms->data[n].index ;
                    if ( ! qcAtoms->data[n].QBOUNDARY ) self->QMM[i] = False ; /* . Pure QC atom which is not MM. */
                }

                /* . Define MM atom charges. */
                self->mmCharges = QCAtomContainer_GetMMCharges ( qcAtoms, mmAtoms, ( qcmmcoupling == QCMMLinkAtomCoupling_RC ) ) ;
                if ( self->mmCharges != NULL )
                {
                    self->mmcoordinates3   = Coordinates3_Allocate ( self->mmCharges->length ) ;
                    self->mmgradients3     = Coordinates3_Allocate   ( self->mmCharges->length ) ;
                    self->mmexclusions     = QCAtomContainer_MakeMMPairList ( qcAtoms, mmAtoms, exclusions     ) ;
                    self->mminteractions14 = QCAtomContainer_MakeMMPairList ( qcAtoms, mmAtoms, interactions14 ) ;
                    SelfPairList_MakeConnections ( self->mmexclusions,     self->mmCharges->length ) ;
                    SelfPairList_MakeConnections ( self->mminteractions14, self->mmCharges->length ) ;
                }
                QOK = QOK && ( self->mmCharges != NULL ) && ( self->mmcoordinates3 != NULL ) && ( self->mmgradients3 != NULL ) &&
                                                          ( ( self->mmexclusions     != NULL ) || ( exclusions     == NULL ) ) &&
                                                          ( ( self->mminteractions14 != NULL ) || ( interactions14 == NULL ) ) ;
            }

            /* . Assign aliases. */
            self->exclusions     = exclusions     ;
            self->interactions14 = interactions14 ;
            self->ljParameters   = ljParameters   ;
            self->ljParameters14 = ljParameters14 ;
            self->mmAtoms        = mmAtoms        ;
            self->qcAtoms        = qcAtoms        ;

            /* . Check that everything is OK. */
            if ( ! QOK ) NBModelORCAState_Deallocate ( &self ) ;
        }
    }
    return self ;
}
