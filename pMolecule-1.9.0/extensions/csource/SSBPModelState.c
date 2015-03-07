/*------------------------------------------------------------------------------
! . File      : SSBPModelState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . State procedures for the SSBP boundary model.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "Memory.h"
# include "SSBPModelState.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
SSBPModelState *SSBPModelState_Allocate ( void )
{
    SSBPModelState *self = NULL ;
    self = ( SSBPModelState * ) Memory_Allocate ( sizeof ( SSBPModelState ) ) ;
    if ( self != NULL )
    {
        self->atomOutsideCavity    = False   ;
        self->numberOfAtoms        =  0      ;
        self->particleIndex        = -1      ;
        self->eAngular             = 0.0e+00 ;
        self->eCavity              = 0.0e+00 ;
        self->eEmpiricalCorrection = 0.0e+00 ;
        self->eHardSphere          = 0.0e+00 ;
        self->eKirkwood            = 0.0e+00 ;
        self->eKirkwoodCheck       = 0.0e+00 ;
        self->eTotal               = 0.0e+00 ;
        self->qTotal               = 0.0e+00 ;
        self->radius               = 0.0e+00 ;
        /* . Aliases. */
        self->coordinates3         = NULL ;
        self->gradients3           = NULL ;
        self->mmAtoms              = NULL ;
	self->qcAtoms              = NULL ;
	self->qcCharges            = NULL ;
	self->qcmmPotentials       = NULL ;
	self->qcqcPotentials       = NULL ;
        /* . Arrays. */
        self->waterAtomIndices     = NULL ;
        self->radii                = NULL ;
        self->cavitySelection      = NULL ;
        self->radiusSelection      = NULL ;
        /* . Kirkwood - MM. */
        self->dQLMi                = NULL ;
        self->dQLMr                = NULL ;
        self->dQLMs                = NULL ;
        self->a                    = NULL ;
        self->dA                   = NULL ;
        self->factorial            = NULL ;
        self->tvar12               = NULL ;
        self->rlr3                 = NULL ;
        self->uS                   = NULL ;
        self->vS                   = NULL ;
        self->xr2                  = NULL ;
        self->xrl                  = NULL ;
        self->yr2                  = NULL ;
        self->yrl                  = NULL ;
        self->zr2                  = NULL ;
        self->zxr3                 = NULL ;
        self->zyr3                 = NULL ;
        self->comr                 = NULL ;
        self->comi                 = NULL ;
        self->rQs                  = NULL ;
        self->tvar                 = NULL ;
        /* . Kirkwood - QC. */
        self->dQLMiQ               = NULL ;
        self->dQLMrQ               = NULL ;
        self->dQLMsQ               = NULL ;
	self->qcGradients3         = NULL ;
        self->tvar12Q              = NULL ;
        self->rlr3Q                = NULL ;
        self->uSQ                  = NULL ;
        self->vSQ                  = NULL ;
        self->xr2Q                 = NULL ;
        self->xrlQ                 = NULL ;
        self->yr2Q                 = NULL ;
        self->yrlQ                 = NULL ;
        self->zr2Q                 = NULL ;
        self->zxr3Q                = NULL ;
        self->zyr3Q                = NULL ;
        self->comrQ                = NULL ;
        self->comiQ                = NULL ;
        self->rQsQ                 = NULL ;
        self->tvarQ                = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SSBPModelState_Deallocate ( SSBPModelState **self )
{
    if ( (*self) != NULL )
    {
# ifdef SSBPTESTQC
{
    QCAtomContainer_Deallocate ( &((*self)->qcAtoms       ) ) ;
    Real1DArray_Deallocate     ( &((*self)->qcCharges     ) ) ;
    Real1DArray_Deallocate     ( &((*self)->qcmmPotentials) ) ;
    SymmetricMatrix_Deallocate ( &((*self)->qcqcPotentials) ) ;
}
# endif
        Integer2DArray_Deallocate ( &((*self)->waterAtomIndices) ) ;
        Real1DArray_Deallocate    ( &((*self)->radii           ) ) ;
        Selection_Deallocate      ( &((*self)->cavitySelection ) ) ;
        Selection_Deallocate      ( &((*self)->radiusSelection ) ) ;
        /* . Kirkwood - MM. */
        Coordinates3_Deallocate       ( &((*self)->dQLMi           ) ) ;
        Coordinates3_Deallocate       ( &((*self)->dQLMr           ) ) ;
        Coordinates3_Deallocate       ( &((*self)->dQLMs           ) ) ;
        Real1DArray_Deallocate    ( &((*self)->a               ) ) ;
        Real1DArray_Deallocate    ( &((*self)->dA              ) ) ;
        Real1DArray_Deallocate    ( &((*self)->factorial       ) ) ;
        Real1DArray_Deallocate    ( &((*self)->tvar12          ) ) ;
        Real1DArray_Deallocate    ( &((*self)->rlr3            ) ) ;
        Real1DArray_Deallocate    ( &((*self)->uS              ) ) ;
        Real1DArray_Deallocate    ( &((*self)->vS              ) ) ;
        Real1DArray_Deallocate    ( &((*self)->xr2             ) ) ;
        Real1DArray_Deallocate    ( &((*self)->xrl             ) ) ;
        Real1DArray_Deallocate    ( &((*self)->yr2             ) ) ;
        Real1DArray_Deallocate    ( &((*self)->yrl             ) ) ;
        Real1DArray_Deallocate    ( &((*self)->zr2             ) ) ;
        Real1DArray_Deallocate    ( &((*self)->zxr3            ) ) ;
        Real1DArray_Deallocate    ( &((*self)->zyr3            ) ) ;
        Real2DArray_Deallocate    ( &((*self)->comr            ) ) ;
        Real2DArray_Deallocate    ( &((*self)->comi            ) ) ;
        Real2DArray_Deallocate    ( &((*self)->rQs             ) ) ;
        Real2DArray_Deallocate    ( &((*self)->tvar            ) ) ;
        /* . Kirkwood - QC. */
        Coordinates3_Deallocate       ( &((*self)->dQLMiQ          ) ) ;
        Coordinates3_Deallocate       ( &((*self)->dQLMrQ          ) ) ;
        Coordinates3_Deallocate       ( &((*self)->dQLMsQ          ) ) ;
        Coordinates3_Deallocate       ( &((*self)->qcGradients3    ) ) ;
        Real1DArray_Deallocate    ( &((*self)->tvar12Q         ) ) ;
        Real1DArray_Deallocate    ( &((*self)->rlr3Q           ) ) ;
        Real1DArray_Deallocate    ( &((*self)->uSQ             ) ) ;
        Real1DArray_Deallocate    ( &((*self)->vSQ             ) ) ;
        Real1DArray_Deallocate    ( &((*self)->xr2Q            ) ) ;
        Real1DArray_Deallocate    ( &((*self)->xrlQ            ) ) ;
        Real1DArray_Deallocate    ( &((*self)->yr2Q            ) ) ;
        Real1DArray_Deallocate    ( &((*self)->yrlQ            ) ) ;
        Real1DArray_Deallocate    ( &((*self)->zr2Q            ) ) ;
        Real1DArray_Deallocate    ( &((*self)->zxr3Q           ) ) ;
        Real1DArray_Deallocate    ( &((*self)->zyr3Q           ) ) ;
        Real2DArray_Deallocate    ( &((*self)->comrQ           ) ) ;
        Real2DArray_Deallocate    ( &((*self)->comiQ           ) ) ;
        Real2DArray_Deallocate    ( &((*self)->rQsQ            ) ) ;
        Real2DArray_Deallocate    ( &((*self)->tvarQ           ) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! .  Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SSBPModelState_Initialize ( SSBPModelState *self, Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    if ( self != NULL )
    {
        /* . Scalars. */
        self->atomOutsideCavity    = False   ;
        self->particleIndex        = -1      ;
        self->eAngular             = 0.0e+00 ;
        self->eCavity              = 0.0e+00 ;
        self->eEmpiricalCorrection = 0.0e+00 ;
        self->eHardSphere          = 0.0e+00 ;
        self->eKirkwood            = 0.0e+00 ;
	self->eKirkwoodCheck       = 0.0e+00 ;
        self->eTotal               = 0.0e+00 ;
        self->radius               = 0.0e+00 ;
        /* . Arrays. */
        self->coordinates3         = coordinates3 ;
        self->gradients3           = gradients3   ;
	/* . QC. */
	Coordinates3_Set        ( self->qcGradients3  , 0.0e+00 ) ;
	Real1DArray_Set     ( self->qcCharges     , 0.0e+00 ) ;
	Real1DArray_Set     ( self->qcmmPotentials, 0.0e+00 ) ;
        SymmetricMatrix_Set ( self->qcqcPotentials, 0.0e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! .  Set up the state.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . The total charge may not be correct if there are QC atoms and this should be corrected.
! . However, the relevant term is very small if the radius is large and so the effects of any
! . error should be minor.
*/

SSBPModelState *SSBPModelState_Setup ( const Boolean doKirkwood, const Integer maximumL, MMAtomContainer *mmAtoms,
                                                             QCAtomContainer *qcAtoms, Selection *cavitySelection,
                                                     Selection *radiusSelection, Integer2DArray *waterAtomIndices,
                                                              Real1DArray *qcCharges, Real1DArray *qcmmPotentials,
                                                                 SymmetricMatrix *qcqcPotentials, Status *status )
{
    SSBPModelState *self = NULL ;
    if ( mmAtoms != NULL )
    {
        auto Integer numberOfAtoms ;
        auto Status  localStatus ;

        /* . Initialization. */
        Status_Set ( &localStatus, Status_Continue ) ;

        /* . Allocation. */
        self = SSBPModelState_Allocate ( ) ;

        /* . Aliases, general data and arrays. */
        numberOfAtoms          = mmAtoms->natoms  ;
	self->numberOfAtoms    = numberOfAtoms    ;
        self->qTotal           = fabs ( MMAtomContainer_TotalCharge ( mmAtoms ) ) ;
        self->mmAtoms          = mmAtoms        ;
	self->qcAtoms          = qcAtoms        ;
	self->qcCharges        = qcCharges      ;
	self->qcmmPotentials   = qcmmPotentials ;
	self->qcqcPotentials   = qcqcPotentials ;
        self->radii            = Real1DArray_Allocate ( numberOfAtoms, &localStatus ) ;
        self->cavitySelection  = cavitySelection  ;
        self->radiusSelection  = radiusSelection  ;
        self->waterAtomIndices = waterAtomIndices ;

# ifdef SSBPTESTQC
{
if ( ( self->cavitySelection != NULL ) && ( self->cavitySelection->nindices > 0 ) )
{
    auto Integer    i, n, s ;
    auto Selection *qcSelection ;
    qcSelection = self->cavitySelection ;
    n = qcSelection->nindices ;
    self->qcAtoms        = QCAtomContainer_Allocate ( n ) ;
    self->qcCharges      = Real1DArray_Allocate     ( n, NULL ) ;
    self->qcmmPotentials = Real1DArray_Allocate     ( n, NULL ) ; Real1DArray_Set     ( self->qcmmPotentials, 0.0e+00 ) ;
    self->qcqcPotentials = SymmetricMatrix_Allocate ( n )       ; SymmetricMatrix_Set ( self->qcqcPotentials, 0.0e+00 ) ;
    for ( s = 0 ; s < qcSelection->nindices ; s++ )
    {
        i = qcSelection->indices[s] ;
	self->qcAtoms->data[s].index = i ;
        Real1DArray_Item ( self->qcCharges, s ) = mmAtoms->data[i].charge ;
    }
}
# endif

        /* . Kirkwood terms. */
        if ( doKirkwood )
        {
	    auto Integer n ;

# ifdef SSBPTESTQC
{
if ( self->qcAtoms != NULL )
{
    auto Integer i, index ;
    for ( i = 0 ; i < self->qcAtoms->natoms ; i++ )
    {
        index = self->qcAtoms->data[i].index ;
        self->mmAtoms->data[index].QACTIVE = False ;
    }
}
# endif

            /* . Find number of active MM atoms. */
	    n = MMAtomContainer_NumberOfActiveAtoms ( self->mmAtoms ) ;

# ifdef SSBPTESTQC
if ( self->qcAtoms != NULL )
{
    auto Integer i, index ;
    for ( i = 0 ; i < self->qcAtoms->natoms ; i++ )
    {
        index = self->qcAtoms->data[i].index ;
        self->mmAtoms->data[index].QACTIVE = True ;
    }
}
# endif
            /* . Allocate arrays. */
            self->a  = Real1DArray_Allocate ( maximumL+1, &localStatus ) ;
            self->dA = Real1DArray_Allocate ( maximumL+1, &localStatus ) ;
 	    if ( n > 0 )
	    {
        	self->dQLMi  = Coordinates3_Allocate    ( n ) ;
        	self->dQLMr  = Coordinates3_Allocate    ( n ) ;
        	self->dQLMs  = Coordinates3_Allocate    ( n ) ;
        	self->tvar12 = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->rlr3   = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->uS     = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->vS     = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->xr2    = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->xrl    = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->yr2    = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->yrl    = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->zr2    = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->zxr3   = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->zyr3   = Real1DArray_Allocate (             n, &localStatus ) ;
        	self->comr   = Real2DArray_Allocate ( maximumL+1, n, &localStatus ) ;
        	self->comi   = Real2DArray_Allocate ( maximumL+1, n, &localStatus ) ;
        	self->rQs    = Real2DArray_Allocate ( maximumL+1, n, &localStatus ) ;
        	self->tvar   = Real2DArray_Allocate ( maximumL+1, n, &localStatus ) ;
	    }

            /* . Factorial. */
            self->factorial = Real1DArray_Allocate ( 2*maximumL+1, &localStatus ) ;
            if ( self->factorial != NULL )
            {
                auto Integer i ;
                Real1DArray_Item ( self->factorial, 0 ) = 1.0e+00 ;
                for ( i = 1 ; i < 2*maximumL+1 ; i++ )
		{
		    Real1DArray_Item ( self->factorial, i ) = ( ( Real ) i ) * Real1DArray_Item ( self->factorial, i-1 ) ;
		}
            }

            /* . QC. */
	    n = QCAtomContainer_Size ( self->qcAtoms ) ;
	    if ( n > 0 )
	    {
		self->dQLMiQ       = Coordinates3_Allocate    ( n ) ;
        	self->dQLMrQ       = Coordinates3_Allocate    ( n ) ;
        	self->dQLMsQ       = Coordinates3_Allocate    ( n ) ;
        	self->qcGradients3 = Coordinates3_Allocate    ( n ) ;
        	self->tvar12Q      = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->rlr3Q        = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->uSQ          = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->vSQ          = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->xr2Q         = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->xrlQ         = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->yr2Q         = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->yrlQ         = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->zr2Q         = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->zxr3Q        = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->zyr3Q        = Real1DArray_Allocate (		n, &localStatus ) ;
        	self->comrQ        = Real2DArray_Allocate ( maximumL+1, n, &localStatus ) ;
        	self->comiQ        = Real2DArray_Allocate ( maximumL+1, n, &localStatus ) ;
        	self->rQsQ         = Real2DArray_Allocate ( maximumL+1, n, &localStatus ) ;
        	self->tvarQ        = Real2DArray_Allocate ( maximumL+1, n, &localStatus ) ;
            }
        }

        /* . Finish up. */
        Status_Set ( status, localStatus ) ;
        if ( ! Status_OK ( &localStatus ) ) SSBPModelState_Deallocate ( &self ) ;
    }
    return self ;
}
