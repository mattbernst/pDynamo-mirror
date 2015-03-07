/*------------------------------------------------------------------------------
! . File      : NBModelMonteCarloState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . This module implements Monte Carlo non-bonding interaction state procedures.
!=============================================================================*/

# include "math.h"
# include "Memory.h"
# include "NBModelMonteCarloState.h"

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
NBModelMonteCarloState *NBModelMonteCarloState_Allocate ( const Integer n )
{
    NBModelMonteCarloState *self = NULL ;
    if ( n > 0 )
    {
        self = ( NBModelMonteCarloState * ) Memory_Allocate ( sizeof ( NBModelMonteCarloState ) ) ;
        if ( self != NULL )
        {
            auto Boolean QOK ;
            /* . Dimension. */
            self->nisolates          = n       ;
            /* . Scaling. */
            self->isolatescale       =     - 1 ;
            self->chargeScale        = 0.0e+00 ;
            self->epsilonScale       = 0.0e+00 ;
            self->sigmaScale         = 0.0e+00 ;
            /* . Energies. */
            self->efmmel             = 0.0e+00 ;
            self->e1mmel             = 0.0e+00 ;
            self->efmmlj             = 0.0e+00 ;
            self->e1mmlj             = 0.0e+00 ;
            /* . Aliases. */
            self->coordinates3       = NULL    ;
            self->isolates           = NULL    ;
            self->ljParameters       = NULL    ;
            self->mmAtoms            = NULL    ;
            self->symmetryParameters = NULL    ;
            /* . Arrays to be allocated. */
            self->centers = ( Vector3 ** ) Memory_Allocate_Array ( n, sizeof ( Vector3 * ) ) ;
            self->QFREE   = Memory_Allocate_Array_Boolean_Initialize ( n, True  ) ;
            QOK = ( self->centers != NULL ) && ( self->QFREE != NULL ) ;
            if ( QOK )
            {
                auto Integer i ;
                for ( i = 0 ; i < n ; i++ )
                {
                    self->centers[i] = Vector3_Allocate ( ) ;
                    QOK = QOK && ( self->centers[i] != NULL ) ;
                }
            }
            /* . Deallocate if there is not enough memory. */
            if ( ! QOK ) NBModelMonteCarloState_Deallocate ( &self ) ;
        }
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void NBModelMonteCarloState_Deallocate ( NBModelMonteCarloState **self )
{
    if ( (*self) != NULL )
    {
        /* . Allocated arrays. */
        if ( (*self)->centers != NULL )
        {
            auto Integer i ;
	    for ( i = 0 ; i < (*self)->nisolates ; i++ ) Vector3_Deallocate ( &((*self)->centers[i]) ) ;
        }
        free ( (*self)->centers ) ;
        Memory_Deallocate_Boolean ( &((*self)->QFREE) ) ;
        /* . Object. */
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Initialization for an energy/gradient calculation.
!-----------------------------------------------------------------------------*/
void NBModelMonteCarloState_Initialize ( NBModelMonteCarloState  *self, Coordinates3 *coordinates3, SymmetryParameters *symmetryParameters )
{
    if ( self != NULL )
    {
        /* . Energies. */
        self->efmmel = 0.0e+00 ;
        self->e1mmel = 0.0e+00 ;
        self->efmmlj = 0.0e+00 ;
        self->e1mmlj = 0.0e+00 ;
        /* . Coordinates. */
        self->coordinates3       = coordinates3       ;
        self->symmetryParameters = symmetryParameters ;
    }
}

/*------------------------------------------------------------------------------
! . Scale the interaction parameters for an isolate.
! . Only one isolate is allowed to be scaled at a time.
!-----------------------------------------------------------------------------*/
void NBModelMonteCarloState_ScaleIsolateInteractionParameters ( NBModelMonteCarloState *self, const Integer isolate, const Real chargeScale, const Real epsilonScale, const Real sigmaScale )
{
    if ( self != NULL )
    {
        self->isolatescale = isolate     ;
        self->chargeScale  = chargeScale                    ; /* . Scaling for qi * qj.           */
        self->epsilonScale = sqrt ( fabs ( epsilonScale ) ) ; /* . Scaling for sqrt ( ei * ej ).  */
        self->sigmaScale   = pow ( sigmaScale, 3.0 )        ; /* . Scaling for ( si**3 * sj**3 ). */
    }
}

/*------------------------------------------------------------------------------
! . Setup the state.
!-----------------------------------------------------------------------------*/
NBModelMonteCarloState *NBModelMonteCarloState_Setup ( SelectionContainer *isolates, MMAtomContainer *mmAtoms, LJParameterContainer *ljParameters, Selection *fixedatoms )
{
    NBModelMonteCarloState *self = NULL ;
    if ( ( isolates != NULL ) && ( mmAtoms != NULL ) && ( ljParameters != NULL ) )
    {
        auto Integer  n ;

        /* . Basic allocation. */
        n = isolates->nitems ;
        self = NBModelMonteCarloState_Allocate ( n ) ;

        /* . Remaining allocation. */
        if ( self != NULL )
        {
            /* . Fixed atoms - if one atom in an isolate is fixed the whole isolate is fixed! */
            if ( fixedatoms != NULL )
            {
                auto Integer             i, s ;
                auto Selection *iselection ;
                Selection_MakeFlags ( fixedatoms, mmAtoms->natoms ) ;
                for ( i = 0 ; i < n ; i++ )
                {
                    iselection = isolates->items[i] ;
                    for ( s = 0 ; s < iselection->nindices ; s++ )
                    {
                        if ( fixedatoms->flags[iselection->indices[s]] )
                        {
                            self->QFREE[i] = False ;
                            break ;
                        }
                    }
                }
            }

            /* . Assign aliases. */
            self->isolates     = isolates     ;
            self->ljParameters = ljParameters ;
            self->mmAtoms      = mmAtoms      ;
        }
    }
    return self ;
}
