/*------------------------------------------------------------------------------
! . File      : MNDOCIModel.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . A basic MNDO CI module.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "DefineStatements.h"
# include "Integer2DArray.h"
# include "Memory.h"
# include "MNDOCIModel.h"
# include "MNDOIntegrals.h"
# include "QCOnePDM.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Default_DegeneracyTolerance          1.0e-3 /* . Hartrees. */
# define Default_FractionalOccupancyTolerance 1.0e-6 /* . Dimensionless. */
# define Default_NumberOfStates               100
# define Default_SpinTolerance                1.0e-2 /* . Dimensionless. */

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOCIModel *MNDOCIModel_Allocate ( Status *status )
{
    MNDOCIModel *self = NULL ;
    MEMORY_ALLOCATE ( self, MNDOCIModel ) ;
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    else
    {
        self->checkAlgorithm               = False ;
        self->doAllStates                  = False ;
        self->identifyRootSpin             = False ;
        self->localizeOrbitals             = False ;
        self->activeElectrons              = 0     ;
        self->activeOrbitals               = 0     ;
        self->localizeStart                = -1    ;
        self->localizeStop                 = -1    ;
        self->minimalMultiplicity          = -1    ;
        self->numberOfStates               = Default_NumberOfStates ;
        self->requiredRoot                 = 0     ;
        self->rootMultiplicity             = -1    ;
        self->degeneracyTolerance          = Default_DegeneracyTolerance          ;
        self->fractionalOccupancyTolerance = Default_FractionalOccupancyTolerance ;
        self->spinTolerance                = Default_SpinTolerance                ;
        self->algorithm                    = MNDOCIAlgorithm_Full                 ;
        self->method                       = MNDOCIMethod_Full                    ;
        self->microstates                  = NULL  ;
        /* . Solver. */
        JDEigenvalueSolver_Initialize ( &(self->eigenvalueSolver) ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOCIModel *MNDOCIModel_Clone ( const MNDOCIModel *self, Status *status )
{
    MNDOCIModel *new = NULL ;
    if ( self != NULL )
    {
        new = MNDOCIModel_Allocate ( status ) ;
        if ( new != NULL )
        {
            new->checkAlgorithm               = self->checkAlgorithm               ;
            new->doAllStates                  = self->doAllStates                  ;
            new->identifyRootSpin             = self->identifyRootSpin             ;
            new->localizeOrbitals             = self->localizeOrbitals             ;
            new->activeElectrons              = self->activeElectrons              ;
            new->activeOrbitals               = self->activeOrbitals               ;
            new->localizeStart                = self->localizeStart                ;
            new->localizeStop                 = self->localizeStop                 ;
            new->minimalMultiplicity          = self->minimalMultiplicity          ;
            new->numberOfStates               = self->numberOfStates               ;
            new->requiredRoot                 = self->requiredRoot                 ;
            new->rootMultiplicity             = self->rootMultiplicity             ;
            new->degeneracyTolerance          = self->degeneracyTolerance          ;
            new->fractionalOccupancyTolerance = self->fractionalOccupancyTolerance ;
            new->spinTolerance                = self->spinTolerance                ;
            new->algorithm                    = self->algorithm                    ;
            new->method                       = self->method                       ;
            new->microstates                  = Integer2DArray_Clone ( self->microstates, status ) ;
            /* . Solver. */
            JDEigenvalueSolver_CopyTo ( &(self->eigenvalueSolver), &(new->eigenvalueSolver) ) ;
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIModel_Deallocate ( MNDOCIModel **self )
{
    if ( (*self) != NULL )
    {
        Integer2DArray_Deallocate ( &((*self)->microstates) ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
        (*self) = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the CI energy.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIModel_Energy ( const MNDOCIModel *self, MNDOCIState *state, const Real eelectronic, const Real enuclear, QCOnePDM *densityp, QCOnePDM *densityq, Status *status )
{
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        auto Boolean     isOK ;
        auto Integer     astart, i, iroot, j, k, m, root, vstart ;
        auto Real        sfactor ;
        auto Real1DArray cslice ;

# ifdef DEBUGMNDOCI
{
    printf ( "\nMNDOCIModel_Energy:\n" ) ;
    printf ( "\nCore occupancies:\n" ) ;
    Real1DArray_Print ( state->ocore ) ;
    printf ( "\nOrbital energies:\n" ) ;
    Real1DArray_Print ( state->energies ) ;
    printf ( "\nOrbitals:\n" ) ;
    Real2DArray_Print ( state->orbitals ) ;
    printf ( "\nActive Orbitals:\n" ) ;
    Real2DArray_Print ( &(state->activemos) ) ;
}
# endif

        /* . Get some indices. */
        astart = state->ncore ;
        vstart = state->ncore + state->nactive ;

        /* . Check for degenerate orbital energies. */
        isOK = True ;
        if ( astart > 0                ) isOK = isOK && ( fabs ( Real1DArray_Item ( state->energies, astart ) - Real1DArray_Item ( state->energies, astart-1 ) ) > self->degeneracyTolerance ) ;
        if ( vstart < state->norbitals ) isOK = isOK && ( fabs ( Real1DArray_Item ( state->energies, vstart ) - Real1DArray_Item ( state->energies, vstart-1 ) ) > self->degeneracyTolerance ) ;
        state->orbitalDegeneracies = ! isOK ;

        /* . Check for inactive orbitals with fractional occupancies. */
        isOK = True ;
        if ( astart > 0                ) isOK = isOK && ( fabs ( Real1DArray_Item ( state->occupancies, astart - 1 ) - 2.0e+00 ) < self->fractionalOccupancyTolerance ) ;
        if ( vstart < state->norbitals ) isOK = isOK && ( fabs ( Real1DArray_Item ( state->occupancies, vstart     )           ) < self->fractionalOccupancyTolerance ) ;
        state->fractionallyOccupiedInactive = ! isOK ;

        /* . Get the transformed integrals. */
        MNDOCIState_FourIndexTransformation ( state ) ;

# ifdef DEBUGMNDOCI
{
    printf ( "\nHalf-transformed integrals:\n" ) ;
    Real2DArray_Print ( state->motei34 ) ;
    printf ( "\nTransformed integrals:\n" ) ;
    DoubleSymmetricMatrix_Print ( state->moteis ) ;
}
# endif

        /* . Make the frozen-core density matrix. */
        QCOnePDM_MakeElementary ( state->pcore, state->ncore, state->ocore, state->orbitals, status ) ;

        /* . Calculate the core Fock matrix. */
        MNDOIntegrals_MakeFockG   ( state->twoelectronintegrals, state->pcore, state->fcore, NULL, NULL ) ;
        SymmetricMatrix_Increment ( state->fcore, state->oneelectronmatrix ) ;

        /* . Calculate the core energy. */
        state->ecore = 0.5e+00 * ( SymmetricMatrix_Multiply2_Trace ( state->pcore, state->oneelectronmatrix ) +
                                   SymmetricMatrix_Multiply2_Trace ( state->pcore, state->fcore             ) ) ;

        /* . Transform the core Fock matrix to the active MO basis. */
        SymmetricMatrix_Transform ( state->fcore, &(state->activemos), False, state->fcoreMO ) ;

        /* . Calculate the CI matrix. */
        MNDOCIState_MakeCIMatrix ( state ) ;

# ifdef DEBUGMNDOCI
{
    printf ( "\nCore energy = %.6f\n", state->ecore ) ;
    printf ( "\nPcore:\n" ) ;
    SymmetricMatrix_Print ( state->pcore ) ;
    printf ( "\nFcore:\n" ) ;
    SymmetricMatrix_Print ( state->fcore ) ;
    printf ( "\nFcore in MO basis:\n" ) ;
    SymmetricMatrix_Print ( state->fcoreMO ) ;
    printf ( "\nCI matrix:\n" ) ;
    SymmetricMatrix_Print ( state->ciMatrixFull ) ;
}
# endif

        /* . Diagonalization. */
        if ( self->algorithm == MNDOCIAlgorithm_Full ) MNDOCIState_DiagonalizeCIMatrix ( state, NULL                     , status ) ;
        else                                           MNDOCIState_DiagonalizeCIMatrix ( state, &(self->eigenvalueSolver), status ) ;

# ifdef DEBUGMNDOCI
{
    printf ( "\nCI energies:\n" ) ;
    Real1DArray_Print ( state->ciEnergies ) ;
    printf ( "\nCI vectors:\n" ) ;
    Real2DArray_Print ( state->ciVectors ) ;
    printf ( "\nSPQRs:\n" ) ;
    for ( j = 0 ; j < state->nconfigurations ; j++ )
    {
        if ( state->configurations[j].nspqr > 0 )
        {
            printf ( "State %d %d :", j, state->configurations[j].nspqr ) ;
            for ( k = 0 ; k < state->configurations[j].nspqr ; k++ ) printf ( " %d", Integer1DArray_Item ( state->configurations[j].spqr, k ) - 1 ) ;
            printf ( "\n" ) ;
        }
    }
}
# endif

        /* . Determine the spins of the configurations and find the required root (either of specified or unspecified multiplicity). */
        iroot = -1 ;
        root  = -1 ;
        Real1DArray_Set ( state->spins, 0.0e+00 ) ;
        for ( i = 0 ; i < state->numberOfStates ; i++ )
        {
            sfactor = 0.5e+00 * ( Real ) state->nelectrons ;
            for ( j = 0 ; j < state->nconfigurations ; j++ )
            {
                sfactor -= Real2DArray_Item ( state->ciVectors, i, j ) * Real2DArray_Item ( state->ciVectors, i, j ) * state->configurations[j].spin / 4.0e+00 ;
                for ( k = 0 ; k < state->configurations[j].nspqr ; k++ )
                {
                    m = Integer1DArray_Item ( state->configurations[j].spqr, k ) ;
                    if ( m < 0 ) sfactor -= Real2DArray_Item ( state->ciVectors, i, j ) * Real2DArray_Item ( state->ciVectors, i, -m-1 ) * 2.0e+00 ;
                    else         sfactor += Real2DArray_Item ( state->ciVectors, i, j ) * Real2DArray_Item ( state->ciVectors, i,  m-1 ) * 2.0e+00 ;
                }
            }
            if ( self->identifyRootSpin && ( fabs ( state->requiredSpin - sfactor ) < self->spinTolerance ) )
            {
                iroot++ ;
                if ( iroot == self->requiredRoot ) root = i ;
            }
            Real1DArray_Item ( state->spins, i ) = sfactor ;
        }
        if ( ! self->identifyRootSpin ) root = self->requiredRoot ;

        /* . Check that a root has been found. */
        state->rootNotFound = ( root < 0 ) ;
        if ( state->rootNotFound ) root = 0 ;

        /* . Save the energies, etc. */
        state->root       = root ;
        state->rootenergy = Real1DArray_Item ( state->ciEnergies, root ) ;
        Real2DArray_RowSlice ( state->ciVectors, root, &cslice, NULL ) ;
        Real1DArray_CopyTo   ( &cslice, state->ciVector, status ) ;

# ifdef DEBUGMNDOCI
{
    if ( ! state->rootNotFound )
    {
        printf ( "\nCI Root, Energy, Spin = %d %.8f %.6f\n", state->root, state->rootenergy, Real1DArray_Item ( state->spins, state->root ) ) ;
        printf ( "\nCI vector:\n" ) ;
        Real1DArray_Print ( state->ciVector ) ;
    }
}
# endif

        /* . Calculate the final energy. */
        state->baseline = state->ecore - eelectronic ;
        state->ciEnergy = state->ecore + enuclear + state->rootenergy ;

        /* . Make the CI densities. */
        MNDOCIState_MakeDensities ( state, densityp, densityq ) ;

        /* . Calculate the Z-matrix for a gradient calculation. */
        MNDOCIState_CalculateZMatrix ( state ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a CI state.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOCIState *MNDOCIModel_MakeState ( const MNDOCIModel *self, const Integer multiplicity, const Integer nalpha, const Integer nbeta, const Integer norbitals, Status *status )
{
    Integer      localizeStart = -1, localizeStop = 0, mlocal, nactive, nclosed, ncore, ncoreelectrons, nelectrons, nopen, nvirtual ;
    MNDOCIState *state = NULL  ;
    Real         requiredSpin ;
    Status       localstatus  ;

    /* . Initialization. */
    Status_Set ( &localstatus, Status_Continue ) ;

    /* . Get counters that are consistent with the model. */
    ncoreelectrons = ( nalpha + nbeta ) - self->activeElectrons ;
    ncore          = ncoreelectrons / 2 ;
    nelectrons     = ( nalpha + nbeta ) - ( 2 * ncore ) ;
    nopen          = nalpha - nbeta ;
    nclosed        = ( self->activeElectrons - nopen ) / 2 ;
    nvirtual       = Minimum ( self->activeOrbitals, norbitals ) - ( nclosed + nopen ) ;
    nactive        = nclosed + nopen + nvirtual ;

    /* . Localization options - orbital range. */
    if ( self->localizeOrbitals )
    {
        if ( self->localizeStart < 0 ) localizeStart = ncore ;
        else                           localizeStart = Minimum ( self->localizeStart, ncore + nactive - 1 ) ;
        if ( self->localizeStop  < 0 ) localizeStop  = ncore + nactive ;
        else                           localizeStop  = Minimum ( self->localizeStop , ncore + nactive     ) ;
    }

    /* . Check for incompatibilities. */
    /* . The number of electrons must be as specified but the number of active orbitals can be less than as specified if there are not enough virtual orbitals. */
    if ( ( self->activeElectrons            > 2 * self->activeOrbitals ) ||
         ( self->activeElectrons            > ( nalpha + nbeta )       ) ||
         ( self->activeElectrons            != nelectrons              ) ||
         ( self->activeOrbitals             < ( nclosed + nopen      ) ) ||
         ( ncoreelectrons                   <  0                       ) ||
         ( ncoreelectrons                   != 2 * ncore               ) ||
         ( nclosed                          <  0                       ) ||
         ( nvirtual                         <  0                       ) ||
         ( ( localizeStop - localizeStart ) <= 0                       ) )
    {
        printf ( "\nMNDOCIModel_MakeState> CI model and system incompatible.\n" ) ;
        Status_Set ( &localstatus, Status_InvalidArgument ) ;
        goto FinishUp ;
    }

    /* . Spin quantities . */
    mlocal = self->rootMultiplicity ;
    if ( mlocal < 1 ) mlocal = multiplicity ;

    /* . 0, 0.75, 2, 3.75, 6, 8.75, 12, 15.75, 20. */
    requiredSpin = 0.0e+00 ;
    if ( self->identifyRootSpin ) requiredSpin = 0.25e+00 * ( Real ) ( ( mlocal - 1 ) * ( mlocal + 1 ) ) ;

    /* . Generate the configurations and associated data. */
    if ( self->method == MNDOCIMethod_Full )
    {
        auto Integer m = 0, ndown, nup ;
        if ( self->minimalMultiplicity > 0 ) m = ( Minimum ( mlocal, self->minimalMultiplicity ) - 1 ) / 2 ;
        nup   = ( nelectrons + 1 ) / 2 + m ;
        ndown = nelectrons - nup ;
        if ( ( nup < 0 ) || ( ndown < 0 ) ) { Status_Set ( &localstatus, Status_InvalidArgument ) ; goto FinishUp ; }
        state = MNDOCIState_MakeFull ( nactive, nup, ndown, &localstatus ) ;
    }
    else if ( self->method == MNDOCIMethod_UserSpecified )
    {
        state = MNDOCIState_MakeUserSpecified ( self->microstates, nactive, nelectrons, &localstatus ) ;
    }
    else
    {
        state = MNDOCIState_MakeSinglesDoubles ( ( self->method == MNDOCIMethod_Singles ) || ( self->method == MNDOCIMethod_SinglesDoubles ),
                                                 ( self->method == MNDOCIMethod_Doubles ) || ( self->method == MNDOCIMethod_SinglesDoubles ),
                                                                                 self->doAllStates, nactive, nclosed, nopen, &localstatus ) ;
    }

    /* . Finish making the state.*/
    if ( Status_OK ( &localstatus ) )
    {
        /* . Set the sparsity. */
        MNDOCIState_GetCIMatrixSparsity ( state ) ;

        /* . Set the number of states. */
        if ( self->numberOfStates <= 0 ) state->numberOfStates =                                 state->nconfigurations   ;
        else                             state->numberOfStates = Minimum ( self->numberOfStates, state->nconfigurations ) ;

        /* . Check that the requested root can exist. */
        if ( state->numberOfStates <= self->requiredRoot ) { Status_Set ( &localstatus, Status_InvalidArgument ) ; goto FinishUp ; }

        /* . State data. */
        state->doFull             = ( self->algorithm == MNDOCIAlgorithm_Full ) || ( self->checkAlgorithm ) ;
        state->doSparse           = ( self->algorithm == MNDOCIAlgorithm_Sparse ) ;
        state->usePreconditioning = self->eigenvalueSolver.usePreconditioning ;
        state->ncore              = ncore        ;
        state->nelectrons         = nelectrons   ;
        state->norbitals          = norbitals    ;
        state->requiredSpin       = requiredSpin ;

        /* . Localization options. */
        if ( self->localizeOrbitals )
        {
            state->localizeStart = localizeStart ;
            state->localizeStop  = localizeStop  ;
        }
    }

# ifdef DEBUGMNDOCI
if ( state == NULL ) printf ( "\nMNDOCIModel_MakeState> No state defined.\n" ) ;
else
{
    printf ( "\nMNDOCIModel_MakeState:\n" ) ;
    printf ( "Number of Active Orbitals  = %d\n", state->nactive         ) ;
    printf ( "Number of Configurations   = %d\n", state->nconfigurations ) ;
    printf ( "Number of States           = %d\n", state->numberOfStates  ) ;
    printf ( "Number of Active Electrons = %d\n", state->nelectrons      ) ;
    printf ( "Number of Alpha  Electrons = %d\n", nalpha                 ) ;
    printf ( "Number of Beta   Electrons = %d\n", nbeta                  ) ;
    printf ( "Number of Orbitals         = %d\n", norbitals              ) ;
    printf ( "Root Multiplicity          = %d\n", mlocal                 ) ;
}
# endif

    /* . Finish up. */
FinishUp:

    if ( ! Status_OK ( &localstatus ) ) MNDOCIState_Deallocate ( &state ) ;
    Status_Set ( status, localstatus ) ;

    return state ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
