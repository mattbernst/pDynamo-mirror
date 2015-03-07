/*------------------------------------------------------------------------------
! . File      : QCModelMNDO.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Procedures for calculating an MNDO quantum chemical energy.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "ChargeConstraintContainer.h"
# include "DefineStatements.h"
# include "GaussianBasisCore.h"
# include "GaussianBasisIntegrals.h"
# include "Memory.h"
# include "MNDODefinitions.h"
# include "MNDOIntegrals.h"
# include "MNDOParameters.h"
# include "QCModelMNDO.h"
# include "Real1DArray.h"
# include "Real2DArray.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Defaults for some variables.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DefaultIsSpinRestricted          True
# define DefaultKeepOrbitalData           False
# define DefaultLinkAtomRatio             True
# define DefaultFermiBroadening           1.0e+3
# define DefaultNumberFractionalAlphaHOOs 1
# define DefaultNumberFractionalAlphaLUOs 0
# define DefaultNumberFractionalBetaHOOs  1
# define DefaultNumberFractionalBetaLUOs  0
# define DefaultNumberFractionalHOOs      1
# define DefaultNumberFractionalLUOs      0
# define DefaultQCChargeModel             QCChargeModel_Lowdin
# define DefaultOccupancyType             QCOnePDMOccupancyType_Cardinal

# define DOCIQCMM

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifdef MNDOCI
static Real CIStateCharacterDeterminant ( const Integer activeElectrons, const Boolean includeCoreOrbitals, const Integer numberCoreOrbitals,
                                                                  const Integer1DArray *iActiveIndices, const Integer1DArray *jActiveIndices,
                                                                const Real2DArray *orbitalTransformation, Real2DArray *work, Status *status ) ;
# endif
static void MakeQsPcTransformation      ( const QCModelMNDO *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3 ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCModelMNDO *QCModelMNDO_Allocate ( void )
{
    QCModelMNDO *self = NULL ;
    self = ( QCModelMNDO * ) Memory_Allocate ( sizeof ( QCModelMNDO ) ) ;
    if ( self != NULL )
    {
        self->isSpinRestricted          = DefaultIsSpinRestricted          ;
        self->keepOrbitalData           = DefaultKeepOrbitalData           ;
        self->linkAtomRatio             = DefaultLinkAtomRatio             ;
        self->fermiBroadening           = DefaultFermiBroadening           ;
        self->numberFractionalAlphaHOOs = DefaultNumberFractionalAlphaHOOs ;
        self->numberFractionalAlphaLUOs = DefaultNumberFractionalAlphaLUOs ;
        self->numberFractionalBetaHOOs  = DefaultNumberFractionalBetaHOOs  ;
        self->numberFractionalBetaLUOs  = DefaultNumberFractionalBetaLUOs  ;
        self->numberFractionalHOOs      = DefaultNumberFractionalHOOs      ;
        self->numberFractionalLUOs      = DefaultNumberFractionalLUOs      ;
# ifdef MNDOCI
        self->ciModel                   = NULL                             ;
# endif
        self->qcChargeModel             = DefaultQCChargeModel             ;
        self->occupancyType             = DefaultOccupancyType             ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fock contributions - charge constraints using Lowdin charges.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void QCModelMNDO_CCFockLowdinWork ( const QCAtomContainer           *qcAtoms           ,
                                           const ChargeConstraintContainer *chargeConstraints ,
                                           const Real1DArray               *lambdas           ,
                                           const Real                       spinFactor        ,
                                                 SymmetricMatrix           *fock              ,
                                                 Real1DArray               *fCC               )
{
    auto Integer a, b, bStart, bStop, c, t ;
    auto Real    lambda, w ;
    auto ChargeConstraint *constraint ;

    /* . Initialization. */
    Real1DArray_Set ( fCC, 0.0e+00 ) ;
/*
printf ( "\nFOCK ON ENTRY:\n" ) ;
SymmetricMatrix_Print ( fock ) ;
printf ( "\nLAMBDAS:\n" ) ;
Real1DArray_Print ( lambdas ) ;
*/
    /* . Loop over constraints. */
    for ( c = 0 ; c < chargeConstraints->numberOfConstraints ; c++ )
    {
        constraint = chargeConstraints->constraints[c] ;
        lambda     = Real1DArray_Item ( lambdas, c ) ;
        for ( t = 0 ; t < constraint->numberOfCharges ; t++ )
        {
            a =   Integer1DArray_Item ( constraint->chargeIndices, t ) ;
            w = - Integer1DArray_Item ( constraint->chargeWeights, t ) * lambda ;
            bStart = qcAtoms->data[a].ostart ;
            bStop  = bStart + qcAtoms->data[a].nobasis ;
            for ( b = bStart ; b < bStop ; b++ ) Real1DArray_Item ( fCC, b ) += w ;
        }
        for ( t = 0 ; t < constraint->numberOfSpins ; t++ )
        {
            a =   Integer1DArray_Item ( constraint->spinIndices, t ) ;
            w = - Integer1DArray_Item ( constraint->spinWeights, t ) * lambda * spinFactor ;
            bStart = qcAtoms->data[a].ostart ;
            bStop  = bStart + qcAtoms->data[a].nobasis ;
            for ( b = bStart ; b < bStop ; b++ ) Real1DArray_Item ( fCC, b ) += w ;
        }
    }
/*
printf ( "\nFCC:\n" ) ;
Real1DArray_Print ( fCC ) ;
*/
    /* . Finish up. */
    SymmetricMatrix_IncrementDiagonalFromArray ( fock, fCC, NULL ) ;
/*
printf ( "\nFOCK ON EXIT:\n" ) ;
SymmetricMatrix_Print ( fock ) ;
*/
}

void QCModelMNDO_CCFockLowdin ( const QCModelMNDO *self, QCModelMNDOState *qcState, const ChargeConstraintContainer *chargeConstraints, const Real1DArray *lambdas, Status *status )
{
    if ( ( self != NULL ) && ( qcState != NULL ) && ( chargeConstraints != NULL ) && ( lambdas != NULL ) )
    {
        auto Real1DArray *fCC ;
        fCC = Real1DArray_Allocate ( qcState->qcAtoms->nobasis, status ) ;
        if ( fCC != NULL )
        {
                                             QCModelMNDO_CCFockLowdinWork ( qcState->qcAtoms, chargeConstraints, lambdas, -1.0e+00, qcState->densityp->fock, fCC ) ;
            if ( qcState->densityq != NULL ) QCModelMNDO_CCFockLowdinWork ( qcState->qcAtoms, chargeConstraints, lambdas,  1.0e+00, qcState->densityq->fock, fCC ) ;
        }
        Real1DArray_Deallocate ( &fCC ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge constraint function, gradient and hessian.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real QCModelMNDO_CCFunction ( const QCModelMNDO               *self              ,
                                    QCModelMNDOState          *qcState           ,
                              const ChargeConstraintContainer *chargeConstraints ,
                              const Real1DArray               *lambdas           ,
                              const Boolean                    buildFock         ,
                                    Real1DArray               *gradients         , 
                                    SymmetricMatrix           *hessian           ,
                                    Status *status                               )
{
    Real f = 0.0e+00 ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( chargeConstraints != NULL ) && ( lambdas != NULL ) )
    {
        auto Integer      a, c ;
        auto Real1DArray *charges = NULL, *deviations = NULL, *spins = NULL ;
        c          = ChargeConstraintContainer_Size ( chargeConstraints ) ;
        a          = QCAtomContainer_Size           ( qcState->qcAtoms  ) ;
        charges    = Real1DArray_Allocate ( a, status ) ;
        deviations = Real1DArray_Allocate ( c, status ) ;
        spins      = Real1DArray_Allocate ( a, status ) ;
        if ( ( charges != NULL ) && ( deviations != NULL ) && ( spins != NULL ) )
        {
            if ( buildFock )
            {
                QCModelMNDOState_RestoreFock ( qcState ) ;
                QCModelMNDO_CCFockLowdin ( self, qcState, chargeConstraints, lambdas, status ) ;
            }
            QCModelMNDOState_MakeDensities ( qcState ) ;
            Real1DArray_Set ( charges, 0.0e+00 ) ;
            Real1DArray_Set ( spins  , 0.0e+00 ) ;
            QCModelMNDO_LowdinCharges ( self, qcState, False, charges ) ;
            QCModelMNDO_LowdinCharges ( self, qcState, True , spins   ) ;
            ChargeConstraintContainer_Deviations ( chargeConstraints, charges, spins, deviations, status ) ;
            QCModelMNDOState_RestoreFock ( qcState ) ;
            f = - ( QCModelMNDOState_EnergyFromFock ( qcState ) + Real1DArray_Dot ( deviations, lambdas, status ) ) ;
            if ( gradients != NULL )
            {
                Real1DArray_CopyTo ( deviations, gradients, status ) ;
                Real1DArray_Scale  ( gradients, -1.0e+00 ) ;
            }
            if ( hessian != NULL ) QCModelMNDO_CCHessian ( self, qcState, chargeConstraints, hessian, status ) ;
        }
        Real1DArray_Deallocate ( &charges    ) ;
        Real1DArray_Deallocate ( &deviations ) ;
        Real1DArray_Deallocate ( &spins      ) ;        
    }
    return f ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . CC hessian.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Tolerance 1.0e-10
static void QCModelMNDO_CCHessianWork ( const QCAtomContainer           *qcAtoms           ,
                                        const ChargeConstraintContainer *chargeConstraints ,
                                        const Real                       spinFactor        ,
                                        const Integer                    occupied          ,
                                        const Integer                    virtual           ,
                                        const Real1DArray               *energies          ,
                                        const Real2DArray               *orbitals          ,
                                              SymmetricMatrix           *hessian           ,
                                              Real2DArray               *eFactors          ,
                                              Real2DArray               *eW2               ,
                                              Real2DArray              **wMatrices         )
{
    auto Integer           a, b, bStart, bStop, c, c1, c2, o, v, t ;
    auto Real              eO, eV, f, w ;
    auto ChargeConstraint *constraint ;
    auto Real2DArray      *w1, *wMatrix ;

    /* . Create the constraint coupling matrices. */
    for ( c = 0 ; c < chargeConstraints->numberOfConstraints ; c++ )
    {
        constraint = chargeConstraints->constraints[c] ;
        wMatrix    = wMatrices[c] ;
        Real2DArray_Set ( wMatrix, 0.0e+00 ) ;
        for ( t = 0 ; t < constraint->numberOfCharges ; t++ )
        {
            a =   Integer1DArray_Item ( constraint->chargeIndices, t ) ;
            w = - Integer1DArray_Item ( constraint->chargeWeights, t ) ;
            bStart = qcAtoms->data[a].ostart ;
            bStop  = bStart + qcAtoms->data[a].nobasis ;
            for ( o = 0 ; o < occupied ; o++ )
            {
                for ( v = 0 ; v < virtual ; v++ )
                {
                    for ( b = bStart ; b < bStop ; b++ ) Real2DArray_Item ( wMatrix, o, v ) += w * Real2DArray_Item ( orbitals, b, o ) * Real2DArray_Item ( orbitals, b, v + occupied ) ;
                }
            }
        }
        for ( t = 0 ; t < constraint->numberOfSpins ; t++ )
        {
            a =   Integer1DArray_Item ( constraint->spinIndices, t ) ;
            w = - Integer1DArray_Item ( constraint->spinWeights, t ) * spinFactor ;
            bStart = qcAtoms->data[a].ostart ;
            bStop  = bStart + qcAtoms->data[a].nobasis ;
            for ( o = 0 ; o < occupied ; o++ )
            {
                for ( v = 0 ; v < virtual ; v++ )
                {
                    for ( b = bStart ; b < bStop ; b++ ) Real2DArray_Item ( wMatrix, o, v ) += w * Real2DArray_Item ( orbitals, b, o ) * Real2DArray_Item ( orbitals, b, v + occupied ) ;
                }
            }
        }
    }

    /* . Create the matrix of eigenvalue factors. */
    Real2DArray_Set ( eFactors, 0.0e+00 ) ;
    for ( o = 0 ; o < occupied ; o++ )
    {
        eO = Real1DArray_Item ( energies, o ) ;
        for ( v = 0 ; v < virtual ; v++ )
        {
            eV = Real1DArray_Item ( energies, v + occupied ) ;
            f  = ( eO - eV ) ;
            if ( fabs ( f ) > _Tolerance ) Real2DArray_Item ( eFactors, o, v ) = 1.0e+00 / f ;
        }
    }

    /* . Form the elements of the matrix. */
    for ( c1 = 0 ; c1 < chargeConstraints->numberOfConstraints ; c1++ )
    {
        w1 = wMatrices[c1] ;
        for ( c2 = 0 ; c2 <= c1 ; c2++ )
        {
            Real2DArray_CopyTo   ( wMatrices[c2], eW2, NULL ) ;
            Real2DArray_Multiply ( eW2, eFactors, NULL ) ;
            SymmetricMatrix_Item ( hessian, c1, c2 ) += Real2DArray_Trace ( w1, eW2, NULL ) ;
        }
    }
}
# undef _Tolerance

void QCModelMNDO_CCHessian ( const QCModelMNDO *self, QCModelMNDOState *qcState, const ChargeConstraintContainer *chargeConstraints, SymmetricMatrix *hessian, Status *status )
{
    SymmetricMatrix_Set ( hessian, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( chargeConstraints != NULL ) && ( hessian != NULL ) )
    {
        auto Integer occupied, virtual ;
        occupied = qcState->densityp->numberOccupied ;
        virtual  = qcState->densityp->numberOrbitals - occupied ;
        if ( qcState->densityq != NULL )
        {
            occupied = Maximum ( occupied, qcState->densityq->numberOccupied ) ;
            virtual  = Maximum ( virtual , qcState->densityq->numberOrbitals - qcState->densityq->numberOccupied ) ;
        }
        if ( ( occupied > 0 ) && ( virtual > 0 ) )
        {
            auto Boolean isOK = False ;
            auto Integer c     ;
            auto Real    scale = -4.0e+00 ;
            auto Real2DArray *eFactors = NULL, *eW2 = NULL, **wMatrices = NULL ;
            eFactors = Real2DArray_Allocate ( occupied, virtual, status ) ;
            eW2      = Real2DArray_Allocate ( occupied, virtual, status ) ;
            MEMORY_ALLOCATEARRAY ( wMatrices, chargeConstraints->numberOfConstraints, Real2DArray * ) ;
            if ( wMatrices != NULL )
            {
                isOK = True ;
                for ( c = 0 ; c < chargeConstraints->numberOfConstraints ; c++ )
                {
                    wMatrices[c] = Real2DArray_Allocate ( occupied, virtual, status ) ;
                    isOK = isOK && ( wMatrices[c] != NULL ) ;
                }
            }
            if ( ( eFactors != NULL ) && ( eW2 != NULL ) && isOK )
            {
                occupied = qcState->densityp->numberOccupied ;
                virtual  = qcState->densityp->numberOrbitals - occupied ;
                QCModelMNDO_CCHessianWork (  qcState->qcAtoms            ,
                                             chargeConstraints           ,
                                            -1.0e+00                     ,
                                             occupied                    ,
                                             virtual                     ,
                                             qcState->densityp->energies ,
                                             qcState->densityp->orbitals ,
                                             hessian                     ,
                                             eFactors                    ,
                                             eW2                         ,
                                             wMatrices                   ) ;
                if ( qcState->densityq != NULL )
                {
                    occupied = qcState->densityq->numberOccupied ;
                    virtual  = qcState->densityq->numberOrbitals - occupied ;
                    QCModelMNDO_CCHessianWork ( qcState->qcAtoms            ,
                                                chargeConstraints           ,
                                                1.0e+00                     ,
                                                occupied                    ,
                                                virtual                     ,
                                                qcState->densityq->energies ,
                                                qcState->densityq->orbitals ,
                                                hessian                     ,
                                                eFactors                    ,
                                                eW2                         ,
                                                wMatrices                   ) ;
                    scale = -2.0e+00 ;
                }
                SymmetricMatrix_Scale ( hessian, scale ) ;
            }
            Real2DArray_Deallocate ( &eFactors ) ;
            Real2DArray_Deallocate ( &eW2      ) ;
            for ( c = 0 ; c < chargeConstraints->numberOfConstraints ; c++ ) Real2DArray_Deallocate ( &wMatrices[c] ) ;
            MEMORY_DEALLOCATE ( wMatrices ) ;
        }
    }
}

# ifdef MNDOCI
/*----------------------------------------------------------------------------------------------------------------------------------
! . CI energy.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_CIEnergy ( const QCModelMNDO *self, QCModelMNDOState *state )
{
    if ( ( self != NULL ) && ( state != NULL ) && ( self->ciModel != NULL ) && ( state->cistate != NULL ) )
    {
        Boolean      doQCMM  ;
        MNDOCIState *cistate ;
        Real         eelectronic, enuclear ;

        /* . Initialization. */
        cistate     = state->cistate ;
        doQCMM      = ( ( state->qcmmstate != NULL ) && ( state->qcmmstate->qcmmPotentials != NULL ) ) ;
        eelectronic = state->eelectronic ;
        enuclear    = state->enuclear    ;

        /* . Set up some data structures from the SCF calculation. */
        cistate = state->cistate ;
        cistate->oneelectronmatrix    = state->oneelectronmatrix     ;
        cistate->twoelectronintegrals = state->twoelectronintegrals  ;

# ifdef DOCIQCMM
/* . QC/MM. */
/* . Add the potentials into the one-electronmatrix. */
if ( doQCMM )
{
    auto Real       e, p ;
    auto Integer          iatom, ibf, istart, istop ;
    auto Real1DArray *work ;

    /* . QC/MM only. */
    work = state->qcmmstate->qcmmPotentials ;

    /* . Add in the contributions to the one-electron matrix. */
    for ( iatom = 0 ; iatom < state->qcAtoms->natoms ; iatom++ )
    {
        p = Real1DArray_Item ( work, iatom ) ;
        istart = state->qcAtoms->data[iatom].ostart ;
        istop  = istart + state->qcAtoms->data[iatom].nobasis ;
        for ( ibf = istart ; ibf < istop ; ibf++ ) SymmetricMatrix_IncrementComponent ( cistate->oneelectronmatrix, ibf, ibf, - p ) ;
    }

    /* . Modify the energies. */
    QCAtomContainer_GetNuclearCharges ( state->qcAtoms, state->qcParameters, state->qcmmstate->qcCharges ) ;
    e = Real1DArray_Dot ( state->qcmmstate->qcCharges, work, NULL ) ;
    eelectronic += ( state->eqcmm - e ) ;
    enuclear    += e ;
}
# endif

        /* . Localize the orbitals if requested. */
        if ( self->ciModel->localizeOrbitals ) QCModelMNDO_LocalizeOrbitals ( self, state, False, cistate->localizeStart, cistate->localizeStop, NULL ) ;

        /* . Calculate the energy. */
        MNDOCIModel_Energy ( self->ciModel, cistate, eelectronic, enuclear, state->densityp, state->densityq, NULL ) ;

# ifdef DOCIQCMM
/* . QC/MM. */
/* . Recalculate the charges with the z-matrix for the gradient calculation. */
/* . This means that these charges cannot be relied upon as charges! */
if ( ( doQCMM ) && ( cistate->doGradients ) )
{
    auto SymmetricMatrix *dtotal ;
    dtotal = state->densityp->density ;
    SymmetricMatrix_AddScaledMatrix ( dtotal,  1.0e+00, state->cistate->zMatrix ) ;
    QCModelMNDO_LowdinCharges ( self, state, False, state->qcmmstate->qcCharges ) ;
    SymmetricMatrix_AddScaledMatrix ( dtotal, -1.0e+00, state->cistate->zMatrix ) ;
}
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . CI state characters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . The configurations are stored as ( core with alpha and beta alternating ) ( alpha active ) ( beta active ).
! . Phases are unnecessary when cores are not included because alphas and betas are already ordered.
! . Likewise when cores are included the phase to get all core betas to the active betas is determined by whether
! . Nc * Na + ( Nc * ( Nc - 1 ) ) / 2 is odd or even. Here Nc is the number of core orbitals and Na is the number
! . of active alpha. However, these phases are always the same for non-zero interactions between states and so
! . can be ignored as they multiply to 1.
*/
Real1DArray *QCModelMNDO_CIStateCharacters ( const QCModelMNDO *self, QCModelMNDOState *state, const Matrix33 *rotation, const Integer1DArray *mapping,
                                                                               const Selection *stateIndices, const Boolean includeCoreOrbitals, Status *status )
{
    Real1DArray *characters = NULL ;
    if ( ( self != NULL ) && ( state != NULL ) && ( self->ciModel != NULL ) && ( state->cistate != NULL ) && ( rotation != NULL ) && ( mapping != NULL ) && ( stateIndices != NULL ) )
    {
        if ( ( state->densityp != NULL ) && ( state->densityp->orbitals != NULL ) )
        {
            auto Integer         i, j, nai, naj, nbi, s, start, stop ;
            auto Integer1DArray *ialphas, *ibetas, *jalphas, *jbetas ;
            auto Real1DArray     inOrbital, inState, outOrbital, *outState ;
            auto Real2DArray     activeOrbitals, *orbitals, *orbitalTransformation, *rotatedOrbitals, *stateTransformation, *work ;

            /* . Get the orbital transformation. */
            orbitals = state->densityp->orbitals ;

            /* . Find the orbital range. */
            if ( includeCoreOrbitals ) start = 0 ;
            else                       start = state->cistate->ncore ;
            stop = ( state->cistate->ncore + state->cistate->nactive ) ;

            /* . Allocate all space - including for operations below. */
            characters            = Real1DArray_Allocate ( stateIndices->nindices         , status ) ;
            outState              = Real1DArray_Allocate ( state->cistate->nconfigurations, status ) ;
            orbitalTransformation = Real2DArray_Allocate ( stop - start     , stop - start, status ) ;
            rotatedOrbitals       = Real2DArray_Allocate ( orbitals->length0, stop - start, status ) ;
            stateTransformation   = Real2DArray_Allocate ( state->cistate->nconfigurations, state->cistate->nconfigurations, status ) ;
            work                  = Real2DArray_Allocate ( stop - start     , stop - start, status ) ;
            if ( ( characters != NULL ) && ( outState != NULL ) && ( orbitalTransformation != NULL ) && ( rotatedOrbitals != NULL ) && ( stateTransformation != NULL ) )
            {
                /* . Rotate orbitals. */
                for ( i = start ; i < stop ; i++ )
                {
                    Real2DArray_ColumnSlice ( orbitals       , i        ,  &inOrbital, status ) ;
                    Real2DArray_ColumnSlice ( rotatedOrbitals, i - start, &outOrbital, status ) ;
                    QCModelMNDO_RotateOrbital ( self, state, rotation, mapping, &inOrbital, &outOrbital ) ;
                }

                /* . Determine the transformation. */
                Real2DArray_Slice ( orbitals, 0, orbitals->length0, 1, start, stop, 1, &activeOrbitals, status ) ;
                Real2DArray_MatrixMultiply ( True, False, 1.0e+00, &activeOrbitals, rotatedOrbitals, 0.0e+00, orbitalTransformation, status ) ;

                /* . Get the state transformation. */
                Real2DArray_Set ( stateTransformation, 0.0e+00 ) ;

                /* . Double loop over configurations. */
                for ( i = 0 ; i < state->cistate->nconfigurations ; i++ )
                {
                    nai     = state->cistate->configurations[i].nalphas ;
                    nbi     = state->cistate->nelectrons - nai ;
                    ialphas = state->cistate->configurations[i].alphas  ;
                    ibetas  = state->cistate->configurations[i].betas   ;
                    for ( j = 0 ; j < state->cistate->nconfigurations ; j++ )
                    {
                        naj     = state->cistate->configurations[j].nalphas ;
                        jalphas = state->cistate->configurations[j].alphas  ;
                        jbetas  = state->cistate->configurations[j].betas   ;

                        /* . Skip if there are different numbers of alpha orbitals in the two configurations. */
                        if ( nai != naj ) continue ;
# ifdef DEBUGMNDOCISYMMETRY
if ( ( Integer1DArray_Sum ( ialphas ) != Integer1DArray_Sum ( jalphas ) ) ||
     ( Integer1DArray_Sum ( ibetas  ) != Integer1DArray_Sum ( jbetas  ) ) ||
     ( Integer1DArray_Sum ( ibetas  ) != nbi ) ) printf ( "\nInvalid number of electrons in CI state character determination.\n" ) ;
# endif
                        /* . Set the value. */
                        Real2DArray_Item ( stateTransformation, i, j ) = CIStateCharacterDeterminant ( nai, includeCoreOrbitals, state->cistate->ncore, ialphas, jalphas, orbitalTransformation, work, status ) *
                                                                         CIStateCharacterDeterminant ( nbi, includeCoreOrbitals, state->cistate->ncore, ibetas , jbetas , orbitalTransformation, work, status ) ;
                    }
                }

# ifdef DEBUGMNDOCISYMMETRY
printf ( "\nCI State Character Determination: orbital range = %d %d\n", start, stop ) ;
printf ( "Input Orbitals\n" ) ;
Real2DArray_Print ( orbitals ) ;
printf ( "Rotated Orbitals\n" ) ;
Real2DArray_Print ( rotatedOrbitals ) ;
printf ( "Orbital Transformation\n" ) ;
Real2DArray_Print ( orbitalTransformation ) ;
printf ( "State Transformation\n" ) ;
Real2DArray_Print ( stateTransformation ) ;
# endif

                /* . Get the characters. */
                for ( s = 0 ; s < stateIndices->nindices ; s++ )
                {
                    i = stateIndices->indices[s] ;
                    Real2DArray_ColumnSlice    ( state->cistate->ciVectors, i, &inState, status ) ;
                    Real2DArray_VectorMultiply ( False, 1.0e+00, stateTransformation, &inState, 0.0e+00, outState, status ) ;
                    Real1DArray_Item ( characters, s ) = Real1DArray_Dot ( &inState, outState, status ) ;
                }
            }
            /* . Out of space. */
            else Status_Set ( status, Status_MemoryAllocationFailure ) ;

            /* . Finish up. */
            Real1DArray_Deallocate ( &outState              ) ;
            Real2DArray_Deallocate ( &orbitalTransformation ) ;
            Real2DArray_Deallocate ( &rotatedOrbitals       ) ;
            Real2DArray_Deallocate ( &stateTransformation   ) ;
            Real2DArray_Deallocate ( &work                  ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return characters ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCModelMNDO *QCModelMNDO_Clone ( const QCModelMNDO *self )
{
    QCModelMNDO *new = NULL ;
    if ( self != NULL )
    {
        new = QCModelMNDO_Allocate ( ) ;
        new->isSpinRestricted          = self->isSpinRestricted          ;
        new->keepOrbitalData           = self->keepOrbitalData           ;
        new->linkAtomRatio             = self->linkAtomRatio             ;
        new->fermiBroadening           = self->fermiBroadening           ;
        new->numberFractionalAlphaHOOs = self->numberFractionalAlphaHOOs ;
        new->numberFractionalAlphaLUOs = self->numberFractionalAlphaLUOs ;
        new->numberFractionalBetaHOOs  = self->numberFractionalBetaHOOs  ;
        new->numberFractionalBetaLUOs  = self->numberFractionalBetaLUOs  ;
        new->numberFractionalHOOs      = self->numberFractionalHOOs      ;
        new->numberFractionalLUOs      = self->numberFractionalLUOs      ;
        new->qcChargeModel             = self->qcChargeModel             ;
        new->occupancyType             = self->occupancyType             ;
# ifdef MNDOCI
        new->ciModel = MNDOCIModel_Clone ( self->ciModel, NULL ) ;
# endif
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_Deallocate ( QCModelMNDO **self )
{
    if ( (*self) != NULL )
    {
# ifdef MNDOCI
        MNDOCIModel_Deallocate ( &((*self)->ciModel) ) ;
# endif
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole moment.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Notes: DIPFC1 = AU_to_Debyes * ANGSTROMS_TO_BOHRS and DIPFC2 = 2 * AU_to_Debyes. */
# define DIPFC1 4.8030e+00
# define DIPFC2 5.0832e+00

Vector3 *QCModelMNDO_DipoleMoment ( const QCModelMNDO *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3, const Vector3 *center )
{
    Vector3 *dipole = NULL;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( coordinates3 != NULL ) )
    {
# ifdef MNDODORBITALS
        auto Real dx, dy, dz, factor, hyfpd, hyfsp, q, xc = 0.0e+00, xq, yc = 0.0e+00, yq, zc = 0.0e+00, zq ;
# else
        auto Real dx, dy, dz, factor, hyfsp, q, xc = 0.0e+00, xq, yc = 0.0e+00, yq, zc = 0.0e+00, zq ;
# endif
        auto Integer    ic, iqm, istart ;
        auto MNDOParameters *idata     ;
        auto Real1DArray    *qcCharges ;

        /* . Initialization. */
        dipole = Vector3_Allocate ( )   ;
        Vector3_Set ( dipole, 0.0e+00 ) ;

        /* . Get the center. */
        if ( center != NULL )
        {
            xc = Vector3_Item ( center, 0 ) ;
            yc = Vector3_Item ( center, 1 ) ;
            zc = Vector3_Item ( center, 2 ) ;
        }

        /* . Get the charges. */
        qcCharges = Real1DArray_Allocate ( qcState->qcAtoms->natoms, NULL ) ;
        Real1DArray_Set ( qcCharges, 0.0e+00 ) ;
        QCModelMNDO_LowdinCharges ( self, qcState, False, qcCharges ) ;

        /* . Get the total density. */
        QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Initialize some constants. */
        factor = 1.0e+00 / sqrt ( 3.0e+00 ) ;

        /* . Loop over the atoms. */
        for ( iqm = 0, dx = dy = dz = 0.0e+00 ; iqm < qcState->qcAtoms->natoms ; iqm++ )
        {
            /* . Get some information about the atom. */
            ic     = qcState->qcAtoms->data[iqm].center ;
            idata  = qcState->qcParameters->centers[ic].mndoparameters ;
            istart = qcState->qcAtoms->data[iqm].ostart ;
            /* . Charge part. */
            q   = DIPFC1 * Real1DArray_Item ( qcCharges, iqm ) ;
            QCAtomContainer_GetAtomCoordinates3 ( qcState->qcAtoms, iqm, coordinates3, &xq, &yq, &zq ) ;
            dx += q * ( xq - xc ) ;
            dy += q * ( yq - yc ) ;
            dz += q * ( zq - zc ) ;
            /* . Electronic part. */
            /* . sp. */
            if ( idata->norbitals >= 4 )
            {
                hyfsp  = DIPFC2 * idata->ddp[1] ; /* . originally ddp[2] or dd; */
                dx -= hyfsp * SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + PX, istart ) ;
                dy -= hyfsp * SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + PY, istart ) ;
                dz -= hyfsp * SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + PZ, istart ) ;
            }
# ifdef MNDODORBITALS
            /*. pd. */
            if ( idata->norbitals >= 9 )
            {
                hyfpd = DIPFC2 * idata->ddp[4] ; /* . originally ddp[5] ; */
                dx -= hyfpd * (          SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DXZ,   istart + PZ ) +
                                         SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DX2Y2, istart + PX ) +
                                         SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DXY,   istart + PY ) -
                                factor * SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DZ2,   istart + PX ) ) ;
                dy -= hyfpd * (          SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DYZ,   istart + PZ ) -
                                         SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DX2Y2, istart + PY ) +
                                         SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DXY,   istart + PX ) -
                                factor * SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DZ2,   istart + PY ) ) ;
                dz -= hyfpd * (          SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DXZ,   istart + PX ) +
                                         SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DYZ,   istart + PY ) +
                      2.0e+00 * factor * SymmetricMatrix_Get_Component ( qcState->densityp->density, istart + DZ2,   istart + PZ ) ) ;
            }
# endif
        }

        /* . Fill the dipole moment vector. */
        Vector3_Item ( dipole, 0 ) = dx ;
        Vector3_Item ( dipole, 1 ) = dy ;
        Vector3_Item ( dipole, 2 ) = dz ;

        /* . Finish up. */
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;
        Real1DArray_Deallocate ( &qcCharges ) ;
    }
    return dipole ;
}
# undef DIPFC1
# undef DIPFC2

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the Fock matrices and the electronic energy.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_Fock ( const QCModelMNDO *self, QCModelMNDOState *qcState, Real *eelectronic )
{
    /* . Initialization.*/
    if ( eelectronic != NULL ) (*eelectronic) = 0.0e+00 ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( qcState->densityp != NULL ) && ( qcState->oneelectronmatrix != NULL ) && ( qcState->twoelectronintegrals != NULL ) )
    {
        /* . Build the two-electron part of the Fock matrices. */
        if ( qcState->densityq == NULL ) MNDOIntegrals_MakeFockG ( qcState->twoelectronintegrals, qcState->densityp->density, qcState->densityp->fock, NULL, NULL ) ;
        else                             MNDOIntegrals_MakeFockG ( qcState->twoelectronintegrals, qcState->densityp->density, qcState->densityp->fock, qcState->densityq->density, qcState->densityq->fock ) ;

        /* . Determine the energies and finish construction of the matrices. */
        qcState->eocc = qcState->densityp->occupancyEnergy ;
        qcState->eoei =           SymmetricMatrix_Multiply2_Trace ( qcState->densityp->density, qcState->oneelectronmatrix ) ;
        qcState->etei = 0.5e+00 * SymmetricMatrix_Multiply2_Trace ( qcState->densityp->density, qcState->densityp->fock    ) ;
        SymmetricMatrix_Increment ( qcState->densityp->fock, qcState->oneelectronmatrix ) ;
        if ( qcState->densityq != NULL )
        {
            qcState->eocc += qcState->densityq->occupancyEnergy ;
            qcState->eoei +=           SymmetricMatrix_Multiply2_Trace ( qcState->densityq->density, qcState->oneelectronmatrix ) ;
            qcState->etei += 0.5e+00 * SymmetricMatrix_Multiply2_Trace ( qcState->densityq->density, qcState->densityq->fock    ) ;
            SymmetricMatrix_Increment ( qcState->densityq->fock, qcState->oneelectronmatrix ) ;
        }

        /* . Sum the total electronic energy. */
        qcState->eelectronic = qcState->eocc + qcState->eoei + qcState->etei ;
        if ( eelectronic != NULL ) (*eelectronic) = qcState->eelectronic ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the gradients for a MNDO method.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_Gradients ( const QCModelMNDO *self, QCModelMNDOState *qcState )
{
    if ( ( self != NULL ) && ( qcState != NULL ) && ( qcState->gradients3 != NULL ) )
    {
        auto Coordinates3      *qcgradients3 ;
        auto Real                bpij, eng[3], *xi, *xj ;
        auto Integer                   i, ic, ifirst, ij, iorbs, iqm, j, jc, jfirst, jorbs, jqm ;
        auto Real2DArray          *dTwoIJ = NULL , *sijx, *sijy, *sijz ;
        auto MNDOParameters  *idata, *jdata ;
        auto SymmetricMatrix *dspin = NULL, *dtotal = NULL ;
        auto Real1DArray          *dOneI = NULL, *dOneJ = NULL ;

        /* . Get a gradient array - automatically initialized. */
        qcgradients3 = Coordinates3_Allocate ( qcState->qcAtoms->natoms ) ;
        Coordinates3_Set ( qcgradients3, 0.0e+00 ) ;

# ifdef MNDOCI
        /* . Set up densities (non-CI and CI cases). */
        if ( qcState->cistate == NULL )
        {
# endif
            QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;
            if ( qcState->densityp != NULL ) dtotal = qcState->densityp->density ;
            if ( qcState->densityq != NULL ) dspin  = qcState->densityq->density ;
# ifdef MNDOCI
        }
        else
        {
            dtotal = qcState->cistate->work1 ;
            dspin  = NULL ;
            SymmetricMatrix_CopyTo          ( qcState->densityp->density, dtotal ) ;
            SymmetricMatrix_AddScaledMatrix ( dtotal, 1.0e+00, qcState->cistate->zMatrix ) ;
# ifdef DEBUGMNDOCIGRADIENTS
printf ( "\nD Total:\n" ) ;
SymmetricMatrix_Print ( dtotal ) ;
printf ( "\npCore:\n" ) ;
SymmetricMatrix_Print ( qcState->cistate->pcore ) ;
printf ( "\nOnePDM:\n" ) ;
SymmetricMatrix_Print ( qcState->cistate->onepdm ) ;
printf ( "\npHF:\n" ) ;
SymmetricMatrix_Print ( qcState->cistate->onepdmHF ) ;
printf ( "\nz-Matrix:\n" ) ;
SymmetricMatrix_Print ( qcState->cistate->zMatrix ) ;
# endif
        }
# endif

        /* . Outer loop over atoms. */
        for ( iqm = 0 ; iqm < qcState->qcAtoms->natoms ; iqm++ )
        {
           /* . Get data for the atom. */
           ic     = qcState->qcAtoms->data[iqm].center ;
           idata  = qcState->qcParameters->centers[ic].mndoparameters ;
           ifirst = qcState->qcAtoms->data[iqm].ostart ;
           iorbs  = qcState->qcParameters->centers[ic].mndoparameters->norbitals ;
           xi     = Coordinates3_RowPointer ( qcState->qccoordinates3, iqm ) ;

           /* . Inner loop over atoms. */
           for ( jqm = 0 ; jqm < iqm ; jqm++ )
           {
              /* . Get data for the atom. */
              jc     = qcState->qcAtoms->data[jqm].center ;
              jdata  = qcState->qcParameters->centers[jc].mndoparameters ;
              jfirst = qcState->qcAtoms->data[jqm].ostart ;
              jorbs  = qcState->qcParameters->centers[jc].mndoparameters->norbitals ;
              xj     = Coordinates3_RowPointer ( qcState->qccoordinates3, jqm ) ;

              /* . The two-center core, oei and tei integral derivatives. */
# ifdef MNDOCI
# ifdef DEBUGMNDOCIGRADIENTS
printf ( "\nAtoms: %d, %d\n", iqm, jqm ) ;
# endif
# endif
              QCModelMNDOState_GetGradientDensityTerms ( qcState, idata, ifirst, jdata, jfirst, dtotal, dspin, &dOneI, &dOneJ, &dTwoIJ ) ;
              MNDOIntegrals_MolecularFrame2CIntegralsD ( idata, ifirst, xi, jdata, jfirst, xj, dOneI, dOneJ, dTwoIJ, eng ) ;
              Real1DArray_Deallocate ( &dOneI  ) ;
              Real1DArray_Deallocate ( &dOneJ  ) ;
              Real2DArray_Deallocate ( &dTwoIJ ) ;

              /* . The resonance integral derivatives. */
              GaussianBasis_2OverlapD ( qcState->qcParameters->centers[ic].orbitalbasis, xi, qcState->qcParameters->centers[jc].orbitalbasis, xj, &sijx, &sijy, &sijz ) ;

              /* . Transform the integrals. */
              QCModelMNDO_TransformIntegrals ( &sijx, qcState->qcParameters->centers[ic].orbitalbasis->c2o, qcState->qcParameters->centers[jc].orbitalbasis->c2o ) ;
              QCModelMNDO_TransformIntegrals ( &sijy, qcState->qcParameters->centers[ic].orbitalbasis->c2o, qcState->qcParameters->centers[jc].orbitalbasis->c2o ) ;
              QCModelMNDO_TransformIntegrals ( &sijz, qcState->qcParameters->centers[ic].orbitalbasis->c2o, qcState->qcParameters->centers[jc].orbitalbasis->c2o ) ;
              if ( ( sijx != NULL ) && ( sijy != NULL ) && ( sijz != NULL ) )
              {
                 for ( i = 0, ij = 0 ; i < iorbs ; i++ )
                 {
                    ij = ( ( i + ifirst ) * ( ( i + ifirst ) + 1 ) ) / 2 + jfirst ;
                    for ( j = 0 ; j < jorbs ; ij++, j++ )
                    {
                        bpij    = ( idata->beta[i] + jdata->beta[j] ) * dtotal->data[ij] ; /* . Note the implicit factor of 2 here. */
                        eng[0] += bpij * Real2DArray_Item ( sijx, i, j ) ;
                        eng[1] += bpij * Real2DArray_Item ( sijy, i, j ) ;
                        eng[2] += bpij * Real2DArray_Item ( sijz, i, j ) ;
                    }
                 }
                 Real2DArray_Deallocate ( &sijx ) ;
                 Real2DArray_Deallocate ( &sijy ) ;
                 Real2DArray_Deallocate ( &sijz ) ;
              }

              /* . Add in the contributions to the gradients. */
              Coordinates3_IncrementRow ( qcgradients3, iqm, eng[0], eng[1], eng[2] ) ;
              Coordinates3_DecrementRow ( qcgradients3, jqm, eng[0], eng[1], eng[2] ) ;
           }
        }

        /* . Finish up. */
# ifdef MNDOCI
        if ( qcState->cistate == NULL ) QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;
# else
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;
# endif

        /* . Copy the gradients into their appropriate places. */
        QCAtomContainer_SetGradients3 ( qcState->qcAtoms, qcState->coordinates3, qcgradients3, True, &(qcState->gradients3) ) ;
        Coordinates3_Deallocate ( &qcgradients3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate grid point densities in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define BASISFUNCTIONTOLERANCE 1.0e-10
void QCModelMNDO_GridPointDensities ( const QCModelMNDO *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data, const Boolean QSPIN )
{
    Real1DArray_Set ( data, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( coordinates3 != NULL ) && ( points != NULL ) && ( data != NULL ) && ( points->length0 <= data->length ) )
    {
        auto Real                   tolerance      = BASISFUNCTIONTOLERANCE ;
        auto Coordinates3          *qcCoordinates3 = NULL ;
        auto GridFunctionDataBlock *basisData      = NULL ;
        auto QCOnePDM              *density        = NULL ;

        /* . Get the Q->P transformation. */
        MakeQsPcTransformation ( self, qcState, coordinates3 ) ;

        /* . Get the appropriately transformed density. */
        QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;
        density = QCOnePDM_Allocate ( qcState->qptransformation->length0, NULL ) ;
        if ( QSPIN ) SymmetricMatrix_Transform ( qcState->densityq->density, qcState->qptransformation, True, density->density ) ;
        else         SymmetricMatrix_Transform ( qcState->densityp->density, qcState->qptransformation, True, density->density ) ;
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Get the QC coordinates. */
        QCAtomContainer_GetCoordinates3 ( qcState->qcAtoms, coordinates3, True, &qcCoordinates3 ) ;

        /* . Generate data. */
        basisData = GridFunctionDataBlock_Allocate ( qcState->qcAtoms->nobasisw, points->length0, 0, NULL ) ;
        if ( basisData != NULL )
        {
            QCAtomContainer_OrbitalBasisGridPointValues ( qcState->qcAtoms, qcState->qcParameters, qcCoordinates3, points, True, &tolerance, basisData, NULL ) ;
            QCOnePDM_DensityGridValues ( density, basisData, data, NULL ) ;
        }

        /* . Deallocate space. */
        Coordinates3_Deallocate          ( &qcCoordinates3 ) ;
        GridFunctionDataBlock_Deallocate ( &basisData      ) ;
        QCOnePDM_Deallocate              ( &density        ) ;
    }
}
# undef BASISFUNCTIONTOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate grid point orbitals in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define BASISFUNCTIONTOLERANCE 1.0e-10
void QCModelMNDO_GridPointOrbitals ( const QCModelMNDO *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data, const Integer norbitals, Integer *orbitalIndices, const Boolean QALPHA )
{
    Real1DArray_Set ( data, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( coordinates3 != NULL ) && ( points != NULL ) && ( data != NULL ) && ( points->length0 <= data->length ) )
    {
        auto Boolean      QOK = True ;
        auto Integer      i, *indices, n ;
        auto QCOnePDM    *density   = NULL ;
        auto Real2DArray *cOrbitals = NULL, *sOrbitals = NULL ;

        /* . Get the Q->P transformation. */
        MakeQsPcTransformation ( self, qcState, coordinates3 ) ;

        /* . Get the correct orbital set. */
        if ( QALPHA ) density = qcState->densityp ;
        else          density = qcState->densityq ;
        if ( density   != NULL ) sOrbitals = density->orbitals ;
        if ( sOrbitals != NULL )
        {
            /* . Transform the orbitals. */
            cOrbitals = Real2DArray_Allocate ( qcState->qptransformation->length0, sOrbitals->length1, NULL ) ;
            Real2DArray_MatrixMultiply ( False, False, 1.0e+00, qcState->qptransformation, sOrbitals, 0.0e+00, cOrbitals, NULL ) ;

            /* . Check the orbital indices. */
            if ( ( norbitals > 0 ) && ( orbitalIndices != NULL ) )
            {
                n       = norbitals ;
                indices = orbitalIndices ;
                for ( i = 0 ; i < n ; i++ )
                {
                    if ( ( indices[i] < 0 ) || ( indices[i] > ( cOrbitals->length1 - 1 ) ) ) { QOK = False ; break ; }
                }
            }
            else
            {
                n          = 1 ;
                indices    = Memory_Allocate_Array_Integer ( n ) ;
                indices[0] = QCOnePDM_HOMOIndex ( density ) ;
                QOK        = ( ( indices[0] >= 0 ) && ( indices[0] < cOrbitals->length1 ) ) ;
            }
            if ( QOK )
            {
                auto Real                   tolerance      = BASISFUNCTIONTOLERANCE ;
                auto Coordinates3          *qcCoordinates3 = NULL ;
                auto GridFunctionDataBlock *basisData      = NULL ;
                auto Real2DArray           *orbitals       = NULL, view ;

                /* . Get the QC coordinates. */
                QCAtomContainer_GetCoordinates3 ( qcState->qcAtoms, coordinates3, True, &qcCoordinates3 ) ;

                /* . Generate basis data. */
                basisData = GridFunctionDataBlock_Allocate ( qcState->qcAtoms->nobasisw, points->length0, 0, NULL ) ;
                if ( basisData != NULL )
                {
                    QCAtomContainer_OrbitalBasisGridPointValues ( qcState->qcAtoms, qcState->qcParameters, qcCoordinates3, points, True, &tolerance, basisData, NULL ) ;

                    /* . Extract orbital data. */
                    orbitals = Real2DArray_Allocate ( basisData->indices->length, n, NULL ) ;
                    if ( orbitals != NULL )
                    {
                        auto Integer i, m, o ;
                        for ( i = 0 ; i < basisData->indices->length ; i++ )
                        {
                            m = Integer1DArray_Item ( basisData->indices, i ) ;
                            for ( o = 0 ; o < n ; o++ ) Real2DArray_Item ( orbitals, i, o ) = Real2DArray_Item ( cOrbitals, m, indices[o] ) ;
                        }
                        /* . Generate data - this needs to be better. */
                        Real2DArray_ViewOfRaw      ( &view, 0, points->length0, n, n, 1, Real1DArray_Data ( data ), data->size, NULL ) ;
                        Real2DArray_MatrixMultiply ( True, False, 1.0e+00, basisData->f, orbitals, 0.0e+00, &view, NULL ) ;
                    }
                }

                /* . Deallocate space. */
                Coordinates3_Deallocate          ( &qcCoordinates3 ) ;
                GridFunctionDataBlock_Deallocate ( &basisData      ) ;
                Real2DArray_Deallocate           ( &orbitals       ) ;
            }
            if ( orbitalIndices == NULL ) Memory_Deallocate_Integer ( &indices ) ;
        }
        Real2DArray_Deallocate ( &cOrbitals ) ;
    }
}
# undef BASISFUNCTIONTOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate grid point potentials in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_GridPointPotentials ( const QCModelMNDO *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *potentials )
{
    Real1DArray_Set ( potentials, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( coordinates3 != NULL ) && ( points != NULL ) && ( potentials != NULL ) && ( points->length0 <= potentials->length ) )
    {
        auto Coordinates3    *qcCoordinates3 = NULL ;
        auto SymmetricMatrix *density  = NULL ;
        auto Real1DArray     *znuclear = NULL ;

        /* . Get the Q->P transformation. */
        MakeQsPcTransformation ( self, qcState, coordinates3 ) ;

        /* . Get the appropriately transformed density. */
        QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;
        density = SymmetricMatrix_Allocate ( qcState->qptransformation->length0 ) ;
        SymmetricMatrix_Transform ( qcState->densityp->density, qcState->qptransformation, True, density ) ;
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Get the QC coordinates. */
        QCAtomContainer_GetCoordinates3 ( qcState->qcAtoms, coordinates3, True, &qcCoordinates3 ) ;

        /* . Electronic contribution. */
        GaussianBasis_Point_Electron ( qcState->qcAtoms, qcState->qcParameters, qcCoordinates3, points, density, potentials ) ;

        /* . Nuclear contribution. */
        znuclear = Real1DArray_Allocate   ( qcState->qcAtoms->natoms, NULL ) ;
        QCAtomContainer_GetNuclearCharges ( qcState->qcAtoms, qcState->qcParameters, znuclear ) ;
        GaussianBasis_Point_Nuclear       ( qcState->qcAtoms, znuclear, qcCoordinates3, points, potentials ) ;

        /* . Deallocate space. */
        Coordinates3_Deallocate    ( &qcCoordinates3 ) ;
        Real1DArray_Deallocate     ( &znuclear       ) ;
        SymmetricMatrix_Deallocate ( &density        ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the integrals for a MNDO method.
! . The procedure returns enuclear and the one- and two-electron integrals in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define MNDO_BLOCKSIZE 1024
# define MNDO_UNDERFLOW 1.0e-12
# define N1COEIS 45 /*. Equivalent to ( n * ( n + 1 ) ) / 2 with n = 9 for spd basis. */
void QCModelMNDO_Integrals ( const QCModelMNDO *self, QCModelMNDOState *qcState )
{
    if ( ( self != NULL ) && ( qcState != NULL ) )
    {
        auto Real  enuc, *xi, *xj ;
        auto Integer     i, ic, iqm, j, jc, jqm, ni, nj ;
        auto Real2DArray         *sij ;
        auto MNDOParameters *idata, *jdata ;
        auto Real1DArray         *e1b, *e2a ;

        /* . One-electron matrix.*/
        qcState->oneelectronmatrix = SymmetricMatrix_Allocate ( qcState->qcAtoms->nobasis ) ;
        SymmetricMatrix_Set_Zero ( qcState->oneelectronmatrix ) ;

        /* . Two-electron integrals. */
        qcState->twoelectronintegrals = BlockStorage_Allocate ( ) ;
        qcState->twoelectronintegrals->blocksize  = MNDO_BLOCKSIZE ;
        qcState->twoelectronintegrals->nindices16 = 4 ;
        qcState->twoelectronintegrals->QUNDERFLOW = True ;
        qcState->twoelectronintegrals->underflow  = MNDO_UNDERFLOW ;

        /* . Outer loop over atoms. */
        for ( iqm = 0 ; iqm < qcState->qcAtoms->natoms ; iqm++ )
        {
            /* . Get data for the atom. */
            ic    = qcState->qcAtoms->data[iqm].center ;
            idata = qcState->qcParameters->centers[ic].mndoparameters ;
            ni    = idata->norbitals ;
            xi    = Coordinates3_RowPointer ( qcState->qccoordinates3, iqm ) ;

            /* . The one-center kinetic energy terms. */
            SymmetricMatrix_Set_Diagonal_V ( qcState->oneelectronmatrix, qcState->qcAtoms->data[iqm].ostart, qcState->qcAtoms->data[iqm].nobasis, idata->uspd ) ;

            /* . The one-center TEIs. */
            MNDOIntegrals_AddInOneCenterTEIs ( idata, qcState->qcAtoms->data[iqm].ostart, qcState->twoelectronintegrals ) ;

            /* . Inner loop over atoms. */
            for ( jqm = 0 ; jqm < iqm ; jqm++ )
            {
               /* . Get data for the atom. */
               jc    = qcState->qcAtoms->data[jqm].center ;
               jdata = qcState->qcParameters->centers[jc].mndoparameters ;
               nj    = jdata->norbitals ;
               xj    = Coordinates3_RowPointer ( qcState->qccoordinates3, jqm ) ;

               /* . The resonance integrals. */
               sij = GaussianBasis_2Overlap ( qcState->qcParameters->centers[ic].orbitalbasis, xi, qcState->qcParameters->centers[jc].orbitalbasis, xj ) ;

               /* . Transform the integrals. */
               QCModelMNDO_TransformIntegrals ( &sij, qcState->qcParameters->centers[ic].orbitalbasis->c2o, qcState->qcParameters->centers[jc].orbitalbasis->c2o ) ;
               if ( sij != NULL )
               {
                  for ( i = 0 ; i < ni ; i++ )
                  {
                     for ( j = 0 ; j < nj ; j++ ) Real2DArray_Item ( sij, i, j ) *= 0.5e+00 * ( idata->beta[i] + jdata->beta[j] ) ;
                  }
/*
printf ( "\nOverlaps:\n" ) ;
Matrix_Print ( sij ) ;
*/
                  SymmetricMatrix_Increment_OB ( qcState->oneelectronmatrix, qcState->qcAtoms->data[iqm].ostart, qcState->qcAtoms->data[iqm].nobasis, qcState->qcAtoms->data[jqm].ostart, qcState->qcAtoms->data[jqm].nobasis, sij->data ) ;
/*
printf ( "\nPartial one-electron matrix:\n" ) ;
SymmetricMatrix_Print ( qcState->oneelectronmatrix ) ;
*/
                  Real2DArray_Deallocate ( &sij ) ;
               }

               /* . Two-center quantities. */
               e1b = Real1DArray_Allocate ( ( ni * ( ni + 1 ) ) / 2, NULL ) ;
               e2a = Real1DArray_Allocate ( ( nj * ( nj + 1 ) ) / 2, NULL ) ;
               MNDOIntegrals_MolecularFrame2CIntegrals ( idata, qcState->qcAtoms->data[iqm].ostart, xi, jdata, qcState->qcAtoms->data[jqm].ostart, xj, &enuc, e1b, e2a, qcState->twoelectronintegrals ) ;

               /* . Increment enuclear. */
               qcState->enuclear += enuc ;

               /* . Add in the electron-nuclear attraction terms to the one-electron matrix. */
               /* . Stop gap until use vector directly. */
               SymmetricMatrix_Increment_DB ( qcState->oneelectronmatrix, qcState->qcAtoms->data[iqm].ostart, qcState->qcAtoms->data[iqm].nobasis, e1b->data ) ;
               SymmetricMatrix_Increment_DB ( qcState->oneelectronmatrix, qcState->qcAtoms->data[jqm].ostart, qcState->qcAtoms->data[jqm].nobasis, e2a->data ) ;

               /* . Finish up. */
               Real1DArray_Deallocate ( &e1b ) ;
               Real1DArray_Deallocate ( &e2a ) ;
            }
        }

# ifdef PRINTMOPACOEIS
{
auto Integer i, ii, j, jj, n, order[1000] ;
auto Real v ;
auto SymmetricMatrix *comparison = NULL ;
for ( i = 0 ; i < 1000 ; i++ ) order[i] = -1000000 ;
for ( iqm = n = 0 ; iqm < qcState->qcAtoms->natoms ; iqm++ )
{
    ic    = qcState->qcAtoms->data[iqm].center ;
    idata = qcState->qcParameters->centers[ic].mndoparameters ;
    ni    = idata->norbitals ;
    order[n] = n ;
    if ( ni > 1 )
    {
        /* . p. */
        order[n+1] = n + 3 ; order[n+2] = n + 1 ; order[n+3] = n + 2 ;
        /* . d. */
        if ( ni > 4 ) { order[n+4] = n + 6 ; order[n+5] = n + 5 ; order[n+6] = n + 7 ; order[n+7] = n + 4 ; order[n+8] = n + 8 ; }
    }
    n += ni ;
/* . Parameters. */
    printf ( "\nParameters for element %d:\n", idata->atomicNumber ) ;
    for ( i = 0 ; i < 6  ; i++ ) printf ( "DD %2d = %20.15f\n", i+1, idata->ddp[i] ) ;
    for ( i = 0 ; i < 9  ; i++ ) printf ( "PO %2d = %20.15f\n", i+1, idata->po [i] ) ;
    for ( i = 0 ; i < ni ; i++ ) printf ( "BETA %2d %20.15f %20.15f\n", i, idata->beta[i], idata->beta0[i] ) ;
}
comparison = SymmetricMatrix_Allocate ( qcState->oneelectronmatrix->dimension ) ;
SymmetricMatrix_Set ( comparison, 0.0e+00 ) ;
for ( i = 0 ; i < n ; i++ )
{
    ii = order[i] ;
    for ( j = 0 ; j <= i ; j++ )
    {
        jj = order[j] ;
        v = SymmetricMatrix_Get_Component ( qcState->oneelectronmatrix, i, j ) ;
        SymmetricMatrix_Set_Component ( comparison, ii, jj, v ) ;
    }
}
SymmetricMatrix_Scale ( comparison, 27.21e+00 ) ;
printf ( "\nFull one-electron matrix:\n" ) ;
SymmetricMatrix_Print ( comparison ) ;
}
# endif

    }
}

# undef MNDO_BLOCKSIZE
# undef MNDO_UNDERFLOW

/*# undef __RANGE_CHECKING*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Localize the occupied orbitals using the approach in Mopac.
! . Fractional occupancies are ignored.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . The procedure maximizes (psi)**4 for each atom.
! . See "A new rapid method for orbital localisation", P.G. Perkins and J.J.P. Stewart, J.C.S. Faraday (II) 78, 285, (1982).
*/

# define IllDefinitionTolerance 1.0e-3
# define LocalizationTolerance  1.0e-7
# define MaximumIterations      100
# define RotationTolerance      1.0e-14

void QCModelMNDO_LocalizeOrbitals ( const QCModelMNDO *self, QCModelMNDOState *state, const Boolean doQ, const Integer startOrbital, const Integer stopOrbital, Status *status )
{
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        auto Boolean          isCandidate ;
        auto Integer          i, ibasis, ibf, ifirst, ii, ilast, index, iqm, iterations, j, jj, kbf, kfirst, klast, kqm, l, loop, n, numberAtoms, numberOccupied = 0, numberOrbitals = 0, numberToLocalize = 0, start = 0, stop = 0 ;
        auto Real             aij, bij, ca, co, coi, coj, dii, dij, djj, pi, pj, sa, sum, sumi, suml, x, xiiii, xiijj, xijij, xijjj, xjiii, xjjjj ;
        auto Integer1DArray  *indices ;
        auto QCOnePDM        *onepdm  ;
        auto Real1DArray     *eigenValues, ivector, jvector, lvector, *newEnergies, *oldEnergies = NULL, *psi1, *psi2 ;
        auto Real2DArray     *eigenVectors, *newOrbitals, *newVectors, *oldOrbitals = NULL, *oldVectors ;
        auto SymmetricMatrix *interaction ;

        /* . Get initial data. */
        if ( doQ ) onepdm = state->densityq ;
        else       onepdm = state->densityp ;
        if ( onepdm != NULL ) {  oldEnergies = onepdm->energies ; oldOrbitals = onepdm->orbitals ; numberOccupied = onepdm->numberOccupied ; numberOrbitals = onepdm->numberOrbitals ; }
        numberAtoms = state->qcAtoms->natoms ;

        /* . Get the range of orbitals to localize - defaults are all occupied orbitals. */
        if ( startOrbital < 0 ) start = 0 ;
        else                    start = Minimum ( startOrbital, numberOrbitals - 1 ) ;
        if ( stopOrbital  < 0 ) stop  = numberOccupied ;
        else                    stop  = Minimum ( stopOrbital , numberOrbitals     ) ;
        numberToLocalize = stop - start ;

        /* . Allocate space. */
        indices     = Integer1DArray_Allocate ( numberOrbitals, status ) ;
        psi1        = Real1DArray_Allocate    ( numberOrbitals, status ) ;
        psi2        = Real1DArray_Allocate    ( numberOrbitals, status ) ;
        newEnergies = Real1DArray_Clone ( oldEnergies, status ) ;
        newOrbitals = Real2DArray_Clone ( oldOrbitals, status ) ;

        /* . Check that everything is OK. */
        if  ( ( oldEnergies == NULL ) || ( oldOrbitals == NULL ) || ( numberOccupied <= 0 ) || ( numberOrbitals <= 0 ) || ( numberToLocalize <= 0 ) ) Status_Set ( status, Status_InvalidArgument ) ;
        else if ( ( indices == NULL ) || ( psi1 == NULL ) || ( psi2 == NULL ) || ( newEnergies == NULL ) || ( newOrbitals == NULL ) ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
        else
        {
            /* . Initial localization. */
            iterations = 0 ;
            do
            {
                iterations ++ ;
                sum = 0.0e+00 ;
                for ( i = start ; i < stop ; i++ )
                {
                    Real2DArray_ColumnSlice ( newOrbitals, i, &ivector, status ) ;
                    for ( j = start ; j < stop ; j++ )
                    {
                        Real2DArray_ColumnSlice ( newOrbitals, j, &jvector, status ) ;
                        if ( j != i )
                        {
                            /* . Initialization. */
                            xijjj = 0.0e+00 ;
                            xjiii = 0.0e+00 ;
                            xiiii = 0.0e+00 ;
                            xjjjj = 0.0e+00 ;
                            xijij = 0.0e+00 ;
                            xiijj = 0.0e+00 ;
                            Real1DArray_CopyTo ( &ivector, psi1, status ) ;
                            Real1DArray_CopyTo ( &jvector, psi2, status ) ;

                            /* . Loop over atoms. */
                           for ( kqm = 0 ; kqm < numberAtoms ; kqm++ )
                           {
                               dii    = 0.0e+00 ;
                               dij    = 0.0e+00 ;
                               djj    = 0.0e+00 ;
                               kfirst =          state->qcAtoms->data[kqm].ostart  ;
                               klast  = kfirst + state->qcAtoms->data[kqm].nobasis ;
                               for ( kbf = kfirst ; kbf < klast ; kbf++ )
                               {
                                   pi  = Real1DArray_Item ( psi1, kbf ) ;
                                   pj  = Real1DArray_Item ( psi2, kbf ) ;
                                   dij += pi * pj ;
                                   dii += pi * pi ;
                                   djj += pj * pj ;
                                }
                                xijjj += dij * djj ;
                                xjiii += dij * dii ;
                                xiiii += dii * dii ;
                                xjjjj += djj * djj ;
                                xijij += dij * dij ;
                                xiijj += dii * djj ;
                            }

                            /* . Determine the rotation. */
                            aij = xijij - ( xiiii + xjjjj - 2.e+00 * xiijj ) / 4.0e+00 ;
                            bij = xjiii - xijjj ;
                            ca  = sqrt ( aij * aij + bij * bij ) ;
                            sa  = aij + ca ;
                            if ( sa > RotationTolerance )
                            {
                                sum += sa ;
                                ca   = - aij / ca ;
                                ca   = ( 1.e+00 + sqrt ( ( 1.e+00 + ca ) / 2.e+00 ) ) / 2.0e+00 ;
                                if ( ( 2.e+00 * ca - 1.e+00 ) * bij < 0.e+00 ) ca = 1.e+00 - ca ;
                                sa = sqrt ( 1.e+00 - ca ) ;
                                ca = sqrt ( ca ) ;

                                /* . Update the orbitals. */
                                Real1DArray_Set ( &ivector, 0.0e+00 ) ;
                                Real1DArray_AddScaledArray ( &ivector,  ca, psi1, status ) ;
                                Real1DArray_AddScaledArray ( &ivector,  sa, psi2, status ) ;
                                Real1DArray_Set ( &jvector, 0.0e+00 ) ;
                                Real1DArray_AddScaledArray ( &jvector, -sa, psi1, status ) ;
                                Real1DArray_AddScaledArray ( &jvector,  ca, psi2, status ) ;
                            }
                        }
                    }
                }
            }
            while ( ( sum > LocalizationTolerance ) && ( iterations < MaximumIterations ) ) ;

            /* . Determine the final localization value (might be useful). */
            sum = 0.0e+00 ;
            for ( i = start ; i < stop ; i++ )
            {
                for ( iqm = 0 ; iqm < numberAtoms ; iqm++ )
                {
                    ifirst =          state->qcAtoms->data[iqm].ostart  ;
                    ilast  = ifirst + state->qcAtoms->data[iqm].nobasis ;
                    x      = 0.0e+00 ;
                    for ( ibf = ifirst ; ibf < ilast ; ibf++ ) x += pow ( Real2DArray_Item ( newOrbitals, ibf, i ), 2 ) ;
                    sum += x * x ;
                }
            }

            /* . Resolve any ill-definition of LMOs that involve the same atoms by ensuring that their interaction is zero. */
            for ( loop = start ; loop < stop ; loop++ )
            {
                /* . Is the LMO a candidate for degeneracy? */
                /* . Only LMOs that are 50% or 100% on an atom are potential candidates. */
                isCandidate = True ;
                for ( iqm = 0 ; iqm < numberAtoms ; iqm++ )
                {
                    ibasis = state->qcAtoms->data[iqm].nobasis ;
                    if ( ibasis > 1 )
                    {
                        ifirst = state->qcAtoms->data[iqm].ostart  ;
                        ilast  = ifirst + ibasis ;
                        sum    = 0.0e+00 ;
                        for ( ibf = ifirst ; ibf < ilast ; ibf++ ) sum += pow ( Real2DArray_Item ( newOrbitals, ibf, loop ), 2 ) ;
                        if ( ( sum <= IllDefinitionTolerance ) || ( sum >= ( 0.5e+00 - IllDefinitionTolerance ) ) ) continue ;
                        isCandidate = False ;
                        break ;
                    }
                }
                if ( ! isCandidate ) continue ;

                /* . The LMO is a candidate so identify any related LMOs. */
                Integer1DArray_Item ( indices, 0 ) = loop ; n = 1 ;
                for ( i = loop+1 ; i < stop ; i++ )
                {
                    for ( iqm = 0 ; iqm < numberAtoms ; iqm++ )
                    {
                        ifirst      =          state->qcAtoms->data[iqm].ostart  ;
                        ilast       = ifirst + state->qcAtoms->data[iqm].nobasis ;
                        isCandidate = True ;
                        sumi        = 0.0e+00 ;
                        suml        = 0.0e+00 ;
                        for ( ibf = ifirst ; ibf < ilast ; ibf++ )
                        {
                            sumi += pow ( Real2DArray_Item ( newOrbitals, ibf, i    ), 2 ) ;
                            suml += pow ( Real2DArray_Item ( newOrbitals, ibf, loop ), 2 ) ;
                        }
                        if ( fabs ( sumi - suml ) > IllDefinitionTolerance ) { isCandidate = False; break ; }
                    }
                    if ( isCandidate ) { Integer1DArray_Item ( indices, n ) = i ; n += 1 ; }
                }

                /* . There are similar orbitals. */
                if ( n > 1 )
                {
                    /* . Allocate space. */
                    eigenValues  = Real1DArray_Allocate      ( n,    status ) ;
                    eigenVectors = Real2DArray_Allocate      ( n, n, status ) ;
                    interaction  = SymmetricMatrix_AllocateN ( n,    status ) ;
                    newVectors   = Real2DArray_Allocate      ( numberOrbitals, n, status ) ;
                    oldVectors   = Real2DArray_Allocate      ( numberOrbitals, n, status ) ;

                    /* . Everything is OK. */
                    if ( ( eigenValues != NULL ) && ( eigenVectors != NULL ) && ( interaction != NULL ) && ( newVectors != NULL ) && ( oldVectors != NULL ) )
                    {
                        /* . Build small secular determinant. */
                        for ( ii = 0 ; ii < n ; ii++ )
                        {
                            i = Integer1DArray_Item ( indices, ii ) ;
                            Real2DArray_ColumnSlice ( newOrbitals,  i, &ivector, status ) ;
                            Real2DArray_ColumnSlice ( oldVectors,  ii, &jvector, status ) ;
                            Real1DArray_CopyTo ( &ivector, &jvector, status ) ;
                            for ( jj = 0 ; jj <= ii ; jj++ )
                            {
                                j = Integer1DArray_Item ( indices, jj ) ;
                                Real2DArray_ColumnSlice ( newOrbitals, j, &jvector, status ) ;
                                sum = 0.0e+00 ;
                                for ( l = start ; l < stop ; l++ )
                                {
                                    Real2DArray_ColumnSlice ( oldOrbitals, l, &lvector, status ) ;
                                    coi = Real1DArray_Dot ( &ivector, &lvector, status ) ;
                                    coj = Real1DArray_Dot ( &jvector, &lvector, status ) ;
                                    sum += coi * Real1DArray_Item ( oldEnergies, l ) * coj ;
                                }
                                SymmetricMatrix_Item ( interaction, ii, jj ) = sum ;
                            }
                        }

                        /* . Diagonalize. */
                        SymmetricMatrix_Diagonalize ( interaction, eigenValues, eigenVectors, NULL ) ;

                        /* . Rotate LMOs and put back in correct place. */
                        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, oldVectors, eigenVectors, 0.0e+00, newVectors, status ) ;
                        for ( ii = 0 ; ii < n ; ii++ )
                        {
                            i = Integer1DArray_Item ( indices, ii ) ;
                            Real2DArray_ColumnSlice ( newOrbitals,  i, &ivector, status ) ;
                            Real2DArray_ColumnSlice ( newVectors,  ii, &jvector, status ) ;
                            Real1DArray_CopyTo ( &jvector, &ivector, status ) ;
                        }
                    }
                    else Status_Set ( status, Status_MemoryAllocationFailure ) ;

                    /* . Finish up. */
                    Real1DArray_Deallocate     ( &eigenValues  ) ;
                    Real2DArray_Deallocate     ( &eigenVectors ) ;
                    Real2DArray_Deallocate     ( &newVectors   ) ;
                    Real2DArray_Deallocate     ( &oldVectors   ) ;
                    SymmetricMatrix_Deallocate ( &interaction  ) ;
                }
            }

            /* . Determine new energies. */
            for ( i = start ; i < stop ; i++ )
            {
                Real2DArray_ColumnSlice ( newOrbitals, i, &ivector, status ) ;
                sum = 0.0e+00 ;
                for ( j = start ; j < stop ; j++ )
                {
                    Real2DArray_ColumnSlice ( oldOrbitals, j, &jvector, status ) ;
                    co   = Real1DArray_Dot ( &ivector, &jvector, NULL ) ;
                    sum += co * co * Real1DArray_Item ( oldEnergies, j ) ;
                }
                Real1DArray_Item ( newEnergies, i ) = sum ;
            }

            /* . Reorder the new energies and orbitals and put back in the old arrays. */
            Real1DArray_SortIndex ( newEnergies, indices, status ) ;
            for ( i = start ; i < stop ; i++ )
            {
                index = Integer1DArray_Item ( indices, i ) ;
                Real1DArray_Item ( oldEnergies, i ) = Real1DArray_Item ( newEnergies, index ) ;
                Real2DArray_ColumnSlice ( newOrbitals, index, &ivector, status ) ;
                Real2DArray_ColumnSlice ( oldOrbitals, i    , &jvector, status ) ;
                Real1DArray_CopyTo ( &ivector, &jvector, status ) ;
            }
        }

        /* . Finish up. */
        Integer1DArray_Deallocate ( &indices     ) ;
        Real1DArray_Deallocate    ( &newEnergies ) ;
        Real1DArray_Deallocate    ( &psi1        ) ;
        Real1DArray_Deallocate    ( &psi2        ) ;
        Real2DArray_Deallocate    ( &newOrbitals ) ;
    }
}

# undef IllDefinitionTolerance
# undef LocalizationTolerance
# undef MaximumIterations
# undef RotationTolerance

/*----------------------------------------------------------------------------------------------------------------------------------
! . Lowdin charge calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_LowdinCharges ( const QCModelMNDO *self, QCModelMNDOState *qcState, const Boolean QSPIN, Real1DArray *qcCharges )
{
    /* . Initialization. */
    Real1DArray_Set ( qcCharges, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( qcCharges != NULL ) && ( qcState->qcAtoms != NULL ) && ( qcState->qcAtoms->natoms <= qcCharges->length ) )
    {
        auto Real ps, scale = 1.0e+00 ;
        auto Integer    iatom, ibf, istart, istop ;

# ifdef MNDOCI
        /* . Check for CI spin density. */
        if ( ( QSPIN ) && ( self->ciModel != NULL ) && ( qcState->cistate != NULL ) && ( qcState->cistate->spinDensity != NULL ) )
        {
            for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
            {
                istart = qcState->qcAtoms->data[iatom].ostart ;
                istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
                for ( ibf = istart, ps = 0.0e+00 ; ibf < istop ; ibf++ ) ps += SymmetricMatrix_Get_Component ( qcState->cistate->spinDensity, ibf, ibf ) ;
                Real1DArray_Item ( qcCharges, iatom ) = ps ;
            }
/*
if ( ( qcState->cistate->onepdma != NULL ) && ( qcState->cistate->onepdmb != NULL ) )
{
auto Real a, b ;
for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
{
istart = qcState->qcAtoms->data[iatom].ostart ;
istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
a = 0.0e+00 ;
b = 0.0e+00 ;
for ( ibf = istart ; ibf < istop ; ibf++ ) { a += SymmetricMatrix_Get_Component ( qcState->cistate->onepdma, ibf, ibf ) ;
b += SymmetricMatrix_Get_Component ( qcState->cistate->onepdmb, ibf, ibf ) ; }
printf ( "Charges : %10d %20.5f %20.5f\n", iatom, a, b ) ;
}
}
*/
        }
        /* . Other cases. */
        else if ( ! QSPIN || ( QSPIN && ( qcState->densityq != NULL ) ) )
        {
# endif
            /* . Initialization. */
            if ( QSPIN ) scale = -1.0e+00 ;
            else QCAtomContainer_GetNuclearCharges ( qcState->qcAtoms, qcState->qcParameters, qcCharges ) ;
            /* . densityp. */
            for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
            {
                istart = qcState->qcAtoms->data[iatom].ostart ;
                istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
                for ( ibf = istart, ps = 0.0e+00 ; ibf < istop ; ibf++ ) ps += SymmetricMatrix_Get_Component ( qcState->densityp->density, ibf, ibf ) ;
                Real1DArray_Item ( qcCharges, iatom ) -= ( scale * ps ) ;
            }
            /* . densityq. */
            if ( qcState->densityq != NULL )
            {
                for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
                {
                    istart = qcState->qcAtoms->data[iatom].ostart ;
                    istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
                    for ( ibf = istart, ps = 0.0e+00 ; ibf < istop ; ibf++ ) ps += SymmetricMatrix_Get_Component ( qcState->densityq->density, ibf, ibf ) ;
                    Real1DArray_Item ( qcCharges, iatom ) -= ps ;
                }
            }
# ifdef MNDOCI
        }
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Mayer bond orders.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_MayerBondOrders ( const QCModelMNDO *self, QCModelMNDOState *qcState, SymmetricMatrix *bondorders, Real1DArray *charges, Real1DArray *freevalence, Real1DArray *totalvalence )
{
    Real1DArray_Set     ( charges     , 0.0e+00 ) ;
    Real1DArray_Set     ( freevalence , 0.0e+00 ) ;
    Real1DArray_Set     ( totalvalence, 0.0e+00 ) ;
    SymmetricMatrix_Set ( bondorders  , 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( bondorders != NULL ) && ( charges != NULL ) && ( freevalence != NULL ) && ( totalvalence != NULL ) &&
                ( qcState->qcAtoms != NULL ) && ( bondorders->dimension == qcState->qcAtoms->natoms ) && ( charges->length == qcState->qcAtoms->natoms ) &&
                                               ( freevalence->length == qcState->qcAtoms->natoms ) && ( totalvalence->length == qcState->qcAtoms->natoms ) )
    {
        auto Real sum ;
        auto Integer    iatom, ibf, istart, istop, jatom, jbf, jstart, jstop ;

        /* . Convert to total and spin densities. */
        QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Total density. */
        for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
        {
            istart = qcState->qcAtoms->data[iatom].ostart ;
            istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
            for ( jatom = 0 ; jatom <= iatom ; jatom++ )
            {
                jstart = qcState->qcAtoms->data[jatom].ostart ;
                jstop  = jstart + qcState->qcAtoms->data[jatom].nobasis ;
                for ( ibf = istart, sum = 0.0e+00 ; ibf < istop ; ibf++ )
                {
                    for ( jbf = jstart ; jbf < jstop ; jbf++ ) sum += pow ( SymmetricMatrix_Get_Component ( qcState->densityp->density, ibf, jbf ), 2 ) ;
                }
                SymmetricMatrix_Set_Component ( bondorders, iatom, jatom, sum ) ;
            }
        }

        /* . Total valence - first contribution. */
        for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ ) Real1DArray_Item ( totalvalence, iatom ) = - SymmetricMatrix_Get_Component ( bondorders, iatom, iatom ) ;

        /* . Spin density. */
        if ( qcState->densityq != NULL )
        {
            for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
            {
                istart = qcState->qcAtoms->data[iatom].ostart ;
                istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
                for ( jatom = 0 ; jatom <= iatom ; jatom++ )
                {
                    jstart = qcState->qcAtoms->data[jatom].ostart ;
                    jstop  = jstart + qcState->qcAtoms->data[jatom].nobasis ;
                    for ( ibf = istart, sum = 0.0e+00 ; ibf < istop ; ibf++ )
                    {
                        for ( jbf = jstart ; jbf < jstop ; jbf++ ) sum += pow ( SymmetricMatrix_Get_Component ( qcState->densityq->density, ibf, jbf ), 2 ) ;
                    }
                    SymmetricMatrix_IncrementComponent ( bondorders, iatom, jatom, sum ) ;
                }
            }
        }

        /* . Convert densities back. */
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Calculate the sum of the bond orders for each atom and start the calculation of the free valence. */
        for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
        {
            for ( jatom = 0, sum = 0.0e+00 ; jatom < qcState->qcAtoms->natoms ; jatom++ ) sum += SymmetricMatrix_Get_Component ( bondorders, iatom, jatom ) ;
            Real1DArray_Item ( charges,     iatom ) = sum ;
            Real1DArray_Item ( freevalence, iatom ) = - sum + SymmetricMatrix_Get_Component ( bondorders, iatom, iatom ) ;
        }

        /* . Complete the calculation of atomic quantities. */
        Real1DArray_AddScaledArray ( totalvalence, 1.0e+00, charges, NULL ) ;
        Real1DArray_Scale ( charges, 0.5e+00 ) ;
        Real1DArray_AddScaledArray ( freevalence, 1.0e+00, totalvalence, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the characters of selected orbitals under a rotation (proper or improper).
!---------------------------------------------------------------------------------------------------------------------------------*/
Real1DArray *QCModelMNDO_OrbitalCharacters ( const QCModelMNDO *self, QCModelMNDOState *state, const Matrix33 *rotation, const Integer1DArray *mapping,
                                                                           const Boolean useDensityP, const Selection *orbitalIndices, Status *status )
{
    Real1DArray *characters = NULL ;
    if ( ( self != NULL ) && ( state != NULL ) && ( rotation != NULL ) && ( mapping != NULL ) && ( orbitalIndices != NULL ) )
    {
        auto QCOnePDM *qcOnePDM = NULL ;

        /* . Decide on the orbitals to use. */
        if ( useDensityP ) qcOnePDM = state->densityp ;
        else               qcOnePDM = state->densityq ;
        if ( ( qcOnePDM != NULL ) && ( qcOnePDM->orbitals != NULL ) )
        {
            auto Integer i, s ;
            auto Real1DArray  inOrbital, *outOrbital ;
            auto Real2DArray *orbitals   ;

            /* . Initialization. */
            orbitals   = qcOnePDM->orbitals ;
            characters = Real1DArray_Allocate ( orbitalIndices->nindices, NULL ) ;
            outOrbital = Real1DArray_Allocate ( Real2DArray_Length ( orbitals, 0 ), NULL ) ;

            /* . Get characters. */
            for ( s = 0 ; s < orbitalIndices->nindices ; s++ )
            {
                i = orbitalIndices->indices[s] ;
                Real2DArray_ColumnSlice ( orbitals, i, &inOrbital, status ) ;
                QCModelMNDO_RotateOrbital ( self, state, rotation, mapping, &inOrbital, outOrbital ) ;
                Real1DArray_Item ( characters, s ) = Real1DArray_Dot ( &inOrbital, outOrbital, status ) ;
            }

            /* . Finish up. */
            Real1DArray_Deallocate ( &outOrbital ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return characters ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the orbital energies and HOMO and LUMO indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real1DArray *QCModelMNDO_OrbitalEnergies ( const QCModelMNDO *self, QCModelMNDOState *qcState, const Boolean QALPHA, Integer *homo, Integer *lumo )
{
    Real1DArray *data = NULL ;
    if ( homo != NULL ) (*homo) = -1 ;
    if ( lumo != NULL ) (*lumo) = -1 ;
    if ( ( self != NULL ) && ( qcState != NULL ) )
    {
        auto QCOnePDM *density  = NULL ;

        /* . Get the correct orbital set. */
        if ( QALPHA ) density = qcState->densityp ;
        else          density = qcState->densityq ;

        /* . Get the data. */
        data = Real1DArray_Clone ( density->energies, NULL ) ;
        if ( homo != NULL ) (*homo) = QCOnePDM_HOMOIndex ( density ) ;
        if ( lumo != NULL ) (*lumo) = QCOnePDM_LUMOIndex ( density ) ;
    }
    return data ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the orbitals.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real2DArray *QCModelMNDO_Orbitals ( const QCModelMNDO *self, QCModelMNDOState *qcState, const Boolean QALPHA )
{
    Real2DArray *data = NULL ;
    if ( ( self != NULL ) && ( qcState != NULL ) )
    {
        auto QCOnePDM *density  = NULL ;

        /* . Get the correct orbital set. */
        if ( QALPHA ) density = qcState->densityp ;
        else          density = qcState->densityq ;

        /* . Get the data. */
        data = Real2DArray_Clone ( density->orbitals, NULL ) ;
    }
    return data ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Rotate an orbital by either a proper or an improper rotation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_RotateOrbital ( const QCModelMNDO *self, const QCModelMNDOState *state, const Matrix33 *rotation, const Integer1DArray *mapping, const Real1DArray *inOrbital, Real1DArray *outOrbital )
{
    if ( ( self != NULL ) && ( state != NULL ) && ( rotation != NULL ) && ( mapping != NULL ) && ( inOrbital != NULL ) && ( outOrbital != NULL ) )
    {
        auto Integer      iAtom, iFirstOrbital, jAtom, jFirstOrbital, numberOrbitals ;
        auto Real1DArray  inSlice, outSlice ;
        auto Real2DArray *dTransformation = NULL, *pTransformation = NULL ;

        /* . Initialization. */
        Real1DArray_Set ( outOrbital, 0.0e+00 ) ;

        /* . Find maximum number of orbitals per atom and calculate the appropriate transformation matrices. */
        numberOrbitals = 0 ;
        for ( iAtom = 0 ; iAtom < state->qcAtoms->natoms ; iAtom++ ) numberOrbitals = Maximum ( numberOrbitals, state->qcAtoms->data[iAtom].nobasis ) ;

        /* . p and d transformations. */
        if ( numberOrbitals > 1 )
        {
            auto Real r00, r0m, r0p, rm0, rmm, rmp, rp0, rpm, rpp ;

            /* . Get matrix elements. */
            r00 = Matrix33_Item ( rotation, 2, 2 ) ;
            r0p = Matrix33_Item ( rotation, 2, 0 ) ;
            r0m = Matrix33_Item ( rotation, 2, 1 ) ;
            rp0 = Matrix33_Item ( rotation, 0, 2 ) ;
            rpp = Matrix33_Item ( rotation, 0, 0 ) ;
            rpm = Matrix33_Item ( rotation, 0, 1 ) ;
            rm0 = Matrix33_Item ( rotation, 1, 2 ) ;
            rmp = Matrix33_Item ( rotation, 1, 0 ) ;
            rmm = Matrix33_Item ( rotation, 1, 1 ) ;

            /* . p transformation - 10, 11, 1-1 = z, x, y. */
            pTransformation = Real2DArray_Allocate ( 3, 3, NULL ) ;
            Real2DArray_Item ( pTransformation, 0, 0 ) = r00 ;
            Real2DArray_Item ( pTransformation, 0, 1 ) = r0p ;
            Real2DArray_Item ( pTransformation, 0, 2 ) = r0m ;
            Real2DArray_Item ( pTransformation, 1, 0 ) = rp0 ;
            Real2DArray_Item ( pTransformation, 1, 1 ) = rpp ;
            Real2DArray_Item ( pTransformation, 1, 2 ) = rpm ;
            Real2DArray_Item ( pTransformation, 2, 0 ) = rm0 ;
            Real2DArray_Item ( pTransformation, 2, 1 ) = rmp ;
            Real2DArray_Item ( pTransformation, 2, 2 ) = rmm ;

            /* . d transformation - 20, 21, 2-1, 22, 2-2. */
            if ( numberOrbitals > 4 )
            {
                auto Real sqrt3 ;
                sqrt3 = sqrt ( 3.0e+00 ) ;
                dTransformation = Real2DArray_Allocate ( 5, 5, NULL ) ;
                Real2DArray_Item ( dTransformation, 0, 0 ) = ( 3.0e+00 * r00 * r00 - 1.0e+00 ) / 2.0e+00 ;
                Real2DArray_Item ( dTransformation, 0, 1 ) =  sqrt3 * r00 * r0p ;
                Real2DArray_Item ( dTransformation, 0, 2 ) =  sqrt3 * r00 * r0m ;
                Real2DArray_Item ( dTransformation, 0, 3 ) =  sqrt3 * ( r0p * r0p - r0m * r0m ) / 2.0e+00 ;
                Real2DArray_Item ( dTransformation, 0, 4 ) =  sqrt3 * r0p * r0m ;
                Real2DArray_Item ( dTransformation, 1, 0 ) =  sqrt3 * rp0 * r00 ;
                Real2DArray_Item ( dTransformation, 1, 1 ) =  rpp * r00 + rp0 * r0p ;
                Real2DArray_Item ( dTransformation, 1, 2 ) =  rpm * r00 + rp0 * r0m ;
                Real2DArray_Item ( dTransformation, 1, 3 ) =  rpp * r0p - rpm * r0m ;
                Real2DArray_Item ( dTransformation, 1, 4 ) =  rpp * r0m + r0p * rpm ;
                Real2DArray_Item ( dTransformation, 2, 0 ) =  sqrt3 * rm0 * r00 ;
                Real2DArray_Item ( dTransformation, 2, 1 ) =  rmp * r00 + r0p * rm0 ;
                Real2DArray_Item ( dTransformation, 2, 2 ) =  rmm * r00 + r0m * rm0 ;
                Real2DArray_Item ( dTransformation, 2, 3 ) =  rmp * r0p - rmm * r0m ;
                Real2DArray_Item ( dTransformation, 2, 4 ) =  rmp * r0m + r0p * rmm ;
                Real2DArray_Item ( dTransformation, 3, 0 ) =  sqrt3 * ( rp0 * rp0 - rm0 * rm0 ) / 2.0e+00 ;
                Real2DArray_Item ( dTransformation, 3, 1 ) =  rpp * rp0 - rmp * rm0 ;
                Real2DArray_Item ( dTransformation, 3, 2 ) =  rpm * rp0 - rmm * rm0 ;
                Real2DArray_Item ( dTransformation, 3, 3 ) =  ( rpp * rpp + rmm * rmm - rmp * rmp - rpm * rpm ) / 2.0e+00 ;
                Real2DArray_Item ( dTransformation, 3, 4 ) =  rpp * rpm - rmp * rmm ;
                Real2DArray_Item ( dTransformation, 4, 0 ) =  sqrt3 * rp0 * rm0 ;
                Real2DArray_Item ( dTransformation, 4, 1 ) =  rpp * rm0 + rmp * rp0 ;
                Real2DArray_Item ( dTransformation, 4, 2 ) =  rpm * rm0 + rmm * rp0 ;
                Real2DArray_Item ( dTransformation, 4, 3 ) =  rpp * rmp - rpm * rmm ;
                Real2DArray_Item ( dTransformation, 4, 4 ) =  rpp * rmm + rmp * rpm ;
            }
        }

        /* . Loop over atoms and rotate each block of orbitals separately. */
        for ( iAtom = 0 ; iAtom < state->qcAtoms->natoms ; iAtom++ )
        {
            numberOrbitals = state->qcAtoms->data[iAtom].nobasis ;
            if ( numberOrbitals > 0 )
            {
                jAtom          = Integer1DArray_Item ( mapping, iAtom ) ;
                iFirstOrbital  = state->qcAtoms->data[iAtom].ostart  ;
                jFirstOrbital  = state->qcAtoms->data[jAtom].ostart  ;
                /* . s. */
                Real1DArray_Item ( outOrbital, jFirstOrbital ) = Real1DArray_Item ( inOrbital, iFirstOrbital ) ;
                /* . p. */
                if ( numberOrbitals > 1 )
                {
                    Real1DArray_Slice (  inOrbital, iFirstOrbital+1, iFirstOrbital+4, 1,  &inSlice, NULL ) ;
                    Real1DArray_Slice ( outOrbital, jFirstOrbital+1, jFirstOrbital+4, 1, &outSlice, NULL ) ;
                    Real2DArray_VectorMultiply ( False, 1.0e+00, pTransformation, &inSlice, 0.0e+00, &outSlice, NULL ) ;
                }
                /* . d. */
                if ( numberOrbitals > 4 )
                {
                    Real1DArray_Slice (  inOrbital, iFirstOrbital+4, iFirstOrbital+9, 1,  &inSlice, NULL ) ;
                    Real1DArray_Slice ( outOrbital, jFirstOrbital+4, jFirstOrbital+9, 1, &outSlice, NULL ) ;
                    Real2DArray_VectorMultiply ( False, 1.0e+00, dTransformation, &inSlice, 0.0e+00, &outSlice, NULL ) ;
                }
            }
        }

        /* . Finish up. */
        Real2DArray_Deallocate ( &dTransformation ) ;
        Real2DArray_Deallocate ( &pTransformation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM Fock contributions - Lowdin charges.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_QCMMFockLowdin ( const QCModelMNDO *self, QCModelMNDOState *qcState, Real *eelectronic )
{
    if ( ( self != NULL ) && ( qcState != NULL ) && ( qcState->qcmmstate != NULL ) )
    {
        auto Real       p ;
        auto Integer          iatom, ibf, istart, istop ;
        auto Real1DArray *work ;

        /* . Get the Lowdin charges. */
        QCModelMNDO_LowdinCharges ( self, qcState, False, qcState->qcmmstate->qcCharges ) ;

        /* . Get the QC/MM energy. */
        qcState->eqcmm = Real1DArray_Dot ( qcState->qcmmstate->qcCharges, qcState->qcmmstate->qcmmPotentials, NULL ) ;

        /* . QC/QC. */
        if ( qcState->qcmmstate->qcqcPotentials != NULL )
        {
            /* . Get F * q. */
            work = Real1DArray_Allocate ( qcState->qcmmstate->qcCharges->length, NULL ) ;
            SymmetricMatrix_VectorMultiply ( qcState->qcmmstate->qcqcPotentials, qcState->qcmmstate->qcCharges, work, NULL ) ;

            /* . Get the QC/QC energy. */
            qcState->eqcqc = Real1DArray_Dot ( qcState->qcmmstate->qcCharges, work, NULL ) ;

            /* . Scale work and add the QC/MM potentials. */
            Real1DArray_Scale ( work, 2.0e+00 ) ;
            Real1DArray_Add   ( work, qcState->qcmmstate->qcmmPotentials, NULL ) ;
        }
        /* . QC/MM only. */
        else work = qcState->qcmmstate->qcmmPotentials ;

        /* . Add in the contributions to the Fock matrices. */
        /* . densityp. */
        for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
        {
            p = Real1DArray_Item ( work, iatom ) ;
            istart = qcState->qcAtoms->data[iatom].ostart ;
            istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
            for ( ibf = istart ; ibf < istop ; ibf++ ) SymmetricMatrix_IncrementComponent ( qcState->densityp->fock, ibf, ibf, - p ) ;
        }
        /* . densityq. */
        if ( qcState->densityq != NULL )
        {
            for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
            {
                p = Real1DArray_Item ( work, iatom ) ;
                istart = qcState->qcAtoms->data[iatom].ostart ;
                istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
                for ( ibf = istart ; ibf < istop ; ibf++ ) SymmetricMatrix_IncrementComponent ( qcState->densityq->fock, ibf, ibf, - p ) ;
            }
        }

        /* . Clear up. */
        if ( qcState->qcmmstate->qcqcPotentials != NULL ) Real1DArray_Deallocate ( &work ) ;

        /* . Add in the energy contributions. */
        if ( eelectronic != NULL ) (*eelectronic) += ( qcState->eqcmm + qcState->eqcqc ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform a matrix of integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDO_TransformIntegrals ( Real2DArray **integrals, const Real2DArray *ic2o, const Real2DArray *jc2o )
{

    if ( ( (*integrals) != NULL ) && ( ic2o != NULL ) && ( jc2o != NULL ) )
    {
        auto Integer nci, noi, noj ;
        auto Real2DArray *temp1 = NULL, *temp2 = NULL ;
        /* . Get dimensions. */
        nci = ic2o->length0 ;
        noi = ic2o->length1 ;
/*        ncj = jc2o->length0 ; */
        noj = jc2o->length1 ;
        /* . Do the multiplication. */
        temp1 = Real2DArray_Allocate ( nci, noj, NULL ) ;
        temp2 = Real2DArray_Allocate ( noi, noj, NULL ) ;
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, (*integrals), jc2o, 0.0e+00, temp1, NULL ) ;
        Real2DArray_MatrixMultiply ( True,  False, 1.0e+00, ic2o,        temp1, 0.0e+00, temp2, NULL ) ;
/*
printf ( "\nDEBUGGING:\n" ) ;
Matrix_Print ( (*integrals) ) ;
Matrix_Print ( ic2o ) ;
Matrix_Print ( jc2o ) ;
Matrix_Print ( temp1 ) ;
Matrix_Print ( temp2 ) ;
*/
        /* . Clear up. */
        Real2DArray_Deallocate ( integrals ) ;
        Real2DArray_Deallocate ( &temp1    ) ;
        (*integrals) = temp2 ;
    }
}

/*==================================================================================================================================
! . Local procedures.
!=================================================================================================================================*/
# ifdef MNDOCI
/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the alpha or beta contribution to an element of the CI state transformation matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CIStateCharacterDeterminant ( const Integer activeElectrons, const Boolean includeCoreOrbitals, const Integer numberCoreOrbitals,
                                                                  const Integer1DArray *iActiveIndices, const Integer1DArray *jActiveIndices,
                                                                const Real2DArray *orbitalTransformation, Real2DArray *work, Status *status )
{
    Real determinant = 1.0e+00 ;
    if ( ( activeElectrons > 0 ) && ( work != NULL ) )
    {
        auto Integer     coreIncrement, i, iActive, iIndex, j, jActive, jIndex, numberActive, totalElectrons ;
        auto Real2DArray matrix ;

        /* . Get space. */
        totalElectrons = activeElectrons ;
        if ( includeCoreOrbitals ) totalElectrons += numberCoreOrbitals ;
        Real2DArray_Slice ( work, 0, totalElectrons, 1, 0, totalElectrons, 1, &matrix, status ) ;
        Real2DArray_Set ( &matrix, 0.0e+00 ) ;

        /* . Cores. */
        if ( includeCoreOrbitals )
        {
            coreIncrement = numberCoreOrbitals ;
            for ( i = 0 ; i < numberCoreOrbitals ; i++ )
            {
                for ( j = 0 ; j < numberCoreOrbitals ; j++ ) Real2DArray_Item ( &matrix, i, j ) = Real2DArray_Item ( orbitalTransformation, i, j ) ;
            }
        }
        else coreIncrement = 0 ;

        /* . Active space. */
        numberActive = Integer1DArray_Length ( iActiveIndices ) ;
        for ( i = coreIncrement, iActive = 0 ; iActive < numberActive ; iActive++ )
        {
            iIndex = Integer1DArray_Item ( iActiveIndices, iActive ) ;
            if ( iIndex > 0 )
            {
                for ( j = coreIncrement, jActive = 0 ; jActive < numberActive ; jActive++ )
                {
                    jIndex = Integer1DArray_Item ( jActiveIndices, jActive ) ;
                    if ( jIndex > 0 )
                    {
                        Real2DArray_Item ( &matrix, i, j ) = Real2DArray_Item ( orbitalTransformation, iActive + coreIncrement, jActive + coreIncrement ) ;
                        j++ ;
                    }
                }
                i++ ;
            }
        }

        /* . Get the determinant. */
        determinant = Real2DArray_Determinant ( &matrix, status ) ;
    }
    return determinant ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the Qs -> Pc transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void MakeQsPcTransformation ( const QCModelMNDO *self, QCModelMNDOState *qcState, const Coordinates3 *coordinates3 )
{
    if ( ( self != NULL ) && ( qcState != NULL ) && ( coordinates3 != NULL ) )
    {
        auto Integer          n ;
        auto Coordinates3    *qccoordinates3 = NULL ;
        auto Real1DArray     *overlapeigenvalues = NULL ;
        auto Real2DArray     *overlapeigenvectors = NULL ;
        auto SymmetricMatrix *scratch = NULL, *overlap = NULL ;

        /* . Deallocate. */
        Real2DArray_Deallocate ( &(qcState->qptransformation) ) ;

        /* . Get the QC coordinates. */
        QCAtomContainer_GetCoordinates3 ( qcState->qcAtoms, coordinates3, True, &qccoordinates3 ) ;

        /* . Get the overlap as C*C and transform it to S*S. */
        GaussianBasis_Kinetic_2Overlap ( qcState->qcAtoms, qcState->qcParameters, qccoordinates3, &scratch, &overlap ) ;
        SymmetricMatrix_Transform_In_Place ( overlap, qcState->c2o ) ;

        /* . Get the inverse square root overlap. */
        n = overlap->dimension ;
        overlapeigenvectors = Real2DArray_Allocate ( n, n, NULL ) ;
        overlapeigenvalues  = Real1DArray_Allocate ( n   , NULL ) ;
        SymmetricMatrix_Diagonalize ( overlap, overlapeigenvalues, overlapeigenvectors, NULL ) ;
        SymmetricMatrix_InverseSquareRoot ( overlap, overlapeigenvalues, overlapeigenvectors ) ;

        /* . Get the transformation. */
        qcState->qptransformation = Real2DArray_Allocate ( qcState->c2o->length0, n, NULL ) ;
        SymmetricMatrix_PreMatrixMultiply ( overlap, qcState->c2o, False, qcState->qptransformation, NULL ) ;

        /* . Finish up. */
        Coordinates3_Deallocate    ( &qccoordinates3      ) ;
        Real1DArray_Deallocate     ( &overlapeigenvalues  ) ;
        Real2DArray_Deallocate     ( &overlapeigenvectors ) ;
        SymmetricMatrix_Deallocate ( &overlap             ) ;
        SymmetricMatrix_Deallocate ( &scratch             ) ;
    }
}
