/*------------------------------------------------------------------------------
! . File      : MNDOCIState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . The state for a MNDO CI calculation.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "CGLinearEquationSolver.h"
# include "DefineStatements.h"
# include "Integer2DArray.h"
# include "Macros.h"
# include "Memory.h"
# include "MNDOCIState.h"
# include "MNDOIntegrals.h"

/*==================================================================================================================================
! . Local procedures declarations.
!=================================================================================================================================*/
static void CIDensity_Diagonal        ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Real aiaj, SymmetricMatrix *onepdma, SymmetricMatrix *onepdmb, DoubleSymmetricMatrix *twopdm ) ;
static void CIDensity_OneAlphaOneBeta ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Integer1DArray *jalphas, const Integer1DArray *jbetas, const Real aiaj, DoubleSymmetricMatrix *twopdm ) ;
static void CIDensity_OneOrbital      ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Integer1DArray *jalphas, const Real aiaj, SymmetricMatrix *onepdm, DoubleSymmetricMatrix *twopdm ) ;
static void CIDensity_TwoOrbitals     ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *jalphas, const Real aiaj, DoubleSymmetricMatrix *twopdm ) ;

static void CIFIT_TransformIndex1    ( const Real2DArray *mos, const RealNDArray *tei234, DoubleSymmetricMatrix *moteis ) ;
static void CIFIT_TransformIndex2    ( const Real2DArray *mos, const Real2DArray *tei34, RealNDArray *tei234 ) ;
static void CIFIT_TransformIndices34 ( const Real2DArray *mos, BlockStorage *twoelectronintegrals, Real2DArray *tei34 ) ;

static void    CIGradient_Allocate             ( MNDOCIState *self, const Boolean doGradients, Status *status ) ;
static void    CIGradient_ApplyCPHFMatrix      ( MNDOCIState *self, const Boolean addDiagonal, const Integer n1, const Integer2DArray *in1, const Integer n2, const Integer2DArray *in2, const Real1DArray *b, Real1DArray *x ) ;
static void    CIGradient_CalculateCPHFVectors ( MNDOCIState *self ) ;
static void    CIGradient_CPHFTransform        ( const Integer n1, const Integer2DArray *in1, const Real1DArray *x1, const Integer n2, const Integer2DArray *in2, const Real1DArray *x2,
                                                                                    const Real2DArray *orbitals, const Boolean doScale, SymmetricMatrix *work, SymmetricMatrix *z ) ;
static void    CIGradient_Finalize             ( MNDOCIState *self ) ;
static void    CIGradient_Initialize           ( MNDOCIState *self ) ;
static void    CIGradient_SolveCPHFEquations   ( MNDOCIState *self ) ;
static void    CIGradient_ApplyMatrix          ( void *object, Real1DArray *x, Real1DArray *y ) ;
static void    CIGradient_ApplyPreconditioner  ( void *object, Real1DArray *x, Real1DArray *y ) ;

static Real CIMatrix_Diagonal        ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const SymmetricMatrix *fcore, const DoubleSymmetricMatrix *moteis  ) ;
static Real CIMatrix_OneAlphaOneBeta ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Integer1DArray *jalphas, const Integer1DArray *jbetas, const DoubleSymmetricMatrix *moteis ) ;
static Real CIMatrix_OneOrbital      ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Integer1DArray *jalphas, const SymmetricMatrix *fcore, const DoubleSymmetricMatrix *moteis ) ;
static Real CIMatrix_TwoOrbitals     ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *jalphas, const DoubleSymmetricMatrix *moteis ) ;

static Integer2DArray *CISetup_MakePermutations ( const Integer m, const Integer n, Status *status ) ;
static void            CISetup_MakeSPQR         ( MNDOCIState *self, Status *status ) ;
static Integer         Factorial                ( const Integer n ) ;

static void CITemporary_Allocate   ( MNDOCIState *self, Status *status ) ;
static void CITemporary_Deallocate ( MNDOCIState *self, const Boolean keepWavefunction ) ;
static void CITemporary_Initialize ( MNDOCIState *self ) ;

/* . Davidson methods. */
static void JDApplyMatrixSparse   ( void *x, void *y, Integer *blockSize, struct primme_params *primme ) ;
static void JDApplyPreconditioner ( void *x, void *y, Integer *blockSize, struct primme_params *primme ) ;

static void MNDOCIConfiguration_AllocateAlphasBetas ( MNDOCIConfiguration *self, const Integer nactive, Status *status ) ;
static void MNDOCIConfiguration_Deallocate          ( MNDOCIConfiguration *self ) ;
static void MNDOCIConfiguration_Initialize          ( MNDOCIConfiguration *self ) ;

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation - basic.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOCIState *MNDOCIState_Allocate ( const Integer nactive, const Integer nconfigurations, Status *status )
{
    MNDOCIState *self = NULL ;
    Status       localstatus ;
    Status_Set ( &localstatus, Status_Continue ) ;
    if ( ( nactive < 0 ) || ( nconfigurations < 0 ) ) Status_Set ( &localstatus, Status_InvalidArgument ) ;
    else
    {
        MEMORY_ALLOCATE ( self, MNDOCIState ) ;
        if ( self != NULL )
        {
            self->nactive         = nactive         ;
            self->nconfigurations = nconfigurations ;
            if ( nconfigurations > 0 )
            {
                MEMORY_ALLOCATEARRAY ( self->configurations, nconfigurations, MNDOCIConfiguration ) ;
                if ( self->configurations == NULL ) Status_Set ( &localstatus, Status_MemoryAllocationFailure ) ;
                else
                {
                    auto Integer i ;
                    for ( i = 0 ; i < nconfigurations ; i++ )   MNDOCIConfiguration_Initialize          ( &(self->configurations[i]) ) ;
                    for ( i = 0 ; i < nconfigurations ; i++ ) { MNDOCIConfiguration_AllocateAlphasBetas ( &(self->configurations[i]), self->nactive, &localstatus ) ; if ( ! Status_OK ( &localstatus ) ) break ; }
                }
            }
            /* . Remaining data. */
            CIGradient_Initialize  ( self ) ;
            CITemporary_Initialize ( self ) ;
        }
    }
    if ( ! Status_OK ( &localstatus ) ) MNDOCIState_Deallocate ( &self ) ;
    Status_Set ( status, localstatus ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a quantity of the form Kpa = Sum_qrs Gamma_pqrs TEI234_aqrs: a runs over all AOs.
! . This is then transformed to the MO basis and scaled by 2.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIState_CalculateKPA ( MNDOCIState *self )
{
    if ( self != NULL )
    {
        auto Integer i, nactive, nbasis, p, pq, q, r, rs, s ;
        auto Real    t ;

        /* . Initialization. */
        nactive = self->nactive ;
        nbasis  = self->kpa->length0 ;
        Real2DArray_Set ( self->kpa, 0.0e+00 ) ;

        /* . Loop over MOs. */
        for ( p = 0 ; p < nactive ; p++ )
        {
            for ( q = 0 ; q < nactive ; q++, pq++ )
            {
                for ( r = rs = 0 ; r < nactive ; r++ )
                {
                    for ( s = 0 ; s <= r ; s++, rs++ )
                    {
                        if ( r == s ) t = DoubleSymmetricMatrix_GetItem ( self->twopdm, p, q, r, r, NULL ) ;
                        else          t = DoubleSymmetricMatrix_GetItem ( self->twopdm, p, q, r, s, NULL ) + DoubleSymmetricMatrix_GetItem ( self->twopdm, p, q, s, r, NULL ) ;
                        /* . This is an increment so ultimately should use slices. */
                        for ( i = 0 ; i < nbasis ; i++ ) Real2DArray_Item ( self->kpa, i, p ) += t * RealNDArray_Item3D ( self->motei234, i, q, rs ) ;
                    }
                }
            }
        }

        /* . Transform to the full MO basis. */
        Real2DArray_MatrixMultiply ( True, False, 1.0e+00, self->orbitals, self->kpa, 0.0e+00, self->kpaMO, NULL ) ;
        Real2DArray_Scale ( self->kpaMO, 2.0e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the Z-matrix for a gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIState_CalculateZMatrix ( MNDOCIState *self )
{
    if ( ( self != NULL ) && ( self->doGradients ) )
    {
        /* . Calculate the CPHF vectors. */
        CIGradient_CalculateCPHFVectors ( self ) ;

        /* . Solve for zNR - put back in qNR. */
        CIGradient_SolveCPHFEquations ( self ) ;

        /* . Extract Z and convert to A.O. basis. */
        CIGradient_CPHFTransform ( self->numberNonRedundant, self->indicesNR, self->qNR, self->numberRedundant, self->indicesR, self->qR, self->orbitals, False, self->work1, self->zMatrix ) ;

# ifdef DEBUGMNDOCIGRADIENTS
printf ( "\nZ-matrix:\n" ) ;
SymmetricMatrix_Print ( self->zMatrix ) ;
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIState_Deallocate ( MNDOCIState **self )
{
    if ( (*self) != NULL )
    {
        if ( (*self)->configurations != NULL )
        {
            auto Integer i ;
            for ( i = 0 ; i < (*self)->nconfigurations ; i++ ) MNDOCIConfiguration_Deallocate ( &((*self)->configurations[i]) ) ;
            MEMORY_DEALLOCATE ( (*self)->configurations ) ;
        }
        CITemporary_Deallocate ( (*self), True ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Diagonalize CI matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIState_DiagonalizeCIMatrix ( MNDOCIState *self, const JDEigenvalueSolver *eigenvalueSolver, Status *status )
{
    if ( self != NULL )
    {
        auto Boolean doCheck ;

        /* . Options. */
        doCheck = ( self->doFull && self->doSparse ) ;

        /* . Sparse case. */
        if ( self->doSparse ) JDEigenvalueSolver_Solve ( eigenvalueSolver, self->eigenvalueSolverState, &(self->eigenvalueSolverReport) ) ;

        /* . Full case. */
        if ( self->doFull )
        {
            auto Real1DArray *eigenvalues  ;
            auto Real2DArray *eigenvectors ;
            if ( doCheck ) { eigenvalues = self->checkEnergies ; eigenvectors = self->checkVectors ; }
            else           { eigenvalues = self->ciEnergies    ; eigenvectors = self->ciVectors    ; }
            SymmetricMatrix_DiagonalizePartial ( self->ciMatrixFull, 0, self->numberOfStates, eigenvalues, eigenvectors, status ) ;
        }

        /* . Check solution. */
        if ( doCheck ) JDEigenvalueSolver_CheckSolution ( eigenvalueSolver, self->eigenvalueSolverState, self->checkEnergies, &(self->eigenvalueSolverReport) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Finalization after an energy and gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIState_Finalize ( MNDOCIState *self, const Boolean keepWavefunction )
{
    if ( self != NULL )
    {
        /* . Deallocate space. */
        CITemporary_Deallocate ( self, keepWavefunction ) ;
        CIGradient_Finalize    ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . A no-nonsense four index transformation for small numbers of MOs only.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIState_FourIndexTransformation ( MNDOCIState *self )
{
    if ( self != NULL )
    {
        CIFIT_TransformIndices34 ( &(self->activemos), self->twoelectronintegrals, self->motei34  ) ;
        CIFIT_TransformIndex2    ( &(self->activemos), self->motei34             , self->motei234 ) ;
        CIFIT_TransformIndex1    ( &(self->activemos), self->motei234            , self->moteis   ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an estimate of the sparsity of the CI matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIState_GetCIMatrixSparsity ( MNDOCIState *self )
{
    if ( self != NULL )
    {
        auto Integer         i, j, k, na, nai, naj, nactive, nb, nOff = 0 ;
        auto Integer1DArray *ialphas, *ibetas, *jalphas, *jbetas ;

        /* . Double loop over configurations. */
        nactive = self->nactive ;
        for ( i = 0 ; i < self->nconfigurations ; i++ )
        {
            nai     = self->configurations[i].nalphas ;
            ialphas = self->configurations[i].alphas  ;
            ibetas  = self->configurations[i].betas   ;
            for ( j = 0 ; j < i ; j++ )
            {
                naj     = self->configurations[j].nalphas ;
                jalphas = self->configurations[j].alphas  ;
                jbetas  = self->configurations[j].betas   ;

                /* . Skip if there are different numbers of alpha orbitals in the two configurations. */
                if ( nai != naj ) continue ;

                /* . Find the differences in the numbers of alpha and beta orbitals including positional information. */
                for ( k = na = nb = 0 ; k < nactive ; k++ )
                {
                    na += abs ( Integer1DArray_Item ( ialphas, k ) - Integer1DArray_Item ( jalphas, k ) ) ;
                    nb += abs ( Integer1DArray_Item ( ibetas,  k ) - Integer1DArray_Item ( jbetas,  k ) ) ;
                }

                /* . Skip if more than two orbitals are different. */
                if ( ( na + nb ) <= 4 ) nOff++ ;
            }
        }

        /* . Finish up. */
        self->ciMatrixNonZero  = ( self->nconfigurations + nOff ) ;
        self->ciMatrixSparsity = 100.0e+00 * ( 1.0e+00 - ( ( Real ) ( self->nconfigurations + 2 * nOff ) ) / ( ( Real ) ( self->nconfigurations * self->nconfigurations ) ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization for an energy and gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean MNDOCIState_Initialize ( MNDOCIState *self, const Boolean doGradients, QCOnePDM *densityp )
{
    Boolean isOK = True ;
    if ( self != NULL )
    {
        auto Integer i ;
        auto Status localstatus ;

        /* . Allocate space - energy. */
        Status_Set ( &localstatus, Status_Continue ) ;
        CITemporary_Allocate ( self, &localstatus ) ;
        isOK = Status_OK ( &localstatus ) ;

        /* . Set up some data structures from the SCF calculation. */
        self->energies    = densityp->energies    ;
        self->occupancies = densityp->occupancies ;
        self->orbitals    = densityp->orbitals    ;
        Real2DArray_Slice ( self->orbitals, 0, self->norbitals, 1, self->ncore, (self->ncore+self->nactive), 1, &(self->activemos), NULL ) ;

        /* . Allocate space - gradients. */
        CIGradient_Allocate ( self, doGradients, &localstatus ) ;
        isOK = ( isOK && Status_OK ( &localstatus ) ) ;

        /* . Various assignments. */
        Real1DArray_Set ( self->ocore, 0.0e+00 ) ;
        if ( self->ocore != NULL ) { for ( i = 0 ; i < self->ncore ; i++ ) Real1DArray_Item ( self->ocore, i ) = 2.0e+00 ; }
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the CI matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIState_MakeCIMatrix ( MNDOCIState *self )
{
    if ( self != NULL )
    {
        auto Boolean         doFull, doSparse ;
        auto Integer         i, j, k, na, nai, naj, nactive, nb ;
        auto Integer1DArray *ialphas, *ibetas, *jalphas, *jbetas ;
        auto Real            value = 0.0e+00 ;

# ifdef DEBUGMNDOCI
auto Real test = 0.0e+00 ;
# endif

        /* . Options. */
        doFull   = self->doFull   ;
        doSparse = self->doSparse ;

        /* . Initialization. */
        if ( doFull   ) SymmetricMatrix_Set ( self->ciMatrixFull, 0.0e+00 ) ;
        if ( doSparse ) SparseSymmetricMatrix_Clear ( self->ciMatrixSparse ) ;

        /* . Double loop over configurations. */
        nactive = self->nactive ;
        for ( i = 0 ; i < self->nconfigurations ; i++ )
        {
            nai     = self->configurations[i].nalphas ;
            ialphas = self->configurations[i].alphas  ;
            ibetas  = self->configurations[i].betas   ;
            for ( j = 0 ; j < i ; j++ )
            {
                naj     = self->configurations[j].nalphas ;
                jalphas = self->configurations[j].alphas  ;
                jbetas  = self->configurations[j].betas   ;

                /* . Skip if there are different numbers of alpha orbitals in the two configurations. */
                if ( nai != naj ) continue ;

                /* . Find the differences in the numbers of alpha and beta orbitals including positional information. */
                for ( k = na = nb = 0 ; k < nactive ; k++ )
                {
                    na += abs ( Integer1DArray_Item ( ialphas, k ) - Integer1DArray_Item ( jalphas, k ) ) ;
                    nb += abs ( Integer1DArray_Item ( ibetas,  k ) - Integer1DArray_Item ( jbetas,  k ) ) ;
                }

                /* . Skip if more than two orbitals are different. */
                if ( ( na + nb ) > 4 ) continue ;

                /* . Two orbitals different. */
                if ( ( na + nb ) == 4 )
                {
                    /* . Two beta orbitals. */
                    if      ( na == 0 ) value = CIMatrix_TwoOrbitals ( nactive, ibetas, jbetas, self->moteis ) ;
                    /* . One alpha and one beta orbital. */
                    else if ( na == 2 ) value = CIMatrix_OneAlphaOneBeta ( nactive, ialphas, ibetas, jalphas, jbetas, self->moteis ) ;
                    /* . Two alpha orbitals. */
                    else                value = CIMatrix_TwoOrbitals ( nactive, ialphas, jalphas, self->moteis ) ;
                }
                /* . One alpha orbital different. */
                else if ( na == 2 )     value = CIMatrix_OneOrbital ( nactive, ialphas, ibetas,  jalphas, self->fcoreMO, self->moteis ) ;
                /* . One beta orbital different. */
                else if ( nb == 2 )     value = CIMatrix_OneOrbital ( nactive, ibetas,  ialphas, jbetas,  self->fcoreMO, self->moteis ) ;

# ifdef DEBUGMNDOCI
{
    /* . Two orbitals different. */
    if ( ( na + nb ) == 4 )
    {
        /* . Two beta orbitals. */
        if      ( na == 0 ) test = CIMatrix_TwoOrbitals ( nactive, jbetas, ibetas, self->moteis ) ;
        /* . One alpha and one beta orbital. */
        else if ( na == 2 ) test = CIMatrix_OneAlphaOneBeta ( nactive, jalphas, jbetas, ialphas, ibetas, self->moteis ) ;
        /* . Two alpha orbitals. */
        else                test = CIMatrix_TwoOrbitals ( nactive, jalphas, ialphas, self->moteis ) ;
    }
    /* . One alpha orbital different. */
    else if ( na == 2 )     test = CIMatrix_OneOrbital ( nactive, jalphas, jbetas,  ialphas, self->fcoreMO, self->moteis ) ;
    /* . One beta orbital different. */
    else if ( nb == 2 )     test = CIMatrix_OneOrbital ( nactive, jbetas,  jalphas, ibetas,  self->fcoreMO, self->moteis ) ;

    if ( fabs ( test - value ) > 1.0e-10 ) printf ( "CI Real2DArray Element Mismatch> %d %d %d %d %25.15f %25.15f\n", i, j, na, nb, value, test ) ;
}
# endif
                /* . Save the value. */
                if ( doFull   ) SymmetricMatrix_Item ( self->ciMatrixFull, i, j ) = value ;
                if ( doSparse ) SparseSymmetricMatrix_AppendItem ( self->ciMatrixSparse, i, j, value, NULL ) ;
            }

            /* . Diagonal elements. */
            value = CIMatrix_Diagonal ( nactive, ialphas, ibetas, self->fcoreMO, self->moteis ) ;
            if ( doFull   ) SymmetricMatrix_Item ( self->ciMatrixFull, i, i ) = value ;
            if ( doSparse ) SparseSymmetricMatrix_AppendItem ( self->ciMatrixSparse, i, i, value, NULL ) ;
        }

        /* . Finalization. */
        if ( doSparse )
        {
            SparseSymmetricMatrix_Canonicalize ( self->ciMatrixSparse, NULL ) ;
            if ( self->usePreconditioning ) SparseSymmetricMatrix_MakeDiagonalPreconditioner ( self->ciMatrixSparse, self->ciMatrixPreconditioner, NULL, NULL ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the CI densities.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOCIState_MakeDensities ( MNDOCIState *self, QCOnePDM *densityp, QCOnePDM *densityq )
{
    if ( self != NULL )
    {
        auto Integer         i, j, k, na, nactive, nai, naj, nb ;
        auto Integer1DArray *ialphas, *ibetas, *jalphas, *jbetas ;
        auto Real            aiaj ;

        /* . Double loop over configurations. */
        nactive = self->nactive ;
        SymmetricMatrix_Set       ( self->onepdmMOa, 0.0e+00 ) ;
        SymmetricMatrix_Set       ( self->onepdmMOb, 0.0e+00 ) ;
        DoubleSymmetricMatrix_Set ( self->twopdm   , 0.0e+00 ) ;
        for ( i = 0 ; i < self->nconfigurations ; i++ )
        {
            nai     = self->configurations[i].nalphas ;
            ialphas = self->configurations[i].alphas  ;
            ibetas  = self->configurations[i].betas   ;
            for ( j = 0 ; j < i ; j++ )
            {
                naj     = self->configurations[j].nalphas ;
                jalphas = self->configurations[j].alphas  ;
                jbetas  = self->configurations[j].betas   ;

                /* . Skip if there are different numbers of alpha orbitals in the two configurations. */
                if ( nai != naj ) continue ;

                /* . Find the differences in the numbers of alpha and beta orbitals including positional information. */
                for ( k = na = nb = 0 ; k < nactive ; k++ )
                {
                    na += abs ( Integer1DArray_Item ( ialphas, k ) - Integer1DArray_Item ( jalphas, k ) ) ;
                    nb += abs ( Integer1DArray_Item ( ibetas,  k ) - Integer1DArray_Item ( jbetas,  k ) ) ;
                }

                /* . Skip if more than two orbitals are different. */
                if ( ( na + nb ) > 4 ) continue ;

                /* . Get the coefficient factor. */
                aiaj = 2.0e+00 * Real1DArray_Item ( self->ciVector, i ) * Real1DArray_Item ( self->ciVector, j ) ;

                /* . Two orbitals different. */
                if ( ( na + nb ) == 4 )
                {
                    /* . Two beta orbitals. */
                    if      ( na == 0 ) CIDensity_TwoOrbitals ( nactive, ibetas, jbetas, aiaj, self->twopdm ) ;
                    /* . One alpha and one beta orbital. */
                    else if ( na == 2 ) CIDensity_OneAlphaOneBeta ( nactive, ialphas, ibetas, jalphas, jbetas, aiaj, self->twopdm ) ;
                    /* . Two alpha orbitals. */
                    else                CIDensity_TwoOrbitals ( nactive, ialphas, jalphas, aiaj, self->twopdm ) ;
                }
                /* . One alpha orbital different. */
                else if ( na == 2 )     CIDensity_OneOrbital ( nactive, ialphas, ibetas,  jalphas, aiaj, self->onepdmMOa, self->twopdm ) ;
                /* . One beta orbital different. */
                else if ( nb == 2 )     CIDensity_OneOrbital ( nactive, ibetas,  ialphas, jbetas,  aiaj, self->onepdmMOb, self->twopdm ) ;
            }

            /* . Diagonal elements. */
            aiaj = Real1DArray_Item ( self->ciVector, i ) * Real1DArray_Item ( self->ciVector, i ) ;
            CIDensity_Diagonal ( nactive, ialphas, ibetas, aiaj, self->onepdmMOa, self->onepdmMOb, self->twopdm ) ;
        }

        /* . Unweight the matrices. */
        SymmetricMatrix_Unweight       ( self->onepdmMOa ) ;
        SymmetricMatrix_Unweight       ( self->onepdmMOb ) ;
        DoubleSymmetricMatrix_Unweight ( self->twopdm    ) ;

        /* . Transform to the A.O. basis. */
        SymmetricMatrix_Transform ( self->onepdmMOa, &(self->activemos), True, self->onepdma ) ;
        SymmetricMatrix_Transform ( self->onepdmMOb, &(self->activemos), True, self->onepdmb ) ;

        /* . Fill the densities for a gradient calculation. */
        if ( self->doGradients )
        {
            /* . HF density. */
            SymmetricMatrix_CopyTo ( densityp->density, self->onepdmHF ) ;
            if ( densityq != NULL ) SymmetricMatrix_AddScaledMatrix ( self->onepdmHF, 1.0e+00, densityq->density ) ;
            /* . Active densities in A.O. and M.O. bases. */
            SymmetricMatrix_CopyTo          ( self->onepdma  ,          self->onepdm    ) ;
            SymmetricMatrix_AddScaledMatrix ( self->onepdm   , 1.0e+00, self->onepdmb   ) ;
            SymmetricMatrix_CopyTo          ( self->onepdmMOa,          self->onepdmMO  ) ;
            SymmetricMatrix_AddScaledMatrix ( self->onepdmMO , 1.0e+00, self->onepdmMOb ) ;
        }

        /* . Move the CI one particle density matrices to the proper places. */
        /* . Core terms. */
        SymmetricMatrix_CopyTo ( self->pcore, densityp->density ) ;
        if ( densityq != NULL )
        {
            SymmetricMatrix_Scale  ( densityp->density, 0.5e+00 ) ;
            SymmetricMatrix_CopyTo ( densityp->density, densityq->density ) ;
        }
        /* . Active terms. */
        SymmetricMatrix_AddScaledMatrix ( densityp->density, 1.0e+00, self->onepdma ) ;
        if ( densityq == NULL ) SymmetricMatrix_AddScaledMatrix ( densityp->density, 1.0e+00, self->onepdmb ) ;
        else                    SymmetricMatrix_AddScaledMatrix ( densityq->density, 1.0e+00, self->onepdmb ) ;
/*
printf ( "\nOne PDMs:\n" ) ;
SymmetricMatrix_Print ( self->onepdma ) ;
SymmetricMatrix_Print ( self->onepdmb ) ;
*/
        /* . Save the spin densities. */
        if ( self->spinDensity != NULL )
        {
            SymmetricMatrix_CopyTo ( self->onepdma, self->spinDensity ) ;
            SymmetricMatrix_AddScaledMatrix ( self->spinDensity, -1.0e+00, self->onepdmb ) ;
        }
/*SymmetricMatrix_Print ( self->spinDensity ) ;*/
# ifdef DEBUGMNDOCI
{
    printf ( "\nMNDOCIModel_Densities:\n" ) ;
    printf ( "\nOne particle density matrix:\n" ) ;
    SymmetricMatrix_Increment ( self->onepdma, self->onepdmb ) ;
    SymmetricMatrix_Print ( self->onepdma ) ;
    printf ( "\nOne particle density matrix - MO basis:\n" ) ;
    SymmetricMatrix_Increment ( self->onepdmMOa, self->onepdmMOb ) ;
    SymmetricMatrix_Print ( self->onepdmMOa ) ;
    printf ( "\nTwo particle density matrix - MO basis:\n" ) ;
    DoubleSymmetricMatrix_Print ( self->twopdm ) ;
}
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate all configurations consistent with a given number of up and down electrons.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOCIState *MNDOCIState_MakeFull ( const Integer nactive, const Integer nup, const Integer ndown, Status *status )
{
    Integer2DArray *apermutations = NULL, *bpermutations = NULL ;
    MNDOCIState    *self = NULL ;

    /* . Create the permutations. */
    apermutations = CISetup_MakePermutations ( nup  , nactive, status ) ;
    bpermutations = CISetup_MakePermutations ( ndown, nactive, status ) ;
    if ( ( apermutations != NULL ) && ( bpermutations != NULL ) )
    {
        auto Integer a, b, n, na, nb, nconfigurations ;

        /* . Set some counters. */
        na = apermutations->length0 ;
        nb = bpermutations->length0 ;
        nconfigurations = na * nb ;

        /* . Set up the configurations. */
        self = MNDOCIState_Allocate ( nactive, nconfigurations, status ) ;
        if ( self != NULL )
        {
            auto Integer1DArray arow, brow ;
            for ( a = n = 0 ; a < na ; a++ )
            {
                Integer2DArray_RowSlice ( apermutations, a, &arow, NULL ) ;
                for ( b = 0 ; b < nb ; b++, n++ )
                {
                    Integer2DArray_RowSlice ( bpermutations, b, &brow , NULL ) ;
                    Integer1DArray_CopyTo   ( &arow, self->configurations[n].alphas, NULL ) ;
                    Integer1DArray_CopyTo   ( &brow, self->configurations[n].betas , NULL ) ;
                }
            }

            /* . Make remaining configuration data.*/
            CISetup_MakeSPQR ( self, status ) ;
        }
        else MNDOCIState_Deallocate ( &self ) ;
    }

    /* . Finish up. */
    Integer2DArray_Deallocate ( &apermutations ) ;
    Integer2DArray_Deallocate ( &bpermutations ) ;

    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate a combination of singles and doubles configurations.
!-----------------------------------------------------------------------------------------------------------------------------------
!
! . Easier to treat in terms of numbers of alpha and betas. Use A = (C+O), B = C, VA = V, VB = (O+V).
!
! . Singles - all:
!
!   A changes: A * VA + A * VB
!   B changes: B * VB + B * VA
!
! . Singles - preserve numbers of alpha and beta:
!
!   A changes: A * VA
!   B changes: B * VB
!
! . Doubles - all:
!
!   2A  -> 2A : A(A-1)/2 * VA*(VA-1)/2
!   2B  -> 2B : B(B-1)/2 * VB*(VB-1)/2
!   A,B -> A,B: A * B * VA * VB
!   2A  -> 2B : A(A-1)/2 * VB*(VB-1)/2
!   2B  -> 2A : B(B-1)/2 * VA*(VA-1)/2
!   2A  -> A,B: A(A-1)/2 * VA*VB
!   2B  -> A,B: B(B-1)/2 * VA*VB
!   A,B -> 2A : A * B * VA(VA-1)/2
!   A.B -> 2B : A * B * VB(VB-1)/2
!
! . Doubles - preserve numbers of alpha and beta:
!
!   2A  -> 2A : A(A-1)/2 * VA*(VA-1)/2
!   2B  -> 2B : B(B-1)/2 * VB*(VB-1)/2
!   A,B -> A,B: A * B * VA * VB
!
! . Extra configurations are needed for open shell cases to ensure that the spin wavefunction is correct.
!
! . "doAllStates" is not currently implemented as it would be nice to have a simpler scheme for generation of configurations
!   in the open shell case.
!
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOCIState *MNDOCIState_MakeSinglesDoubles ( const Boolean doSingles, const Boolean doDoubles, const Boolean doAllStates, const Integer nactive, const Integer nclosed, const Integer nopen, Status *status )
{
    auto Integer      na, naa, nab, nb, nbb, nconfigurations, nvirtual, va, vaa, vab, vb, vbb ;
    auto MNDOCIState *self = NULL ;

    /* . Set up some counters. */
    nvirtual = nactive - ( nclosed + nopen ) ;
    na  = nclosed + nopen ;
    nb  = nclosed ;
    va  = nactive - na ;
    vb  = nactive - nb ;
    naa = ( na * ( na - 1 ) ) / 2 ;
    nab = na * nb ;
    nbb = ( nb * ( nb - 1 ) ) / 2 ;
    vaa = ( va * ( va - 1 ) ) / 2 ;
    vab = va * vb ;
    vbb = ( vb * ( vb - 1 ) ) / 2 ;

    /* . Determine the number of configurations - including the ground state! */
    nconfigurations = 1 ;
    if ( doSingles ) nconfigurations += ( na * va + nb * vb ) ;
    if ( doDoubles ) nconfigurations += ( naa * vaa + nbb * vbb + nab * vab ) ;
/*
    if ( doAllStates )
    {
        if ( doSingles ) nconfigurations += ( na * vb + nb * va ) ;
        if ( doDoubles ) nconfigurations += ( naa * ( vab + vbb ) + nbb * ( vab + vaa ) + nab * ( vaa + vbb ) ) ;
    }
*/

    /* . Extra configurations needed for open-shell systems. */
    if ( nopen > 0 )
    {
        if ( doSingles && ( ! doDoubles ) ) nconfigurations += ( nopen * nclosed * nvirtual ) ;
        if ( doDoubles )
        {
            nconfigurations += nopen * ( 4 * nbb * vaa + nopen * ( nbb * nvirtual + vaa * nclosed ) ) + ( nopen * ( nopen - 1 ) ) / 2 * nbb * vaa ;
            if ( ! doSingles ) nconfigurations += ( 2 * nclosed * nvirtual ) ;
        }
    }

    /* . Allocate space. */
    self = MNDOCIState_Allocate ( nactive, nconfigurations, status ) ;

    /* . Set up the configurations. */
    if ( self != NULL )
    {
        auto Integer        i, j, n = 1, o, p, u, v ;
        auto Integer1DArray slice ;

        /* . Initialize all configurations to the ground state. */
        for ( i = 0 ; i < nconfigurations ; i++ )
        {
            Integer1DArray_Set ( self->configurations[i].alphas, 0 ) ; Integer1DArray_Slice ( self->configurations[i].alphas, 0, na, 1, &slice, NULL ) ; Integer1DArray_Set ( &slice, 1 ) ;
            Integer1DArray_Set ( self->configurations[i].betas , 0 ) ; Integer1DArray_Slice ( self->configurations[i].betas , 0, nb, 1, &slice, NULL ) ; Integer1DArray_Set ( &slice, 1 ) ;
        }

        /* . Singles. */
        if ( doSingles )
        {
            for ( i = 0 ; i < na ; i++ )
            {
                for ( u = na ; u < nactive ; n++, u++ ) { Integer1DArray_Item ( self->configurations[n].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; }
            }
            for ( i = 0 ; i < nb ; i++ )
            {
                for ( u = nb ; u < nactive ; n++, u++ ) { Integer1DArray_Item ( self->configurations[n].betas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas, u ) = 1 ; }
            }
/*
            if ( doAllStates )
            {
                for ( i = 0 ; i < na ; i++ )
                {
                    for ( u = nb ; u < nactive ; n++, u++ ) { Integer1DArray_Item ( self->configurations[n].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas, u ) = 1 ; }
                }
                for ( i = 0 ; i < nb ; i++ )
                {
                    for ( u = na ; u < nactive ; n++, u++ ) { Integer1DArray_Item ( self->configurations[n].betas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; }
                }
            }
*/

            /* . Add in the extra doubles that are necessary to have a consistent spin wavefunction. */
            if ( ! doDoubles && ( nopen > 0 ) )
            {
                for ( i = 0 ; i < nb ; i++ )
                {
                    for ( o = nb ; o < na ; o++ )
                    {
                        for ( u = na ; u < nactive ; n++, u++ ) { Integer1DArray_Item ( self->configurations[n].betas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].alphas, o ) = 0 ;
                                                                  Integer1DArray_Item ( self->configurations[n].betas, o ) = 1 ; Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; }
                    }
                }
            }
        }

        /* . Doubles. */
        if ( doDoubles )
        {
            for ( i = 1 ; i < na ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    for ( u = (na+1) ; u < nactive ; u++ )
                    {
                        for ( v = na ; v < u ; n++, v++ ) { Integer1DArray_Item ( self->configurations[n].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].alphas, j ) = 0 ;
                                                            Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].alphas, v ) = 1 ; }
                    }
                }
            }
            for ( i = 1 ; i < nb ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    for ( u = (nb+1) ; u < nactive ; u++ )
                    {
                        for ( v = nb ; v < u ; n++, v++ ) { Integer1DArray_Item ( self->configurations[n].betas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas, j ) = 0 ;
                                                            Integer1DArray_Item ( self->configurations[n].betas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].betas, v ) = 1 ; }
                    }
                }
            }
            for ( i = 0 ; i < na ; i++ )
            {
                for ( j = 0 ; j < nb ; j++ )
                {
                    for ( u = na ; u < nactive ; u++ )
                    {
                        for ( v = nb ; v < nactive ; n++, v++ ) { Integer1DArray_Item ( self->configurations[n].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas, j ) = 0 ;
                                                                  Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].betas, v ) = 1 ; }
                    }
                }
            }
/*
            if ( doAllStates )
            {
                for ( i = 1 ; i < na ; i++ )
                {
                    for ( j = 0 ; j < i ; j++ )
                    {
                        for ( u = na ; u < nactive ; u++ )
                        {
                            for ( v = nb ; v < nactive ; n++, v++ ) { Integer1DArray_Item ( self->configurations[n].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].alphas, j ) = 0 ;
                                                                      Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].betas , v ) = 1 ; }
                        }
                        for ( u = (nb+1) ; u < nactive ; u++ )
                        {
                            for ( v = nb ; v < u ; n++, v++ ) { Integer1DArray_Item ( self->configurations[n].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].alphas, j ) = 0 ;
                                                                Integer1DArray_Item ( self->configurations[n].betas , u ) = 1 ; Integer1DArray_Item ( self->configurations[n].betas , v ) = 1 ; }
                        }
                    }
                }
                for ( i = 1 ; i < nb ; i++ )
                {
                    for ( j = 0 ; j < i ; j++ )
                    {
                        for ( u = na ; u < nactive ; u++ )
                        {
                            for ( v = nb ; v < nactive ; n++, v++ ) { Integer1DArray_Item ( self->configurations[n].betas , i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas, j ) = 0 ;
                                                                      Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].betas, v ) = 1 ; }
                        }
                        for ( u = (na+1) ; u < nactive ; u++ )
                        {
                            for ( v = na ; v < u ; n++, v++ ) { Integer1DArray_Item ( self->configurations[n].betas , i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , j ) = 0 ;
                                                                Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].alphas, v ) = 1 ; }
                        }
                    }
                }
                for ( i = 0 ; i < na ; i++ )
                {
                    for ( j = 0 ; j < nb ; j++ )
                    {
                        for ( u = (na+1) ; u < nactive ; u++ )
                        {
                            for ( v = na ; v < u ; n++, v++ ) { Integer1DArray_Item ( self->configurations[n].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , j ) = 0 ;
                                                                Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].alphas, v ) = 1 ; }
                        }
                        for ( u = (nb+1) ; u < nactive ; u++ )
                        {
                            for ( v = nb ; v < u ; n++, v++ ) { Integer1DArray_Item ( self->configurations[n].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas, j ) = 0 ;
                                                                Integer1DArray_Item ( self->configurations[n].betas , u ) = 1 ; Integer1DArray_Item ( self->configurations[n].betas, v ) = 1 ; }
                        }
                    }
                }
            }
*/

            /* . Add in the extra configurations that are necessary to have a consistent spin wavefunction. */
            if ( nopen > 0 )
            {
                for ( i = 1 ; i < nb ; i++ )
                {
                    for ( j = 0 ; j < i ; j++ )
                    {
                        for ( u = (na+1) ; u < nactive ; u++ )
                        {
                            for ( v = na ; v < u ; v++ )
                            {
                                for ( o = nb ; o < na ; n += 4, o++ )
                                {
                                     /* . 2 closed alpha -> 2 virtual alpha. */
                                     Integer1DArray_Item ( self->configurations[n  ].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n  ].betas , j ) = 0 ;
                                     Integer1DArray_Item ( self->configurations[n  ].alphas, o ) = 0 ; Integer1DArray_Item ( self->configurations[n  ].betas , o ) = 1 ;
                                     Integer1DArray_Item ( self->configurations[n  ].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n  ].alphas, v ) = 1 ;
                                     Integer1DArray_Item ( self->configurations[n+1].alphas, j ) = 0 ; Integer1DArray_Item ( self->configurations[n+1].betas , i ) = 0 ;
                                     Integer1DArray_Item ( self->configurations[n+1].alphas, o ) = 0 ; Integer1DArray_Item ( self->configurations[n+1].betas , o ) = 1 ;
                                     Integer1DArray_Item ( self->configurations[n+1].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n+1].alphas, v ) = 1 ;
                                     /* . 2 closed beta -> 2 virtual beta. */
                                     Integer1DArray_Item ( self->configurations[n+2].betas , i ) = 0 ; Integer1DArray_Item ( self->configurations[n+2].betas , j ) = 0 ;
                                     Integer1DArray_Item ( self->configurations[n+2].alphas, o ) = 0 ; Integer1DArray_Item ( self->configurations[n+2].betas , o ) = 1 ;
                                     Integer1DArray_Item ( self->configurations[n+2].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n+2].betas , v ) = 1 ;
                                     Integer1DArray_Item ( self->configurations[n+3].betas , i ) = 0 ; Integer1DArray_Item ( self->configurations[n+3].betas , j ) = 0 ;
                                     Integer1DArray_Item ( self->configurations[n+3].alphas, o ) = 0 ; Integer1DArray_Item ( self->configurations[n+3].betas , o ) = 1 ;
                                     Integer1DArray_Item ( self->configurations[n+3].alphas, v ) = 1 ; Integer1DArray_Item ( self->configurations[n+3].betas , u ) = 1 ;
                                }
                            }
                        }
                    }
                }
                for ( o = nb ; o < na ; o++ )
                {
                    for ( p = nb ; p < na ; p++ )
                    {
                        if ( p != o )
                        {
                            /* . 1 closed alpha, 1 open alpha -> 2 virtual alpha. */
                            for ( i = 0 ; i < nb ; i++ )
                            {
                                for ( u = (na+1) ; u < nactive ; u++ )
                                {
                                    for ( v = na ; v < u ; n++, v++ )
                                    {
                                         Integer1DArray_Item ( self->configurations[n].betas , i ) = 0 ; Integer1DArray_Item ( self->configurations[n].alphas, o ) = 0 ;
                                         Integer1DArray_Item ( self->configurations[n].alphas, p ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , p ) = 1 ;
                                         Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].alphas, v ) = 1 ;
                                    }
                                }
                            }
                            /* . 2 closed beta -> 1 open beta, 1 virtual beta (or 1 open alpha, 1 closed beta -> 1 virtual alpha, 1 virtual beta). */
                            for ( i = 1 ; i < nb ; i++ )
                            {
                                for ( j = 0 ; j < i ; j++ )
                                {
                                    for ( u = na ; u < nactive ; n++, u++ )
                                    {
                                         Integer1DArray_Item ( self->configurations[n].betas , i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , j ) = 0 ;
                                         Integer1DArray_Item ( self->configurations[n].betas , o ) = 1 ; Integer1DArray_Item ( self->configurations[n].alphas, p ) = 0 ;
                                         Integer1DArray_Item ( self->configurations[n].betas , p ) = 1 ; Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ;
                                    }
                                }
                            }
                        }
                    }
                }
                for ( o = nb ; o < na ; o++ )
                {
                    /* . 1 closed alpha, 1 closed beta (same closed) -> 1 virtual alpha, 1 virtual beta (different virtual). */
                    for ( i = 0 ; i < nb ; i++ )
                    {
                        for ( u = (na+1) ; u < nactive ; u++ )
                        {
                            for ( v = na ; v < u ; n++, v++ )
                            {
                                 Integer1DArray_Item ( self->configurations[n].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , i ) = 0 ;
                                 Integer1DArray_Item ( self->configurations[n].alphas, o ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , o ) = 1 ;
                                 Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].alphas, v ) = 1 ;
                            }
                        }
                    }
                    /* . 1 closed alpha, 1 closed beta (different closed) -> 1 virtual alpha, 1 virtual beta (same virtual). */
                    for ( i = 1 ; i < nb ; i++ )
                    {
                        for ( j = 0 ; j < i ; j++ )
                        {
                            for ( u = na ; u < nactive ; n++, u++ )
                            {
                                 Integer1DArray_Item ( self->configurations[n].betas , i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , j ) = 0 ;
                                 Integer1DArray_Item ( self->configurations[n].alphas, o ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , o ) = 1 ;
                                 Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].betas , u ) = 1 ;
                            }
                        }
                    }
                }
                /* . Quadruples arising from when two open shell orbitals have beta spin. */
                for ( o = (nb+1) ; o < na ; o++ )
                {
                    for ( p = nb ; p < o ; p++ )
                    {
                        for ( i = 1 ; i < nb ; i++ )
                        {
                            for ( j = 0 ; j < i ; j++ )
                            {
                                for ( u = (na+1) ; u < nactive ; u++ )
                                {
                                    for ( v = na ; v < u ; n++, v++ )
                                    {
                                         Integer1DArray_Item ( self->configurations[n].betas , i ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , j ) = 0 ;
                                         Integer1DArray_Item ( self->configurations[n].alphas, o ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , o ) = 1 ;
                                         Integer1DArray_Item ( self->configurations[n].alphas, p ) = 0 ; Integer1DArray_Item ( self->configurations[n].betas , p ) = 1 ;
                                         Integer1DArray_Item ( self->configurations[n].alphas, u ) = 1 ; Integer1DArray_Item ( self->configurations[n].alphas, v ) = 1 ;
                                    }
                                }
                            }
                        }
                    }
                }
                if ( ! doSingles )
                {
                    for ( i = 0 ; i < nb ; i++ )
                    {
                        for ( u = na ; u < nactive ; n += 2, u++ )
                        {
                            Integer1DArray_Item ( self->configurations[n  ].alphas, i ) = 0 ; Integer1DArray_Item ( self->configurations[n  ].alphas, u ) = 1 ;
                            Integer1DArray_Item ( self->configurations[n+1].betas , i ) = 0 ; Integer1DArray_Item ( self->configurations[n+1].betas , u ) = 1 ;
                        }
                    }
                }
            }
        }

        /* . Make remaining configuration data.*/
        CISetup_MakeSPQR ( self, status ) ;
    }

    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate a state given a set of user-specified microstates.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOCIState *MNDOCIState_MakeUserSpecified ( const Integer2DArray *microstates, const Integer activeOrbitals, const Integer activeElectrons, Status *status )
{
    MNDOCIState *self = NULL ;
    if ( microstates != NULL )
    {
        auto Boolean isOK ;
        auto Integer i, nconfigurations ;
        auto Integer1DArray state ;

        /* . Basic checks. */
        isOK            = True ;
        nconfigurations = microstates->length0 ;
        for ( i = 0 ; i < nconfigurations ; i++ )
        {
            Integer2DArray_RowSlice ( microstates, i, &state, NULL ) ;
            if ( Integer1DArray_Sum ( &state ) != activeElectrons ) { isOK = False ; break ; }
        }

        /* . Basic checks. */
        if ( isOK && ( microstates->length1 == ( 2 * activeOrbitals ) ) )
        {
            /* . Set up the configurations. */
            self = MNDOCIState_Allocate ( activeOrbitals, nconfigurations, status ) ;
            if ( self != NULL )
            {
                auto Integer1DArray alphas, betas ;
                for ( i = 0 ; i < nconfigurations ; i++ )
                {
                    Integer2DArray_1DSlice ( microstates, i, i+1, 1, 0             ,   activeOrbitals, 1, &alphas, status ) ;
                    Integer2DArray_1DSlice ( microstates, i, i+1, 1, activeOrbitals, 2*activeOrbitals, 1, &betas , status ) ;
                    Integer1DArray_CopyTo  ( &alphas, self->configurations[i].alphas, status ) ;
                    Integer1DArray_CopyTo  ( &betas , self->configurations[i].betas , status ) ;
                }

                /* . Make remaining configuration data.*/
                CISetup_MakeSPQR ( self, status ) ;
            }
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return self ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/

/*==================================================================================================================================
! . CI density procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIDensity_Diagonal ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Real aiaj, SymmetricMatrix *onepdma, SymmetricMatrix *onepdmb, DoubleSymmetricMatrix *twopdm )
{
    Integer i, j ;

    /* . Loop over active alpha and beta orbitals. */
    for ( i = 0 ; i < nactive ; i++ )
    {
        if ( Integer1DArray_Item ( ialphas, i ) != 0 )
        {
            SymmetricMatrix_Item ( onepdma, i, i ) += aiaj ;

            /* . Alpha/alpha terms. */
            for ( j = 0 ; j < i ; j++ )
            {
                if ( Integer1DArray_Item ( ialphas, j ) != 0 )
                {
                    DoubleSymmetricMatrix_IncrementItem ( twopdm, i, i, j, j,  aiaj, NULL ) ;
                    DoubleSymmetricMatrix_IncrementItem ( twopdm, i, j, i, j, -aiaj, NULL ) ;
                }
            }

            /* . Alpha/beta terms. */
            for ( j = 0 ; j < nactive ; j++ )
            {
                if ( Integer1DArray_Item ( ibetas, j ) != 0 ) DoubleSymmetricMatrix_IncrementItem ( twopdm, i, i, j, j, aiaj, NULL ) ;
            }
        }
        if ( Integer1DArray_Item ( ibetas, i ) != 0 )
        {
            SymmetricMatrix_Item ( onepdmb, i, i ) += aiaj ;

            /* . Beta/beta terms. */
            for ( j = 0 ; j < i ; j++ )
            {
                if ( Integer1DArray_Item ( ibetas, j ) != 0 )
                {
                    DoubleSymmetricMatrix_IncrementItem ( twopdm, i, i, j, j,  aiaj, NULL ) ;
                    DoubleSymmetricMatrix_IncrementItem ( twopdm, i, j, i, j, -aiaj, NULL ) ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by one alpha and one beta orbital.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIDensity_OneAlphaOneBeta ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Integer1DArray *jalphas, const Integer1DArray *jbetas, const Real aiaj, DoubleSymmetricMatrix *twopdm )
{
    Integer i, j, k, l, n, p ;
    Real    factor ;

    /* . Find the alpha orbitals that differ (i and j with j > i). */
    for ( i = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { i = n ; break ; } }
    for ( j = -1, n = i+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { j = n ; break ; } }

    /* . Find the beta orbitals that differ (k and l with l > k). */
    for ( k = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ibetas , n ) != Integer1DArray_Item ( jbetas , n ) ) { k = n ; break ; } }
    for ( l = -1, n = k+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ibetas , n ) != Integer1DArray_Item ( jbetas , n ) ) { l = n ; break ; } }

    /* . Check the parity. */
    p = 0 ;
    /* . i in state 2, j in state 1. */
    if ( Integer1DArray_Item ( ialphas, i ) == 0 )
    {
        for ( n = 0 ; n <= j ; n++ ) p += Integer1DArray_Item ( ialphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Integer1DArray_Item ( jalphas, n ) ;
    }
    /* . i in state 1, j in state 2. */
    else
    {
        for ( n = 0 ; n <= j ; n++ ) p += Integer1DArray_Item ( jalphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Integer1DArray_Item ( ialphas, n ) ;
    }
    /* . k in state 2, l in state 1. */
    if ( Integer1DArray_Item ( ibetas, k ) == 0 )
    {
        for ( n = 0 ; n <= l ; n++ ) p += Integer1DArray_Item ( ibetas , n ) ;
        for ( n = 0 ; n <= k ; n++ ) p -= Integer1DArray_Item ( jbetas , n ) ;
    }
    /* . k in state 1, l in state 2. */
    else
    {
        for ( n = 0 ; n <= l ; n++ ) p += Integer1DArray_Item ( jbetas , n ) ;
        for ( n = 0 ; n <= k ; n++ ) p -= Integer1DArray_Item ( ibetas , n ) ;
    }
    if ( IsOdd ( p ) ) factor = -aiaj ;
    else               factor =  aiaj ;

    /* . Calculate P2. */
    DoubleSymmetricMatrix_IncrementItem ( twopdm, i, j, k, l, factor, NULL ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by one alpha or one beta orbital.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIDensity_OneOrbital ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Integer1DArray *jalphas, const Real aiaj, SymmetricMatrix *onepdm, DoubleSymmetricMatrix *twopdm )
{
    Integer i, j, n, p ;
    Real    factor ;

    /* . Find the orbitals that differ (j > i). */
    for ( i = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { i = n ; break ; } }
    for ( j = -1, n = i+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { j = n ; break ; } }

    /* . Check the parity. */
    p = 0 ;
    /* . i in state 2, j in state 1. */
    if ( Integer1DArray_Item ( ialphas, i ) == 0 )
    {
        for ( n = 0 ; n <= j ; n++ ) p += Integer1DArray_Item ( ialphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Integer1DArray_Item ( jalphas, n ) ;
    }
    /* . i in state 1, j in state 2. */
    else
    {
        for ( n = 0 ; n <= j ; n++ ) p += Integer1DArray_Item ( jalphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Integer1DArray_Item ( ialphas, n ) ;
    }
    if ( IsOdd ( p ) ) factor = -aiaj ;
    else               factor =  aiaj ;

    /* . Calculate P1 and P2. */
    SymmetricMatrix_IncrementItem ( onepdm, i, j, factor, NULL ) ;
    for ( n = 0 ; n < nactive ; n++ )
    {
        /* . Common alpha. */
        if ( ( Integer1DArray_Item ( ialphas, n ) != 0 ) && ( Integer1DArray_Item ( jalphas, n ) != 0 ) )
        {
            DoubleSymmetricMatrix_IncrementItem ( twopdm, i, j, n, n,  factor, NULL ) ;
            DoubleSymmetricMatrix_IncrementItem ( twopdm, i, n, j, n, -factor, NULL ) ;
        }
        /* . Common beta. */
        if ( Integer1DArray_Item ( ibetas , n ) != 0 ) DoubleSymmetricMatrix_IncrementItem ( twopdm, i, j, n, n,  factor, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by two alpha or two beta orbitals.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIDensity_TwoOrbitals ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *jalphas, const Real aiaj, DoubleSymmetricMatrix *twopdm )
{
    Integer i, j, k, l, n, p ;
    Real    factor ;

    /* . Find the orbitals that differ (i and j in state 2 with j > i and k and l in state 1 with l > k). */
    for ( i = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) < Integer1DArray_Item ( jalphas, n ) ) { i = n ; break ; } }
    for ( j = -1, n = i+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) < Integer1DArray_Item ( jalphas, n ) ) { j = n ; break ; } }
    for ( k = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) > Integer1DArray_Item ( jalphas, n ) ) { k = n ; break ; } }
    for ( l = -1, n = k+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) > Integer1DArray_Item ( jalphas, n ) ) { l = n ; break ; } }

    /* . Check the parity. */
    p = 0 ;
    for ( n = k+1 ; n <= l ; n++ ) p += Integer1DArray_Item ( ialphas, n ) ;
    for ( n = i+1 ; n <= j ; n++ ) p -= Integer1DArray_Item ( jalphas, n ) ;
    if ( IsOdd ( p ) ) factor = -aiaj ;
    else               factor =  aiaj ;

    /* . Calculate P2. */
    DoubleSymmetricMatrix_IncrementItem ( twopdm, i, k, j, l,  factor, NULL ) ;
    DoubleSymmetricMatrix_IncrementItem ( twopdm, i, l, k, j, -factor, NULL ) ;
}

/*==================================================================================================================================
! . CI four index transformation procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform index 1 by reading the hybrid integrals already with indices 2, 3 and 4 transformed.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIFIT_TransformIndex1 ( const Real2DArray *mos, const RealNDArray *tei234, DoubleSymmetricMatrix *moteis )
{
    if  ( ( mos != NULL ) && ( tei234 != NULL ) && ( moteis != NULL ) )
    {
        auto Integer i, nactive, nbasis, p, pq, q, r, rs, s, upper ;
        auto Real    sum ;

        /* . Initialization. */
        nactive = mos->length1 ;
        nbasis  = mos->length0 ;
        DoubleSymmetricMatrix_Set ( moteis, 0.0e+00 ) ;

        /* . Loop over MOs. */
        for ( p = pq = 0 ; p < nactive ; p++ )
        {
            for ( q = 0 ; q <= p ; q++, pq++ )
            {
                for ( r = rs = 0 ; r <= p ; r++ )
                {
                    if ( r == p ) upper = q ;
                    else          upper = r ;
                    for ( s = 0 ; s <= upper ; s++, rs++ )
                    {
                        /* . This is a dot-product so ultimately should use slices. */
                        for ( i = 0, sum = 0.0e+00 ; i < nbasis ; i++ ) sum += Real2DArray_Item ( mos, i, p ) * RealNDArray_Item3D ( tei234, i, q, rs ) ;
                        DoubleSymmetricMatrix_SetItem ( moteis, p, q, r, s, sum, NULL ) ;
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform index 2 by reading the hybrid integrals already with indices 3 and 4 transformed.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIFIT_TransformIndex2 ( const Real2DArray *mos, const Real2DArray *tei34, RealNDArray *tei234 )
{
    if  ( ( mos != NULL ) && ( tei34 != NULL ) && ( tei234 != NULL ) )
    {
        auto Integer i, ij, j, nactive, nbasis,q, r, rs, s ;
        auto Real    t ;

        /* . Initialization. */
        nactive = mos->length1 ;
        nbasis  = mos->length0 ;
        RealNDArray_Set ( tei234, 0.0e+00 ) ;

        /* . Loop over MO pairs. */
        for ( r = rs = 0 ; r < nactive ; r++ )
        {
            for ( s = 0 ; s <= r ; s++, rs++ )
            {
                /* . Loop over AOs. */
                for ( i = ij = 0 ; i < nbasis ; i++ )
                {
                    for ( j = 0 ; j <= i ; ij++, j++ )
                    {
                        t  = Real2DArray_Item ( tei34, ij, rs ) ;
                        if ( i == j ) t *= 0.5e+00 ;
                        for ( q = 0 ; q < nactive ; q++ )
                        {
                            RealNDArray_Item3D ( tei234, i, q, rs ) += t * Real2DArray_Item ( mos, j, q ) ;
                            RealNDArray_Item3D ( tei234, j, q, rs ) += t * Real2DArray_Item ( mos, i, q ) ;
                        }
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform indices 3 and 4 together by reading the A.O. integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIFIT_TransformIndices34 ( const Real2DArray *mos, BlockStorage *twoelectronintegrals, Real2DArray *tei34 )
{
    if  ( ( mos != NULL ) && ( twoelectronintegrals != NULL ) && ( tei34 != NULL ) )
    {
        auto Integer     i, ii, ij, j, k, kl, l, m, n, nactive, r, rs, s ;
        auto Real        t, wij, wkl ;
        auto Block *block ;

        /* . Initialization. */
        nactive = mos->length1 ;
        Real2DArray_Set ( tei34, 0.0e+00 ) ;

        /* . Loop over the integral blocks. */
        List_Iterate_Initialize ( twoelectronintegrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( twoelectronintegrals ) ) != NULL )
        {
            /* . Loop over the integrals. */
            for ( ii = 0, n = 0 ; ii < block->ndata ; ii++, n += 4 )
            {
                /* . Get the data. */
                i = block->indices16[n  ] ;
                j = block->indices16[n+1] ;
                k = block->indices16[n+2] ;
                l = block->indices16[n+3] ;
                t = block->data[ii] ;

                /* . Shuffle the index pairs. */
	        if ( i < j ) { m = i ; i = j ; j = m ; }
                if ( k < l ) { m = k ; k = l ; l = m ; }

                /* . Shuffle the indices of both pairs. */
                if ( ( i < k ) || ( ( i == k ) && ( j < l  ) ) ) { m = i ; i = k ; k = m ; m = j ; j = l ; l = m ; }

                /* . Real2DArray indices. */
                ij = SymmetricMatrix_DiagonalIndex ( i ) + j ;
                kl = SymmetricMatrix_DiagonalIndex ( k ) + l ;

                /* . Get the appropriate scaling factors. */
                wij = 1.0e+00 ;
                wkl = 1.0e+00 ;
	        if ( i == j ) wij *= 0.5e+00 ;
	        if ( k == l ) wkl *= 0.5e+00 ;
                if ( ij == kl ) { wij *= 0.5e+00 ; wkl *= 0.5e+00 ; }

                /* . Loop over MO pairs. */
                for ( r = rs = 0 ; r < nactive ; r++ )
                {
                    for ( s = 0 ; s <= r ; rs++, s++ )
                    {
                        Real2DArray_Item ( tei34, ij, rs ) += t * ( Real2DArray_Item ( mos, k, r ) * Real2DArray_Item ( mos, l, s ) + Real2DArray_Item ( mos, l, r ) * Real2DArray_Item ( mos, k, s ) ) * wkl ;
                        Real2DArray_Item ( tei34, kl, rs ) += t * ( Real2DArray_Item ( mos, i, r ) * Real2DArray_Item ( mos, j, s ) + Real2DArray_Item ( mos, j, r ) * Real2DArray_Item ( mos, i, s ) ) * wij ;
                    }
                }
            }
        }
    }
}

/*==================================================================================================================================
! . CI gradient procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation for a CI gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------------------------------
!
! . There are five groups of orbitals: inactive double, active double, fractional, active virtual, inactive virtual.
!   Their interactions make up 15 blocks.
!
! . Variables:
!
!   Non-redundant - all orbital pairs, xy, which have different occupancies (eight blocks).
!   Redundant     - all orbital pairs involving at least one active orbital that have the same occupancy (five blocks).
!
!   The two inactive blocks (double/double and virtual/virtual) do not enter.
!
!---------------------------------------------------------------------------------------------------------------------------------*/
# define OccupancyTolerance 1.0e-06

static void CIGradient_Allocate ( MNDOCIState *self, const Boolean doGradients, Status *status )
{
    if ( self != NULL )
    {
        /* . Initialization. */
        CIGradient_Initialize ( self ) ;
        self->doGradients = doGradients ;

        /* . Set up the calculation. */
        if ( doGradients )
        {
            auto Boolean isOK ;
            auto Integer i, lengths[3], numberActive, numberActiveDouble, numberActiveVirtual, numberFractional, numberCore, numberInactiveVirtual, numberNonRedundant, numberOrbitals, numberRedundant ;
            auto Real    fractionalOccupancy, o ;
            auto Status  localstatus ;

            /* . Initialization. */
            numberActive   = self->nactive   ;
            numberCore     = self->ncore     ;
            numberOrbitals = self->norbitals ;
            Status_Set ( &localstatus, Status_Continue ) ;

            /* . Determine some counters. */
            fractionalOccupancy   = 0.0e+00 ;
            numberActiveDouble    = 0 ;
            numberActiveVirtual   = 0 ;
            numberFractional      = 0 ;
            isOK = True ;
            for ( i = numberCore ; i < ( numberCore + numberActive ) ; i++ )
            {
                o = Real1DArray_Item ( self->occupancies, i ) ;
                if       ( fabs ( o - 2.0e+00 ) < OccupancyTolerance ) numberActiveDouble  += 1 ;
                else if  ( fabs ( o           ) < OccupancyTolerance ) numberActiveVirtual += 1 ;
                else
                {
                    if ( ( numberFractional > 0 ) && ( fabs ( o - fractionalOccupancy ) > OccupancyTolerance ) ) { isOK = False ; }
                    fractionalOccupancy  = o ;
                    numberFractional    += 1 ;
                }
            }
            self->ciGradientError = ( ! isOK ) ;
            numberInactiveVirtual = numberOrbitals - ( numberCore + numberActive ) ;

            /* . Number of variables. */
            numberNonRedundant =  numberCore * ( numberFractional + numberActiveVirtual + numberInactiveVirtual ) + \
                                  numberActiveDouble   * ( numberFractional + numberActiveVirtual + numberInactiveVirtual ) + \
                                  numberFractional     * (                    numberActiveVirtual + numberInactiveVirtual ) ;
            numberRedundant    =  numberCore * numberActiveDouble + \
                                ( numberActiveDouble   * ( numberActiveDouble  - 1 ) ) / 2 + \
                                ( numberFractional     * ( numberFractional    - 1 ) ) / 2 + \
                                ( numberActiveVirtual  * ( numberActiveVirtual - 1 ) ) / 2 + \
                                  numberActiveVirtual  * numberInactiveVirtual ;

            /* . Set some counters. */
            self->numberNonRedundant = numberNonRedundant ;
            self->numberRedundant    = numberRedundant    ;

            /* . Allocation. */
            lengths[0] = numberActive ; lengths[1] = numberActive ; lengths[2] = numberActive ;
            self->indicesNR      = Integer2DArray_Allocate   ( numberNonRedundant, 2, &localstatus ) ;
            self->indicesR       = Integer2DArray_Allocate   ( numberRedundant   , 2, &localstatus ) ;
            self->aDiagonal      = Real1DArray_Allocate      ( numberNonRedundant,    &localstatus ) ;
            self->qNR            = Real1DArray_Allocate      ( numberNonRedundant,    &localstatus ) ;
            self->qR             = Real1DArray_Allocate      ( numberRedundant   ,    &localstatus ) ;
            self->preconditioner = Real1DArray_Allocate      ( numberNonRedundant,    &localstatus ) ;
            self->tpdm1          = Real1DArray_Allocate      ( numberActive      ,    &localstatus ) ;
            self->tpdm2          = Real2DArray_Allocate      ( numberOrbitals    , numberActive , &localstatus ) ;
            self->tpdm3          = RealNDArray_Allocate      ( 3, lengths        ,    &localstatus ) ;
            self->onepdm         = SymmetricMatrix_AllocateN ( numberOrbitals    ,    &localstatus ) ;
            self->onepdmHF       = SymmetricMatrix_AllocateN ( numberOrbitals    ,    &localstatus ) ;
            self->onepdmMO       = SymmetricMatrix_AllocateN ( numberActive      ,    &localstatus ) ;
            self->work1          = SymmetricMatrix_AllocateN ( numberOrbitals    ,    &localstatus ) ;
            self->work2          = SymmetricMatrix_AllocateN ( numberOrbitals    ,    &localstatus ) ;
            self->zMatrix        = SymmetricMatrix_AllocateN ( numberOrbitals    ,    &localstatus ) ;

            /* . Finish up. */
            if ( ! Status_OK ( &localstatus ) ) CIGradient_Finalize ( self ) ;
            Status_Set ( status, localstatus ) ;
        }
    }
}
# undef OccupancyTolerance

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate A * B where A is the CPHF TEI matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIGradient_ApplyCPHFMatrix ( MNDOCIState *self, const Boolean addDiagonal, const Integer n1, const Integer2DArray *in1, const Integer n2, const Integer2DArray *in2, const Real1DArray *b, Real1DArray *x )
{
    if ( self != NULL )
    {
        auto Integer i, j, n ;
        auto SymmetricMatrix *work1, *work2 ;

        /* . Initialization. */
        Real1DArray_Set ( x, 0.0e+00 ) ;
        work1 = self->work1 ;
        work2 = self->work2 ;

        /* . Transform B to the A.O. basis - in work2. */
        CIGradient_CPHFTransform ( n2, in2, b, 0, in2, b, self->orbitals, True, work1, work2 ) ;

        /* . Build Y in the A.O. basis in work1. */
        MNDOIntegrals_MakeFockG ( self->twoelectronintegrals, work2, work1, NULL, NULL ) ;

        /* . Transform Y to the M.O. basis - in work2.*/
        SymmetricMatrix_Transform ( work1, self->orbitals, False, work2 ) ;

        /* . Fill X and scale. */
        for ( n = 0 ; n < n1 ; n++ )
        {
            i = Integer2DArray_Item ( in1, n, 0 ) ;
            j = Integer2DArray_Item ( in1, n, 1 ) ;
            Real1DArray_Item ( x, n ) = SymmetricMatrix_Item ( work2, j, i ) ;
        }
        Real1DArray_Scale ( x, 4.0e+00 ) ;

        /* . Add in the diagonal terms. */
        if ( addDiagonal )
        {
            for ( i = 0 ; i < n1 ; i++ ) Real1DArray_Item ( x, i ) += ( Real1DArray_Item ( self->aDiagonal, i ) * Real1DArray_Item ( b, i ) ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the vectors required for solution of the CPHF equations.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define OccupancyTolerance         1.0e-06
# define OrbitalDegeneracyTolerance 1.0e-06
# define PreconditionerTolerance    1.0e-06
# define ZeroTolerance              1.0e-12

static void CIGradient_CalculateCPHFVectors ( MNDOCIState *self )
{
    if ( self != NULL )
    {
        auto Integer          a, i, numberActive, numberCore, numberCoreAndActive, numberDegenerateRedundant, numberNonRedundant, numberOrbitals, numberRedundant, p, q, r ;
        auto Real             energyDifference, f, iEnergy, iOccupancy, occupancyDifference, pEnergy, pOccupancy ;
        auto Integer2DArray  *indicesNR, *indicesR ;
        auto Real1DArray     *qNR, *qR, *workv ;
        auto Real2DArray     *twoXY ;
        auto SymmetricMatrix *fTransformed, *gGamma ;

        /* . Initialization. */
        numberActive        = self->nactive ;
        numberCore          = self->ncore   ;
        numberCoreAndActive = numberCore + numberActive ;
        numberOrbitals      = self->norbitals ;
        fTransformed        = self->work1     ;
        gGamma              = self->work2     ;
        indicesNR           = self->indicesNR ;
        indicesR            = self->indicesR  ;
        qNR                 = self->qNR       ;
        qR                  = self->qR        ;

        /* . Calculate Ggamma - use FTRANSFORMED as scratch. */
        if ( numberCore > 0 )
        {
            MNDOIntegrals_MakeFockG   ( self->twoelectronintegrals, self->onepdm, fTransformed, NULL, NULL ) ;
            SymmetricMatrix_Transform ( fTransformed, self->orbitals, False, gGamma ) ;
            SymmetricMatrix_Scale     ( gGamma, 2.0e+00 ) ;
        }

        /* . Transform fCore to the M.O. basis. */
        SymmetricMatrix_Transform ( self->fcore, self->orbitals, False, fTransformed ) ;

        /* . Gamma terms. */
        MNDOCIState_CalculateKPA ( self ) ;
        twoXY = self->kpaMO ;

# ifdef DEBUGMNDOCIGRADIENTS
if ( numberCore > 0 )
{
printf ( "\ngGamma:\n" ) ;
SymmetricMatrix_Print ( gGamma ) ;
}
printf ( "\nfTransformed:\n" ) ;
SymmetricMatrix_Print ( fTransformed ) ;
printf ( "\nTwoXY:\n" ) ;
Real2DArray_Print ( twoXY ) ;
# endif

        /* . Fill the elements in order - non-redundant and redundant at the same time. */
        numberDegenerateRedundant = 0 ;
        numberNonRedundant        = 0 ;
        numberRedundant           = 0 ;
        Real1DArray_Set ( qNR, 0.0e+00 ) ;
        Real1DArray_Set ( qR , 0.0e+00 ) ;

        /* . Core-active. */
        for ( i = 0 ; i < numberCore ; i++ )
        {
            iEnergy    = Real1DArray_Item ( self->energies   , i ) ;
            iOccupancy = Real1DArray_Item ( self->occupancies, i ) ;
            for ( p = numberCore ; p < numberCoreAndActive ; p++ )
            {
                energyDifference    = Real1DArray_Item ( self->energies   , p ) - iEnergy    ;
                occupancyDifference = 0.5e+00 * ( iOccupancy - Real1DArray_Item ( self->occupancies, p ) ) ;
                for ( f = 0.0e+00, r = 0 ; r < numberActive ; r++ ) f += ( SymmetricMatrix_Get_Component ( self->onepdmMO, p - numberCore, r ) * SymmetricMatrix_Item ( fTransformed, r + numberCore, i ) ) ;
                f += ( Real2DArray_Item ( twoXY, i, p - numberCore ) - 2.0e+00 * SymmetricMatrix_Item ( fTransformed, p, i ) - SymmetricMatrix_Item ( gGamma, p, i ) ) ;
                if ( fabs ( occupancyDifference ) > OccupancyTolerance )
                {
                    Integer2DArray_Item ( indicesNR, numberNonRedundant, 0 ) = i ;
                    Integer2DArray_Item ( indicesNR, numberNonRedundant, 1 ) = p ;
                    Real1DArray_Item    ( self->aDiagonal, numberNonRedundant ) = energyDifference / occupancyDifference ;
                    Real1DArray_Item    ( qNR            , numberNonRedundant ) = f                / occupancyDifference ;
                    numberNonRedundant += 1 ;
                }
                else if ( fabs ( f ) > ZeroTolerance )
                {
                    if ( fabs ( energyDifference ) > OrbitalDegeneracyTolerance )
                    {
                        Integer2DArray_Item ( indicesR, numberRedundant, 0 ) = i ;
                        Integer2DArray_Item ( indicesR, numberRedundant, 1 ) = p ;
                        Real1DArray_Item    ( qR      , numberRedundant    ) = f / energyDifference ;
                        numberRedundant += 1 ;
                    }
                    else numberDegenerateRedundant += 1 ;
                }
            }
        }

        /* . Core-virtual (only non-redundant). */
        for ( i = 0 ; i < numberCore ; i++ )
        {
            iEnergy    = Real1DArray_Item ( self->energies   , i ) ;
            iOccupancy = Real1DArray_Item ( self->occupancies, i ) ;
            for ( a = numberCoreAndActive ; a < numberOrbitals ; a++ )
            {
                energyDifference    = Real1DArray_Item ( self->energies   , a ) - iEnergy    ;
                occupancyDifference = 0.5e+00 * ( iOccupancy - Real1DArray_Item ( self->occupancies, a ) ) ;
                f = - 2.0e+00 * SymmetricMatrix_Item ( fTransformed, a, i ) - SymmetricMatrix_Item ( gGamma, a, i ) ;
                Integer2DArray_Item ( indicesNR, numberNonRedundant, 0 ) = i ;
                Integer2DArray_Item ( indicesNR, numberNonRedundant, 1 ) = a ;
                Real1DArray_Item    ( self->aDiagonal, numberNonRedundant ) = energyDifference / occupancyDifference ;
                Real1DArray_Item    ( qNR            , numberNonRedundant ) = f                / occupancyDifference ;
                numberNonRedundant += 1 ;
            }
        }

        /* . Active-active. */
        for ( p = numberCore ; p < numberCoreAndActive ; p++ )
        {
            pEnergy    = Real1DArray_Item ( self->energies   , p ) ;
            pOccupancy = Real1DArray_Item ( self->occupancies, p ) ;
            for ( q = ( p + 1 ) ; q < numberCoreAndActive ; q++ )
            {
                energyDifference    = Real1DArray_Item ( self->energies   , q ) - pEnergy    ;
                occupancyDifference = 0.5e+00 * ( pOccupancy - Real1DArray_Item ( self->occupancies, q ) ) ;
                for ( f = 0.0e+00, r = 0 ; r < numberActive ; r++ ) f += ( SymmetricMatrix_Get_Component ( self->onepdmMO, q - numberCore, r ) * SymmetricMatrix_Get_Component ( fTransformed, p, r + numberCore ) -
                                                                           SymmetricMatrix_Get_Component ( self->onepdmMO, p - numberCore, r ) * SymmetricMatrix_Get_Component ( fTransformed, q, r + numberCore ) ) ;
                f += ( Real2DArray_Item ( twoXY, p, q - numberCore ) - Real2DArray_Item ( twoXY, q, p - numberCore ) ) ;
                if ( fabs ( occupancyDifference ) > OccupancyTolerance )
                {
                    Integer2DArray_Item ( indicesNR, numberNonRedundant, 0 ) = p ;
                    Integer2DArray_Item ( indicesNR, numberNonRedundant, 1 ) = q ;
                    Real1DArray_Item    ( self->aDiagonal, numberNonRedundant ) = energyDifference / occupancyDifference ;
                    Real1DArray_Item    ( qNR            , numberNonRedundant ) = f                / occupancyDifference ;
                    numberNonRedundant += 1 ;
                }
                else if ( fabs ( f ) > ZeroTolerance )
                {
                    if ( fabs ( energyDifference ) > OrbitalDegeneracyTolerance )
                    {
                        Integer2DArray_Item ( indicesR, numberRedundant, 0 ) = p ;
                        Integer2DArray_Item ( indicesR, numberRedundant, 1 ) = q ;
                        Real1DArray_Item    ( qR      , numberRedundant    ) = f / energyDifference ;
                        numberRedundant += 1 ;
                    }
                    else numberDegenerateRedundant += 1 ;
                }
            }
        }

        /* . Active-virtual. */
        for ( p = numberCore ; p < numberCoreAndActive ; p++ )
        {
            pEnergy    = Real1DArray_Item ( self->energies   , p ) ;
            pOccupancy = Real1DArray_Item ( self->occupancies, p ) ;
            for ( a = numberCoreAndActive ; a < numberOrbitals ; a++ )
            {
                energyDifference    = Real1DArray_Item ( self->energies   , a ) - pEnergy    ;
                occupancyDifference = 0.5e+00 * ( pOccupancy - Real1DArray_Item ( self->occupancies, a ) ) ;
                for ( f = 0.0e+00, r = 0 ; r < numberActive ; r++ ) f -= ( SymmetricMatrix_Get_Component ( self->onepdmMO, p - numberCore, r ) * SymmetricMatrix_Item ( fTransformed, a, r + numberCore ) ) ;
                f -= Real2DArray_Item ( twoXY, a, p - numberCore ) ;
                if ( fabs ( occupancyDifference ) > OccupancyTolerance )
                {
                    Integer2DArray_Item ( indicesNR, numberNonRedundant, 0 ) = p ;
                    Integer2DArray_Item ( indicesNR, numberNonRedundant, 1 ) = a ;
                    Real1DArray_Item    ( self->aDiagonal, numberNonRedundant ) = energyDifference / occupancyDifference ;
                    Real1DArray_Item    ( qNR      , numberNonRedundant       ) = f                / occupancyDifference ;
                    numberNonRedundant += 1 ;
                }
                else if ( fabs ( f ) > ZeroTolerance )
                {
                    if ( fabs ( energyDifference ) > OrbitalDegeneracyTolerance )
                    {
                        Integer2DArray_Item ( indicesR, numberRedundant, 0 ) = p ;
                        Integer2DArray_Item ( indicesR, numberRedundant, 1 ) = a ;
                        Real1DArray_Item    ( qR      , numberRedundant    ) = f / energyDifference ;
                        numberRedundant += 1 ;
                    }
                    else numberDegenerateRedundant += 1 ;
                }
            }
        }

# ifdef DEBUGMNDOCIGRADIENTS
printf ( "\nQNR before modification:\n" ) ;
Real1DArray_Print ( qNR ) ;
if ( numberRedundant > 0 )
{
    printf ( "\nQR (first %d only valid):\n", numberRedundant ) ;
    Real1DArray_Print ( qR ) ;
}
# endif

        /* . Remove redundant terms from qNR by using the redundant A matrix and qR. */
        if ( numberRedundant > 0 )
        {
            workv = self->preconditioner ;
            CIGradient_ApplyCPHFMatrix ( self, False, numberNonRedundant, indicesNR, numberRedundant, indicesR, qR, workv ) ;
            Real1DArray_AddScaledArray ( qNR, -1.0e+00, workv, NULL ) ;
        }

        /* . Determine the preconditioner.*/
        /* . Should there be a square root here? */
        for ( i = 0 ; i < numberNonRedundant ; i++ )
        {
            f = Real1DArray_Item ( self->aDiagonal, i ) ;
            if ( fabs ( f ) > PreconditionerTolerance ) Real1DArray_Item ( self->preconditioner, i ) = 1.0e+00 / sqrt ( fabs ( f ) ) ;
            else                                        Real1DArray_Item ( self->preconditioner, i ) = 1.0e+00 / PreconditionerTolerance ;
        }

# ifdef DEBUGMNDOCIGRADIENTS
printf ( "\nQNR after modification:\n" ) ;
Real1DArray_Print ( qNR ) ;
printf ( "\nA diagonal:\n" ) ;
Real1DArray_Print ( self->aDiagonal ) ;
printf ( "\nPreconditioner:\n" ) ;
Real1DArray_Print ( self->preconditioner ) ;
# endif

        /* . Finish up. */
        self->numberDegenerateRedundant = numberDegenerateRedundant ;
        self->numberRedundant           = numberRedundant           ;
    }
}
# undef OccupancyTolerance

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert packed vectors indexed by M.O.s to the A.O. basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . Note X1 and X2 are not symmetric although the output vector, Z, will be.
! . No index pair occurs more than once and there are no diagonal pairs.
! . The process of symmetrizing introduces a factor of 2 which should be corrected for.
*/
static void CIGradient_CPHFTransform ( const Integer n1, const Integer2DArray *in1, const Real1DArray *x1, const Integer n2, const Integer2DArray *in2, const Real1DArray *x2,
                                                                                const Real2DArray *orbitals, const Boolean doScale, SymmetricMatrix *work, SymmetricMatrix *z )
{
    auto Integer i, j, n ;

    /* . Initialization. */
    SymmetricMatrix_Set ( work, 0.0e+00 ) ;
    SymmetricMatrix_Set ( z   , 0.0e+00 ) ;

    /* . Unpack the elements (first index always less than second). */
    if ( n1 > 0 )
    {
        for ( n = 0 ; n < n1 ; n++ )
        {
            i = Integer2DArray_Item ( in1, n, 0 ) ;
            j = Integer2DArray_Item ( in1, n, 1 ) ;
            SymmetricMatrix_Item ( work, j, i ) = Real1DArray_Item ( x1, n ) ;
        }
    }
    if ( n2 > 0 )
    {
        for ( n = 0 ; n < n2 ; n++ )
        {
            i = Integer2DArray_Item ( in2, n, 0 ) ;
            j = Integer2DArray_Item ( in2, n, 1 ) ;
            SymmetricMatrix_Item ( work, j, i ) = Real1DArray_Item ( x2, n ) ;
        }
    }

    /* . Transform the Z-matrix to the A.O. basis. */
    SymmetricMatrix_Transform ( work, orbitals, True, z ) ;

    /* . Scale if necessary. */
    if ( doScale ) SymmetricMatrix_Scale ( z, 0.5e+00 ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Finalization for a CI gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIGradient_Finalize ( MNDOCIState *self )
{
    if ( self != NULL )
    {
        /* . Initialization. */
        self->ciGradientError           = False ;
        self->doGradients               = False ;
        self->cphfErrorFlag             = 0 ;
        self->cphfIterations            = 0 ;
        self->numberDegenerateRedundant = 0 ;
        self->numberNonRedundant        = 0 ;
        self->numberRedundant           = 0 ;
        Integer2DArray_Deallocate  ( &(self->indicesNR     ) ) ;
        Integer2DArray_Deallocate  ( &(self->indicesR      ) ) ;
        Real1DArray_Deallocate     ( &(self->aDiagonal     ) ) ;
        Real1DArray_Deallocate     ( &(self->qNR           ) ) ;
        Real1DArray_Deallocate     ( &(self->qR            ) ) ;
        Real1DArray_Deallocate     ( &(self->preconditioner) ) ;
        Real1DArray_Deallocate     ( &(self->tpdm1         ) ) ;
        Real2DArray_Deallocate     ( &(self->tpdm2         ) ) ;
        RealNDArray_Deallocate     ( &(self->tpdm3         ) ) ;
        SymmetricMatrix_Deallocate ( &(self->onepdm        ) ) ;
        SymmetricMatrix_Deallocate ( &(self->onepdmHF      ) ) ;
        SymmetricMatrix_Deallocate ( &(self->onepdmMO      ) ) ;
        SymmetricMatrix_Deallocate ( &(self->work1         ) ) ;
        SymmetricMatrix_Deallocate ( &(self->work2         ) ) ;
        SymmetricMatrix_Deallocate ( &(self->zMatrix       ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization for a CI gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIGradient_Initialize ( MNDOCIState *self )
{
    if ( self != NULL )
    {
        self->ciGradientError           = False ;
        self->doGradients               = False ;
        self->cphfErrorFlag             = 0 ;
        self->cphfIterations            = 0 ;
        self->numberDegenerateRedundant = 0 ;
        self->numberNonRedundant        = 0 ;
        self->numberRedundant           = 0 ;
        self->indicesNR                 = NULL ;
        self->indicesR                  = NULL ;
        self->aDiagonal                 = NULL ;
        self->qNR                       = NULL ;
        self->qR                        = NULL ;
        self->preconditioner            = NULL ;
        self->tpdm1                     = NULL ;
        self->tpdm2                     = NULL ;
        self->tpdm3                     = NULL ;
        self->onepdm                    = NULL ;
        self->onepdmHF                  = NULL ;
        self->onepdmMO                  = NULL ;
        self->work1                     = NULL ;
        self->work2                     = NULL ;
        self->zMatrix                   = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Solve the CPHF equations for zNR.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIGradient_SolveCPHFEquations ( MNDOCIState *self )
{
    Integer n ;
    CGLinearEquationSolver        solver ;
    CGLinearEquationSolverReport  report ;
    CGLinearEquationSolverState  *state  ;
    CGLinearEquationSolverTarget  target ;
    Real1DArray                  *rhs    ;
# ifdef DEBUGMNDOCICPHF
    Real1DArray                  *save   ;
# endif

    /* . Initialize target arrays. */
    n    = self->numberNonRedundant ;
    rhs  = Real1DArray_Allocate ( n, NULL ) ;
# ifdef DEBUGMNDOCICPHF
    save = Real1DArray_Allocate ( n, NULL ) ;
# endif
    Real1DArray_CopyTo ( self->qNR, rhs , NULL ) ; /* . RHS. */
# ifdef DEBUGMNDOCICPHF
    Real1DArray_CopyTo ( self->qNR, save, NULL ) ; /* . RHS for check. */
# endif
    Real1DArray_Set    ( self->qNR, 0.0e+00 )   ; /* . Initial guess at solution. */

    /* . Set up the solver (defaults). */
    solver.convergenceMode    = 1       ;
    solver.errorTolerance     = 1.0e-10 ;
    solver.maximumIterations  = 500     ;

    /* . Set up the target. */
    CGLinearEquationSolverTarget_Initialize ( &target ) ;
    target.object              = ( void * ) self ;
    target.ApplyMatrix         = ( * CIGradient_ApplyMatrix         ) ;
    target.ApplyPreconditioner = ( * CIGradient_ApplyPreconditioner ) ;
    target.rhs                 = rhs       ;
    target.solution            = self->qNR ;

    /* . Set up state. */
    state = CGLinearEquationSolverState_SetupFromTarget ( &target, NULL ) ;

    /* . Solve. */
    CGLinearEquationSolver_PCGSolver ( &solver, state, &report ) ;

    /* . Finish up. */
    CGLinearEquationSolverState_Deallocate ( &state ) ;
    Real1DArray_Deallocate                 ( &rhs   ) ;
    if ( report.isConverged ) self->cphfErrorFlag = 0 ;
    else                      self->cphfErrorFlag = 1 ;
    self->cphfIterations = report.iterations  ;

# ifdef DEBUGMNDOCICPHF
    /* . Check solution. */
    Real1DArray_Set ( rhs, 0.0e+00 ) ;
    CIGradient_ApplyCPHFMatrix ( self, True, self->numberNonRedundant, self->indicesNR, self->numberNonRedundant, self->indicesNR, self->qNR, rhs ) ;
    Real1DArray_AddScaledArray ( rhs, -1.0e+00, save, NULL ) ;
    printf ( "\nMaximum error = %20.10f\n", Real1DArray_AbsoluteMaximum ( rhs ) ) ;
    printf (   "RMS     error = %20.10f\n", Real1DArray_RootMeanSquare  ( rhs ) ) ;
    Real1DArray_Deallocate ( &save ) ;
# endif

# ifdef DEBUGMNDOCIGRADIENTS
    /* . Printing. */
    printf ( "\nSolver Report:\n" ) ;
    printf ( "Is Converged     = %10u\n"  , report.isConverged     ) ;
    printf ( "Iterations       = %10d\n"  , report.iterations      ) ;
    printf ( "Final   Residual = %10.4e\n", report.finalResidual   ) ;
    printf ( "Initial Residual = %10.4e\n", report.initialResidual ) ;
    printf ( "RHS Norm2        = %10.4e\n", report.rhsNorm2        ) ;
    printf ( "\nCPHF Solution:\n" ) ;
    Real1DArray_Print ( self->qNR ) ;
# endif

}

/* . Procedures needed by the solver. */
/* . Apply the matrix. */
static void CIGradient_ApplyMatrix ( void *object, Real1DArray *x, Real1DArray *y )
{
    if ( ( object != NULL ) && ( x != NULL ) && ( y != NULL ) )
    {
        MNDOCIState *self ;
        self = ( MNDOCIState * ) object ;
        CIGradient_ApplyCPHFMatrix ( self, True, self->numberNonRedundant, self->indicesNR, self->numberNonRedundant, self->indicesNR, x, y ) ;
# ifdef DEBUGMNDOCIGRADIENTS
        printf ( "\nApplyMatrix - input vector:\n" ) ;
        Real1DArray_Print ( x ) ;
        printf ( "\nApplyMatrix - output vector:\n" ) ;
        Real1DArray_Print ( y ) ;
# endif
    }
}

/* . Apply the preconditioner. */
static void CIGradient_ApplyPreconditioner ( void *object, Real1DArray *x, Real1DArray *y )
{
    if ( ( object != NULL ) && ( x != NULL ) && ( y != NULL ) )
    {
        MNDOCIState *self ;
        self = ( MNDOCIState * ) object ;
        /* . Diagonal preconditioner. */
        Real1DArray_CopyTo   ( x, y                   , NULL ) ;
        Real1DArray_Multiply ( y, self->preconditioner, NULL ) ;
# ifdef DEBUGMNDOCIGRADIENTS
        printf ( "\nApplyPreconditioner - input vector:\n" ) ;
        Real1DArray_Print ( x ) ;
        printf ( "\nApplyPreconditioner - output vector:\n" ) ;
        Real1DArray_Print ( y ) ;
# endif
    }
}

/*==================================================================================================================================
! . CI matrix element procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CIMatrix_Diagonal ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const SymmetricMatrix *fcore, const DoubleSymmetricMatrix *moteis  )
{
    Integer i, j ;
    Real    hij = 0.0e+00 ;

    /* . Loop over active orbitals. */
    for ( i = 0 ; i < nactive ; i++ )
    {
        /* . Alpha orbital. */
        if ( Integer1DArray_Item ( ialphas, i ) != 0 )
        {
            hij += SymmetricMatrix_Item ( fcore, i, i ) ;

            /* . Alpha/alpha terms. */
            for ( j = 0 ; j < i ; j++ )
            {
                if ( Integer1DArray_Item ( ialphas, j ) != 0 ) hij += ( DoubleSymmetricMatrix_GetItem ( moteis, i, i, j, j, NULL ) - DoubleSymmetricMatrix_GetItem ( moteis, i, j, i, j, NULL ) ) ;
            }

            /* . Alpha/beta terms. */
            for ( j = 0 ; j < nactive ; j++ )
            {
                if ( Integer1DArray_Item ( ibetas, j ) != 0 ) hij += DoubleSymmetricMatrix_GetItem ( moteis, i, i, j, j, NULL ) ;
            }
        }

        /* . Beta orbital. */
        if ( Integer1DArray_Item ( ibetas, i ) != 0 )
        {
            hij += SymmetricMatrix_Item ( fcore, i, i ) ;

            /* . Beta/beta terms. */
            for ( j = 0 ; j < i ; j++ )
            {
                if ( Integer1DArray_Item ( ibetas, j ) != 0 ) hij += ( DoubleSymmetricMatrix_GetItem ( moteis, i, i, j, j, NULL ) - DoubleSymmetricMatrix_GetItem ( moteis, i, j, i, j, NULL ) ) ;
            }
        }
    }
    return hij ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by one alpha and one beta orbital.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CIMatrix_OneAlphaOneBeta ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Integer1DArray *jalphas, const Integer1DArray *jbetas, const DoubleSymmetricMatrix *moteis )
{
    Integer i, j, k, l, n, p ;
    Real    hij = 0.0e+00 ;

    /* . Find the alpha orbitals that differ (i and j with j > i). */
    for ( i = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { i = n ; break ; } }
    for ( j = -1, n = i+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { j = n ; break ; } }

    /* . Find the beta orbitals that differ (k and l with l > k). */
    for ( k = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ibetas , n ) != Integer1DArray_Item ( jbetas , n ) ) { k = n ; break ; } }
    for ( l = -1, n = k+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ibetas , n ) != Integer1DArray_Item ( jbetas , n ) ) { l = n ; break ; } }

    /* . Calculate the matrix element. */
    hij = DoubleSymmetricMatrix_GetItem ( moteis, i, j, k, l, NULL ) ;

    /* . Check the parity. */
    p = 0 ;
    /* . i in state 2, j in state 1. */
    if ( Integer1DArray_Item ( ialphas, i ) == 0 )
    {
        for ( n = 0 ; n <= j ; n++ ) p += Integer1DArray_Item ( ialphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Integer1DArray_Item ( jalphas, n ) ;
    }
    /* . i in state 1, j in state 2. */
    else
    {
        for ( n = 0 ; n <= j ; n++ ) p += Integer1DArray_Item ( jalphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Integer1DArray_Item ( ialphas, n ) ;
    }
    /* . k in state 2, l in state 1. */
    if ( Integer1DArray_Item ( ibetas, k ) == 0 )
    {
        for ( n = 0 ; n <= l ; n++ ) p += Integer1DArray_Item ( ibetas , n ) ;
        for ( n = 0 ; n <= k ; n++ ) p -= Integer1DArray_Item ( jbetas , n ) ;
    }
    /* . k in state 1, l in state 2. */
    else
    {
        for ( n = 0 ; n <= l ; n++ ) p += Integer1DArray_Item ( jbetas , n ) ;
        for ( n = 0 ; n <= k ; n++ ) p -= Integer1DArray_Item ( ibetas , n ) ;
    }
    if ( IsOdd ( p ) ) hij *= -1.0e+00 ;

    return hij ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by one alpha or one beta orbital.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CIMatrix_OneOrbital ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *ibetas, const Integer1DArray *jalphas, const SymmetricMatrix *fcore, const DoubleSymmetricMatrix *moteis )
{
    Integer i, j, n, p ;
    Real    hij ;

    /* . Find the orbitals that differ (j > i). */
    for ( i = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { i = n ; break ; } }
    for ( j = -1, n = i+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { j = n ; break ; } }

    /* . Calculate the matrix element. */
    hij = SymmetricMatrix_Item ( fcore, j, i ) ;
    for ( n = 0 ; n < nactive ; n++ )
    {
        /* . Common alpha. */
        if ( ( Integer1DArray_Item ( ialphas, n ) != 0 ) && ( Integer1DArray_Item ( jalphas, n ) != 0 ) ) hij += ( DoubleSymmetricMatrix_GetItem ( moteis, i, j, n, n, NULL ) - DoubleSymmetricMatrix_GetItem ( moteis, i, n, j, n, NULL ) ) ;
        /* . Common beta. */
        if (   Integer1DArray_Item ( ibetas , n ) != 0                                                  ) hij +=   DoubleSymmetricMatrix_GetItem ( moteis, i, j, n, n, NULL ) ;
    }

    /* . Check the parity. */
    p = 0 ;
    /* . i in state 2, j in state 1. */
    if ( Integer1DArray_Item ( ialphas, i ) == 0 )
    {
        for ( n = 0 ; n <= j ; n++ ) p += Integer1DArray_Item ( ialphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Integer1DArray_Item ( jalphas, n ) ;
    }
    /* . i in state 1, j in state 2. */
    else
    {
        for ( n = 0 ; n <= j ; n++ ) p += Integer1DArray_Item ( jalphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Integer1DArray_Item ( ialphas, n ) ;
    }
    if ( IsOdd ( p ) ) hij *= -1.0e+00 ;

    return hij ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by two alpha or two beta orbitals.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CIMatrix_TwoOrbitals ( const Integer nactive, const Integer1DArray *ialphas, const Integer1DArray *jalphas, const DoubleSymmetricMatrix *moteis )
{
    Integer i, j, k, l, n, p ;
    Real    hij = 0.0e+00 ;

    /* . Find the orbitals that differ (i and j in state 2 with j > i and k and l in state 1 with l > k). */
    for ( i = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) < Integer1DArray_Item ( jalphas, n ) ) { i = n ; break ; } }
    for ( j = -1, n = i+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) < Integer1DArray_Item ( jalphas, n ) ) { j = n ; break ; } }
    for ( k = -1, n = 0   ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) > Integer1DArray_Item ( jalphas, n ) ) { k = n ; break ; } }
    for ( l = -1, n = k+1 ; n < nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) > Integer1DArray_Item ( jalphas, n ) ) { l = n ; break ; } }

    /* . Calculate the matrix element. */
    hij = DoubleSymmetricMatrix_GetItem ( moteis, i, k, j, l, NULL ) - DoubleSymmetricMatrix_GetItem ( moteis, i, l, k, j, NULL ) ;

    /* . Check the parity. */
    p = 0 ;
    for ( n = k+1 ; n <= l ; n++ ) p += Integer1DArray_Item ( ialphas, n ) ;
    for ( n = i+1 ; n <= j ; n++ ) p -= Integer1DArray_Item ( jalphas, n ) ;
    if ( IsOdd ( p ) ) hij *= -1.0e+00 ;

    return hij ;
}

/*==================================================================================================================================
! . CI setup procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate all possible permutations (m in n) - small numbers only!
!---------------------------------------------------------------------------------------------------------------------------------*/
# define MAXIMUMFACTORIAL 12
static Integer Factorial ( const Integer n )
{
    auto Integer f = -1, i ;
    if ( n <= MAXIMUMFACTORIAL )
    {
        f = 1 ;
        for ( i = 2 ; i <= n ; i++ ) f *= i ;
    }
    return f ;
}

static Integer2DArray *CISetup_MakePermutations ( const Integer m, const Integer n, Status *status )
{
    Integer npermutations ;
    Integer2DArray *permutations = NULL ;

    /* . Find the number of permutations. */
    npermutations = Factorial ( n ) / ( Factorial ( m ) * Factorial ( n - m ) ) ;

    /* . Allocate space. */
    permutations = Integer2DArray_Allocate ( npermutations, n, status ) ;

    /* . Generate the permutations. */
    if ( permutations != NULL )
    {
        auto Boolean isOK = True ;
        auto Integer i, j, k ;
        auto Integer1DArray *indices = NULL ;

        /* . Initialization. */
        Integer2DArray_Set ( permutations, 0 ) ;

        /* . Define the initial set of indices. */
        indices = Integer1DArray_Allocate ( m, NULL ) ;
        for ( i = 0 ; i < m ; i++ ) Integer1DArray_Item ( indices, i ) = i ;

        /* . First combination. */
        for ( i = 0 ; i < m ; i++ ) Integer2DArray_Item ( permutations, 0, i ) = 1 ;

        /* . Subsequent combinations. */
        for ( k = 1 ; k < npermutations ; k++ )
        {
            i = m - 1 ;
            while ( ( i > 0 ) && ( Integer1DArray_Item ( indices, i ) == n - m + i ) ) i -= 1 ;
            if ( ( i == 0 ) && ( Integer1DArray_Item ( indices, i ) == n - m + i ) ) { isOK = False ; break ; }
            Integer1DArray_Item ( indices, i ) += 1 ;
            for ( j = i ; j < m-1 ; j++ ) Integer1DArray_Item ( indices, j+1 ) = Integer1DArray_Item ( indices, j ) + 1 ;
            for ( i = 0 ; i < m ; i++ ) Integer2DArray_Item ( permutations, k, Integer1DArray_Item ( indices, i ) ) = 1 ;
        }
        if ( isOK ) Status_Set ( status, Status_Continue ) ;
        else
        {
            Integer2DArray_Deallocate ( &permutations ) ;
            Status_Set ( status, Status_LogicError ) ;
        }

        /* . Finish up. */
        Integer1DArray_Deallocate ( &indices ) ;
    }
    return permutations ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the SPQR data which is necessary for the calculation of state spins.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . 1 is added to the spqr values to avoid the problem of having a configuration of index
!  . "0" with an odd parity. It would be better to have a separate parity array. */

static void CISetup_MakeSPQR ( MNDOCIState *self, Status *status )
{
    if ( ( self != NULL ) && ( Status_OK ( status ) ) )
    {
        auto Integer i, na, nab, nb ;

        /* . Make basic spin data. */
        for ( i = 0 ; i < self->nconfigurations ; i++ )
        {
            na  = Integer1DArray_Sum ( self->configurations[i].alphas ) ;
            nb  = Integer1DArray_Sum ( self->configurations[i].betas  ) ;
            nab = Integer1DArray_Dot ( self->configurations[i].alphas, self->configurations[i].betas, NULL ) ;
            self->configurations[i].nalphas = na ;
            self->configurations[i].spin    = 4.0e+00 * ( ( Real ) nab ) - pow ( ( Real ) ( na - nb ), 2 ) ;
        }

        /* . Make SPQR. */
        if ( self->nconfigurations > 1 )
        {
            auto Integer         j, k, n, nai, naj, nspqr = 0, p, q, r, s, t ;
            auto Integer1DArray *ialphas, *ibetas, *jalphas, *jbetas, slice, *spqr ;

            /* . Temporary space. */
            spqr = Integer1DArray_Allocate ( self->nconfigurations - 1, status ) ;

            /* . Double loop over configurations. */
            if ( spqr != NULL )
            {
                for ( i = 1 ; i < self->nconfigurations ; i++ )
                {
                    nai     = self->configurations[i].nalphas ;
                    ialphas = self->configurations[i].alphas  ;
                    ibetas  = self->configurations[i].betas   ;
                    for ( j = nspqr = 0 ; j < i ; j++ )
                    {
                        naj     = self->configurations[j].nalphas ;
                        jalphas = self->configurations[j].alphas  ;
                        jbetas  = self->configurations[j].betas   ;

                        /* . Skip if there are different numbers of alpha orbitals in the two configurations. */
                        if ( nai != naj ) continue ;

                        /* . Find the differences in the numbers of alpha and beta orbitals including positional information. */
                        for ( k = na = nb = 0 ; k < self->nactive ; k++ )
                        {
                            na += abs ( Integer1DArray_Item ( ialphas, k ) - Integer1DArray_Item ( jalphas, k ) ) ;
                            nb += abs ( Integer1DArray_Item ( ibetas,  k ) - Integer1DArray_Item ( jbetas,  k ) ) ;
                        }

                        /* . States differing by one alpha and one beta orbital. */
                        if ( ( na == 2 ) && ( nb == 2 ) )
                        {
                            /* . Find the alpha orbitals that differ (i and j with j > i). */
                            for ( p = -1, n = 0   ; n < self->nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { p = n ; break ; } }
                            for ( q = -1, n = p+1 ; n < self->nactive ; n++ ) { if ( Integer1DArray_Item ( ialphas, n ) != Integer1DArray_Item ( jalphas, n ) ) { q = n ; break ; } }

                            /* . Find the beta orbitals that differ (k and l with l > k). */
                            for ( r = -1, n = 0   ; n < self->nactive ; n++ ) { if ( Integer1DArray_Item ( ibetas , n ) != Integer1DArray_Item ( jbetas , n ) ) { r = n ; break ; } }
                            for ( s = -1, n = r+1 ; n < self->nactive ; n++ ) { if ( Integer1DArray_Item ( ibetas , n ) != Integer1DArray_Item ( jbetas , n ) ) { s = n ; break ; } }

                            /* . Save configurations of type (paqb and pbqa) and calculate the parity. */
                            if ( ( p == r ) && ( q == s ) && ( Integer1DArray_Item ( ialphas, p ) != Integer1DArray_Item ( ibetas, p ) ) )
                            {
                                t = -1 ;
                                for ( n = p+1 ; n < self->nactive ; n++ ) t += Integer1DArray_Item ( ialphas, n ) ;
                                for ( n = q+1 ; n < self->nactive ; n++ ) t += Integer1DArray_Item ( jalphas, n ) ;
                                for ( n = 0   ; n <= q            ; n++ ) t += Integer1DArray_Item ( ibetas , n ) ;
                                for ( n = 0   ; n <= p            ; n++ ) t += Integer1DArray_Item ( jbetas , n ) ;
                                if ( IsOdd ( t ) ) Integer1DArray_Item ( spqr, nspqr ) = -(j+1) ;
                                else               Integer1DArray_Item ( spqr, nspqr ) =   j+1 ;
                                nspqr++ ;
                            }
                        }
                    }
                    self->configurations[i].nspqr = nspqr ;
                    if ( nspqr > 0 )
                    {
                        Integer1DArray_Slice ( spqr, 0, nspqr, 1, &slice, NULL ) ;
                        self->configurations[i].spqr = Integer1DArray_Clone ( &slice, status ) ;
                        if ( ! Status_OK ( status ) ) break ;
                    }
                }
            }

            /* . Finish up.*/
            Integer1DArray_Deallocate ( &spqr ) ;
        }
    }
}

/*==================================================================================================================================
! . CI temporary data handling.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CITemporary_Allocate ( MNDOCIState *self, Status *status )
{
    if ( self != NULL )
    {
        auto Integer lengths[3], natr, notr ;
        natr = ( self->nactive   * ( self->nactive   + 1 ) ) / 2 ;
        notr = ( self->norbitals * ( self->norbitals + 1 ) ) / 2 ;
        lengths[0] = self->norbitals ; lengths[1] = self->nactive ; lengths[2] = natr ;
        self->moteis       = DoubleSymmetricMatrix_Allocate ( self->nactive         , status ) ;
        self->twopdm       = DoubleSymmetricMatrix_Allocate ( self->nactive         , status ) ;
        self->ocore        = Real1DArray_Allocate           ( self->norbitals       , status ) ;
        self->kpa          = Real2DArray_Allocate           ( self->norbitals       , self->nactive , status ) ;
        self->kpaMO        = Real2DArray_Allocate           ( self->norbitals       , self->nactive , status ) ;
        self->motei34      = Real2DArray_Allocate           ( notr                  , natr          , status ) ;
        self->motei234     = RealNDArray_Allocate           ( 3                     , lengths       , status ) ;
        self->fcore        = SymmetricMatrix_AllocateN      ( self->norbitals       , status ) ;
        self->fcoreMO      = SymmetricMatrix_AllocateN      ( self->nactive         , status ) ;
        self->onepdma      = SymmetricMatrix_AllocateN      ( self->norbitals       , status ) ;
        self->onepdmb      = SymmetricMatrix_AllocateN      ( self->norbitals       , status ) ;
        self->onepdmMOa    = SymmetricMatrix_AllocateN      ( self->nactive         , status ) ;
        self->onepdmMOb    = SymmetricMatrix_AllocateN      ( self->nactive         , status ) ;
        self->pcore        = SymmetricMatrix_AllocateN      ( self->norbitals       , status ) ;
        /* . Wavefunction data that may already exist. */
        if ( self->ciEnergies  == NULL ) self->ciEnergies  = Real1DArray_Allocate      ( self->numberOfStates , status ) ;
        if ( self->ciVector    == NULL ) self->ciVector    = Real1DArray_Allocate      ( self->nconfigurations, status ) ;
        if ( self->spins       == NULL ) self->spins       = Real1DArray_Allocate      ( self->nconfigurations, status ) ;
        if ( self->ciVectors   == NULL ) self->ciVectors   = Real2DArray_Allocate      ( self->numberOfStates, self->nconfigurations, status ) ;
        if ( self->spinDensity == NULL ) self->spinDensity = SymmetricMatrix_AllocateN ( self->norbitals      , status ) ;
        /* . Solver options. */
        /* . Full solver. */
        if ( self->doFull )
        {
            self->ciMatrixFull = SymmetricMatrix_AllocateN ( self->nconfigurations, status ) ;
            if ( self->doSparse )
            {
                self->checkEnergies = Real1DArray_Allocate (                       self->nconfigurations, status ) ;
                self->checkVectors  = Real2DArray_Allocate ( self->numberOfStates, self->nconfigurations, status ) ;
            }
        }
        /* . Sparse solver. */
        if ( self->doSparse )
        {
            if ( self->usePreconditioning ) self->ciMatrixPreconditioner = Real1DArray_Allocate ( self->nconfigurations, status ) ;
            self->ciMatrixSparse = SparseSymmetricMatrix_Allocate ( self->nconfigurations, self->ciMatrixNonZero, status ) ;

            /* . Solver target. */
            JDEigenvalueSolverTarget_Initialize ( &(self->eigenvalueSolverTarget) ) ;
            self->eigenvalueSolverTarget.eigenvalues         = self->ciEnergies      ;
            self->eigenvalueSolverTarget.eigenvectors        = self->ciVectors       ;
            self->eigenvalueSolverTarget.object              = self                  ;
            self->eigenvalueSolverTarget.ApplyMatrix         = JDApplyMatrixSparse   ;
            if ( self->usePreconditioning ) self->eigenvalueSolverTarget.ApplyPreconditioner = JDApplyPreconditioner ;

            /* . Solver state. */
            self->eigenvalueSolverState = JDEigenvalueSolverState_SetupFromTarget ( &(self->eigenvalueSolverTarget), status ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CITemporary_Deallocate ( MNDOCIState *self, const Boolean keepWavefunction )
{
    if ( self != NULL )
    {
        self->twoelectronintegrals = NULL ;
        self->oneelectronmatrix    = NULL ;
        DoubleSymmetricMatrix_Deallocate ( &(self->moteis                ) ) ;
        DoubleSymmetricMatrix_Deallocate ( &(self->twopdm                ) ) ;
        Real1DArray_Deallocate           ( &(self->checkEnergies         ) ) ;
        Real1DArray_Deallocate           ( &(self->ciMatrixPreconditioner) ) ;
        Real1DArray_Deallocate           ( &(self->ocore                 ) ) ;
        Real2DArray_Deallocate           ( &(self->checkVectors          ) ) ;
        Real2DArray_Deallocate           ( &(self->kpa                   ) ) ;
        Real2DArray_Deallocate           ( &(self->kpaMO                 ) ) ;
        Real2DArray_Deallocate           ( &(self->motei34               ) ) ;
        RealNDArray_Deallocate           ( &(self->motei234              ) ) ;
        SparseSymmetricMatrix_Deallocate ( &(self->ciMatrixSparse        ) ) ;
        SymmetricMatrix_Deallocate       ( &(self->ciMatrixFull          ) ) ;
        SymmetricMatrix_Deallocate       ( &(self->fcore                 ) ) ;
        SymmetricMatrix_Deallocate       ( &(self->fcoreMO               ) ) ;
        SymmetricMatrix_Deallocate       ( &(self->onepdma               ) ) ;
        SymmetricMatrix_Deallocate       ( &(self->onepdmb               ) ) ;
        SymmetricMatrix_Deallocate       ( &(self->onepdmMOa             ) ) ;
        SymmetricMatrix_Deallocate       ( &(self->onepdmMOb             ) ) ;
        SymmetricMatrix_Deallocate       ( &(self->pcore                 ) ) ;
        /* . Solver. */
        JDEigenvalueSolverState_Deallocate ( &(self->eigenvalueSolverState) ) ;
        /* . Wavefunction data that may be kept. */
        if ( ! keepWavefunction )
        {
            Real1DArray_Deallocate     ( &(self->ciEnergies ) ) ;
            Real1DArray_Deallocate     ( &(self->ciVector   ) ) ;
            Real1DArray_Deallocate     ( &(self->spins      ) ) ;
            Real2DArray_Deallocate     ( &(self->ciVectors  ) ) ;
            SymmetricMatrix_Deallocate ( &(self->spinDensity) ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CITemporary_Initialize ( MNDOCIState *self )
{
    if ( self != NULL )
    {
        self->doFull                       = False ;
        self->doSparse                     = False ;
        self->fractionallyOccupiedInactive = False ;
        self->orbitalDegeneracies          = False ;
        self->rootNotFound                 = False ;
        self->usePreconditioning           = False ;
        self->ciMatrixNonZero              = -1 ;
        self->localizeStart                = -1 ;
        self->localizeStop                 = -1 ;
        self->ncore                        =  0 ;
        self->nelectrons                   =  0 ;
        self->norbitals                    =  0 ;
        self->numberOfStates               = -1 ;
        self->root                         =  0 ;
        self->baseline                     = 0.0e+00 ;
        self->ciEnergy                     = 0.0e+00 ;
        self->ciMatrixSparsity             = 0.0e+00 ;
        self->ecore                        = 0.0e+00 ;
        self->requiredSpin                 = 0.0e+00 ;
        self->rootenergy                   = 0.0e+00 ;
        self->twoelectronintegrals         = NULL ;
        self->oneelectronmatrix            = NULL ;
        self->moteis                       = NULL ;
        self->twopdm                       = NULL ;
        self->checkEnergies                = NULL ;
        self->ciEnergies                   = NULL ;
        self->ciMatrixPreconditioner       = NULL ;
        self->ciVector                     = NULL ;
        self->energies                     = NULL ;
        self->occupancies                  = NULL ;
        self->ocore                        = NULL ;
        self->spins                        = NULL ;
        self->checkVectors                 = NULL ;
        self->ciVectors                    = NULL ;
        self->kpa                          = NULL ;
        self->kpaMO                        = NULL ;
        self->motei34                      = NULL ;
        self->orbitals                     = NULL ;
        self->motei234                     = NULL ;
        self->ciMatrixSparse               = NULL ;
        self->ciMatrixFull                 = NULL ;
        self->fcore                        = NULL ;
        self->fcoreMO                      = NULL ;
        self->onepdma                      = NULL ;
        self->onepdmb                      = NULL ;
        self->onepdmMOa                    = NULL ;
        self->onepdmMOb                    = NULL ;
        self->pcore                        = NULL ;
        self->spinDensity                  = NULL ;
        /* . Solver. */
        self->eigenvalueSolverState        = NULL ;
        JDEigenvalueSolverReport_Initialize ( &(self->eigenvalueSolverReport) ) ;
        JDEigenvalueSolverTarget_Initialize ( &(self->eigenvalueSolverTarget) ) ;
    }
}

/*==================================================================================================================================
! . Davidson procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply the CI matrix to a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void JDApplyMatrixSparse ( void *xvoid, void *yvoid, Integer *blockSize, primme_params *primme )
{
    Real                    *x, *y   ;
    Real1DArray              hV, v   ;
    JDEigenvalueSolverState *jdState ;
    MNDOCIState             *self    ;

    /* . Casts. */
    jdState = ( JDEigenvalueSolverState * ) primme->matrix ;
    self    = ( MNDOCIState * ) jdState->target->object ;
    x       = ( Real * ) xvoid ;
    y       = ( Real * ) yvoid ;

    /* . Views. */
    Real1DArray_ViewOfRaw (  &v, 0, primme->n, 1, x, primme->n, NULL ) ;
    Real1DArray_ViewOfRaw ( &hV, 0, primme->n, 1, y, primme->n, NULL ) ;

    /* . Calculate H * v. */
    SparseSymmetricMatrix_VectorMultiply ( self->ciMatrixSparse, &v, &hV, NULL ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply the CI matrix preconditioner to a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void JDApplyPreconditioner ( void *xvoid, void *yvoid, Integer *blockSize, primme_params *primme )
{
    Real                    *x, *y ;
    Real1DArray              hV, v ;
    JDEigenvalueSolverState *jdState ;
    MNDOCIState             *self    ;

    /* . Casts. */
    jdState = ( JDEigenvalueSolverState * ) primme->matrix ;
    self    = ( MNDOCIState * ) jdState->target->object ;
    x       = ( Real * ) xvoid ;
    y       = ( Real * ) yvoid ;

    /* . Views. */
    Real1DArray_ViewOfRaw (  &v, 0, primme->n, 1, x, primme->n, NULL ) ;
    Real1DArray_ViewOfRaw ( &hV, 0, primme->n, 1, y, primme->n, NULL ) ;

    /* . Calculate v / Diagonal ( h ). */
    Real1DArray_CopyTo   ( &v , &hV                         , NULL ) ;
    Real1DArray_Multiply ( &hV, self->ciMatrixPreconditioner, NULL ) ;
}

/*==================================================================================================================================
! . MNDO CI configuration procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate alphas and betas.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void MNDOCIConfiguration_AllocateAlphasBetas ( MNDOCIConfiguration *self, const Integer nactive, Status *status )
{
    if ( self != NULL )
    {
        self->alphas  = Integer1DArray_Allocate ( nactive, status ) ;
        self->betas   = Integer1DArray_Allocate ( nactive, status ) ;
        Integer1DArray_Set ( self->alphas, 0 ) ;
        Integer1DArray_Set ( self->betas , 0 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void MNDOCIConfiguration_Deallocate ( MNDOCIConfiguration *self )
{
    if ( self != NULL )
    {
        Integer1DArray_Deallocate ( &(self->alphas) ) ;
        Integer1DArray_Deallocate ( &(self->betas ) ) ;
        Integer1DArray_Deallocate ( &(self->spqr  ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void MNDOCIConfiguration_Initialize ( MNDOCIConfiguration *self )
{
    if ( self != NULL )
    {
        self->nalphas = 0 ;
        self->nspqr   = 0 ;
        self->spin    = 0.0e+00 ;
        self->alphas  = NULL ;
        self->betas   = NULL ;
        self->spqr    = NULL ;
    }
}
