/*------------------------------------------------------------------------------
! . File      : QCModelDFT.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Procedures for calculating the DFT quantum chemical energy.
!=================================================================================================================================*/

# include "DFTIntegrator.h"
# include "GaussianBasisDerivatives.h"
# include "GaussianBasisIntegrals.h"
# include "Memory.h"
# include "Real.h"
# include "QCModelDFT.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Defaults for some variables.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DefaultInCore                    False
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
# define DefaultAccuracy                  DFTGridAccuracy_Medium
# define DefaultQCChargeModel             QCChargeModel_Lowdin
# define DefaultOccupancyType             QCOnePDMOccupancyType_FractionalVariable

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void QCModelDFT_MakeLowdinTransformation ( const QCModelDFT *self, QCModelDFTState *qcState ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCModelDFT *QCModelDFT_Allocate ( void )
{
    QCModelDFT *self = NULL ;
    self = ( QCModelDFT * ) Memory_Allocate ( sizeof ( QCModelDFT ) ) ;
    if ( self != NULL )
    {
        self->inCore                    = DefaultInCore                    ;
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
        self->accuracy                  = DefaultAccuracy                  ;
        self->functionalModel           = NULL                             ;
        self->qcChargeModel             = DefaultQCChargeModel             ;
        self->occupancyType             = DefaultOccupancyType             ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Atomic charge calculation.
! . The input charge model, if specified, overrides the one in self.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFT_AtomicCharges ( const QCModelDFT *self, QCModelDFTState *qcState, const QCChargeModel *qcChargeModel, const Boolean QSPIN, Real1DArray *qcCharges )
{
    /* . Initialization. */
    Real1DArray_Set ( qcCharges, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( qcCharges != NULL ) && ( qcState->qcAtoms != NULL ) && ( qcState->qcAtoms->natoms <= qcCharges->length ) && ! ( QSPIN && ( qcState->densityq == NULL ) ) )
    {
        auto QCChargeModel  model   ;
        auto QCOnePDM      *density ;
        /* . Get the model. */
        if ( qcChargeModel == NULL ) model = self->qcChargeModel ;
        else                         model = (*qcChargeModel)    ;
        /* . Initialize the charges to the nuclear charges. */
        if ( ! QSPIN ) QCAtomContainer_GetNuclearCharges ( qcState->qcAtoms, qcState->qcParameters, qcCharges ) ;
        /* . Get the denisty to use (Lowdin and Mulliken only for the moment). */
        if ( qcState->densityq != NULL ) QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;
        /* . Get the density pointer. */
        if ( QSPIN ) density = qcState->densityq ;
        else         density = qcState->densityp ;
        /* . Switch on the type of analysis to do. */
        switch ( model )
        {
             /* . Coulomb fitting. */
             /* . Temporarily we will assume that the coefficients exist. This means that QSPIN analyses cannot be done. */
            case QCChargeModel_CoulombFitting:
            {
                auto Real q ;
                auto Integer    iatom, ibf, istart, istop ;
                for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
                {
                    istart = qcState->qcAtoms->data[iatom].fstartw ;
                    istop  = istart + qcState->qcAtoms->data[iatom].nfbasisw ;
                    for ( ibf = istart, q = 0.0e+00 ; ibf < istop ; ibf++ ) q += Real1DArray_Item ( qcState->fpotential, ibf ) * Real1DArray_Item ( qcState->fselfoverlap, ibf ) ;
                    Real1DArray_Item ( qcCharges, iatom ) -= q ;
                }
            }
            break ;
            /* . Lowdin. */
            case QCChargeModel_Lowdin:
            {
                auto Real f, g ;
                auto Integer    iatom, i, istart, istop, j, k, n ;

                /* . Make the Lowdin transformation.*/
                QCModelDFT_MakeLowdinTransformation ( self, qcState ) ;

                /* . Electronic contribution - charge indices are proper indices, others are work indices. */
                n = qcState->lowdintransformation->length0 ;
                for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
                {
                    istart = qcState->qcAtoms->data[iatom].ostart ;
                    istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
                    for ( i = istart, f = 0.0e+00 ; i < istop ; i++ )
                    {
                        for ( j = 0 ; j < n ; j++ )
                        {
                            for ( k = 0, g = 0.0e+00 ; k < n ; k++ )
                            {
                                g += SymmetricMatrix_Get_Component ( density->density, j, k ) * Real2DArray_Item ( qcState->lowdintransformation, k, i ) ;
                            }
                            f += Real2DArray_Item ( qcState->lowdintransformation, j, i ) * g ;
                        }
                    }
                    Real1DArray_Item ( qcCharges, iatom ) -= f ;
                }
            }
            break ;
            /* . Mulliken. */
            case QCChargeModel_Mulliken:
            {
                auto Real ps ;
                auto Integer    iatom, ibf, istart, istop, jbf ;
                /* . Electronic contribution. */
                for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
                {
                    istart = qcState->qcAtoms->data[iatom].ostartw ;
                    istop  = istart + qcState->qcAtoms->data[iatom].nobasisw ;
                    for ( ibf = istart, ps = 0.0e+00 ; ibf < istop ; ibf++ )
                    {
                        for ( jbf = 0 ; jbf < density->density->dimension ; jbf++ )
                        {
                            ps += SymmetricMatrix_Get_Component ( density->density, ibf, jbf ) * SymmetricMatrix_Get_Component ( qcState->overlap, ibf, jbf ) ;
                        }
                    }
                    Real1DArray_Item ( qcCharges, iatom ) -= ps ;
                }
            }
            break ;
        }
        /* . Reset the densities. */
        if ( qcState->densityq != NULL ) QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCModelDFT *QCModelDFT_Clone ( const QCModelDFT *self )
{
    QCModelDFT *new = NULL ;
    if ( self != NULL )
    {
        new = QCModelDFT_Allocate ( ) ;
        new->inCore                    = self->inCore                    ;
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
        new->accuracy                  = self->accuracy                  ;
        new->qcChargeModel             = self->qcChargeModel             ;
        new->occupancyType             = self->occupancyType             ;
        new->functionalModel           = DFTFunctionalModel_Clone ( self->functionalModel, NULL ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFT_Deallocate ( QCModelDFT **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole moment.
!---------------------------------------------------------------------------------------------------------------------------------*/
Vector3 *QCModelDFT_DipoleMoment ( const QCModelDFT *self, QCModelDFTState *qcState, const Coordinates3 *coordinates3, const Vector3 *center )
{
    Vector3 *dipole = NULL;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( coordinates3 != NULL ) )
    {
        auto Real dx, dy, dz, q, rc[3], xq, yq, zq ;
        auto Integer    iqm ;
        auto Coordinates3    *qccoordinates3 = NULL ;
        auto SymmetricMatrix *dipx = NULL, *dipy = NULL, *dipz = NULL ;
        /* . Get the QC coordinates. */
        QCAtomContainer_GetCoordinates3 ( qcState->qcAtoms, coordinates3, True, &qccoordinates3 ) ;
        /* . Get the center. */
        if ( center == NULL ) { rc[0] = rc[1] = rc[2] = 0.0e+00 ; }
        else
        {
            rc[0] = Vector3_Item ( center, 0 ) ;
            rc[1] = Vector3_Item ( center, 1 ) ;
            rc[2] = Vector3_Item ( center, 2 ) ;
        }
        /* . Nuclear contribution. */
        for ( iqm = 0, dx = dy = dz = 0.0e+00 ; iqm < qcState->qcAtoms->natoms ; iqm++ )
        {
            q  = ( Real ) qcState->qcParameters->centers[qcState->qcAtoms->data[iqm].center].atomicNumber ;
            Coordinates3_GetRow ( qccoordinates3, iqm, xq, yq, zq ) ;
            dx += q * ( xq - rc[0] ) ;
            dy += q * ( yq - rc[1] ) ;
            dz += q * ( zq - rc[2] ) ;
        }

        /* . Electronic contribution. */
        GaussianBasisIntegrals_Dipole ( qcState->qcAtoms, qcState->qcParameters, qccoordinates3, rc, &dipx, &dipy, &dipz ) ;
        QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;
        dx -= SymmetricMatrix_Multiply2_Trace ( qcState->densityp->density, dipx ) ;
        dy -= SymmetricMatrix_Multiply2_Trace ( qcState->densityp->density, dipy ) ;
        dz -= SymmetricMatrix_Multiply2_Trace ( qcState->densityp->density, dipz ) ;
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Deallocate space. */
        Coordinates3_Deallocate ( &qccoordinates3 ) ;
        SymmetricMatrix_Deallocate ( &dipx ) ;
        SymmetricMatrix_Deallocate ( &dipy ) ;
        SymmetricMatrix_Deallocate ( &dipz ) ;

        /* . Set up the dipole vector. */
        dx *= UNITS_DIPOLE_ATOMIC_UNITS_TO_DEBYES ;
        dy *= UNITS_DIPOLE_ATOMIC_UNITS_TO_DEBYES ;
        dz *= UNITS_DIPOLE_ATOMIC_UNITS_TO_DEBYES ;
        dipole = Vector3_Allocate ( )   ;
        Vector3_Item ( dipole, 0 ) = dx ;
        Vector3_Item ( dipole, 1 ) = dy ;
        Vector3_Item ( dipole, 2 ) = dz ;
    }
    return dipole ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the Fock matrices and the electronic energy.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFT_Fock ( const QCModelDFT *self, QCModelDFTState *qcState, Real *eelectronic )
{
    /* . Initialization.*/
    if ( eelectronic != NULL ) (*eelectronic) = 0.0e+00 ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( qcState->qccoordinates3 != NULL ) && ( qcState->densityp != NULL )  && ( qcState->dftgrid != NULL ) && ( qcState->fitintegrals != NULL )
                                                                                                                  && ( qcState->fpotential != NULL ) && ( qcState->oneelectronmatrix != NULL ) )
    {
# ifndef USEOPENMP
        auto SymmetricMatrix *fockq = NULL ;
# endif
        /* . Occupancy energy. */
        qcState->eocc = qcState->densityp->occupancyEnergy ;
        if ( qcState->densityq != NULL ) qcState->eocc += qcState->densityq->occupancyEnergy ;

        /* . Put the total density in densityp. */
        QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . One-electron energy. */
        qcState->eoei = SymmetricMatrix_Multiply2_Trace ( qcState->densityp->density, qcState->oneelectronmatrix ) ;

        /* . Two-electron energy and contributions to the Fock matrix. */
        qcState->etei = GaussianBasis_Fock ( qcState->qcAtoms->nfbasisw, qcState->densityp, qcState->fitintegrals, qcState->fpotential, qcState->inversefitmatrix ) ;

        /* . Convert back densityp. */
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Add in the one-electron terms to the Fock matrices. */
        SymmetricMatrix_Increment ( qcState->densityp->fock, qcState->oneelectronmatrix ) ;
        if ( qcState->densityq != NULL ) SymmetricMatrix_Transfer ( qcState->densityq->fock, qcState->densityp->fock ) ;

        /* . DFT contribution. */
# ifdef USEOPENMP
        QCModelDFTState_InitializeFockArray ( qcState, qcState->densityp, qcState->fockAArray ) ;
        QCModelDFTState_InitializeFockArray ( qcState, qcState->densityq, qcState->fockBArray ) ;
# else
        if ( qcState->densityq != NULL ) fockq = qcState->densityq->fock ;
# endif
        DFTIntegrator_Integrate ( self->functionalModel         ,
                                  qcState->dftgrid              ,
                                  qcState->qcAtoms              ,
                                  qcState->qcParameters         ,
                                  qcState->qccoordinates3       ,
                                  qcState->densityp             ,
                                  qcState->densityq             ,
                                  self->inCore                  ,
                                  ( qcState->densityq != NULL ) ,
                                  &(qcState->equad)             ,
                                  &(qcState->rhoquad)           ,
# ifdef USEOPENMP
                                  qcState->fockAArray           ,
                                  qcState->fockBArray           ,
# else
                                  qcState->densityp->fock       ,
                                  fockq                         ,
# endif
                                  NULL                          ,
                                  NULL                          ) ;

# ifdef USEOPENMP
        /* . Accumulate the Fock matrices. */
        QCModelDFTState_AccumulateFock ( qcState, qcState->fockAArray ) ;
        qcState->fockAArray[0] = NULL ;
        if ( qcState->densityq != NULL )
        {
            QCModelDFTState_AccumulateFock ( qcState, qcState->fockBArray ) ;
            qcState->fockBArray[0] = NULL ;
        }
# endif

        /* . Sum the total electronic energy. */
        qcState->eelectronic = qcState->eocc + qcState->eoei + qcState->equad + qcState->etei ;
        if ( eelectronic != NULL ) (*eelectronic) = qcState->eelectronic ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFT_Gradients ( const QCModelDFT *self, QCModelDFTState *qcState )
{
    if ( ( self != NULL ) && ( qcState != NULL ) && ( qcState->gradients3 != NULL ) )
    {
        auto Coordinates3 *qcGradients3 ;

        /* . Make the weighted density matrix. */
        qcState->wdensity = QCOnePDM_MakeWeighted ( qcState->densityp, qcState->densityq, NULL ) ;
        QCModelDFT_QCMMMakeWeightedDensity ( self, qcState ) ;

        /* . Get a gradient array - automatically initialized. */
        qcGradients3 = Coordinates3_Allocate ( qcState->qcAtoms->natoms ) ;
        Coordinates3_Set ( qcGradients3, 0.0e+00 ) ;
# ifdef USEOPENMP
        qcState->qcGradients3Array[0] = qcGradients3 ;
# endif

        /* . Gaussian basis derivatives - EF/EN/FF/KE + Overlap/NN . */
        QCOnePDM_AlphaBetaToTotalSpin   ( qcState->densityp, qcState->densityq ) ;
        GaussianBasis_Electron_FitD     ( qcState->qcAtoms, qcState->qcParameters, qcState->qccoordinates3, qcState->densityp->density, qcState->fpotential, qcState->wvector, qcGradients3 ) ;
        GaussianBasis_Electron_NuclearD ( qcState->qcAtoms, qcState->qcParameters, qcState->qccoordinates3, qcState->densityp->density, qcGradients3 ) ;
        GaussianBasis_Kinetic_2OverlapD ( qcState->qcAtoms, qcState->qcParameters, qcState->qccoordinates3, qcState->densityp->density, qcState->wdensity, qcGradients3 ) ;
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;
        GaussianBasis_Fit_FitD          ( qcState->qcAtoms, qcState->qcParameters, qcState->qccoordinates3, qcState->fpotential, qcState->wvector, qcGradients3 ) ;
        GaussianBasis_Nuclear_NuclearD  ( qcState->qcAtoms, qcState->qccoordinates3, qcGradients3 ) ;

        /* . DFT derivatives. */
        DFTIntegrator_Integrate ( self->functionalModel         ,
                                  qcState->dftgrid              ,
                                  qcState->qcAtoms              ,
                                  qcState->qcParameters         ,
                                  qcState->qccoordinates3       ,
                                  qcState->densityp             ,
                                  qcState->densityq             ,
                                  False                         ,
                                  ( qcState->densityq != NULL ) ,
                                  NULL                          ,
                                  NULL                          ,
                                  NULL                          ,
                                  NULL                          ,
# ifdef USEOPENMP
                                  qcState->qcGradients3Array    ,
# else
                                  qcGradients3                  ,
# endif
                                  NULL                          ) ;

# ifdef USEOPENMP
        /* . Accumulate gradients. */
        QCModelDFTState_AccumulateGradients3 ( qcState ) ;
        qcState->qcGradients3Array[0] = NULL ;
# endif

        /* . Copy the gradients into their appropriate places. */
        QCAtomContainer_SetGradients3 ( qcState->qcAtoms, qcState->coordinates3, qcGradients3, True, &(qcState->gradients3) ) ;
        Coordinates3_Deallocate ( &qcGradients3 ) ;

# ifdef TESTGRADIENTS
        QCModelDFT_GradientsTest ( qcState->qcAtoms        ,
                                   qcState->qcParameters   ,
                                   qcState->qccoordinates3 ,
                                   self->functionalModel   ,
                                   qcState->dftgrid        ,
                                   qcState->densityp       ,
                                   qcState->densityq       ,
                                   qcState->fpotential     ,
                                   qcState->wdensity       ,
                                   qcState->gradients3     ) ;
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate grid point densities in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define BASISFUNCTIONTOLERANCE 1.0e-10
void QCModelDFT_GridPointDensities ( const QCModelDFT *self, QCModelDFTState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data, const Boolean QSPIN )
{
    Real1DArray_Set ( data, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( coordinates3 != NULL ) && ( points != NULL ) && ( data != NULL ) && ( points->length0 <= data->length ) )
    {
        auto Real                   tolerance      = BASISFUNCTIONTOLERANCE ;
        auto Coordinates3          *qcCoordinates3 = NULL ;
        auto GridFunctionDataBlock *basisData      = NULL ;
        auto QCOnePDM              *density        = NULL ;

        /* . Get the QC coordinates. */
        QCAtomContainer_GetCoordinates3 ( qcState->qcAtoms, coordinates3, True, &qcCoordinates3 ) ;

        /* . Get the correct density. */
        QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;
        if ( QSPIN ) density = qcState->densityq ;
        else         density = qcState->densityp ;

        /* . Generate data. */
        basisData = GridFunctionDataBlock_Allocate ( qcState->qcAtoms->nobasisw, points->length0, 0, NULL ) ;
        if ( basisData != NULL )
        {
            QCAtomContainer_OrbitalBasisGridPointValues ( qcState->qcAtoms, qcState->qcParameters, qcCoordinates3, points, True, &tolerance, basisData, NULL ) ;
            QCOnePDM_DensityGridValues ( density, basisData, data, NULL ) ;
        }

        /* . Finish up. */
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Deallocate space. */
        Coordinates3_Deallocate          ( &qcCoordinates3 ) ;
        GridFunctionDataBlock_Deallocate ( &basisData      ) ;
    }
}
# undef BASISFUNCTIONTOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate grid point orbitals in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define BASISFUNCTIONTOLERANCE 1.0e-10
void QCModelDFT_GridPointOrbitals ( const QCModelDFT *self, QCModelDFTState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *data, const Integer norbitals, Integer *orbitalIndices, const Boolean QALPHA )
{
    Real1DArray_Set ( data, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( coordinates3 != NULL ) && ( points != NULL ) && ( data != NULL ) && ( points->length0 <= data->length ) )
    {
        auto Boolean       QOK = True ;
        auto Integer       i, *indices, n ;
        auto Coordinates3 *qcCoordinates3 = NULL ;
        auto QCOnePDM     *density  = NULL ;
        auto Real2DArray  *orbitals = NULL ;

        /* . Get the QC coordinates. */
        QCAtomContainer_GetCoordinates3 ( qcState->qcAtoms, coordinates3, True, &qcCoordinates3 ) ;

        /* . Get the correct orbital set. */
        if ( QALPHA ) density = qcState->densityp ;
        else          density = qcState->densityq ;
        if ( density  != NULL ) orbitals = density->orbitals ;
        if ( orbitals != NULL )
        {
            /* . Check the orbital indices. */
            if ( ( norbitals > 0 ) && ( orbitalIndices != NULL ) )
            {
                n       = norbitals ;
                indices = orbitalIndices ;
                for ( i = 0 ; i < n ; i++ )
                {
                    if ( ( indices[i] < 0 ) || ( indices[i] > ( orbitals->length1 - 1 ) ) ) { QOK = False ; break ; }
                }
            }
            else
            {
                n          = 1 ;
                indices    = Memory_Allocate_Array_Integer ( n ) ;
                indices[0] = QCOnePDM_HOMOIndex ( density ) ;
                QOK        = ( ( indices[0] >= 0 ) && ( indices[0] < orbitals->length1 ) ) ;
            }
            if ( QOK )
            {
                auto Real                   tolerance      = BASISFUNCTIONTOLERANCE ;
                auto GridFunctionDataBlock *basisData      = NULL ;
                auto Real2DArray           *reduced        = NULL, view ;

                /* . Generate basis data. */
                basisData = GridFunctionDataBlock_Allocate ( qcState->qcAtoms->nobasisw, points->length0, 0, NULL ) ;
                if ( basisData != NULL )
                {
                    QCAtomContainer_OrbitalBasisGridPointValues ( qcState->qcAtoms, qcState->qcParameters, qcCoordinates3, points, True, &tolerance, basisData, NULL ) ;

                    /* . Extract orbital data. */
                    reduced = Real2DArray_Allocate ( basisData->indices->length, n, NULL ) ;
                    if ( reduced != NULL )
                    {
                        auto Integer i, m, o ;
                        for ( i = 0 ; i < basisData->indices->length ; i++ )
                        {
                            m = Integer1DArray_Item ( basisData->indices, i ) ;
                            for ( o = 0 ; o < n ; o++ ) Real2DArray_Item ( reduced, i, o ) = Real2DArray_Item ( orbitals, m, indices[o] ) ;
                        }
                        /* . Generate data - this needs to be better. */
                        Real2DArray_ViewOfRaw      ( &view, 0, points->length0, n, n, 1, Real1DArray_Data ( data ), data->size, NULL ) ;
                        Real2DArray_MatrixMultiply ( True, False, 1.0e+00, basisData->f, reduced, 0.0e+00, &view, NULL ) ;
                    }
                }

                /* . Deallocate space. */
                GridFunctionDataBlock_Deallocate ( &basisData ) ;
                Real2DArray_Deallocate           ( &reduced   ) ;
            }
            if ( orbitalIndices == NULL ) Memory_Deallocate_Integer ( &indices ) ;
        }
        Coordinates3_Deallocate ( &qcCoordinates3 ) ;
    }
}
# undef BASISFUNCTIONTOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate grid point potentials in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFT_GridPointPotentials ( const QCModelDFT *self, QCModelDFTState *qcState, const Coordinates3 *coordinates3, const Coordinates3 *points, Real1DArray *potentials )
{
    Real1DArray_Set ( potentials, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( qcState != NULL ) && ( coordinates3 != NULL ) && ( points != NULL ) && ( potentials != NULL ) && ( points->length0 <= potentials->length ) )
    {
        auto Coordinates3 *qcCoordinates3 = NULL ;
        auto Real1DArray  *znuclear = NULL ;

        /* . Get the QC coordinates. */
        QCAtomContainer_GetCoordinates3 ( qcState->qcAtoms, coordinates3, True, &qcCoordinates3 ) ;

        /* . Electronic contribution. */
        QCOnePDM_AlphaBetaToTotalSpin   ( qcState->densityp, qcState->densityq ) ;
        GaussianBasis_Point_Electron    ( qcState->qcAtoms, qcState->qcParameters, qcCoordinates3, points, qcState->densityp->density, potentials ) ;
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Nuclear contribution. */
        znuclear = Real1DArray_Allocate   ( qcState->qcAtoms->natoms, NULL ) ;
        QCAtomContainer_GetNuclearCharges ( qcState->qcAtoms, qcState->qcParameters, znuclear ) ;
        GaussianBasis_Point_Nuclear       ( qcState->qcAtoms, znuclear, qcCoordinates3, points, potentials ) ;

        /* . Deallocate space. */
        Coordinates3_Deallocate ( &qcCoordinates3 ) ;
        Real1DArray_Deallocate  ( &znuclear       ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the integrals and make the orthogonalizing transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFT_Integrals ( const QCModelDFT *self, QCModelDFTState *qcState )
{
    if ( ( self != NULL ) && ( qcState != NULL ) )
    {
        /* . Gaussian basis integrals - EF/EN/FF/KE + Overlap/NN . */
        GaussianBasis_SelfOverlap      ( qcState->qcAtoms, qcState->qcParameters, &(qcState->fselfoverlap) ) ;
        GaussianBasis_Kinetic_2Overlap ( qcState->qcAtoms, qcState->qcParameters, qcState->qccoordinates3, &(qcState->oneelectronmatrix), &(qcState->overlap) ) ;
        GaussianBasis_Electron_Nuclear ( qcState->qcAtoms, qcState->qcParameters, qcState->qccoordinates3, qcState->oneelectronmatrix ) ;
        GaussianBasis_Fit_Fit          ( qcState->qcAtoms, qcState->qcParameters, qcState->qccoordinates3, qcState->fselfoverlap, &(qcState->inversefitmatrix) ) ;
        GaussianBasis_Electron_Fit     ( qcState->qcAtoms, qcState->qcParameters, qcState->qccoordinates3, &(qcState->fitintegrals) ) ;
        qcState->enuclear = GaussianBasis_Nuclear_Nuclear ( qcState->qcAtoms, qcState->qccoordinates3 ) ;
        /* . Orthogonalizing transformation. */
        QCModelDFT_MakeOrthogonalizingTransformation ( self, qcState ) ;
        /* . Allocate orbital data. */
        QCOnePDM_AllocateOrbitalData ( qcState->densityp, qcState->orthogonalizer, NULL ) ;
        QCOnePDM_AllocateOrbitalData ( qcState->densityq, qcState->orthogonalizer, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the Lowdin transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void QCModelDFT_MakeLowdinTransformation ( const QCModelDFT *self, QCModelDFTState *qcState )
{
    if ( ( self != NULL ) && ( qcState != NULL ) )
    {
        if ( ( qcState->lowdintransformation == NULL ) || ( qcState->overlapeigenvalues == NULL ) || ( qcState->overlapeigenvectors == NULL ) )
        {
            auto Integer n = qcState->c2o->length1 ;
            auto SymmetricMatrix *overlap = NULL ;

            /* . Deallocate. */
            Real2DArray_Deallocate ( &(qcState->lowdintransformation) ) ;
            Real2DArray_Deallocate ( &(qcState->overlapeigenvectors ) ) ;
            Real1DArray_Deallocate ( &(qcState->overlapeigenvalues  ) ) ;

            /* . Reallocate. */
            qcState->lowdintransformation = Real2DArray_Allocate ( qcState->o2c->length0, n, NULL ) ;
            qcState->overlapeigenvectors  = Real2DArray_Allocate ( n, n, NULL ) ;
            qcState->overlapeigenvalues   = Real1DArray_Allocate ( n   , NULL ) ;

            /* . Get the eigenvalues and eigenvectors of the transformed overlap. */
            overlap = SymmetricMatrix_Clone ( qcState->overlap ) ;
            SymmetricMatrix_Transform_In_Place ( overlap, qcState->c2o ) ;
            SymmetricMatrix_Diagonalize ( overlap, qcState->overlapeigenvalues, qcState->overlapeigenvectors, NULL ) ;

            /* . Create the Lowdin overlap. */
            SymmetricMatrix_SquareRoot ( overlap, qcState->overlapeigenvalues, qcState->overlapeigenvectors ) ;

            /* . Form o2c * overlap. */
            SymmetricMatrix_PreMatrixMultiply ( overlap, qcState->o2c, False, qcState->lowdintransformation, NULL ) ;

            /* . Finish up. */
            SymmetricMatrix_Deallocate ( &overlap ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the orthogonalizing transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFT_MakeOrthogonalizingTransformation ( const QCModelDFT *self, QCModelDFTState *qcState )
{
    if ( ( self != NULL ) && ( qcState != NULL ) )
    {
        auto Real2DArray     *transformation ;
        auto SymmetricMatrix *overlap        ;

        /* . Deallocate existing space. */
        Real2DArray_Deallocate ( &(qcState->orthogonalizer) ) ;

        /* . Make the transformation. */
        overlap = SymmetricMatrix_AllocateN ( qcState->c2o->length1, NULL ) ;
        SymmetricMatrix_Transform ( qcState->overlap, qcState->c2o, False, overlap ) ;
        transformation = SymmetricMatrix_OrthogonalizingTransformation ( overlap, NULL, NULL, NULL ) ;
        qcState->orthogonalizer = Real2DArray_Allocate ( qcState->c2o->length0, transformation->length1, NULL ) ;
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, qcState->c2o, transformation, 0.0e+00, qcState->orthogonalizer, NULL ) ;

        /* . Clean up.*/
        Real2DArray_Deallocate     ( &transformation ) ;
        SymmetricMatrix_Deallocate ( &overlap        ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Mayer bond orders.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFT_MayerBondOrders ( const QCModelDFT *self, QCModelDFTState *qcState, const QCChargeModel *qcChargeModel, SymmetricMatrix *bondorders, Real1DArray *charges, Real1DArray *freevalence, Real1DArray *totalvalence )
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
        auto QCChargeModel model ;

        /* . Get the model. */
        if ( qcChargeModel == NULL ) model = self->qcChargeModel ;
        else                         model = (*qcChargeModel)    ;

        /* . Convert to total and spin densities. */
        QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Mulliken bond orders.*/
        if ( model == QCChargeModel_Mulliken )
        {
            auto Real2DArray *ps = NULL ;

            /* . Allocate space. */
            ps = Real2DArray_Allocate ( qcState->overlap->dimension, qcState->overlap->dimension, NULL ) ;
            Real2DArray_Set ( ps, 0.0e+00 ) ;

            /* . Total density. */
            SymmetricMatrix_Multiply2 ( qcState->densityp->density, qcState->overlap, ps ) ;
            for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
            {
                istart = qcState->qcAtoms->data[iatom].ostartw ;
                istop  = istart + qcState->qcAtoms->data[iatom].nobasisw ;
                for ( jatom = 0 ; jatom <= iatom ; jatom++ )
                {
                    jstart = qcState->qcAtoms->data[jatom].ostartw ;
                    jstop  = jstart + qcState->qcAtoms->data[jatom].nobasisw ;
                    for ( ibf = istart, sum = 0.0e+00 ; ibf < istop ; ibf++ )
                    {
                        for ( jbf = jstart ; jbf < jstop ; jbf++ ) sum += Real2DArray_Item ( ps, ibf, jbf ) * Real2DArray_Item ( ps, jbf, ibf ) ;
                    }
                    SymmetricMatrix_Set_Component ( bondorders, iatom, jatom, sum ) ;
                }
            }

            /* . Total valence - first contribution. */
            for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ ) Real1DArray_Item ( totalvalence, iatom ) = - SymmetricMatrix_Get_Component ( bondorders, iatom, iatom ) ;

            /* . Spin density. */
            if ( qcState->densityq != NULL )
            {
                SymmetricMatrix_Multiply2 ( qcState->densityq->density, qcState->overlap, ps ) ;
                for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
                {
                    istart = qcState->qcAtoms->data[iatom].ostartw ;
                    istop  = istart + qcState->qcAtoms->data[iatom].nobasisw ;
                    for ( jatom = 0 ; jatom <= iatom ; jatom++ )
                    {
                        jstart = qcState->qcAtoms->data[jatom].ostartw ;
                        jstop  = jstart + qcState->qcAtoms->data[jatom].nobasisw ;
                        for ( ibf = istart, sum = 0.0e+00 ; ibf < istop ; ibf++ )
                        {
                            for ( jbf = jstart ; jbf < jstop ; jbf++ ) sum += Real2DArray_Item ( ps, ibf, jbf ) * Real2DArray_Item ( ps, jbf, ibf ) ;
                        }
                        SymmetricMatrix_IncrementComponent ( bondorders, iatom, jatom, sum ) ;
                    }
                }
            }

            /* . Finish up. */
            Real2DArray_Deallocate ( &ps ) ;
        }
        /* . Other cases - defaults to Lowdin. */
        else
        {
            auto SymmetricMatrix *ps = NULL ;

            /* . Make the Lowdin transformation.*/
            QCModelDFT_MakeLowdinTransformation ( self, qcState ) ;

            /* . Allocate space. */
            ps = SymmetricMatrix_Allocate ( qcState->lowdintransformation->length1 ) ;
            SymmetricMatrix_Set ( ps, 0.0e+00 ) ;

            /* . Total density. */
            SymmetricMatrix_Transform ( qcState->densityp->density, qcState->lowdintransformation, False, ps ) ;
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
                        for ( jbf = jstart ; jbf < jstop ; jbf++ ) sum += pow ( SymmetricMatrix_Get_Component ( ps, ibf, jbf ), 2 ) ;
                    }
                    SymmetricMatrix_Set_Component ( bondorders, iatom, jatom, sum ) ;
                }
            }

/* . Is there a problem with the definition of totalvalence for a Lowdin orthogonalized density? They are often too large. */

            /* . Total valence - first contribution. */
            for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ ) Real1DArray_Item ( totalvalence, iatom ) = - SymmetricMatrix_Get_Component ( bondorders, iatom, iatom ) ;

            /* . Spin density. */
            if ( qcState->densityq != NULL )
            {
                SymmetricMatrix_Transform ( qcState->densityq->density, qcState->lowdintransformation, False, ps ) ;
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
                            for ( jbf = jstart ; jbf < jstop ; jbf++ ) sum += pow ( SymmetricMatrix_Get_Component ( ps, ibf, jbf ), 2 ) ;
                        }
                        SymmetricMatrix_IncrementComponent ( bondorders, iatom, jatom, sum ) ;
                    }
                }
            }

            /* . Finish up. */
            SymmetricMatrix_Deallocate ( &ps ) ;
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
        Real1DArray_Scale          ( charges, 0.5e+00 ) ;
        Real1DArray_AddScaledArray ( freevalence, 1.0e+00, totalvalence, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the orbital energies and HOMO and LUMO indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real1DArray *QCModelDFT_OrbitalEnergies ( const QCModelDFT *self, QCModelDFTState *qcState, const Boolean QALPHA, Integer *homo, Integer *lumo )
{
    Real1DArray *data = NULL ;
    if ( homo != NULL ) (*homo) = -1 ;
    if ( lumo != NULL ) (*lumo) = -1 ;
    if ( ( self != NULL ) && ( qcState != NULL ) )
    {
        auto QCOnePDM *density = NULL ;

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
Real2DArray *QCModelDFT_Orbitals ( const QCModelDFT *self, QCModelDFTState *qcState, const Boolean QALPHA )
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
! . QC/MM Fock contributions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . This should really be called from QCModelDFT_Fock (after Gaussian_Fock) to avoid allocating extra storage. */
void QCModelDFT_QCMMFock ( const QCModelDFT *self, QCModelDFTState *qcState, Real *eelectronic )
{
    if ( ( self != NULL ) && ( qcState != NULL ) && ( qcState->qcmmstate != NULL ) )
    {
        auto Real1DArray     *work ;
        auto SymmetricMatrix *fock ;

        /* . Get the charges. */
        QCModelDFT_AtomicCharges ( self, qcState, NULL, False, qcState->qcmmstate->qcCharges ) ;

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

        /* . Get the place where to build Fock. */
        if ( qcState->densityq == NULL ) fock = qcState->densityp->fock ;
        else
        {
            fock = SymmetricMatrix_Allocate ( qcState->densityp->fock->dimension ) ;
            SymmetricMatrix_Set ( fock, 0.0e+00 ) ;
        }

        /* . Switch on the type of charge model. */
        switch ( self->qcChargeModel )
        {
             /* . Coulomb fitting. */
            case QCChargeModel_CoulombFitting:
            {
                auto Real  f ;
                auto Integer     iatom, ibf, istart, istop ;
                auto Real1DArray *vvector ;

                /* . Allocate space. */
                vvector = Real1DArray_Allocate ( Real1DArray_Length ( qcState->fpotential ), NULL ) ;
                Real1DArray_Set ( vvector, 0.0e+00 ) ;

                /* . Create dE/df from dE/dq. */
                for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
                {
                    istart = qcState->qcAtoms->data[iatom].fstartw ;
                    istop  = istart + qcState->qcAtoms->data[iatom].nfbasisw ;
                    f      = Real1DArray_Item ( work, iatom ) ;
                    for ( ibf = istart; ibf < istop ; ibf++ ) Real1DArray_Item ( vvector, ibf ) = - f * Real1DArray_Item ( qcState->fselfoverlap, ibf ) ;
                }

                /* . Find the w vector by multiplying by the inverse matrix. */
                SymmetricMatrix_VectorMultiply ( qcState->inversefitmatrix, vvector, qcState->wvector, NULL ) ;

                /* . Clear up. */
                Real1DArray_Deallocate ( &vvector ) ;

                /* . Build Fock. */
                {
                    auto Block *block ;
                    auto Integer         i     ;
                    /* . Loop over the integral blocks. */
                    List_Iterate_Initialize ( qcState->fitintegrals->blocks ) ;
                    while ( ( block = BlockStorage_Iterate ( qcState->fitintegrals ) ) != NULL )
                    {
                        for ( i = 0 ; i < block->ndata ; i++ ) fock->data[block->indices32[i]] += qcState->wvector->data[block->indices16[i]] * block->data[i] ;
                    }
                }
            }
            break ;
            /* . Lowdin. */
            case QCChargeModel_Lowdin:
            {
                auto Real f, p ;
                auto Integer    iatom, i, istart, istop, j, k ;

                /* . Build Fock. */
                for ( j = 0 ; j < fock->dimension ; j++ )
                {
                    for ( k = 0 ; k <= j ; k++ )
                    {
                        for ( iatom = 0, f = 0.0e+00 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
                        {
                            p = Real1DArray_Item ( work, iatom ) ;
                            istart = qcState->qcAtoms->data[iatom].ostart ;
                            istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
                            for ( i = istart ; i < istop ; i++ ) f += p * Real2DArray_Item ( qcState->lowdintransformation, j, i ) * Real2DArray_Item ( qcState->lowdintransformation, k, i );
                        }
                        SymmetricMatrix_IncrementComponent ( fock, j, k, -f ) ;
                    }
                }
            }
            break ;
            /* . Mulliken. */
            case QCChargeModel_Mulliken:
            {
                auto Real p, pq, q ;
                auto Integer    iatom, ibf, istart, istop, jatom, jbf, jstart, jstop ;

                /* . Build Fock. */
                for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
                {
                    p = Real1DArray_Item ( work, iatom ) ;
                    istart = qcState->qcAtoms->data[iatom].ostartw ;
                    istop  = istart + qcState->qcAtoms->data[iatom].nobasisw ;
                    /* . Off-diagonal. */
                    for ( jatom = 0 ; jatom < iatom ; jatom++ )
                    {
                        q  = Real1DArray_Item ( work, jatom ) ;
                        pq = 0.5e+00 * ( p + q ) ;
                        jstart = qcState->qcAtoms->data[jatom].ostartw ;
                        jstop  = jstart + qcState->qcAtoms->data[jatom].nobasisw ;
                        for ( ibf = istart ; ibf < istop ; ibf++ )
                        {
                            for ( jbf = jstart ; jbf < jstop ; jbf++ ) SymmetricMatrix_IncrementComponent ( fock, ibf, jbf, - pq * SymmetricMatrix_Get_Component ( qcState->overlap, ibf, jbf ) ) ;
                        }
                    }
                    /* . Diagonal. */
                    for ( ibf = istart ; ibf < istop ; ibf++ )
                    {
                        for ( jbf = istart ; jbf <= ibf ; jbf++ ) SymmetricMatrix_IncrementComponent ( fock, ibf, jbf, - p * SymmetricMatrix_Get_Component ( qcState->overlap, ibf, jbf ) ) ;
                    }
                }
            }
            break ;
        }

        /* . Make sure the contributions are in the correct place when there are two Fock matrices. */
        if ( qcState->densityq != NULL )
        {
            SymmetricMatrix_AddScaledMatrix ( qcState->densityp->fock, 1.0e+00, fock ) ;
            SymmetricMatrix_AddScaledMatrix ( qcState->densityq->fock, 1.0e+00, fock ) ;
            SymmetricMatrix_Deallocate ( &fock ) ;
        }

        /* . Clear up. */
        if ( qcState->qcmmstate->qcqcPotentials != NULL ) Real1DArray_Deallocate ( &work ) ;

        /* . Add in the energy contributions. */
        if ( eelectronic != NULL ) (*eelectronic) += qcState->eqcmm + qcState->eqcqc ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM weighted density matrix contributions - not Coulomb fitting charges.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define EIGENVALUETOLERANCE 1.0e-10
void QCModelDFT_QCMMMakeWeightedDensity ( const QCModelDFT *self, QCModelDFTState *qcState )
{
    if ( ( self != NULL ) && ( self->qcChargeModel != QCChargeModel_CoulombFitting ) && ( qcState != NULL ) && ( qcState->qcmmstate != NULL ) && ( qcState->wdensity != NULL ) )
    {
        auto Real1DArray *work ;

        /* . QC/QC. */
        if ( qcState->qcmmstate->qcqcPotentials != NULL )
        {
            /* . Get 2 * F * q + phi. */
            work = Real1DArray_Allocate ( qcState->qcmmstate->qcCharges->length, NULL ) ;
            SymmetricMatrix_VectorMultiply ( qcState->qcmmstate->qcqcPotentials, qcState->qcmmstate->qcCharges, work, NULL ) ;
            Real1DArray_Scale ( work, 2.0e+00 ) ;
            Real1DArray_Add   ( work, qcState->qcmmstate->qcmmPotentials, NULL ) ;
        }
        /* . QC/MM only. */
        else work = qcState->qcmmstate->qcmmPotentials ;

        /* . Make total density. */
        QCOnePDM_AlphaBetaToTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Lowdin. */
        if ( self->qcChargeModel == QCChargeModel_Lowdin )
        {
            auto Real a, ab, b, f ;
            auto Integer    i, iatom, ij, istart, istop, j, jatom, jstart, jstop, k, nc, no ;
            auto Real2DArray          *po2c ;
            auto SymmetricMatrix *temp, *temp2 ;

            /* . Get dimensions. */
            nc = qcState->o2c->length0 ;
            no = qcState->o2c->length1 ;

            /* . Allocate space. */
            po2c = Real2DArray_Allocate ( nc, no, NULL ) ;
            temp = SymmetricMatrix_Allocate ( no ) ;
            SymmetricMatrix_Set ( temp, 0.0e+00 ) ;

            /* . Form P * o2c. */
            SymmetricMatrix_PostMatrixMultiply ( qcState->densityp->density, qcState->o2c, False, po2c ) ;

            /* . Calculate the symmetrized core matrix. */
            /* . Outer atom loop. */
            for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
            {
                a = Real1DArray_Item ( work, iatom ) ;
                istart = qcState->qcAtoms->data[iatom].ostart ;
                istop  = istart + qcState->qcAtoms->data[iatom].nobasis ;
                for ( i = istart ; i < istop ; i++ )
                {
                    /* . Inner atom loop. */
                    for ( jatom = 0 ; jatom <= iatom ; jatom++ )
                    {
                        b = Real1DArray_Item ( work, jatom ) ;
                        jstart = qcState->qcAtoms->data[jatom].ostart ;
                        jstop  = Minimum ( i, jstart + qcState->qcAtoms->data[jatom].nobasis - 1 ) ;
                        for ( j = jstart ; j <= jstop ; j++ )
                        {
                            for ( k = 0, f = 0.0e+00 ; k < nc ; k++ )
                            {
                                f += b * Real2DArray_Item ( po2c, k, i ) * Real2DArray_Item ( qcState->lowdintransformation, k, j ) +
                                     a * Real2DArray_Item ( po2c, k, j ) * Real2DArray_Item ( qcState->lowdintransformation, k, i ) ;
                            }
                            SymmetricMatrix_Set_Component ( temp, i, j, f ) ;
                        }
                    }
                }
            }

            /* . First transformation. */
            SymmetricMatrix_Transform_In_Place ( temp, qcState->overlapeigenvectors ) ;

            /* . Scale by the eigenvalue factors. */
/* . Should elements be zeroed for each a and b that are small or for each (a+b) individually? */
            for ( i = ij = 0 ; i < no ; i++ )
            {
                a = Real1DArray_Item ( qcState->overlapeigenvalues, i ) ;
                for ( j = 0 ; j <= i ; ij++, j++ )
                {
                    b  = Real1DArray_Item ( qcState->overlapeigenvalues, j ) ;
                    ab = a + b ;
                    if ( ab > EIGENVALUETOLERANCE ) temp->data[ij] /= ab ;
                    else                            temp->data[ij] = 0.0e+00 ;
                }
            }

            /* . Second transformation. */
            Real2DArray_Transpose ( qcState->overlapeigenvectors, NULL ) ;
            SymmetricMatrix_Transform_In_Place ( temp, qcState->overlapeigenvectors ) ;
            Real2DArray_Transpose ( qcState->overlapeigenvectors, NULL ) ;

            /* . Back transform temp by Z * temp * Z^T. */
            temp2 = SymmetricMatrix_Allocate ( qcState->wdensity->dimension ) ;
            SymmetricMatrix_Transform ( temp, qcState->c2o, True, temp2 ) ;

            /* . Add in the contributions to wdensity. */
            SymmetricMatrix_AddScaledMatrix ( qcState->wdensity, -2.0e+00, temp2 ) ;

            /* . Deallocate space. */
            Real2DArray_Deallocate          ( &po2c  ) ;
            SymmetricMatrix_Deallocate ( &temp  ) ;
            SymmetricMatrix_Deallocate ( &temp2 ) ;
        }
        /* . Mulliken. */
        else if ( self->qcChargeModel == QCChargeModel_Mulliken )
        {
            auto Real p, pq, q ;
            auto Integer    iatom, ibf, istart, istop, jatom, jbf, jstart, jstop ;

            /* . Off-diagonal atom-atom contributions only. */
            for ( iatom = 0 ; iatom < qcState->qcAtoms->natoms ; iatom++ )
            {
                p = Real1DArray_Item ( work, iatom ) ;
                istart = qcState->qcAtoms->data[iatom].ostartw ;
                istop  = istart + qcState->qcAtoms->data[iatom].nobasisw ;
                /* . Off-diagonal only. */
                for ( jatom = 0 ; jatom < iatom ; jatom++ )
                {
                    q  = Real1DArray_Item ( work, jatom ) ;
                    pq = - ( p + q ) ;
                    jstart = qcState->qcAtoms->data[jatom].ostartw ;
                    jstop  = jstart + qcState->qcAtoms->data[jatom].nobasisw ;
                    for ( ibf = istart ; ibf < istop ; ibf++ )
                    {
                        for ( jbf = jstart ; jbf < jstop ; jbf++ )
                        {
                            SymmetricMatrix_IncrementComponent ( qcState->wdensity, ibf, jbf, pq * SymmetricMatrix_Get_Component ( qcState->densityp->density, ibf, jbf ) ) ;
                        }
                    }
                }
            }
        }

        /* . Convert back. */
        QCOnePDM_AlphaBetaFromTotalSpin ( qcState->densityp, qcState->densityq ) ;

        /* . Clear up. */
        if ( qcState->qcmmstate->qcqcPotentials != NULL ) Real1DArray_Deallocate ( &work ) ;
    }
}
# undef EIGENVALUETOLERANCE
