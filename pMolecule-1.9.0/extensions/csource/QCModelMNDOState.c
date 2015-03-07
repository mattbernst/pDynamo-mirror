/*------------------------------------------------------------------------------
! . File      : QCModelMNDOState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module implements MNDO QC model state procedures.
!=================================================================================================================================*/

# include "DefineStatements.h"
# include "Memory.h"
# include "QCModelMNDOState.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCModelMNDOState *QCModelMNDOState_Allocate ( void )
{
    QCModelMNDOState *self = NULL ;
    /* . Allocation. */
    self = ( QCModelMNDOState * ) Memory_Allocate ( sizeof ( QCModelMNDOState ) ) ;
    /* . Initialization. */
    if ( self != NULL )
    {
        /* . Scalars. */
        self->isConverged  = False   ;
        self->eelectronic = 0.0e+00 ;
        self->enuclear    = 0.0e+00 ;
        self->eocc        = 0.0e+00 ;
        self->eoei        = 0.0e+00 ;
        self->eqcmm       = 0.0e+00 ;
        self->eqcqc       = 0.0e+00 ;
        self->etei        = 0.0e+00 ;
        self->ncycles     = -1      ;
        /* . Aliases - persistent. */
        self->qcAtoms              = NULL ;
        self->qcParameters         = NULL ;
        /* . Aliases - temporary. */
        self->coordinates3         = NULL ;
        self->gradients3           = NULL ;
        self->qcmmstate            = NULL ;
        /* . Allocated data - persistent. */
        self->c2o                  = NULL ;
        self->densityp             = NULL ;
        self->densityq             = NULL ;
        /* . Allocated data - temporary. */
        self->qccoordinates3       = NULL ;
        self->qptransformation     = NULL ;
        self->oneelectronmatrix    = NULL ;
        self->twoelectronintegrals = NULL ;
        /* . Backup data. */
        self->fockA                = NULL ;
        self->fockB                = NULL ;
# ifdef MNDOCI
        /* . CI state. */
        self->cistate              = NULL ;
# endif
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Backup the Fock matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDOState_BackupFock ( QCModelMNDOState *self, Status *status )
{
    if ( self != NULL )
    {
        auto Boolean isOK = True ;
        if ( self->densityp != NULL )
        {
            if ( self->fockA == NULL ) self->fockA = SymmetricMatrix_Clone ( self->densityp->fock ) ;
            else                       SymmetricMatrix_CopyTo ( self->densityp->fock, self->fockA ) ;
            isOK = isOK && ( self->fockA != NULL ) && ( self->densityp->fock != NULL ) ;
        }
        if ( self->densityq != NULL )
        {
            if ( self->fockB == NULL ) self->fockB = SymmetricMatrix_Clone ( self->densityq->fock ) ;
            else                       SymmetricMatrix_CopyTo ( self->densityq->fock, self->fockB ) ;
            isOK = isOK && ( self->fockB != NULL ) && ( self->densityq->fock != NULL ) ;
        }
        if ( ! isOK ) Status_Set ( status, Status_InvalidArrayOperation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDOState_Deallocate ( QCModelMNDOState **self )
{
    if ( (*self) != NULL )
    {
        Real2DArray_Deallocate     ( &((*self)->c2o                 ) ) ;
        QCOnePDM_Deallocate        ( &((*self)->densityp            ) ) ;
        QCOnePDM_Deallocate        ( &((*self)->densityq            ) ) ;
        Coordinates3_Deallocate    ( &((*self)->qccoordinates3      ) ) ;
        Real2DArray_Deallocate     ( &((*self)->qptransformation    ) ) ;
        SymmetricMatrix_Deallocate ( &((*self)->oneelectronmatrix   ) ) ;
        BlockStorage_Deallocate    ( &((*self)->twoelectronintegrals) ) ;
        /* . Backup data. */
        SymmetricMatrix_Deallocate ( &((*self)->fockA               ) ) ;
        SymmetricMatrix_Deallocate ( &((*self)->fockB               ) ) ;
# ifdef MNDOCI
        MNDOCIState_Deallocate     ( &((*self)->cistate             ) ) ;
# endif
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the electronic energy given complete Fock matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real QCModelMNDOState_EnergyFromFock ( const QCModelMNDOState *self )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->densityp != NULL ) && ( self->oneelectronmatrix != NULL ) )
    {
/*
        energy = self->densityp->occupancyEnergy + 0.5e+00 * ( SymmetricMatrix_Multiply2_Trace ( self->densityp->density, self->oneelectronmatrix ) +
                                                               SymmetricMatrix_Multiply2_Trace ( self->densityp->density, self->densityp->fock    ) ) ;
        if ( self->densityq != NULL )
        {
            energy += self->densityq->occupancyEnergy + 0.5e+00 * ( SymmetricMatrix_Multiply2_Trace ( self->densityq->density, self->oneelectronmatrix ) +
                                                                    SymmetricMatrix_Multiply2_Trace ( self->densityq->density, self->densityq->fock    ) ) ;
        }
*/
/*
printf ( "\nENERGY DENSITY:\n" ) ;
SymmetricMatrix_Print ( self->densityp->density ) ;
printf ( "\nENERGY FOCK:\n" ) ;
SymmetricMatrix_Print ( self->densityp->fock ) ;
*/
            energy  = ( self->densityp->occupancyEnergy + SymmetricMatrix_Multiply2_Trace ( self->densityp->density, self->densityp->fock ) ) ;
/*
printf ( "\nENERGIES = %20.10f %20.10f %20.10f\n", energy, self->densityp->occupancyEnergy, energy - self->densityp->occupancyEnergy ) ;
*/
        if ( self->densityq != NULL )
        {
            energy += ( self->densityq->occupancyEnergy + SymmetricMatrix_Multiply2_Trace ( self->densityq->density, self->densityq->fock ) ) ;
        }
    }
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Finalization for an energy/gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDOState_Finalize ( QCModelMNDOState *self, const Boolean QKEEPDATA )
{
    if ( self != NULL )
    {
# ifdef MNDOCI
        auto Boolean keepWavefunction ;
# endif
        if ( ! QKEEPDATA )
        {
            /* . Integrals, etc. */
            Coordinates3_Deallocate    ( &(self->qccoordinates3      ) ) ;
            Real2DArray_Deallocate     ( &(self->qptransformation    ) ) ;
            SymmetricMatrix_Deallocate ( &(self->oneelectronmatrix   ) ) ;
            BlockStorage_Deallocate    ( &(self->twoelectronintegrals) ) ;
            /* . Orbital data. */
            QCOnePDM_DeallocateOrbitalData ( self->densityp ) ;
            QCOnePDM_DeallocateOrbitalData ( self->densityq ) ;
        }
# ifdef MNDOCI
        /* . CI. */
        keepWavefunction = QKEEPDATA ;
        MNDOCIState_Finalize ( self->cistate, keepWavefunction ) ;
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the densities for the gradient terms between atoms i and j.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDOState_GetGradientDensityTerms ( QCModelMNDOState *self, const MNDOParameters *idata, const Integer ifirst, const MNDOParameters *jdata, const Integer jfirst,
                                                                                                                 const SymmetricMatrix *dtotal, const SymmetricMatrix *dspin,
                                                                                                                       Real1DArray **dOneI, Real1DArray **dOneJ, Real2DArray **dTwoIJ )
{
    Integer i, ij, j, k, kk, kl, l, ll, m, mk, ml, mn, n, ni, nj, nk, nl ;
    Real    aa, bb, f ;

    Real2DArray *densityij ;
    Real1DArray *densityi, *densityj ;

    /* . Allocate space. */
    ni = ( idata->norbitals * ( idata->norbitals + 1 ) ) / 2 ;
    nj = ( jdata->norbitals * ( jdata->norbitals + 1 ) ) / 2 ;
    densityi  = Real1DArray_Allocate ( ni    , NULL ) ; Real1DArray_Set ( densityi , 0.0e+00 ) ;
    densityj  = Real1DArray_Allocate ( nj    , NULL ) ; Real1DArray_Set ( densityj , 0.0e+00 ) ;
    densityij = Real2DArray_Allocate ( ni, nj, NULL ) ; Real2DArray_Set ( densityij, 0.0e+00 ) ;
    if ( ( densityi != NULL ) && ( densityj != NULL ) && ( densityij != NULL ) )
    {
        /* . One-center terms - same for CI and non-CI cases. */
        for ( i = ij = 0 ; i < idata->norbitals ; i++, ij++ )
        {
            mn = ( ( i + ifirst ) * ( i + ifirst + 1 ) ) / 2 + ifirst ;
            for ( j = 0 ; j < i ; ij++, j++, mn++ ) Real1DArray_Item ( densityi, ij ) = 2.0e+00 * dtotal->data[mn] ;
            Real1DArray_Item ( densityi, ij ) = dtotal->data[mn] ;
        }
        for ( i = ij = 0 ; i < jdata->norbitals ; i++, ij++ )
        {
            mn = ( ( i + jfirst ) * ( i + jfirst + 1 ) ) / 2 + jfirst ;
            for ( j = 0 ; j < i ; ij++, j++, mn++ ) Real1DArray_Item ( densityj, ij ) = 2.0e+00 * dtotal->data[mn] ;
            Real1DArray_Item ( densityj, ij ) = dtotal->data[mn] ;
        }

        /* . Two-center terms. */
# ifdef MNDOCI
        /* . Non-CI case. */
        if ( self->cistate == NULL )
        {
# endif
            /* . Exchange. */
            for ( k = ifirst ; k < ifirst + idata->norbitals ; ++k )
            {
	        aa = 1.0e+00 ;
	        kk = k * ( k + 1 ) / 2 ;
	        for ( l = k ; l < ifirst + idata->norbitals ; ++l )
                {
	            ll = l * ( l + 1 ) / 2;
	            kl = ( k - ifirst ) + ( ( l - ifirst ) * ( l - ifirst + 1 ) ) / 2 ;
	            for ( m = jfirst ; m < jfirst + jdata->norbitals ; ++m )
                    {
		        bb = 1.0e+00 ;
		        for ( n = m ; n < jfirst + jdata->norbitals ; ++n )
                        {
		            mn = ( m - jfirst ) + ( ( n - jfirst ) * ( n - jfirst + 1 ) ) / 2 ;
		            mk = m + kk ;
		            nk = n + kk ;
		            ml = m + ll ;
		            nl = n + ll ;
                            f  = dtotal->data[mk] * dtotal->data[nl] + dtotal->data[nk] * dtotal->data[ml] ;
                            if ( dspin != NULL ) f += ( dspin ->data[mk] * dspin ->data[nl] + dspin ->data[nk] * dspin ->data[ml] ) ;
                            Real2DArray_Item ( densityij, kl, mn ) = 0.25e+00 * aa * bb * f ;
		            bb = 2.0e+00 ;
		        }
	            }
	            aa = 2.0e+00 ;
	        }
            }

            /*. Coulomb. */
            for ( i = 0 ; i < ni ; i++ )
            {
                f = Real1DArray_Item ( densityi, i ) ;
                for ( j = 0 ; j < nj ; j++ ) Real2DArray_Item ( densityij, i, j ) -= ( f * Real1DArray_Item ( densityj, j ) ) ;
            }
# ifdef MNDOCI
        }
        /* . CI case. */
        else
        {
            auto Integer          nActive, nInactive, kl0, mn0, p, q, r, s ;
            auto Real             f1, f2, f3, f4 ;
            auto MNDOCIState     *cistate ;
            auto SymmetricMatrix *onePDM, *pCore, *pHF, *zMatrix ;

            /* . Aliases. */
            cistate   = self->cistate ;
            nActive   = cistate->nactive   ;
            nInactive = cistate->ncore     ;
            onePDM    = cistate->onepdm    ;
            pCore     = cistate->pcore     ;
            pHF       = cistate->onepdmHF  ;
            zMatrix   = cistate->zMatrix   ;

            /* . Loop over basis functions. */
            for ( k = ifirst ; k < ifirst + idata->norbitals ; ++k )
            {
	        aa = 1.0e+00 ;
	        kk = k * ( k + 1 ) / 2 ;

                /* . First orbital transformation. */
                for ( s = 0 ; s < nActive ; s++ )
                {
                    for ( r = 0 ; r < nActive ; r++ )
                    {
                        for ( q = 0 ; q < nActive ; q++ )
                        {
                            for ( f = 0.0e+00, p = 0 ; p < nActive ; p++ ) f += ( Real2DArray_Item ( cistate->orbitals, k, p+nInactive ) * DoubleSymmetricMatrix_GetItem ( cistate->twopdm, p, q, r, s, NULL ) ) ;
                            RealNDArray_Item3D ( cistate->tpdm3, q, r, s ) = f ;
                        }
                    }
                }

	        for ( l = k ; l < ifirst + idata->norbitals ; ++l )
                {
	            ll = l * ( l + 1 ) / 2;

                    /* . Second orbital transformation. */
                    for ( s = 0 ; s < nActive ; s++ )
                    {
                        for ( r = 0 ; r < nActive ; r++ )
                        {
                            for ( f = 0.0e+00, q = 0 ; q < nActive ; q++ ) f += ( Real2DArray_Item ( cistate->orbitals, l, q+nInactive ) * RealNDArray_Item3D ( cistate->tpdm3, q, r, s ) ) ;
                            Real2DArray_Item ( cistate->tpdm2, r, s ) = f ;
                        }
                    }

	            kl0 = ( k - ifirst ) + ( ( l - ifirst ) * ( l - ifirst + 1 ) ) / 2 ;
                    kl  = k + ll ;
	            for ( m = jfirst ; m < jfirst + jdata->norbitals ; ++m )
                    {
		        bb = 1.0e+00 ;

                        /* . Third orbital transformation. */
                        for ( s = 0 ; s < nActive ; s++ )
                        {
                            for ( f = 0.0e+00, r = 0 ; r < nActive ; r++ ) f += ( Real2DArray_Item ( cistate->orbitals, m, r+nInactive ) * Real2DArray_Item ( cistate->tpdm2, r, s ) ) ;
                            Real1DArray_Item ( cistate->tpdm1, s ) = f ;
                        }

		        for ( n = m ; n < jfirst + jdata->norbitals ; ++n )
                        {
		            mn0 = ( m - jfirst ) + ( ( n - jfirst ) * ( n - jfirst + 1 ) ) / 2 ;
                            mn  = m + ( n * ( n + 1 ) ) / 2 ;
		            mk  = m + kk ;
		            nk  = n + kk ;
		            ml  = m + ll ;
		            nl  = n + ll ;
                            /* . PCore/PCore term. */
                            f1 = ( pCore->data[kl] * pCore->data[mn] ) - 0.25e+00 * ( pCore->data[mk] * pCore->data[nl] + pCore->data[nk] * pCore->data[ml] ) ;
                            /* . OnePDM/PCore term. */
                            f2 = 0.5e+00   * ( onePDM->data[kl] * pCore ->data[mn] + pCore ->data[kl] * onePDM->data[mn] ) -
                                 0.125e+00 * ( onePDM->data[mk] * pCore ->data[nl] + onePDM->data[nk] * pCore ->data[ml] +
                                               pCore ->data[mk] * onePDM->data[nl] + pCore ->data[nk] * onePDM->data[ml] ) ;
                            /* . TwoPDM term. */
                            for ( f3 = 0.0e+00, s = 0 ; s < nActive ; s++ ) f3 += ( Real2DArray_Item ( cistate->orbitals, n, s+nInactive ) * Real1DArray_Item ( cistate->tpdm1, s ) ) ;
                            /* . Zmatrix term. */
                            f4 = 0.5e+00   * ( pHF    ->data[kl] * zMatrix->data[mn] + zMatrix->data[kl] * pHF    ->data[mn] ) -
                                 0.125e+00 * ( pHF    ->data[mk] * zMatrix->data[nl] + pHF    ->data[nk] * zMatrix->data[ml] +
                                               zMatrix->data[mk] * pHF    ->data[nl] + zMatrix->data[nk] * pHF    ->data[ml] ) ;
# ifdef DEBUGMNDOCIGRADIENTS
printf ( "TEI> %5d %5d %5d %5d %20.10f %20.10f %20.10f %20.10f\n", k, l, m, n, f1, f2, f3, f4 ) ;
# endif
                            /* . Total contribution. */
                            Real2DArray_Item ( densityij, kl0, mn0 ) = - aa * bb * ( f1 + 2.0e+00 * ( f2 + f3 + f4 ) ) ;
		            bb = 2.0e+00 ;
		        }
	            }
	            aa = 2.0e+00 ;
	        }
            }
        }
# endif
    }
# ifdef MNDOCI
# ifdef DEBUGMNDOCIGRADIENTS
{
auto Integer i ;
printf ( "\nDensity I\n"  ) ;
for ( i = 0 ; i < densityi->length0 ; i++ ) printf ( "%20.10f\n", densityi->data[i] ) ;
printf ( "\nDensity J\n" ) ;
for ( i = 0 ; i < densityj->length0 ; i++ ) printf ( "%20.10f\n", densityj->data[i] ) ;
printf ( "\nDensity IJ\n" ) ;
Matrix_Print ( densityij ) ;
}
# endif
# endif

    /* . Finish up. */
    (*dOneI)  = densityi  ;
    (*dOneJ)  = densityj  ;
    (*dTwoIJ) = densityij ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization for an energy/gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDOState_Initialize ( QCModelMNDOState  *self, Coordinates3 *coordinates3, QCMMInteractionState *qcmmstate, Coordinates3 *gradients3 )
{
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
# ifdef MNDOCI
        auto Boolean isOK ;
# endif
        /* . Scalars. */
        self->isConverged   = False   ;
        self->eelectronic  = 0.0e+00 ;
        self->enuclear     = 0.0e+00 ;
        self->eocc         = 0.0e+00 ;
        self->eoei         = 0.0e+00 ;
        self->eqcmm        = 0.0e+00 ;
        self->eqcqc        = 0.0e+00 ;
        self->etei         = 0.0e+00 ;
        self->ncycles      = -1      ;
        /* . Aliases - temporary. */
        self->coordinates3 = coordinates3 ;
        self->gradients3   = gradients3   ;
        self->qcmmstate    = qcmmstate    ;
        /* . Deallocate everything. */
        QCModelMNDOState_Finalize ( self, True ) ;
        /* . Allocated data - temporary. */
        QCAtomContainer_GetCoordinates3 ( self->qcAtoms, coordinates3, True, &(self->qccoordinates3) ) ;
        self->oneelectronmatrix    = NULL ;
        self->twoelectronintegrals = NULL ;
        /* . Allocate orbital data. */
        QCOnePDM_AllocateOrbitalData ( self->densityp, NULL, NULL ) ;
        QCOnePDM_AllocateOrbitalData ( self->densityq, NULL, NULL ) ;
# ifdef MNDOCI
        /* . CI. */
        isOK = MNDOCIState_Initialize ( self->cistate, ( gradients3 != NULL ), self->densityp ) ;
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make densities.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real QCModelMNDOState_MakeDensities ( QCModelMNDOState *self )
{
    auto Real rmsDifference = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Real rmsDifferenceA = 0.0e+00, rmsDifferenceB = 0.0e+00 ;
        rmsDifferenceA = QCOnePDM_MakeFromFock ( self->densityp, NULL, NULL ) ;
        rmsDifferenceB = QCOnePDM_MakeFromFock ( self->densityq, NULL, NULL ) ;
        rmsDifference = Maximum ( rmsDifferenceA, rmsDifferenceB ) ;
    }
    return rmsDifference ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the orbital transformations - should switch to a sparse representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDOState_MakeOrbitalTransformations ( QCModelMNDOState *self )
{
    if ( ( self != NULL ) && ( self->qcAtoms != NULL ) )
    {
        auto Integer        i, ic, iqm, istart, istartw, j ;
        auto GaussianBasis *ibasis ;
        auto Real2DArray   *m      ;

        /* . Free space. */
        Real2DArray_Deallocate ( &(self->c2o) ) ;

        /* . Reallocate space. */
        self->c2o = Real2DArray_Allocate ( self->qcAtoms->nobasisw, self->qcAtoms->nobasis, NULL ) ; Real2DArray_Set ( self->c2o, 0.0e+00 ) ;

        /* . Loop over the orbital bases for each atom. */
        for ( iqm = 0 ; iqm < self->qcAtoms->natoms ; iqm++ )
        {
            ic      = self->qcAtoms->data[iqm].center ;
            ibasis  = self->qcParameters->centers[ic].orbitalbasis ;
            istart  = self->qcAtoms->data[iqm].ostart  ;
            istartw = self->qcAtoms->data[iqm].ostartw ;
            m = ibasis->c2o ;
            for ( i = 0 ; i < m->length0 ; i++ )
            {
                for ( j = 0 ; j < m->length1 ; j++ ) Real2DArray_Item ( self->c2o, i + istartw, j + istart ) =  Real2DArray_Item ( m, i, j ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Restore the Fock matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelMNDOState_RestoreFock ( QCModelMNDOState *self )
{
    if ( self != NULL )
    {
        auto Boolean doCopy ;
        doCopy = ( self->fockA != NULL ) && ( self->densityp != NULL ) && ( self->densityp->fock != NULL ) ;
        if ( doCopy ) SymmetricMatrix_CopyTo ( self->fockA, self->densityp->fock ) ;
        doCopy = ( self->fockB != NULL ) && ( self->densityq != NULL ) && ( self->densityq->fock != NULL ) ;
        if ( doCopy ) SymmetricMatrix_CopyTo ( self->fockB, self->densityq->fock ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setup the state.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCModelMNDOState *QCModelMNDOState_Setup ( QCAtomContainer *qcAtoms, QCParameter *qcParameters, const Real alphaCharge, const Real betaCharge, const QCOnePDMOccupancyType occupancyType, const Boolean isSpinRestricted,
                                                                                                                                                                   const Integer numberFractionalHOOs     , const Integer numberFractionalLUOs     ,
                                                                                                                                                                   const Integer numberFractionalAlphaHOOs, const Integer numberFractionalAlphaLUOs,
                                                                                                                                                                   const Integer numberFractionalBetaHOOs , const Integer numberFractionalBetaLUOs ,
                                                                                                                                                                                                              const Real fermiBroadening )
{
    QCModelMNDOState *self = NULL ;
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) )
    {
        auto Boolean isOK = True ;
        auto Integer length ;
        length = qcAtoms->nobasis ;
        /* . Allocate space. */
        self = QCModelMNDOState_Allocate ( ) ;
        if ( self != NULL )
        {
            /* . Set persistent aliases. */
            self->qcAtoms      = qcAtoms      ;
            self->qcParameters = qcParameters ;
            /* . Allocate densities - simple diagonal guess. */
            if ( isSpinRestricted )
            {
                self->densityp = QCOnePDM_FromDiagonalGuess ( QCOnePDMDensityType_Total, occupancyType, length, alphaCharge + betaCharge, numberFractionalHOOs, numberFractionalLUOs, NULL ) ;
                if ( self->densityp == NULL ) isOK = False ;
                else                          self->densityp->fermiBroadening = fermiBroadening ;
            }
            else
            {
                self->densityp = QCOnePDM_FromDiagonalGuess ( QCOnePDMDensityType_Alpha, occupancyType, length, alphaCharge, numberFractionalAlphaHOOs, numberFractionalAlphaLUOs, NULL ) ;
                if ( self->densityp == NULL ) isOK = False ;
                else                          self->densityp->fermiBroadening = fermiBroadening ;
                self->densityq = QCOnePDM_FromDiagonalGuess ( QCOnePDMDensityType_Beta , occupancyType, length, betaCharge , numberFractionalBetaHOOs , numberFractionalBetaLUOs , NULL ) ;
                if ( self->densityq == NULL ) isOK = False ;
                else                          self->densityq->fermiBroadening = fermiBroadening ;
            }
            /* . Make the full orbital transformations. */
            QCModelMNDOState_MakeOrbitalTransformations ( self ) ;
        }
        if ( ! isOK ) QCModelMNDOState_Deallocate ( &self ) ;
    }
    return self ;
}
