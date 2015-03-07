/*------------------------------------------------------------------------------
! . File      : QCModelDFTState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module implements DFT QC model state procedures.
!=================================================================================================================================*/

# include "ExecutionEnvironment.h"
# include "Memory.h"
# include "QCModelDFTState.h"

# ifdef USEOPENMP
/*----------------------------------------------------------------------------------------------------------------------------------
! . Accumulation for a Fock matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFTState_AccumulateFock ( QCModelDFTState *self, SymmetricMatrix **fockArray )
{
    if ( ( self->numberOfThreads > 1 ) && ( fockArray != NULL ) )
    {
	auto Integer t ;
        auto SymmetricMatrix *fock = fockArray[0] ;
        for ( t = 1 ; t < self->numberOfThreads ; t++ ) SymmetricMatrix_AddScaledMatrix ( fock, 1.0e+00, fockArray[t] ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Accumulation for gradients - could certainly be better here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFTState_AccumulateGradients3 ( QCModelDFTState *self )
{
    if ( ( self->numberOfThreads > 1 ) && ( self->qcGradients3Array != NULL ) )
    {
        #pragma omp parallel
        {
	    auto Integer       r, t ;
            auto Real          gX, gY, gZ ;
            auto Coordinates3 *qcGradients3 = self->qcGradients3Array[0] ;
            #pragma omp for schedule(dynamic)
            for ( r = 0 ; r < Coordinates3_Rows ( qcGradients3 ) ; r++ )
            {
                for ( t = 1 ; t < self->numberOfThreads ; t++ )
                {
                    Coordinates3_GetRow       ( self->qcGradients3Array[t] , r, gX, gY, gZ ) ;
                    Coordinates3_IncrementRow ( qcGradients3               , r, gX, gY, gZ ) ;
                }
            }
        }
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCModelDFTState *QCModelDFTState_Allocate ( void )
{
    QCModelDFTState *self = NULL ;
    /* . Allocation. */
    self = ( QCModelDFTState * ) Memory_Allocate ( sizeof ( QCModelDFTState ) ) ;
    /* . Initialization. */
    if ( self != NULL )
    {
        /* . Scalars. */
        self->isConverged = False   ;
        self->eelectronic = 0.0e+00 ;
        self->enuclear    = 0.0e+00 ;
        self->eocc        = 0.0e+00 ;
        self->eoei        = 0.0e+00 ;
        self->eqcmm       = 0.0e+00 ;
        self->eqcqc       = 0.0e+00 ;
        self->equad       = 0.0e+00 ;
        self->etei        = 0.0e+00 ;
        self->rhoquad     = 0.0e+00 ;
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
        self->o2c                  = NULL ;
        self->densityp             = NULL ;
        self->densityq             = NULL ;
        /* . Allocated data - temporary. */
        self->fitintegrals         = NULL ;
        self->qccoordinates3       = NULL ;
        self->dftgrid              = NULL ;
        self->lowdintransformation = NULL ;
        self->orthogonalizer       = NULL ;
        self->overlapeigenvectors  = NULL ;
        self->inversefitmatrix     = NULL ;
        self->oneelectronmatrix    = NULL ;
        self->overlap              = NULL ;
        self->wdensity             = NULL ;
        self->fpotential           = NULL ;
        self->fselfoverlap         = NULL ;
        self->overlapeigenvalues   = NULL ;
        self->wvector              = NULL ;
# ifdef USEOPENMP
        self->numberOfThreads      =   -1 ;
        self->fockAArray           = NULL ;
        self->fockBArray           = NULL ;
        self->qcGradients3Array    = NULL ;
# endif
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFTState_Deallocate ( QCModelDFTState **self )
{
    if ( (*self) != NULL )
    {
        auto QCModelDFTState *localSelf = (*self) ;
        Real2DArray_Deallocate     ( &(localSelf->c2o                 ) ) ;
        Real2DArray_Deallocate     ( &(localSelf->o2c                 ) ) ;
        QCOnePDM_Deallocate        ( &(localSelf->densityp            ) ) ;
        QCOnePDM_Deallocate        ( &(localSelf->densityq            ) ) ;
        BlockStorage_Deallocate    ( &(localSelf->fitintegrals        ) ) ;
        Coordinates3_Deallocate    ( &(localSelf->qccoordinates3      ) ) ;
        DFTGrid_Deallocate         ( &(localSelf->dftgrid             ) ) ;
        Real2DArray_Deallocate     ( &(localSelf->lowdintransformation) ) ;
        Real2DArray_Deallocate     ( &(localSelf->overlapeigenvectors ) ) ;
        Real2DArray_Deallocate     ( &(localSelf->orthogonalizer      ) ) ;
        SymmetricMatrix_Deallocate ( &(localSelf->inversefitmatrix    ) ) ;
        SymmetricMatrix_Deallocate ( &(localSelf->oneelectronmatrix   ) ) ;
        SymmetricMatrix_Deallocate ( &(localSelf->overlap             ) ) ;
        SymmetricMatrix_Deallocate ( &(localSelf->wdensity            ) ) ;
        Real1DArray_Deallocate     ( &(localSelf->fpotential          ) ) ;
        Real1DArray_Deallocate     ( &(localSelf->fselfoverlap        ) ) ;
        Real1DArray_Deallocate     ( &(localSelf->overlapeigenvalues  ) ) ;
        Real1DArray_Deallocate     ( &(localSelf->wvector             ) ) ;
# ifdef USEOPENMP
        /* . Thread data. */
        if ( localSelf->numberOfThreads > 0 )
        {
            auto Integer t ;
            if ( localSelf->fockAArray != NULL )
            {
                for ( t = 1 ; t < localSelf->numberOfThreads ; t++ ) SymmetricMatrix_Deallocate ( &(localSelf->fockAArray[t]) ) ;
                MEMORY_DEALLOCATE ( localSelf->fockAArray ) ;
            }
            if ( localSelf->fockBArray != NULL )
            {
                for ( t = 1 ; t < localSelf->numberOfThreads ; t++ ) SymmetricMatrix_Deallocate ( &(localSelf->fockBArray[t]) ) ;
                MEMORY_DEALLOCATE ( localSelf->fockBArray ) ;
            }
            if ( localSelf->qcGradients3Array != NULL )
            {
                for ( t = 0 ; t < localSelf->numberOfThreads ; t++ ) Coordinates3_Deallocate ( &(localSelf->qcGradients3Array[t]) ) ;
                MEMORY_DEALLOCATE ( localSelf->qcGradients3Array ) ;
            }
        }
# endif
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Finalization for an energy/gradient calculation.
! . Keep the overlap and orthogonalizer temporarily.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFTState_Finalize ( QCModelDFTState *self, const Boolean QKEEPDATA )
{
    if ( self != NULL )
    {
        if ( ! QKEEPDATA )
        {
            /* . Integrals, etc. */
            BlockStorage_Deallocate        ( &(self->fitintegrals     ) ) ;
            Coordinates3_Deallocate        ( &(self->qccoordinates3   ) ) ;
            DFTGrid_Deallocate             ( &(self->dftgrid          ) ) ;
            SymmetricMatrix_Deallocate     ( &(self->inversefitmatrix ) ) ;
            SymmetricMatrix_Deallocate     ( &(self->oneelectronmatrix) ) ;
            SymmetricMatrix_Deallocate     ( &(self->wdensity         ) ) ;
/*
            Real1DArray_Deallocate          ( &(self->fpotential       ) ) ;
            Real1DArray_Deallocate          ( &(self->fselfoverlap     ) ) ;
*/
            Real1DArray_Deallocate         ( &(self->wvector          ) ) ;
            /* . Orbital data. */
            QCOnePDM_DeallocateOrbitalData ( self->densityp ) ;
            QCOnePDM_DeallocateOrbitalData ( self->densityq ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization for an energy/gradient calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFTState_Initialize ( QCModelDFTState *self, Coordinates3 *coordinates3, QCMMInteractionState *qcmmstate, Coordinates3 *gradients3 )
{
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
        auto Integer ndim ;
        /* . Scalars. */
        self->isConverged  = False   ;
        self->eelectronic  = 0.0e+00 ;
        self->enuclear     = 0.0e+00 ;
        self->eocc         = 0.0e+00 ;
        self->eoei         = 0.0e+00 ;
        self->eqcmm        = 0.0e+00 ;
        self->eqcqc        = 0.0e+00 ;
        self->equad        = 0.0e+00 ;
        self->etei         = 0.0e+00 ;
        self->rhoquad      = 0.0e+00 ;
        self->ncycles      = -1      ;
        /* . Aliases - temporary. */
        self->coordinates3 = coordinates3 ;
        self->gradients3   = gradients3   ;
        self->qcmmstate    = qcmmstate    ;
        /* . Deallocate various items if they still exist. */
        QCModelDFTState_Finalize   ( self, False ) ;
        Real2DArray_Deallocate     ( &(self->lowdintransformation) ) ;
        Real2DArray_Deallocate     ( &(self->overlapeigenvectors ) ) ;
        Real2DArray_Deallocate     ( &(self->orthogonalizer      ) ) ;
        SymmetricMatrix_Deallocate ( &(self->overlap             ) ) ;
        Real1DArray_Deallocate     ( &(self->fpotential          ) ) ;
        Real1DArray_Deallocate     ( &(self->fselfoverlap        ) ) ;
        Real1DArray_Deallocate     ( &(self->overlapeigenvalues  ) ) ;
        /* . Allocated data - temporary. */
        QCAtomContainer_GetCoordinates3 ( self->qcAtoms, coordinates3, True, &(self->qccoordinates3) ) ;
        /* . Allocate some more space - (nfbasis + 1) because the total charge is always constrained. */
        ndim = self->qcAtoms->nfbasisw + 1 ;
        self->fpotential = Real1DArray_Allocate ( ndim, NULL ) ; Real1DArray_Set ( self->fpotential, 0.0e+00 ) ;
        self->wvector    = Real1DArray_Allocate ( ndim, NULL ) ; Real1DArray_Set ( self->wvector,    0.0e+00 ) ;
# ifdef USEOPENMP
        /* . Thread data. */
        if ( ( self->numberOfThreads > 0 ) && ( self->qcGradients3Array != NULL ) )
        {
            auto Integer t ;
            for ( t = 1 ; t < self->numberOfThreads ; t++ ) Coordinates3_Set ( self->qcGradients3Array[t] , 0.0e+00 ) ;
        }
# endif
    }
}

# ifdef USEOPENMP
/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialize a Fock array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFTState_InitializeFockArray ( QCModelDFTState *self, QCOnePDM *qcOnePDM, SymmetricMatrix **fockArray )
{
    if ( ( self != NULL ) && ( qcOnePDM != NULL ) && ( fockArray != NULL ) )
    {
        auto Integer t ;
        fockArray[0] = qcOnePDM->fock ;
        for ( t = 1 ; t < self->numberOfThreads ; t++ ) SymmetricMatrix_Set ( fockArray[t], 0.0e+00 ) ;
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the orbital transformations - should switch to a sparse representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCModelDFTState_MakeOrbitalTransformations ( QCModelDFTState *self )
{
    if ( ( self != NULL ) && ( self->qcAtoms != NULL ) )
    {
        auto Integer i, ic, iqm, istart, istartw, j ;
        auto GaussianBasis *ibasis ;
        auto Real2DArray *m ;

        /* . Free space. */
        Real2DArray_Deallocate ( &(self->c2o) ) ;
        Real2DArray_Deallocate ( &(self->o2c) ) ;

        /* . Reallocate space. */
        self->c2o = Real2DArray_Allocate ( self->qcAtoms->nobasisw, self->qcAtoms->nobasis, NULL ) ; Real2DArray_Set ( self->c2o, 0.0e+00 ) ;
        self->o2c = Real2DArray_Allocate ( self->qcAtoms->nobasisw, self->qcAtoms->nobasis, NULL ) ; Real2DArray_Set ( self->o2c, 0.0e+00 ) ;

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
            m = ibasis->o2c ;
            for ( i = 0 ; i < m->length0 ; i++ )
            {
                for ( j = 0 ; j < m->length1 ; j++ ) Real2DArray_Item ( self->o2c, i + istartw, j + istart ) =  Real2DArray_Item ( m, i, j ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setup the state.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCModelDFTState *QCModelDFTState_Setup (       QCAtomContainer      *qcAtoms                   ,
                                               QCParameter          *qcParameters              ,
                                         const Real                  alphaCharge               ,
                                         const Real                  betaCharge                ,
                                         const QCOnePDMOccupancyType occupancyType             ,
                                         const Boolean               isSpinRestricted          ,
                                         const Integer               numberFractionalHOOs      ,
                                         const Integer               numberFractionalLUOs      ,
                                         const Integer               numberFractionalAlphaHOOs ,
                                         const Integer               numberFractionalAlphaLUOs ,
                                         const Integer               numberFractionalBetaHOOs  ,
                                         const Integer               numberFractionalBetaLUOs  ,
                                         const Real                  fermiBroadening           )
{
    QCModelDFTState *self = NULL ;
    if ( ( qcAtoms != NULL ) && ( qcParameters != NULL ) )
    {
        auto Boolean isOK = True ;
        auto Integer length ;
        length = qcAtoms->nobasisw ; /* . Could there be a problem here - maybe need both nobasis and nobasisw? */
        /* . Allocate space. */
        self = QCModelDFTState_Allocate ( ) ;
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
            QCModelDFTState_MakeOrbitalTransformations ( self ) ;
# ifdef USEOPENMP
            /* . Thread data. */
            {
                auto Integer n = self->qcAtoms->natoms, t ;
                self->numberOfThreads = omp_get_max_threads ( ) ;
                /* . Fock matrices. */
                MEMORY_ALLOCATEARRAY ( self->fockAArray, self->numberOfThreads, SymmetricMatrix * ) ;
                if ( self->fockAArray == NULL ) isOK = False ;
                else
                {
                    self->fockAArray[0] = NULL ;
                    for ( t = 1 ; t < self->numberOfThreads ; t++ )
                    {
                        self->fockAArray[t] = SymmetricMatrix_Allocate ( length ) ;
                        isOK = isOK && ( self->fockAArray[t] != NULL ) ;
                    }
                }
                if ( ! isSpinRestricted )
                {
                    MEMORY_ALLOCATEARRAY ( self->fockBArray, self->numberOfThreads, SymmetricMatrix * ) ;
                    if ( self->fockBArray == NULL ) isOK = False ;
                    else
                    {
                        self->fockBArray[0] = NULL ;
                        for ( t = 1 ; t < self->numberOfThreads ; t++ )
                        {
                            self->fockBArray[t] = SymmetricMatrix_Allocate ( length ) ;
                            isOK = isOK && ( self->fockBArray[t] != NULL ) ;
                        }
                    }
                }
                /* . Gradients. */
                MEMORY_ALLOCATEARRAY ( self->qcGradients3Array, self->numberOfThreads, Coordinates3 * ) ;
                if ( self->qcGradients3Array == NULL ) isOK = False ;
                else
                {
                    self->qcGradients3Array[0] = NULL ;
                    for ( t = 1 ; t < self->numberOfThreads ; t++ )
                    {
                        self->qcGradients3Array[t] = Coordinates3_Allocate ( n ) ;
                        isOK = isOK && ( self->qcGradients3Array[t] != NULL ) ;
                    }
                }
            }
# endif
        }
        else isOK = False ;
        if ( ! isOK ) QCModelDFTState_Deallocate ( &self ) ;
    }
    return self ;
}
