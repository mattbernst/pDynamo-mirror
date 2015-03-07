/*------------------------------------------------------------------------------
! . File      : PairwiseInteraction.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module contains various procedures for generating or manipulating different forms of pairwise interaction.
!=================================================================================================================================*/

/* # define SPLINEPRINTING */
# ifdef SPLINEPRINTING
# include <stdio.h>
# endif

# include "ExecutionEnvironment.h"
# include "Memory.h"
# include "PairwiseInteraction.h"
# include "Units.h"

/* . Need to add "useAnalyticForm" option to the QCMM and QCQC procedures. */

/*==================================================================================================================================
! . Basic procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionABFS *PairwiseInteractionABFS_Allocate ( void )
{
    PairwiseInteractionABFS *self = NULL ;
    self = ( PairwiseInteractionABFS * ) Memory_Allocate ( sizeof ( PairwiseInteractionABFS ) ) ;
    if ( self != NULL )
    {
        self->useAnalyticForm       = True ;
        self->splinePointDensity    = 50   ;
        self->dampingCutoff         =  0.5e+00 ;
        self->innerCutoff           =  8.0e+00 ;
        self->outerCutoff           = 12.0e+00 ;
        self->electrostaticSpline   = NULL ;
        self->lennardJonesASpline   = NULL ;
        self->lennardJonesBSpline   = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionABFS *PairwiseInteractionABFS_Clone ( const PairwiseInteractionABFS *self, Status *status )
{
    PairwiseInteractionABFS *new = NULL ;
    if ( self != NULL )
    {
        new = PairwiseInteractionABFS_Allocate ( ) ;
        if ( new != NULL )
        {
            auto Status localStatus ;
            new->useAnalyticForm       = self->useAnalyticForm       ;
            new->splinePointDensity    = self->splinePointDensity    ;
            new->dampingCutoff         = self->dampingCutoff         ;
            new->innerCutoff           = self->innerCutoff           ;
            new->outerCutoff           = self->outerCutoff           ;
            localStatus = CubicSpline_Clone ( self->electrostaticSpline, &(new->electrostaticSpline) ) ;
            localStatus = CubicSpline_Clone ( self->lennardJonesASpline, &(new->lennardJonesASpline) ) ;
            localStatus = CubicSpline_Clone ( self->lennardJonesBSpline, &(new->lennardJonesBSpline) ) ;
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_Deallocate ( PairwiseInteractionABFS **self )
{
    if ( (*self) != NULL )
    {
        CubicSpline_Deallocate ( &((*self)->electrostaticSpline) ) ;
        CubicSpline_Deallocate ( &((*self)->lennardJonesASpline) ) ;
        CubicSpline_Deallocate ( &((*self)->lennardJonesBSpline) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a factors object.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_MakeFactors ( const PairwiseInteractionABFS *self, PairwiseInteractionABFSFactors *factors )
{
    if ( ( self != NULL ) && ( factors != NULL ) )
    {
        auto Real f, g, gamma, s12, s6 ;

        /* . Cutoff factors. */
        factors->r2Damp   = self->dampingCutoff * self->dampingCutoff ;
        factors->r2Off    = self->outerCutoff   * self->outerCutoff   ;
        factors->r2On     = self->innerCutoff   * self->innerCutoff   ;

        /* . Electrostatic factors. */
        /* . Interaction. */
        gamma             = pow ( factors->r2Off - factors->r2On, 3 ) ;
        factors->a        = factors->r2Off * factors->r2Off * ( factors->r2Off - 3.0e+00 * factors->r2On ) / gamma ;
        factors->b        = 6.0e+00 * factors->r2Off * factors->r2On / gamma ;
        factors->c        = - ( factors->r2Off + factors->r2On )     / gamma ;
        factors->d        = 0.4e+00                                  / gamma ;
        factors->qShift1  = 8.0e+00 * ( factors->r2Off * factors->r2On * ( self->outerCutoff - self->innerCutoff ) -
                            0.2e+00 * ( self->outerCutoff * factors->r2Off * factors->r2Off - self->innerCutoff * factors->r2On * factors->r2On ) ) / gamma ;
        factors->qShift2  = - ( factors->a / self->outerCutoff ) + factors->b * self->outerCutoff + factors->c * self->outerCutoff * factors->r2Off +
                                                                                           factors->d * self->outerCutoff * factors->r2Off * factors->r2Off ;
        /* . Damping. */
        f                 =   1.0e+00 / self->dampingCutoff + factors->qShift1 ;
        g                 = - 1.0e+00 / factors->r2Damp ;
        factors->qF0      = f - 0.5e+00 * self->dampingCutoff * g ;
        factors->qAlpha   =   - 0.5e+00 * g / self->dampingCutoff ;

        /* . Lennard-Jones A factors. */
        /* . Interaction. */
        factors->aF6      = 1.0e+00 / ( factors->r2Off * factors->r2Off * factors->r2Off ) ;
        factors->aK12     = pow ( factors->r2Off, 3 ) / ( pow ( factors->r2Off, 3 ) - pow ( factors->r2On, 3 ) ) ;
        factors->aShift12 = 1.0e+00 / pow ( self->innerCutoff * self->outerCutoff, 6 ) ;
        /* . Damping. */
        s12               = 1.0e+00 / pow ( self->dampingCutoff, 12 ) ;
        f                 = s12 - factors->aShift12 ;
        g                 = - 12.0e+00 * s12 / self->dampingCutoff ;
        factors->aF0      = f - 0.5e+00 * self->dampingCutoff * g ;
        factors->aAlpha   =   - 0.5e+00 * g / self->dampingCutoff ;

        /* . Lennard-Jones B factors. */
        /* . Interaction. */
        factors->bF3      = 1.0e+00 / ( self->outerCutoff * factors->r2Off ) ;
        factors->bK6      = ( self->outerCutoff * factors->r2Off ) / ( self->outerCutoff * factors->r2Off - self->innerCutoff * factors->r2On ) ;
        factors->bShift6  = 1.0e+00 / pow ( self->innerCutoff * self->outerCutoff, 3 ) ;
        /* . Damping. */
        s6                = 1.0e+00 / pow ( self->dampingCutoff, 6 ) ;
        f                 = - s6 + factors->bShift6 ;
        g                 = 6.0e+00 * s6 / self->dampingCutoff ;
        factors->bF0      = f - 0.5e+00 * self->dampingCutoff * g ;
        factors->bAlpha   =   - 0.5e+00 * g / self->dampingCutoff ;
    }
}

/*==================================================================================================================================
! . Spline procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a spline of the ABFS electrostatic interaction.
!---------------------------------------------------------------------------------------------------------------------------------*/
CubicSpline *PairwiseInteractionABFS_MakeElectrostaticSpline ( const PairwiseInteractionABFS *self, const Boolean useAtomicUnits, Status *status )
{
    CubicSpline *spline = NULL ;
    if ( self != NULL )
    {
        auto Integer numberOfPoints ;
        auto Real1DArray *x = NULL, *y = NULL ;
        /* . Allocate space. */
        PairwiseInteractionABFS_NumberOfSplinePoints ( self, numberOfPoints ) ;
        x = Real1DArray_Allocate ( numberOfPoints, status ) ;
        y = Real1DArray_Allocate ( numberOfPoints, status ) ;
        if ( ( x != NULL ) && ( y != NULL ) )
        {
            auto Integer i ;
            auto Real    dF, dR, f, r, r2, s, s2 ;
            auto PairwiseInteractionABFSFactors factors ;
            /* . Make factors. */
            PairwiseInteractionABFS_MakeFactors ( self, &factors ) ;
            /* . Make x and y. */
            dR = self->outerCutoff / ( Real ) ( numberOfPoints - 1 ) ;
            for ( i = 0 ; i < numberOfPoints - 1 ; i++ )
            {
                /* . x. */
                r  = dR * ( Real ) i ;
                r2 = r * r ;
                Real1DArray_Item ( x, i ) = r2 ;
                /* . y. */
                PairwiseInteractionABFS_CheckDistances ( factors, r2, s, s2 ) ;
                PairwiseInteractionABFS_ElectrostaticTerm ( factors, r2, s, 1.0e+00, f, dF ) ;
                Real1DArray_Item ( y, i ) = f ;
            }
            Real1DArray_Item ( x, numberOfPoints - 1 ) = factors.r2Off ;
            Real1DArray_Item ( y, numberOfPoints - 1 ) = 0.0e+00 ;
# ifdef SPLINEPRINTING
printf ( "\n\nElectroststatic Spline:\n" ) ;
for ( i = 0 ; i < numberOfPoints ; i++ ) printf ( "%30.10f %30.10f\n", sqrt ( Real1DArray_Item ( x, i ) ), Real1DArray_Item ( y, i ) ) ;
# endif
            /* . Ensure y has the correct units. */
            if ( useAtomicUnits ) Real1DArray_Scale ( y, 1.0e+00 / UNITS_LENGTH_ANGSTROMS_TO_BOHRS       ) ;
            else                  Real1DArray_Scale ( y, UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE ) ;
            /* . Make the spline. */
            CubicSpline_MakeFromReal1DArrays ( &spline, &x, &y, 1, 0.0e+00, 1, 0.0e+00 ) ;
        }
    }
    return spline ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a spline of the ABFS Lennard-Jones A interaction.
!---------------------------------------------------------------------------------------------------------------------------------*/
CubicSpline *PairwiseInteractionABFS_MakeLennardJonesASpline ( const PairwiseInteractionABFS *self, Status *status )
{
    CubicSpline *spline = NULL ;
    if ( self != NULL )
    {
        auto Integer numberOfPoints ;
        auto Real1DArray *x = NULL, *y = NULL ;
        /* . Allocate space. */
        PairwiseInteractionABFS_NumberOfSplinePoints ( self, numberOfPoints ) ;
        x = Real1DArray_Allocate ( numberOfPoints, status ) ;
        y = Real1DArray_Allocate ( numberOfPoints, status ) ;
        if ( ( x != NULL ) && ( y != NULL ) )
        {
            auto Integer i ;
            auto Real    dF, dR, f, r, r2, s, s2 ;
            auto PairwiseInteractionABFSFactors factors ;
            /* . Make factors. */
            PairwiseInteractionABFS_MakeFactors ( self, &factors ) ;
            /* . Make x and y. */
            dR = self->outerCutoff / ( Real ) ( numberOfPoints - 1 ) ;
            for ( i = 0 ; i < numberOfPoints - 1 ; i++ )
            {
                /* . x. */
                r  = dR * ( Real ) i ;
                r2 = r * r ;
                Real1DArray_Item ( x, i ) = r2 ;
                /* . y. */
                PairwiseInteractionABFS_CheckDistances ( factors, r2, s, s2 ) ;
                PairwiseInteractionABFS_LennardJonesATerm ( factors, r2, s, s2, 1.0e+00, f, dF ) \
                Real1DArray_Item ( y, i ) = f ;
            }
            Real1DArray_Item ( x, numberOfPoints - 1 ) = factors.r2Off ;
            Real1DArray_Item ( y, numberOfPoints - 1 ) = 0.0e+00 ;
            /* . Make the spline. */
            CubicSpline_MakeFromReal1DArrays ( &spline, &x, &y, 1, 0.0e+00, 1, 0.0e+00 ) ;
        }
    }
    return spline ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a spline of the ABFS Lennard-Jones B interaction.
!---------------------------------------------------------------------------------------------------------------------------------*/
CubicSpline *PairwiseInteractionABFS_MakeLennardJonesBSpline ( const PairwiseInteractionABFS *self, Status *status )
{
    CubicSpline *spline = NULL ;
    if ( self != NULL )
    {
        auto Integer numberOfPoints ;
        auto Real1DArray *x = NULL, *y = NULL ;
        /* . Allocate space. */
        PairwiseInteractionABFS_NumberOfSplinePoints ( self, numberOfPoints ) ;
        x = Real1DArray_Allocate ( numberOfPoints, status ) ;
        y = Real1DArray_Allocate ( numberOfPoints, status ) ;
        if ( ( x != NULL ) && ( y != NULL ) )
        {
            auto Integer i ;
            auto Real    dF, dR, f, r, r2, s, s2 ;
            auto PairwiseInteractionABFSFactors factors ;
            /* . Make factors. */
            PairwiseInteractionABFS_MakeFactors ( self, &factors ) ;
            /* . Make x and y. */
            dR = self->outerCutoff / ( Real ) ( numberOfPoints - 1 ) ;
            for ( i = 0 ; i < numberOfPoints - 1 ; i++ )
            {
                /* . x. */
                r  = dR * ( Real ) i ;
                r2 = r * r ;
                Real1DArray_Item ( x, i ) = r2 ;
                /* . y. */
                PairwiseInteractionABFS_CheckDistances ( factors, r2, s, s2 ) ;
                PairwiseInteractionABFS_LennardJonesBTerm ( factors, r2, s, s2, 1.0e+00, f, dF ) \
                Real1DArray_Item ( y, i ) = f ;
            }
            Real1DArray_Item ( x, numberOfPoints - 1 ) = factors.r2Off ;
            Real1DArray_Item ( y, numberOfPoints - 1 ) = 0.0e+00 ;
            /* . Make the spline. */
            CubicSpline_MakeFromReal1DArrays ( &spline, &x, &y, 1, 0.0e+00, 1, 0.0e+00 ) ;
        }
    }
    return spline ;
}

/*==================================================================================================================================
! . Interaction procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . MM/MM energy and gradients given a pairlist using either analytic (default) or spline evaluation.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . Either electrostatic or Lennard-Jones or both energies can be calculated. crd1/crd2 and grd1/grd2 can refer to the
! . same or different objects. If different, the pairList must be such that i refers to crd1/grd1 and j to crd2/grd2.
*/
void PairwiseInteractionABFS_MMMMEnergy ( const PairwiseInteractionABFS *self               ,
                                          const MMAtomContainer         *mmAtoms            ,
                                          const LJParameterContainer    *ljParameters       ,
                                                PairList                *pairList           ,
                                          const Real                     electrostaticScale ,
                                          const Real                     lennardJonesScale  ,
                                          const Coordinates3            *crd1               ,
                                          const Coordinates3            *crd2               ,
                                                Real                    *eElectrostatic     ,
                                                Real                    *eLennardJones      ,
# ifdef USEOPENMP
                                          const Integer                  numberOfThreads    ,
                                                Coordinates3           **grds1              ,
                                                Coordinates3           **grds2              )
# else
                                                Coordinates3            *grd1               ,
                                                Coordinates3            *grd2               )
# endif
{
    Boolean doElectrostatics, doLennardJones ;
    doElectrostatics = ( eElectrostatic != NULL ) ;
    doLennardJones   = ( eLennardJones  != NULL ) ;
    if ( doElectrostatics ) (*eElectrostatic) = 0.0e+00 ;
    if ( doLennardJones   ) (*eLennardJones ) = 0.0e+00 ;
    if ( ( self     != NULL ) &&
         ( mmAtoms  != NULL ) && 
         ( pairList != NULL ) &&
         ( crd1     != NULL ) &&
         ( crd2     != NULL ) &&
         ( doElectrostatics || doLennardJones ) &&
         ( ( ljParameters != NULL ) || ( ! doLennardJones && ( ljParameters == NULL ) ) ) )
    {
        auto Boolean doGradients ;
        auto Real    eLJ, eQQ ;
        auto Integer numberOfLJTypes ;

        /* . Initialization. */
# ifdef USEOPENMP
        doGradients     = ( grds1 != NULL ) && ( grds2 != NULL ) ;
# else
        doGradients     = ( grd1  != NULL ) && ( grd2  != NULL ) ;
# endif
        eLJ             = 0.0e+00 ;
        eQQ             = 0.0e+00 ;
        numberOfLJTypes = ljParameters->ntypes ;
        PairList_MakeRecords ( pairList ) ;

        /* . Use the analytic form. */
        if ( self->useAnalyticForm )
        {
            auto Real eScale ;
            auto PairwiseInteractionABFSFactors abfs ;

            /* . Initialization. */
            eScale = electrostaticScale * UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE ;
            PairwiseInteractionABFS_MakeFactors ( self, &abfs ) ;

# ifdef USEOPENMP
            # pragma omp parallel num_threads(numberOfThreads) shared(abfs,eScale)
# endif
            {
                auto Integer           i, j, n, r, ti, tij ;
                auto Real              Aij, Bij, dF, eTerm, qi, qij, r2, s, s2, xi, xij, xj, yi, yij, yj, zi, zij, zj ;
                auto Coordinates3     *localGrd1, *localGrd2 ;
                auto IndexedSelection *pair ;

                if ( doGradients )
                {
# ifdef USEOPENMP
                    auto Integer thisThread = omp_get_thread_num ( ) ;
                    localGrd1  = grds1[thisThread] ;
                    localGrd2  = grds2[thisThread] ;
# else
                    localGrd1  = grd1 ;
                    localGrd2  = grd2 ;
# endif
                }

                /* . Loop over the interactions. */
# ifdef USEOPENMP
                #pragma omp for reduction(+:eLJ,eQQ) schedule(dynamic) /* . The reduction must be here. */
# endif
                for ( r = 0 ; r < pairList->numberOfRecords ; r++ )
                {
                    /* . Get the pair data. */
                    pair = pairList->records[r] ;

                    /* . First atom. */
                    i  = pair->index ;
	            qi = eScale * mmAtoms->data[i].charge ;
                    ti = numberOfLJTypes * mmAtoms->data[i].ljtype ;
                    Coordinates3_GetRow ( crd1, i, xi, yi, zi ) ;

       	            for ( n = 0 ; n < pair->nindices ; n++ )
	            {
	                /* . Second atom. */
	                j = pair->indices[n] ;

                        /* . Coordinate displacement. */
	                Coordinates3_GetRow ( crd2, j, xj, yj, zj ) ;
                        xij = xi - xj ;
                        yij = yi - yj ;
                        zij = zi - zj ;
                        r2  = ( xij * xij + yij * yij + zij * zij ) ;

                        /* . Skip interactions outside of cutoff. */
                        PairwiseInteractionABFS_CheckDistances ( abfs, r2, s, s2 ) ;

                        /* . Calculate interactions. */
                        dF = 0.0e+00 ;
                        if ( doElectrostatics )
                        {
                            qij = qi * mmAtoms->data[j].charge ;
                            PairwiseInteractionABFS_ElectrostaticTerm ( abfs, r2, s, qij, eTerm, dF ) ;
                            eQQ += eTerm ;
                        }
                        if ( doLennardJones )
                        {
	                    tij = ljParameters->tableindex[ti+mmAtoms->data[j].ljtype] ;
                            Aij = ljParameters->tableA[tij] * lennardJonesScale ;
                            Bij = ljParameters->tableB[tij] * lennardJonesScale ;
                            PairwiseInteractionABFS_LennardJonesTerm ( abfs, r2, s, s2, Aij, Bij, eTerm, dF ) ;
                            eLJ += eTerm ;
                        }

                        /* . Gradients. */
                        if ( doGradients )
                        {
                            xij *= ( 2.0e+00 * dF ) ;
                            yij *= ( 2.0e+00 * dF ) ;
                            zij *= ( 2.0e+00 * dF ) ;
                            Coordinates3_IncrementRow ( localGrd1, i, xij, yij, zij ) ;
                            Coordinates3_DecrementRow ( localGrd2, j, xij, yij, zij ) ;
                        }
	            }
                }
            }
        }
        /* . Use the spline form. */
        else
        {
            Real         r2Off ;
            CubicSpline *referenceSpline ;

            /* . Initialization. */
            r2Off = self->outerCutoff * self->outerCutoff ;
            if ( doElectrostatics ) referenceSpline = self->electrostaticSpline ;
            else                    referenceSpline = self->lennardJonesASpline ;

# ifdef USEOPENMP
            # pragma omp parallel num_threads(numberOfThreads) shared(referenceSpline,r2Off)
# endif
            {
                auto Integer           i, j, l, n, r, ti, tij, u ;
                auto Real              Aij, Bij, d, dF, dG, f, qi, qij, r2, s, t, xi, xij, xj, yi, yij, yj, zi, zij, zj ;
                auto Coordinates3     *localGrd1, *localGrd2 ;
                auto IndexedSelection *pair ;

                if ( doGradients )
                {
# ifdef USEOPENMP
                    auto Integer thisThread = omp_get_thread_num ( ) ;
                    localGrd1  = grds1[thisThread] ;
                    localGrd2  = grds2[thisThread] ;
# else
                    localGrd1  = grd1 ;
                    localGrd2  = grd2 ;
# endif
                }

                /* . Loop over the interactions. */
# ifdef USEOPENMP
                #pragma omp for reduction(+:eLJ,eQQ) schedule(dynamic) /* . The reduction must be here. */
# endif
                for ( r = 0 ; r < pairList->numberOfRecords ; r++ )
                {
                    /* . Get the pair data. */
                    pair = pairList->records[r] ;

                    /* . First atom. */
                    i  = pair->index ;
	            qi = electrostaticScale * mmAtoms->data[i].charge ;
                    ti = numberOfLJTypes    * mmAtoms->data[i].ljtype ;
                    Coordinates3_GetRow ( crd1, i, xi, yi, zi ) ;

       	            for ( n = 0 ; n < pair->nindices ; n++ )
	            {
	                /* . Second atom. */
	                j = pair->indices[n] ;

                        /* . Coordinate displacement. */
	                Coordinates3_GetRow ( crd2, j, xj, yj, zj ) ;
                        xij = xi - xj ;
                        yij = yi - yj ;
                        zij = zi - zj ;
                        r2  = ( xij * xij + yij * yij + zij * zij ) ;

                        /* . Skip interactions outside of cutoff. */
                        if ( r2 > r2Off ) continue ;

                        /* . Get table look-up quantities (assuming all have the same abscissae). */
                        CubicSpline_EvaluateLUDST ( referenceSpline, r2, &l, &u, &d, &s, &t ) ;

                        /* . Electrostatics. */
                        dG = 0.0e+00 ;
                        if ( doElectrostatics )
                        {
                            qij = qi * mmAtoms->data[j].charge ;
                            CubicSpline_FastEvaluateFG ( self->electrostaticSpline, l, u, d, s, t, f, dF ) ;
                            eQQ += ( qij * f ) ;
                            dG  += ( 2.0e+00 * qij * dF ) ;
                        }

                        /* . Lennard-Jones. */
                        if ( doLennardJones )
                        {
 	                    tij  = ljParameters->tableindex[ti+mmAtoms->data[j].ljtype] ;
                            Aij  = ljParameters->tableA[tij] * lennardJonesScale ;
                            Bij  = ljParameters->tableB[tij] * lennardJonesScale ;
                            CubicSpline_FastEvaluateFG ( self->lennardJonesASpline, l, u, d, s, t, f, dF ) ;
                            eLJ += Aij * f  ;
                            dG  += ( 2.0e+00 * Aij * dF ) ;
                            CubicSpline_FastEvaluateFG ( self->lennardJonesBSpline, l, u, d, s, t, f, dF ) ;
                            eLJ += Bij * f  ;
                            dG  += ( 2.0e+00 * Bij * dF ) ;
                        }

                        /* . Gradients. */
                        if ( doGradients )
                        {
	                    xij *= dG ;
	                    yij *= dG ;
	                    zij *= dG ;
	                    Coordinates3_IncrementRow ( localGrd1, i, xij, yij, zij ) ;
	                    Coordinates3_DecrementRow ( localGrd2, j, xij, yij, zij ) ;
                        }
	            }
                }
            }
        }
        /* . Finish up energies. */
        if ( doElectrostatics ) (*eElectrostatic) = eQQ ;
        if ( doLennardJones   ) (*eLennardJones ) = eLJ ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM electrostatic gradients in regular units.
! . crd1/grd1 refer to QC atoms and crd2/grd2 to MM atoms.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_QCMMGradients ( const PairwiseInteractionABFS *self               ,
                                             const Real1DArray             *mmCharges          ,
                                             const QCAtomContainer         *qcAtoms            ,
                                                   PairList                *pairList           ,
                                             const Real                     electrostaticScale ,
                                             const Coordinates3            *crd1               ,
                                             const Coordinates3            *crd2               ,
                                             const Real1DArray             *qcCharges          ,
                                                   Coordinates3            *grd1               ,
                                                   Coordinates3            *grd2               )
{
    if ( ( self != NULL ) && ( mmCharges != NULL ) && ( qcAtoms != NULL ) && ( pairList != NULL ) && ( electrostaticScale != 0.0e+00 ) &&
                                  ( crd1 != NULL ) && ( crd2 != NULL ) && ( qcCharges != NULL ) && ( grd1 != NULL ) && ( grd2 != NULL ) )
    {
        auto Integer           m, n, q ;
        auto Real              dF, gx, gy, gz, qq, qqm, r2, r2Off, scale, xm, xq, xqm, ym, yq, yqm, zm, zq, zqm ;
        auto CubicSpline      *spline ;
        auto IndexedSelection *pair   ;

        /* . Initialization. */
        r2Off  = self->outerCutoff * self->outerCutoff ;
        scale  = UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * electrostaticScale ;
        spline = self->electrostaticSpline ;

        /* . Loop over the interactions. */
        List_Iterate_Initialize ( pairList->pairs ) ;
        while ( ( pair = PairList_Iterate ( pairList ) ) != NULL )
        {
            /* . QC atom. */
            q = pair->index ;
            QCAtomContainer_GetAtomCoordinates3 ( qcAtoms, q, crd1, &xq, &yq, &zq ) ; /* . In Angstroms. */
            qq = scale * Real1DArray_Item ( qcCharges, q ) ;
       	    for ( n = 0, gx = gy = gz = 0.0e+00 ; n < pair->nindices ; n++ )
	    {
	        /* . MM atom. */
                m   = pair->indices[n] ;
                qqm = qq * Real1DArray_Item ( mmCharges, m ) ;
                Coordinates3_GetRow ( crd2, m, xm, ym, zm ) ;

                /* . Coordinate displacement. */
	        xqm = xq - xm ;
                yqm = yq - ym ;
                zqm = zq - zm ;
                r2  = ( xqm * xqm + yqm * yqm + zqm * zqm ) ;

                /* . Skip interactions outside of cutoff. */
                if ( r2 > r2Off ) continue ;

                /* . Get interaction. */
                CubicSpline_Evaluate ( spline, r2, NULL, &dF, NULL ) ;
                dF *= 2.0e+00 * qqm ;

                /* . Gradients. */
                xqm *= dF ;
                yqm *= dF ;
                zqm *= dF ;
                Coordinates3_DecrementRow ( grd2, m, xqm, yqm, zqm ) ;
                gx += xqm ;
                gy += yqm ;
                gz += zqm ;
	    }
            QCAtomContainer_SetAtomGradients3 ( qcAtoms, q, crd1, gx, gy, gz, grd1 ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM electrostatic potentials in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . crd1 refers to QC atoms and crd2 to MM atoms. Potentials is incremented not reset. */
void PairwiseInteractionABFS_QCMMPotentials ( const PairwiseInteractionABFS *self               ,
                                              const Real1DArray             *mmCharges          ,
                                              const QCAtomContainer         *qcAtoms            ,
                                                    PairList                *pairList           ,
                                              const Real                     electrostaticScale ,
                                              const Coordinates3            *crd1               ,
                                              const Coordinates3            *crd2               ,
                                              const Real1DArray             *potentials         )
{
    if ( ( self != NULL ) && ( mmCharges != NULL ) && ( qcAtoms != NULL ) && ( pairList != NULL ) && ( electrostaticScale != 0.0e+00 ) &&
                                                                         ( crd1 != NULL ) && ( crd2 != NULL ) && ( potentials != NULL ) )
    {
        auto Integer           m, n, q ;
        auto Real              fij, pot, qm, r2, r2Off, xm, xq, xqm, ym, yq, yqm, zm, zq, zqm ;
        auto CubicSpline      *spline ;
        auto IndexedSelection *pair   ;

        /* . Initialization. */
        r2Off  = self->outerCutoff * self->outerCutoff ;
        spline = self->electrostaticSpline ;

        /* . Loop over the interactions. */
        List_Iterate_Initialize ( pairList->pairs ) ;
        while ( ( pair = PairList_Iterate ( pairList ) ) != NULL )
        {
            /* . QC atom. */
            q = pair->index ;
            QCAtomContainer_GetAtomCoordinates3 ( qcAtoms, q, crd1, &xq, &yq, &zq ) ; /* . In Angstroms. */
       	    for ( n = 0, pot = 0.0e+00 ; n < pair->nindices ; n++ )
	    {
	        /* . MM atom. */
                m   = pair->indices[n] ;
                qm  = electrostaticScale * Real1DArray_Item ( mmCharges, m ) ;
                Coordinates3_GetRow ( crd2, m, xm, ym, zm ) ;

                /* . Coordinate displacement. */
	        xqm = xq - xm ;
                yqm = yq - ym ;
                zqm = zq - zm ;
                r2  = ( xqm * xqm + yqm * yqm + zqm * zqm ) ;

                /* . Skip interactions outside of cutoff. */
                if ( r2 > r2Off ) continue ;

                /* . Get interaction. */
                CubicSpline_Evaluate ( spline, r2, &fij, NULL, NULL ) ;
                pot += qm * fij ;
	    }
            Real1DArray_Item ( potentials, q ) += pot ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/QC electrostatic gradients in regular units.
! . Note these are with respect to QC gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_QCQCGradients ( const PairwiseInteractionABFS *self               ,
                                                   PairList                *pairList           ,
                                             const Real                     electrostaticScale ,
                                             const Coordinates3            *crd1               ,
                                             const Coordinates3            *crd2               ,
                                             const Real1DArray             *qcCharges          ,
                                                   Coordinates3            *grd1               ,
                                                   Coordinates3            *grd2               )
{
    if ( ( self != NULL ) && ( pairList != NULL ) && ( electrostaticScale != 0.0e+00 ) && ( crd1 != NULL ) && ( crd2 != NULL ) &&
                                                                  ( qcCharges != NULL ) && ( grd1 != NULL ) && ( grd2 != NULL ) )
    {
        auto Integer           i, j, n ;
        auto Real              dF, gx, gy, gz, qi, qij, r2, r2Off, scale, xi, xij, xj, yi, yij, yj, zi, zij, zj ;
        auto CubicSpline      *spline ;
        auto IndexedSelection *pair   ;

        /* . Initialization. */
        r2Off  = self->outerCutoff * self->outerCutoff ;
        scale  = UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * electrostaticScale ;
        spline = self->electrostaticSpline ;

        /* . Loop over the interactions. */
        List_Iterate_Initialize ( pairList->pairs ) ;
        while ( ( pair = PairList_Iterate ( pairList ) ) != NULL )
        {
            /* . First atom. */
            i  = pair->index ;
            qi = scale * Real1DArray_Item ( qcCharges, i ) ;
            Coordinates3_GetRow ( crd1, i, xi, yi, zi ) ;
       	    for ( n = 0, gx = gy = gz = 0.0e+00 ; n < pair->nindices ; n++ )
	    {
	        /* . Second atom. */
                j   = pair->indices[n] ;
                qij = qi * Real1DArray_Item ( qcCharges, j ) ;
                Coordinates3_GetRow ( crd2, j, xj, yj, zj ) ;

                /* . Coordinate displacement. */
	        xij = xi - xj ;
                yij = yi - yj ;
                zij = zi - zj ;
                r2  = ( xij * xij + yij * yij + zij * zij ) ;

                /* . Skip interactions outside of cutoff. */
                if ( r2 > r2Off ) continue ;

                /* . Get interaction. */
                CubicSpline_Evaluate ( spline, r2, NULL, &dF, NULL ) ;
                dF *= 2.0e+00 * qij ;

                /* . Gradients. */
                xij *= dF ;
                yij *= dF ;
                zij *= dF ;
                Coordinates3_DecrementRow ( grd2, j, xij, yij, zij ) ;
                gx += xij ;
                gy += yij ;
                gz += zij ;
	    }
            Coordinates3_IncrementRow ( grd1, i, gx, gy, gz ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/QC electrostatic potentials in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Potentials is incremented not reset. */
void PairwiseInteractionABFS_QCQCPotentials ( const PairwiseInteractionABFS *self               ,
                                                    PairList                *pairList           ,
                                              const Real                     electrostaticScale ,
                                              const Coordinates3            *crd1               ,
                                              const Coordinates3            *crd2               ,
                                                    SymmetricMatrix         *potentials         )
{
    if ( ( self != NULL ) && ( pairList != NULL ) && ( electrostaticScale != 0.0e+00 ) && ( crd1 != NULL ) && ( crd2 != NULL ) && ( potentials != NULL ) )
    {
        auto Integer           i, j, n ;
        auto Real              fij, r2, r2Off, xi, xij, xj, yi, yij, yj, zi, zij, zj ;
        auto CubicSpline      *spline ;
        auto IndexedSelection *pair   ;

        /* . Initialization. */
        r2Off  = self->outerCutoff * self->outerCutoff ;
        spline = self->electrostaticSpline ;

        /* . Loop over the interactions. */
        List_Iterate_Initialize ( pairList->pairs ) ;
        while ( ( pair = PairList_Iterate ( pairList ) ) != NULL )
        {
            /* . First atom. */
            i = pair->index ;
            Coordinates3_GetRow ( crd1, i, xi, yi, zi ) ;
       	    for ( n = 0 ; n < pair->nindices ; n++ )
	    {
	        /* . Second atom. */
                j = pair->indices[n] ;
                Coordinates3_GetRow ( crd2, j, xj, yj, zj ) ;

                /* . Coordinate displacement. */
	        xij = xi - xj ;
                yij = yi - yj ;
                zij = zi - zj ;
                r2  = ( xij * xij + yij * yij + zij * zij ) ;

                /* . Skip interactions outside of cutoff. */
                if ( r2 > r2Off ) continue ;

                /* . Get interaction. */
                CubicSpline_Evaluate ( spline, r2, &fij, NULL, NULL ) ;

                /* . Set the potential. */
                fij *= electrostaticScale ;
                if ( i != j ) fij *= 0.5e+00 ;
                SymmetricMatrix_IncrementComponent ( potentials, i, j, fij ) ;
	    }
        }
    }
}
