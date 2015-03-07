/*------------------------------------------------------------------------------
! . File      : NBModelFull.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . This module implements full non-bonding interaction procedures.
! . Implementations are straightforward with no optimization.
!=============================================================================*/

# include "Memory.h"
# include "NBModelFull.h"
# include "Units.h"

/* . Defaults for some variables. */
# define DEFAULT_DIELECTRIC            1.0
# define DEFAULT_ELECTROSTATICSCALE14  1.0
# define DEFAULT_QCQCMMLinkAtomCoupling_MM          QCMMLinkAtomCoupling_RC

/*------------------------------------------------------------------------------
! . Compilation options.
!-----------------------------------------------------------------------------*/
/* . The RDIELECTRIC code is for comparison and testing of MM/MM interactions only. */
/*#define RDIELECTRIC */

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
NBModelFull *NBModelFull_Allocate ( void )
{
    NBModelFull *self = NULL ;
    self = ( NBModelFull * ) Memory_Allocate ( sizeof ( NBModelFull ) ) ;
    if ( self != NULL )
    {
        self->dielectric           = DEFAULT_DIELECTRIC           ;
        self->electrostaticscale14 = DEFAULT_ELECTROSTATICSCALE14 ;
        self->qcmmcoupling         = DEFAULT_QCQCMMLinkAtomCoupling_MM         ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
NBModelFull *NBModelFull_Clone ( const NBModelFull *self )
{
    NBModelFull *new = NULL ;
    if ( self != NULL )
    {
        new = NBModelFull_Allocate ( ) ;
        new->dielectric           = self->dielectric           ;
        new->electrostaticscale14 = self->electrostaticscale14 ;
        new->qcmmcoupling         = self->qcmmcoupling         ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void NBModelFull_Deallocate ( NBModelFull **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Electrostatic and LJ MM/MM energy and gradients.
!-----------------------------------------------------------------------------*/
void NBModelFull_MMMMEnergy ( const NBModelFull *self, NBModelFullState *nbState )
{
    if ( ( self != NULL ) && ( nbState != NULL ) )
    {
        /* . Aliases. */
        auto Boolean                      *QE14, *QFREE, *QINCL, *QMM    ;
        auto Coordinates3         *coordinates3                  ;
        auto Coordinates3           *gradients3                    ;
        auto LJParameterContainer *ljParameters, *ljParameters14 ;
        auto MMAtomContainer      *mmAtoms                       ;
        auto PairList             *exclusions, *interactions14   ;

        /* . Variables. */
        auto Boolean   QEXCL, QEX14, QFI, QGRADIENTS ;
        auto Real as6 = 0.0e+00, df, eel = 0.0e+00, elj = 0.0e+00, emmel, emmel14, emmlj, emmlj14, escale14, fact, qi, s, s3, xij, yij, zij ;
        auto Integer    i, ifirst = 0, ilast = 0, j, n, nljtypes, nljtypes14, ti, tij, ti14, xfirst = 0, xlast = 0 ;

        /* . Aliases. */
        QE14           = nbState->QE14              ;
        QFREE          = nbState->QFREE             ;
        QINCL          = nbState->QINCL             ;
        QMM            = nbState->QMM               ;
        coordinates3   = nbState->coordinates3      ;
        escale14       = self->electrostaticscale14 ;
        exclusions     = nbState->exclusions        ;
        gradients3     = nbState->gradients3        ;
        interactions14 = nbState->interactions14    ;
        ljParameters   = nbState->ljParameters      ;
        ljParameters14 = nbState->ljParameters14    ;
        mmAtoms        = nbState->mmAtoms           ;

        /* . Flags. */
        QEXCL      = ( exclusions     != NULL ) ;
        QEX14      = ( interactions14 != NULL ) ;
        QGRADIENTS = ( gradients3     != NULL ) ;

        /* . Initialization. */
        nljtypes   = ljParameters->ntypes   ;
        nljtypes14 = ljParameters14->ntypes ;

        /* . Introduce the conversion factor. */
        fact = UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE / self->dielectric ;

        /* . Double loop over MM atoms. */
        /* . First atom. */
        for ( i = 0, emmel = emmel14 = emmlj = emmlj14 = 0.0e+00 ; i < mmAtoms->natoms ; i++ )
        {
            if ( QMM[i] )
            {
                /* . Flag excluded atoms for i. */
                if ( QEXCL )
                {
                    xfirst = exclusions->connectionsi[i]   ;
                    xlast  = exclusions->connectionsi[i+1] ;
                    for ( n = xfirst ; n < xlast ; n++ )
                    {
                        j = exclusions->connectionsj[n] ;
                        if ( j >= i ) break ;
                        QINCL[j] = False ;
                    }
                }
                /* . Flag 1-4 atoms for i. */
                if ( QEX14 )
                {
                    ifirst = interactions14->connectionsi[i]   ;
                    ilast  = interactions14->connectionsi[i+1] ;
                    for ( n = ifirst ; n < ilast ; n++ )
                    {
                        j = interactions14->connectionsj[n] ;
                        if ( j >= i ) break ;
                        QE14[j] = True ;
                    }
                }

                /* . Get information for atom i. */
                QFI  = QFREE[i] ;
  	        qi   = fact       * mmAtoms->data[i].charge ;
                ti   = nljtypes   * mmAtoms->data[i].ljtype ;
                ti14 = nljtypes14 * mmAtoms->data[i].ljtype ;

	        /* . Second atom. */
                /* . Note that i cannot interact with itself here (only for QC/MM electrostatic interactions). */
                for ( j = 0 ; j < i ; j++ )
	        {
                    /* . Check for an MM atom that is not an exclusion. One atom of the pair must be moving. */
                    if ( QMM[j] && ( QE14[j] || QINCL[j] ) && ( QFI || QFREE[j] ) )
                    {
                        /* . Coordinate displacement. */
	                Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
	                s  = 1.0e+00 / ( xij * xij + yij * yij + zij * zij ) ;
                        s3 = s * s * s ;

                        /* . Electrostatic. */
# ifdef RDIELECTRIC
                        eel = qi * mmAtoms->data[j].charge * s ;
# else
                        eel = qi * mmAtoms->data[j].charge * sqrt ( s ) ;
# endif

                        /* . Energy. */
                        if ( QE14[j] )
                        {
                            /* . Electrostatic. */
                            eel     *= escale14 ;
	                    emmel14 += eel      ;
                            /* . LJ. */
	                    tij      = ljParameters14->tableindex[ti14+mmAtoms->data[j].ljtype] ;
	                    as6      = ljParameters14->tableA[tij] * s3 * s3 ;
	                    elj      = as6 - ljParameters14->tableB[tij] * s3 ;
                            emmlj14 += elj ;
                        }
                        else
                        {
                            /* . Electrostatic. */
	                    emmel   += eel ;
                            /* . LJ. */
	                    tij      = ljParameters->tableindex[ti+mmAtoms->data[j].ljtype] ;
	                    as6      = ljParameters->tableA[tij] * s3 * s3 ;
	                    elj      = as6 - ljParameters->tableB[tij] * s3 ;
                            emmlj   += elj ;
                        }

                        /* . Gradients. */
	                if ( QGRADIENTS )
	                {
# ifdef RDIELECTRIC
                            df   =           eel ;
# else
                            df   = 0.5e+00 * eel ;
# endif
                            df  += 3.0e+00 * ( as6 + elj ) ;
                            df  *= -2.0e+00 * s ;
	                    xij *= df ;
	                    yij *= df ;
	                    zij *= df ;
	                    Coordinates3_IncrementRow ( gradients3, i, xij, yij, zij ) ;
	                    Coordinates3_DecrementRow ( gradients3, j, xij, yij, zij ) ;
	                }
                    }
                }

                /* . Unflag excluded atoms for i. */
                if ( QEXCL )
                {
                    for ( n = xfirst ; n < xlast ; n++ )
                    {
                        j = exclusions->connectionsj[n] ;
                        if ( j >= i ) break ;
                        QINCL[j] = True ;
                    }
                }
                /* . Unflag 1-4 atoms for i. */
                if ( QEX14 )
                {
                    for ( n = ifirst ; n < ilast ; n++ )
                    {
                        j = interactions14->connectionsj[n] ;
                        if ( j >= i ) break ;
                        QE14[j] = False ;
                    }
                }
            }
        }

        /* . Set the energies. */
        nbState->emmel   = emmel   ;
        nbState->emmel14 = emmel14 ;
        nbState->emmlj   = emmlj   ;
        nbState->emmlj14 = emmlj14 ;
    }
}

/*------------------------------------------------------------------------------
! . QC/MM LJ energy and gradients.
!-----------------------------------------------------------------------------*/
void NBModelFull_QCMMEnergyLJ ( const NBModelFull *self, NBModelFullState *nbState )
{
    if ( ( self != NULL ) && ( nbState != NULL ) && ( nbState->qcAtoms != NULL ) )
    {
        /* . Aliases. */
        auto Boolean                 *QE14, *QFREE, *QINCL, *QMM    ;
        auto Coordinates3         *coordinates3                  ;
        auto Coordinates3           *gradients3                    ;
        auto LJParameterContainer *ljParameters, *ljParameters14 ;
        auto MMAtomContainer      *mmAtoms                       ;
        auto PairList             *exclusions, *interactions14   ;
        auto QCAtomContainer      *qcAtoms                       ;

        /* . Variables. */
        auto Boolean   QEXCL, QEX14, QFI, QGRADIENTS ;
        auto Real as6, df, elj, energy, energy14, gx, gy, gz, s, s3, xm, xq, xqm, ym, yq, yqm, zm, zq, zqm ;
        auto Integer    i, ifirst = 0, ilast = 0, j, nljtypes, nljtypes14, q, ti, tij, ti14, xfirst = 0, xlast = 0 ;

        /* . Aliases. */
        QE14           = nbState->QE14           ;
        QFREE          = nbState->QFREE          ;
        QINCL          = nbState->QINCL          ;
        QMM            = nbState->QMM            ;
        coordinates3   = nbState->coordinates3   ;
        exclusions     = nbState->exclusions     ;
        gradients3     = nbState->gradients3     ;
        interactions14 = nbState->interactions14 ;
        ljParameters   = nbState->ljParameters   ;
        ljParameters14 = nbState->ljParameters14 ;
        mmAtoms        = nbState->mmAtoms        ;
        qcAtoms        = nbState->qcAtoms        ;

        /* . Flags. */
        QEXCL      = ( exclusions     != NULL ) ;
        QEX14      = ( interactions14 != NULL ) ;
        QGRADIENTS = ( gradients3     != NULL ) ;

        /* . Initialization. */
        nljtypes   = ljParameters->ntypes   ;
        nljtypes14 = ljParameters14->ntypes ;

        /* . Double loop over atoms. */
        /* . QC atom. */
        for ( q = 0, energy = energy14 = 0.0e+00 ; q < qcAtoms->natoms ; q++ )
        {
            /* . Skip boundary atoms. */
            if ( qcAtoms->data[q].QBOUNDARY ) continue ;

            /* . Atom index for q. */
            i    = qcAtoms->data[q].index ;
            QFI  = QFREE[i] ;
            ti   = nljtypes   * mmAtoms->data[i].ljtype ;
            ti14 = nljtypes14 * mmAtoms->data[i].ljtype ;
            QCAtomContainer_GetAtomCoordinates3 ( qcAtoms, q, coordinates3, &xq, &yq, &zq ) ; /* . In Angstroms. */

            /* . Flag excluded atoms for i. */
            if ( QEXCL )
            {
                xfirst = exclusions->connectionsi[i]   ;
                xlast  = exclusions->connectionsi[i+1] ;
                for ( j = xfirst ; j < xlast ; j++ ) QINCL[exclusions->connectionsj[j]] = False ;
            }
            /* . Flag 1-4 atoms for i. */
            if ( QEX14 )
            {
                ifirst = interactions14->connectionsi[i]   ;
                ilast  = interactions14->connectionsi[i+1] ;
                for ( j = ifirst ; j < ilast ; j++ ) QE14[interactions14->connectionsj[j]] = True ;
            }

	    /* . MM atom. */
            for ( j = 0, gx = gy = gz = 0.0e+00 ; j < mmAtoms->natoms ; j++ )
	    {
                /* . One atom must be moving. */
                if ( QFI || QFREE[j] )
                {
                    /* . Include non-self interactions, inclusions and interactions with MM atoms. */
                    if ( ( i != j ) && QMM[j] && ( QE14[j] || QINCL[j] ) )
                    {
                        /* . Coordinate displacement. */
                        Coordinates3_GetRow ( coordinates3, j, xm, ym, zm ) ;
	                xqm = xq - xm ;
                        yqm = yq - ym ;
                        zqm = zq - zm ;
	                s   = 1.0e+00 / ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
                        s3  = s * s * s ;
                        /* . LJ energy. */
                        if ( QE14[j] )
                        {
	                    tij       = ljParameters14->tableindex[ti14+mmAtoms->data[j].ljtype] ;
	                    as6       = ljParameters14->tableA[tij] * s3 * s3 ;
	                    elj       = as6 - ljParameters14->tableB[tij] * s3 ;
                            energy14 += elj ;
                        }
                        else
                        {
	                    tij       = ljParameters->tableindex[ti+mmAtoms->data[j].ljtype] ;
	                    as6       = ljParameters->tableA[tij] * s3 * s3 ;
	                    elj       = as6 - ljParameters->tableB[tij] * s3 ;
                            energy   += elj ;
                        }
                        /* . Gradients. */
	                if ( QGRADIENTS )
	                {
                            df   = -6.0e+00 * ( as6 + elj ) * s ;
                            xqm *= df ;
                            yqm *= df ;
                            zqm *= df ;
                            Coordinates3_DecrementRow ( gradients3, j, xqm, yqm, zqm ) ;
                            gx += xqm ;
                            gy += yqm ;
                            gz += zqm ;
                        }
	            }
                }
            }

            /* . Save the gradients. */
            if ( QGRADIENTS )
            {
                QCAtomContainer_SetAtomGradients3 ( qcAtoms, q, coordinates3, gx, gy, gz, gradients3 ) ;
            }

            /* . Unflag excluded atoms for i. */
            if ( QEXCL )
            {
                for ( j = xfirst ; j < xlast ; j++ ) QINCL[exclusions->connectionsj[j]] = True ;
            }
            /* . Unflag 1-4 atoms for i. */
            if ( QEX14 )
            {
                for ( j = ifirst ; j < ilast ; j++ ) QE14[interactions14->connectionsj[j]] = False ;
            }
        }

        /* . Set the energies. */
        nbState->eqcmmlj   = energy   ;
        nbState->eqcmmlj14 = energy14 ;
    }
}

/*------------------------------------------------------------------------------
! . QC/MM electrostatic gradients in regular units.
!-----------------------------------------------------------------------------*/
void NBModelFull_QCMMGradients ( const NBModelFull *self, NBModelFullState *nbState )
{
    if ( ( self != NULL ) && ( nbState != NULL ) && ( nbState->qcAtoms != NULL ) && ( nbState->gradients3 != NULL ) )
    {
        /* . Aliases. */
        auto Boolean            *QE14, *QINCL, *QMM          ;
        auto Coordinates3    *coordinates3                ;
        auto Coordinates3      *gradients3                  ;
        auto Integer             *baindex                     ;
        auto MMAtomContainer *mmAtoms                     ;
        auto PairList        *exclusions, *interactions14 ;
        auto QCAtomContainer *qcAtoms                     ;
        auto Real1DArray     *qcCharges                   ;

        /* . Variables. */
        auto Boolean   QEXCL, QEX14, QRD, Q2 ;
        auto Real df, escale14, fact, gx, gy, gz, qq, q0, s, xm, xp, xq, xqm, ym, yp, yq, yqm, zm, zp, zq, zqm ;
        auto Integer    b, i, ifirst = 0, ilast = 0, j, m, n, p, q, xfirst = 0, xlast = 0 ;

        /* . Aliases. */
        QMM            = nbState->QMM               ;
        coordinates3   = nbState->coordinates3      ;
        gradients3     = nbState->gradients3        ;
        mmAtoms        = nbState->mmAtoms           ;
        qcAtoms        = nbState->qcAtoms           ;
        qcCharges      = nbState->qcCharges         ;

        /* . Constants. */
        fact = ( UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE / UNITS_LENGTH_ANGSTROMS_TO_BOHRS ) / self->dielectric ;

        /* . MM coupling. */
        if ( self->qcmmcoupling == QCMMLinkAtomCoupling_MM )
        {
            /* . Aliases. */
            QE14           = nbState->QE14              ;
            QINCL          = nbState->QINCL             ;
            escale14       = self->electrostaticscale14 ;
            exclusions     = nbState->exclusions        ;
            interactions14 = nbState->interactions14    ;

            /* . Flags. */
            QEXCL = ( exclusions     != NULL ) ;
            QEX14 = ( interactions14 != NULL ) ;

            /* . Double loop over atoms. */
            /* . QC atom. */
            for ( q = 0 ; q < qcAtoms->natoms ; q++ )
            {
                /* . Atom index for q. */
                i  = qcAtoms->data[q].index ;
                qq = fact * Real1DArray_Item ( qcCharges, q ) ;
                QCAtomContainer_GetAtomCoordinates3 ( qcAtoms, q, coordinates3, &xq, &yq, &zq ) ; /* . In Angstroms. */

                /* . Flag excluded atoms for i. */
                if ( QEXCL )
                {
                    xfirst = exclusions->connectionsi[i]   ;
                    xlast  = exclusions->connectionsi[i+1] ;
                    for ( j = xfirst ; j < xlast ; j++ ) QINCL[exclusions->connectionsj[j]] = False ;
                }
                /* . Flag 1-4 atoms for i. */
                if ( QEX14 )
                {
                    ifirst = interactions14->connectionsi[i]   ;
                    ilast  = interactions14->connectionsi[i+1] ;
                    for ( j = ifirst ; j < ilast ; j++ ) QE14[interactions14->connectionsj[j]] = True ;
                }

	        /* . MM atom. */
                for ( j = 0, gx = gy = gz = 0.0e+00 ; j < mmAtoms->natoms ; j++ )
	        {
                    /* . Include non-self interactions, inclusions and interactions with MM atoms. */
                    if ( ( i != j ) && QMM[j] && ( QE14[j] || QINCL[j] ) )
                    {
                        /* . Coordinate displacement. */
                        Coordinates3_GetRow ( coordinates3, j, xm, ym, zm ) ;
	                xqm = xq - xm ;
                        yqm = yq - ym ;
                        zqm = zq - zm ;
	                s   = 1.0e+00 / ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
                        /* . Electrostatic. */
                        df  = - ( qq * mmAtoms->data[j].charge ) * sqrt ( s ) * s ;
                        if ( QE14[j] ) df *= escale14 ;
                        xqm *= df ;
                        yqm *= df ;
                        zqm *= df ;
                        Coordinates3_DecrementRow ( gradients3, j, xqm, yqm, zqm ) ;
                        gx += xqm ;
                        gy += yqm ;
                        gz += zqm ;
	            }
                }

                /* . Save the gradients. */
                QCAtomContainer_SetAtomGradients3 ( qcAtoms, q, coordinates3, gx, gy, gz, gradients3 ) ;

                /* . Unflag excluded atoms for i. */
                if ( QEXCL )
                {
                    for ( j = xfirst ; j < xlast ; j++ ) QINCL[exclusions->connectionsj[j]] = True ;
                }
                /* . Unflag 1-4 atoms for i. */
                if ( QEX14 )
                {
                    for ( j = ifirst ; j < ilast ; j++ ) QE14[interactions14->connectionsj[j]] = False ;
                }
            }
        }
        /* . RC and RD coupling. */
        else
        {
            /* . Aliases. */
            baindex = nbState->baindex ;

            /* . Set the RD flag. */
            QRD = ( self->qcmmcoupling == QCMMLinkAtomCoupling_RD ) ;

            /* . Double loop over atoms. */
            /* . QC atom. */
            for ( q = 0 ; q < qcAtoms->natoms ; q++ )
            {
                /* . Atom index for q. */
                i  = qcAtoms->data[q].index ;
                qq = fact * Real1DArray_Item ( qcCharges, q ) ;
                QCAtomContainer_GetAtomCoordinates3 ( qcAtoms, q, coordinates3, &xq, &yq, &zq ) ; /* . In Angstroms. */

	        /* . MM atom. */
                for ( j = 0, gx = gy = gz = 0.0e+00 ; j < mmAtoms->natoms ; j++ )
	        {
                    if ( QMM[j] )
                    {
                        /* . Get the coordinates. */
                        Coordinates3_GetRow ( coordinates3, j, xm, ym, zm ) ;

                        /* . Boundary atom. */
                        b = baindex[j] ;
                        if ( b >= 0 )
                        {
                            /* . Loop over the MM partners - n should never be zero. */
                            n  = qcAtoms->data[b].nmmpartners ;
                            q0 = qq * mmAtoms->data[j].charge / ( Real ) n ;
                            for ( m = 0 ; m < n ; m++ )
                            {
                                /* . Get the partner and check for a second contribution. */
                                p    = qcAtoms->data[b].mmpartners[m] ;
                                Q2   = ( baindex[p] < 0 ) && QRD      ; /* . Only for RD non-boundary atom partners. */
                                /* . First contribution. */
                                Coordinates3_GetRow ( coordinates3, p, xp, yp, zp ) ;
	                        xqm  = xq - 0.5e+00 * ( xm + xp ) ;
                                yqm  = yq - 0.5e+00 * ( ym + yp ) ;
                                zqm  = zq - 0.5e+00 * ( zm + zp ) ;
                                s    = 1.0e+00 / ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
                                df   = - q0 * sqrt ( s ) * s ;
                                if ( Q2 ) df *= 2.0e+00 ;
                                xqm *= df ;
                                yqm *= df ;
                                zqm *= df ;
                                gx  += xqm ;
                                gy  += yqm ;
                                gz  += zqm ;
                                xqm *= 0.5e+00 ;
                                yqm *= 0.5e+00 ;
                                zqm *= 0.5e+00 ;
                                Coordinates3_DecrementRow ( gradients3, j, xqm, yqm, zqm ) ;
                                Coordinates3_DecrementRow ( gradients3, p, xqm, yqm, zqm ) ;
                                /* . Second contribution. */
                                if ( Q2 )
                                {
	                            xqm  = xq - xp ;
                                    yqm  = yq - yp ;
                                    zqm  = zq - zp ;
                                    s    = 1.0e+00 / ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
                                    df   = q0 * sqrt ( s ) * s ;
                                    xqm *= df ;
                                    yqm *= df ;
                                    zqm *= df ;
                                    gx  += xqm ;
                                    gy  += yqm ;
                                    gz  += zqm ;
                                    Coordinates3_DecrementRow ( gradients3, p, xqm, yqm, zqm ) ;
                                }
                            }
                        }
                        /* . Normal atom. */
                        else
                        {
	                    xqm = xq - xm ;
                            yqm = yq - ym ;
                            zqm = zq - zm ;
	                    s   = 1.0e+00 / ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
                            df  = - ( qq * mmAtoms->data[j].charge ) * sqrt ( s ) * s ;
                            xqm *= df ;
                            yqm *= df ;
                            zqm *= df ;
                            Coordinates3_DecrementRow ( gradients3, j, xqm, yqm, zqm ) ;
                            gx  += xqm ;
                            gy  += yqm ;
                            gz  += zqm ;
                        }
	            }
                }

                /* . Save the gradients. */
                QCAtomContainer_SetAtomGradients3 ( qcAtoms, q, coordinates3, gx, gy, gz, gradients3 ) ;
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . QC/MM electrostatic potentials in atomic units.
!-----------------------------------------------------------------------------*/
void NBModelFull_QCMMPotentials  ( const NBModelFull *self, NBModelFullState *nbState )
{
    if ( ( self != NULL ) && ( nbState != NULL ) && ( nbState->qcAtoms != NULL ) )
    {
        /* . Aliases. */
        auto Boolean            *QE14, *QINCL, *QMM            ;
        auto Coordinates3    *coordinates3                  ;
        auto Integer             *baindex                       ;
        auto MMAtomContainer *mmAtoms                       ;
        auto PairList        *exclusions, *interactions14   ;
        auto QCAtomContainer *qcAtoms                       ;
        auto Real1DArray     *qcmmPotentials                ;

        /* . Variables. */
        auto Boolean   QEXCL, QEX14, QRC ;
        auto Real f, pot, pot14, q0, scalingFactor, xm, xp, xq, xqm, ym, yp, yq, yqm, zm, zp, zq, zqm ;
        auto Integer    b, i, ifirst = 0, ilast = 0, j, m, n, p, q, xfirst = 0, xlast = 0    ;

        /* . Aliases. */
        QMM            = nbState->QMM            ;
        coordinates3   = nbState->coordinates3   ;
        mmAtoms        = nbState->mmAtoms        ;
        qcAtoms        = nbState->qcAtoms        ;
        qcmmPotentials = nbState->qcmmPotentials ;

        /* . Initialization. */
        scalingFactor = 1.0e+00 / ( self->dielectric * UNITS_LENGTH_ANGSTROMS_TO_BOHRS ) ;

        /* . MM coupling. */
        if ( self->qcmmcoupling == QCMMLinkAtomCoupling_MM )
        {
            /* . Aliases. */
            QE14           = nbState->QE14           ;
            QINCL          = nbState->QINCL          ;
            exclusions     = nbState->exclusions     ;
            interactions14 = nbState->interactions14 ;

            /* . Flags. */
            QEXCL = ( exclusions     != NULL ) ;
            QEX14 = ( interactions14 != NULL ) ;

            /* . Double loop over atoms. */
            /* . QC atom. */
            for ( q = 0 ; q < qcAtoms->natoms ; q++ )
            {
                /* . Atom index for q. */
                i = qcAtoms->data[q].index ;
                QCAtomContainer_GetAtomCoordinates3 ( qcAtoms, q, coordinates3, &xq, &yq, &zq ) ; /* . In Angstroms. */

                /* . Flag excluded atoms for i. */
                if ( QEXCL )
                {
                    xfirst = exclusions->connectionsi[i]   ;
                    xlast  = exclusions->connectionsi[i+1] ;
                    for ( j = xfirst ; j < xlast ; j++ ) QINCL[exclusions->connectionsj[j]] = False ;
                }
                /* . Flag 1-4 atoms for i. */
                if ( QEX14 )
                {
                    ifirst = interactions14->connectionsi[i]   ;
                    ilast  = interactions14->connectionsi[i+1] ;
                    for ( j = ifirst ; j < ilast ; j++ ) QE14[interactions14->connectionsj[j]] = True ;
                }

	        /* . MM atom. */
                for ( j = 0, pot = pot14 = 0.0e+00 ; j < mmAtoms->natoms ; j++ )
	        {
                    /* . Include non-self interactions, inclusions and interactions with MM atoms. */
                    if ( ( i != j ) && QMM[j] && ( QE14[j] || QINCL[j] ) )
                    {
                        /* . Coordinate displacement. */
                        Coordinates3_GetRow ( coordinates3, j, xm, ym, zm ) ;
	                xqm = xq - xm ;
                        yqm = yq - ym ;
                        zqm = zq - zm ;
                        /* . Electrostatic. */
                        f = mmAtoms->data[j].charge / sqrt ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
                        if ( QE14[j] ) pot14 += f ;
                        else           pot   += f ;
	            }
                }

                /* . Save the potential. */
                Real1DArray_Item ( qcmmPotentials, q ) += ( ( pot + self->electrostaticscale14 * pot14 ) * scalingFactor ) ;

                /* . Unflag excluded atoms for i. */
                if ( QEXCL )
                {
                    for ( j = xfirst ; j < xlast ; j++ ) QINCL[exclusions->connectionsj[j]] = True ;
                }
                /* . Unflag 1-4 atoms for i. */
                if ( QEX14 )
                {
                    for ( j = ifirst ; j < ilast ; j++ ) QE14[interactions14->connectionsj[j]] = False ;
                }
            }
        }
        /* . RC and RD coupling. */
        else
        {
            /* . Aliases. */
            baindex = nbState->baindex ;

            /* . Set the RC flag. */
            QRC = ( self->qcmmcoupling == QCMMLinkAtomCoupling_RC ) ;

            /* . Double loop over atoms. */
            /* . QC atom. */
            for ( q = 0 ; q < qcAtoms->natoms ; q++ )
            {
                /* . Atom index for q. */
                i = qcAtoms->data[q].index ;
                QCAtomContainer_GetAtomCoordinates3 ( qcAtoms, q, coordinates3, &xq, &yq, &zq ) ; /* . In Angstroms. */

	        /* . MM atom. */
                for ( j = 0, pot = 0.0e+00 ; j < mmAtoms->natoms ; j++ )
	        {
                    if ( QMM[j] )
                    {
                        /* . Get the coordinates. */
                        Coordinates3_GetRow ( coordinates3, j, xm, ym, zm ) ;

                        /* . Boundary atom. */
                        b = baindex[j] ;
                        if ( b >= 0 )
                        {
                            /* . Loop over the MM partners - n should never be zero. */
                            n  = qcAtoms->data[b].nmmpartners ;
                            q0 = mmAtoms->data[j].charge / ( Real ) n ;
                            for ( m = 0 ; m < n ; m++ )
                            {
                                p    = qcAtoms->data[b].mmpartners[m] ;
                                Coordinates3_GetRow ( coordinates3, p, xp, yp, zp ) ;
	                        xqm  = xq - 0.5e+00 * ( xm + xp ) ;
                                yqm  = yq - 0.5e+00 * ( ym + yp ) ;
                                zqm  = zq - 0.5e+00 * ( zm + zp ) ;
                                /* . RC or RD boundary atom partners. */
                                if ( ( baindex[p] >= 0 ) || QRC )
                                {
                                    pot += q0 / sqrt ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
                                }
                                /* . RD. */
                                else
                                {
                                    /* . First contribution. */
                                    pot += 2.0e+00 * q0 / sqrt ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
                                    /* . Second contribution. */
	                            xqm  = xq - xp ;
                                    yqm  = yq - yp ;
                                    zqm  = zq - zp ;
                                    pot -= q0 / sqrt ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
                                }
                            }
                        }
                        /* . Normal atom. */
                        else
                        {
	                    xqm = xq - xm ;
                            yqm = yq - ym ;
                            zqm = zq - zm ;
                            pot += mmAtoms->data[j].charge / sqrt ( xqm * xqm + yqm * yqm + zqm * zqm ) ;
	                }
                    }
                }

                /* . Save the potential. */
                Real1DArray_Item ( qcmmPotentials, q ) += ( pot * scalingFactor ) ;
            }
        }
    }
}
