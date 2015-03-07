/*------------------------------------------------------------------------------
! . File      : NBModelORCA.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . This module implements full non-bonding interaction procedures that are
! . appropriate for use with the ORCA program.
!=============================================================================*/

# include "Memory.h"
# include "NBModelORCA.h"
# include "Units.h"

/* . Defaults for some variables. */
# define DEFAULT_ELECTROSTATICSCALE14 1.0
# define DEFAULT_QCMMCOUPLING         QCMMLinkAtomCoupling_RC

/*------------------------------------------------------------------------------
! . Compilation options.
!-----------------------------------------------------------------------------*/
/* . The RDIELECTRIC code is for comparison and testing of MM/MM interactions only. */
/*#define RDIELECTRIC */

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
NBModelOrca *NBModelORCA_Allocate ( void )
{
    NBModelOrca *self = NULL ;
    self = ( NBModelOrca * ) Memory_Allocate ( sizeof ( NBModelOrca ) ) ;
    if ( self != NULL )
    {
        self->electrostaticscale14 = DEFAULT_ELECTROSTATICSCALE14 ;
        self->qcmmcoupling         = DEFAULT_QCMMCOUPLING         ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
NBModelOrca *NBModelORCA_Clone ( const NBModelOrca *self )
{
    NBModelOrca *new = NULL ;
    if ( self != NULL )
    {
        new = NBModelORCA_Allocate ( ) ;
        new->electrostaticscale14 = self->electrostaticscale14 ;
        new->qcmmcoupling         = self->qcmmcoupling         ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void NBModelORCA_Deallocate ( NBModelOrca **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Electrostatic and LJ MM/MM energy and gradients.
!-----------------------------------------------------------------------------*/
void NBModelORCA_MMMMEnergy ( const NBModelOrca *self, NBModelOrcaState *nbState )
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
        fact = UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE ;

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
! . QC/MM LJ energy and gradients - standard.
!-----------------------------------------------------------------------------*/
void NBModelORCA_QCMMEnergyLJ ( const NBModelOrca *self, NBModelOrcaState *nbState )
{
    if ( ( self != NULL ) && ( nbState != NULL ) && ( nbState->qcAtoms != NULL ) )
    {
        /* . Aliases. */
        auto Boolean                      *QE14, *QFREE, *QINCL, *QMM    ;
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
