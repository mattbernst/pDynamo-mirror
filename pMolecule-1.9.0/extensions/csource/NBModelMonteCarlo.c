/*------------------------------------------------------------------------------
! . File      : NBModelMonteCarlo.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . This module implements a simple Monte Carlo non-bonding interaction method.
! . There are no derivatives, no intraisolate interactions and no QC!
!=============================================================================*/

# include "Memory.h"
# include "NBModelMonteCarlo.h"
# include "Units.h"

/* . Defaults for some variables. */
# define DEFAULT_BUFFER      0.5
# define DEFAULT_CUTOFF      8.5
# define DEFAULT_DIELECTRIC  1.0
# define DEFAULT_UNDERFLOWEL 0.5
# define DEFAULT_UNDERFLOWLJ 0.2

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
NBModelMonteCarlo *NBModelMonteCarlo_Allocate ( void )
{
    NBModelMonteCarlo *self = NULL ;
    self = ( NBModelMonteCarlo * ) Memory_Allocate ( sizeof ( NBModelMonteCarlo ) ) ;
    if ( self != NULL )
    {
        self->buffer      = DEFAULT_BUFFER      ;
        self->cutoff      = DEFAULT_CUTOFF      ;
        self->dielectric  = DEFAULT_DIELECTRIC  ;
        self->underflowel = DEFAULT_UNDERFLOWEL ;
        self->underflowlj = DEFAULT_UNDERFLOWLJ ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
NBModelMonteCarlo *NBModelMonteCarlo_Clone ( const NBModelMonteCarlo *self )
{
    NBModelMonteCarlo *new = NULL ;
    if ( self != NULL )
    {
        new = NBModelMonteCarlo_Allocate ( ) ;
        new->buffer      = self->buffer      ;
        new->cutoff      = self->cutoff      ;
        new->dielectric  = self->dielectric  ;
        new->underflowel = self->underflowel ;
        new->underflowlj = self->underflowlj ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void NBModelMonteCarlo_Deallocate ( NBModelMonteCarlo **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Full electrostatic and LJ MM/MM energy.
!-----------------------------------------------------------------------------*/
double NBModelMonteCarlo_MMMMEnergyFull ( const NBModelMonteCarlo *self, NBModelMonteCarloState *nbState )
{
    Real emmel = 0.0e+00, emmlj = 0.0e+00 ;
    if ( ( self != NULL ) && ( nbState != NULL ) )
    {
        /* . Aliases. */
        auto Boolean                      *QFREE                       ;
        auto Coordinates3         *coordinates3                ;
        auto LJParameterContainer *ljParameters                ;
        auto MMAtomContainer      *mmAtoms                     ;
        auto Selection            *iselection, *jselection     ;
        auto SelectionContainer   *isolates                    ;
        auto SymmetryParameters   *symmetryParameters          ;
        auto Vector3             **centers, *icenter, *jcenter ;

        /* . Variables. */
        auto Boolean   QFI ;
        auto Real cutsq, dxt, dyt, dzt, eel = 0.0e+00, elj = 0.0e+00, escale, fact, lowersq, qi, qscale, r2, rel2, rij2, rlj2, scale, sscale, s2, s6, underflowel2, underflowlj2, xij, yij, zij ;
        auto Integer    i, j, m, n, nljtypes, s, t, ti, tij ;

        /* . Aliases. */
        centers            = nbState->centers            ;
        coordinates3       = nbState->coordinates3       ;
        isolates           = nbState->isolates           ;
        ljParameters       = nbState->ljParameters       ;
        mmAtoms            = nbState->mmAtoms            ;
        QFREE              = nbState->QFREE              ;
        symmetryParameters = nbState->symmetryParameters ;

        /* . Initialization. */
        cutsq        = self->cutoff * self->cutoff            ;
        lowersq      = pow ( self->cutoff - self->buffer, 2 ) ;
        nljtypes     = ljParameters->ntypes                   ;
        underflowel2 = pow ( self->underflowel, 2 )           ;
        underflowlj2 = pow ( self->underflowlj, 2 )           ;

        /* . Introduce the conversion factor. */
        fact = UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE / self->dielectric ;

        /* . Double loop over isolates. */
        /* . First isolate. */
        for ( i = 0 ; i < isolates->nitems ; i++ )
        {

            /* . Get information about the isolate. */
            icenter    = centers[i] ;
            iselection = isolates->items[i] ;
            QFI        = QFREE[i] ;
            Coordinates3_Center ( coordinates3, iselection, NULL, &icenter ) ;

            /* . Second isolate. */
            for ( j = 0 ; j < i ; j++ )
            {

                /* . At least one isolate must be free. */
                if ( QFI || QFREE[j] )
                {

                    /* . Get information about the isolate. */
                    jcenter    = centers[j] ;
                    jselection = isolates->items[j] ;

                    /* . Calculate the distance between the isolate centers (applying the minimum image convention). */
                    xij = icenter->data[0] - jcenter->data[0] ;
                    yij = icenter->data[1] - jcenter->data[1] ;
                    zij = icenter->data[2] - jcenter->data[2] ;
                    SymmetryParameters_MakeMinimumImageXYZ ( symmetryParameters, &xij, &yij, &zij, &dxt, &dyt, &dzt ) ;
                    rij2 = xij * xij + yij * yij + zij * zij ;

                    /* . Check to see whether the molecules are within the cutoff distance. */
                    if ( rij2 < cutsq )
                    {

                        /* . Check for scaling. */
                        if ( ( nbState->isolatescale == i ) || ( nbState->isolatescale == j ) )
                        {
                            qscale = nbState->chargeScale  ;
                            escale = nbState->epsilonScale ;
                            sscale = nbState->sigmaScale   ;
                        }
                        else
                        {
                            qscale = escale = sscale = 1.0e+00 ;
                        }

                        /* . Double loop over atoms. */
                        /* . First atom. */
                        for ( s = 0, eel = elj = 0.0e+00 ; s < iselection->nindices ; s++ )
                        {

                            /* . Get information about the atom. */
                            m  = iselection->indices[s] ;
  	                    qi = fact     * mmAtoms->data[m].charge ;
                            ti = nljtypes * mmAtoms->data[m].ljtype ;

    	                    /* . Second atom. */
                            for ( t = 0 ; t < jselection->nindices ; t++ )
	                    {
                                /* . Get information about the atom. */
                                n = jselection->indices[t] ;

                                /* . Coordinate displacement. */
	                        Coordinates3_DifferenceRow ( coordinates3, m, n, xij, yij, zij ) ;
                                xij += dxt ;
                                yij += dyt ;
                                zij += dzt ;
                                r2   = xij * xij + yij * yij + zij * zij ;
                                rel2 = Maximum ( r2, underflowel2 ) ;
                                rlj2 = Maximum ( r2, underflowlj2 ) ;

                                /* . Distance factors. */
	                        s2   = 1.0e+00 / rlj2 ;
                                s6   = s2 * s2 * s2 * sscale ;

                                /* . Electrostatic. */
                                eel += qi * mmAtoms->data[n].charge / sqrt ( rel2 ) ;

                                /* . LJ. */
	                        tij  = ljParameters->tableindex[ti+mmAtoms->data[n].ljtype] ;
	                        elj += ( ljParameters->tableA[tij] * s6 - ljParameters->tableB[tij] ) * s6 ;
                            }
                        }

                        /* . Calculate the scale factor. */
                        if ( rij2 > lowersq ) scale = ( cutsq - rij2 ) / ( cutsq - lowersq ) ;
                        else                  scale = 1.0e+00                                ;

                        /* . Accumulate the total energies. */
                        emmel += eel * qscale * scale ;
                        emmlj += elj * escale * scale ;

                    }
                }
            }
        }

        /* . Set the energies. */
        nbState->efmmel = emmel ;
        nbState->efmmlj = emmlj ;
    }
    return emmel + emmlj ;
}

/*------------------------------------------------------------------------------
! . Single isolate electrostatic and LJ MM/MM energy and gradients.
!-----------------------------------------------------------------------------*/
double NBModelMonteCarlo_MMMMEnergySingle ( const NBModelMonteCarlo *self, const Integer isolate, NBModelMonteCarloState *nbState )
{
    Real emmel = 0.0e+00, emmlj = 0.0e+00 ;
    if ( ( self != NULL ) && ( nbState != NULL ) )
    {
        /* . Aliases. */
        auto Boolean                      *QFREE                       ;
        auto Coordinates3         *coordinates3                ;
        auto LJParameterContainer *ljParameters                ;
        auto MMAtomContainer      *mmAtoms                     ;
        auto Selection            *iselection, *jselection     ;
        auto SelectionContainer   *isolates                    ;
        auto SymmetryParameters   *symmetryParameters          ;
        auto Vector3             **centers, *icenter, *jcenter ;

        /* . Variables. */
        auto Boolean   QFI ;
        auto Real cutsq, dxt, dyt, dzt, eel = 0.0e+00, elj = 0.0e+00, escale, fact, lowersq, qi, qscale, r2, rel2, rij2, rlj2, sscale, scale, s2, s6, underflowel2, underflowlj2, xij, yij, zij ;
        auto Integer    j, m, n, nljtypes, s, t, ti, tij ;

        /* . Aliases. */
        centers            = nbState->centers            ;
        coordinates3       = nbState->coordinates3       ;
        isolates           = nbState->isolates           ;
        ljParameters       = nbState->ljParameters       ;
        mmAtoms            = nbState->mmAtoms            ;
        QFREE              = nbState->QFREE              ;
        symmetryParameters = nbState->symmetryParameters ;

        /* . Initialization. */
        cutsq        = self->cutoff * self->cutoff            ;
        lowersq      = pow ( self->cutoff - self->buffer, 2 ) ;
        nljtypes     = ljParameters->ntypes                   ;
        underflowel2 = pow ( self->underflowel, 2 )           ;
        underflowlj2 = pow ( self->underflowlj, 2 )           ;

        /* . Introduce the conversion factor. */
        fact = UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE / self->dielectric ;

        /* . Get information about the first isolate. */
        icenter    = centers[isolate] ;
        iselection = isolates->items[isolate] ;
        QFI        = QFREE[isolate] ;
        Coordinates3_Center ( coordinates3, iselection, NULL, &icenter ) ;

        /* . Loop over the remaining isolates. */
        for ( j = 0 ; j < isolates->nitems ; j++ )
        {

            /* . Skip self interaction and interactions between two fixed isolates. */
            if ( ( j != isolate ) && ( QFI || QFREE[j] ) )
            {

                /* . Get information about the isolate. */
                jcenter    = centers[j] ;
                jselection = isolates->items[j] ;
                Coordinates3_Center ( coordinates3, jselection, NULL, &jcenter ) ;

                /* . Calculate the distance between the isolate centers (applying the minimum image convention). */
                xij = icenter->data[0] - jcenter->data[0] ;
                yij = icenter->data[1] - jcenter->data[1] ;
                zij = icenter->data[2] - jcenter->data[2] ;
                SymmetryParameters_MakeMinimumImageXYZ ( symmetryParameters, &xij, &yij, &zij, &dxt, &dyt, &dzt ) ;
                rij2 = xij * xij + yij * yij + zij * zij ;

                /* . Check to see whether the molecules are within the cutoff distance. */
                if ( rij2 < cutsq )
                {

                    /* . Check for scaling. */
                    if ( ( nbState->isolatescale == isolate ) || ( nbState->isolatescale == j ) )
                    {
                        qscale = nbState->chargeScale  ;
                        escale = nbState->epsilonScale ;
                        sscale = nbState->sigmaScale   ;
                    }
                    else
                    {
                        qscale = escale = sscale = 1.0e+00 ;
                    }

                    /* . Double loop over atoms. */
                    /* . First atom. */
                    for ( s = 0, eel = elj = 0.0e+00 ; s < iselection->nindices ; s++ )
                    {

                        /* . Get information about the atom. */
                        m  = iselection->indices[s] ;
  	                qi = fact     * mmAtoms->data[m].charge ;
                        ti = nljtypes * mmAtoms->data[m].ljtype ;

    	                /* . Second atom. */
                        for ( t = 0 ; t < jselection->nindices ; t++ )
	                {
                            /* . Get information about the atom. */
                            n = jselection->indices[t] ;

                            /* . Coordinate displacement. */
	                    Coordinates3_DifferenceRow ( coordinates3, m, n, xij, yij, zij ) ;
                            xij += dxt ;
                            yij += dyt ;
                            zij += dzt ;
                            r2   = xij * xij + yij * yij + zij * zij ;
                            rel2 = Maximum ( r2, underflowel2 ) ;
                            rlj2 = Maximum ( r2, underflowlj2 ) ;

                            /* . Distance factors. */
	                    s2   = 1.0e+00 / rlj2 ;
                            s6   = s2 * s2 * s2 * sscale ;

                            /* . Electrostatic. */
                            eel += qi * mmAtoms->data[n].charge / sqrt ( rel2 ) ;

                            /* . LJ. */
	                    tij  = ljParameters->tableindex[ti+mmAtoms->data[n].ljtype] ;
	                    elj += ( ljParameters->tableA[tij] * s6 - ljParameters->tableB[tij] ) * s6 ;
                        }
                    }

                    /* . Calculate the scale factor. */
                    if ( rij2 > lowersq ) scale = ( cutsq - rij2 ) / ( cutsq - lowersq ) ;
                    else                  scale = 1.0e+00                                ;

                    /* . Accumulate the total energies. */
                    emmel += eel * qscale * scale ;
                    emmlj += elj * escale * scale ;

                }
            }
        }

        /* . Set the energies. */
        nbState->e1mmel = emmel ;
        nbState->e1mmlj = emmlj ;
    }
    return emmel + emmlj ;
}
