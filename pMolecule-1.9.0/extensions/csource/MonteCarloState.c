/*------------------------------------------------------------------------------
! . File      : MonteCarloState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module implements Monte Carlo state procedures.
!=================================================================================================================================*/

# include <math.h>
/*# include <stdio.h>*/

# include "Memory.h"
# include "MonteCarloState.h"
# include "Units.h"

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Check to see whether a move has been rejected.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define EXPUNDERFLOW 75.0e+00
static Boolean IsMoveRejected ( Real deltaeb, Real random )
{
    auto Boolean QACCEPT = False ;
    if ( deltaeb < EXPUNDERFLOW )
    {
        QACCEPT = ( deltaeb <= 0.0e+00 ) || ( exp ( - deltaeb ) > random ) ;
    }
    return ! QACCEPT ;
}
# undef EXPUNDERFLOW

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
MonteCarloState *MonteCarloState_Allocate ( const Integer nparticles, const Integer nrandom )
{
    MonteCarloState *self = NULL ;
    if ( ( nparticles > 0 ) && ( nrandom > 0 ) )
    {
        self = ( MonteCarloState * ) Memory_Allocate ( sizeof ( MonteCarloState ) ) ;
        if ( self != NULL )
        {
            auto Boolean QOK ;
            /* . Counters. */
            self->blocks             = 0 ;
            self->moves              = 0 ;
            self->nreject            = 0 ;
            self->nrejectm           = 0 ;
            self->nrejectt           = 0 ;
            self->nrejectv           = 0 ;
            self->ntrym              = 0 ;
            self->ntryv              = 0 ;
            /* . Current values. */
            self->beta               = 0.0e+00 ;
            self->ecurrent           = 0.0e+00 ;
            self->pressure           = 0.0e+00 ;
            self->tfact              = 0.0e+00 ;
            self->volume             = 0.0e+00 ;
            /* . Move sizes. */
            self->acceptanceratio    = 0.0e+00 ;
            self->rmax               = 0.0e+00 ;
            self->tmax               = 0.0e+00 ;
            self->vmax               = 0.0e+00 ;
            /* . Statistics. */
            self->eav                = 0.0e+00 ;
            self->eav2               = 0.0e+00 ;
            self->etot               = 0.0e+00 ;
            self->etot2              = 0.0e+00 ;
            self->etotb              = 0.0e+00 ;
            self->etotb2             = 0.0e+00 ;
            self->hav                = 0.0e+00 ;
            self->hav2               = 0.0e+00 ;
            self->htot               = 0.0e+00 ;
            self->htot2              = 0.0e+00 ;
            self->htotb              = 0.0e+00 ;
            self->htotb2             = 0.0e+00 ;
            self->vav                = 0.0e+00 ;
            self->vav2               = 0.0e+00 ;
            self->vtot               = 0.0e+00 ;
            self->vtot2              = 0.0e+00 ;
            self->vtotb              = 0.0e+00 ;
            self->vtotb2             = 0.0e+00 ;
            /* . Aliases. */
            self->coordinates3       = NULL ;
            self->nbModel            = NULL ;
            self->nbState            = NULL ;
            self->isolates           = NULL ;
            self->symmetryParameters = NULL ;
            /* . Arrays to allocate. */
            self->random                = Memory_Allocate_Array_Real_Initialize ( nrandom, 0.0e+00 ) ;
            self->oldcoordinates3       = Coordinates3_Allocate ( nparticles ) ;
            self->rotation              = Matrix33_Allocate           ( ) ;
            self->oldsymmetryParameters = SymmetryParameters_Allocate ( ) ;
            self->translation           = Vector3_Allocate            ( ) ;
            /* . Deallocate if there is not enough memory. */
            QOK = ( self->random != NULL ) && ( self->oldcoordinates3 != NULL ) && ( self->rotation != NULL ) && ( self->oldsymmetryParameters != NULL ) && ( self->translation != NULL ) ;
            if ( ! QOK ) MonteCarloState_Deallocate ( &self ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloState_Deallocate ( MonteCarloState **self )
{
    if ( (*self) != NULL )
    {
        /* . Allocated arrays. */
        Memory_Deallocate_Real        ( &((*self)->random)                ) ;
        Coordinates3_Deallocate       ( &((*self)->oldcoordinates3)       ) ;
        Matrix33_Deallocate           ( &((*self)->rotation)              ) ;
        SymmetryParameters_Deallocate ( &((*self)->oldsymmetryParameters) ) ;
        Vector3_Deallocate            ( &((*self)->translation)           ) ;
        /* . Object. */
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Adjust the move sizes.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DOWN 0.95
# define UP   1.05
void MonteCarloState_AdjustMoveSizes ( MonteCarloState *self )
{
    if ( self != NULL )
    {
        /* . Adjust the rotation and translation move sizes. */
        if ( self->ntrym > 0 )
        {
            if ( ( ( Real ) ( self->ntrym - self->nrejectm ) / ( Real ) ( self->ntrym ) ) > self->acceptanceratio )
            {
                self->rmax *=   UP ;
                self->tmax *=   UP ;
            }
            else
            {
                self->rmax *= DOWN ;
                self->tmax *= DOWN ;
            }
            self->nrejectm = 0 ; self->ntrym = 0 ;
        }

        /* . Adjust the volume move size. */
        if ( self->ntryv > 0 )
        {
            if ( ( ( Real ) ( self->ntryv - self->nrejectv ) / ( Real ) ( self->ntryv ) ) > self->acceptanceratio )
            {
                self->vmax *=   UP ;
            }
            else
            {
                self->vmax *= DOWN ;
            }
            self->nrejectv = 0 ; self->ntryv = 0 ;
        }
    }
}
# undef DOWN
# undef UP

/*----------------------------------------------------------------------------------------------------------------------------------
! . Do an isolate move.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status MonteCarloState_MoveIsolate ( MonteCarloState *self )
{
    Status status = Status_Null ;
    if ( self != NULL )
    {
        auto Real             angle, eafter, ebefore, oldenergy ;
        auto Integer                axis, chosen, i, rindex = 0 ;
        auto Coordinates3 *coordinates3 ;
        auto Matrix33     *rotation     ;
        auto Selection    *cSelection   ;
        auto Vector3      *translation  ;

        /* . Set some aliases. */
        coordinates3 = self->coordinates3 ;
        rotation     = self->rotation     ;
        translation  = self->translation  ;

        /* . Increment the number of tries. */
        self->ntrym += 1 ;

        /* . Choose a isolate to move. */
        chosen     = ( Integer ) floor ( ( Real ) self->isolates->nitems * self->random[rindex] ) ; rindex += 1 ;
        cSelection = self->isolates->items[chosen] ;

        /* . Calculate the energy of the isolate at the old configuration. */
        ebefore = NBModelMonteCarlo_MMMMEnergySingle ( self->nbModel, chosen, self->nbState ) ;

        /* . Save the old configuration. */
        oldenergy = self->ecurrent ;
        Coordinates3_Gather ( self->oldcoordinates3, coordinates3, cSelection ) ;

        /* . Calculate the center of the isolate and translate to the origin. */
        Coordinates3_Center    ( coordinates3, cSelection, NULL, &translation ) ;
        Vector3_Scale          ( translation, -1.0e+00 ) ;
        Coordinates3_Translate ( coordinates3, translation, cSelection ) ;

        /* . Do a rotation but only if the isolate has more than one particle. */
        if ( cSelection->nindices > 1 )
        {
            angle = 2.0e+00 * self->rmax * ( self->random[rindex] - 0.5e+00 ) * UNITS_ANGLE_DEGREES_TO_RADIANS ; rindex += 1 ;
            axis  = ( Integer ) floor ( 3.0e+00 * self->random[rindex] ) ;                                           rindex += 1 ;
            switch ( axis )
            {
                case 0: Matrix33_RotationAboutAxis ( &rotation, angle, 1.0e+00, 0.0e+00, 0.0e+00 ) ; break ;
                case 1: Matrix33_RotationAboutAxis ( &rotation, angle, 0.0e+00, 1.0e+00, 0.0e+00 ) ; break ;
                case 2: Matrix33_RotationAboutAxis ( &rotation, angle, 0.0e+00, 0.0e+00, 1.0e+00 ) ; break ;
            }
            Coordinates3_Rotate ( coordinates3, rotation, cSelection ) ;
        }

        /* . Calculate the translation for the isolate within the minimum image convention. */
        Vector3_Scale ( translation, -1.0e+00 ) ;
        for ( i = 0 ; i < 3 ; i++, rindex++ ) translation->data[i] += 2.0e+00 * self->tmax * ( self->random[rindex] - 0.5e+00 ) ;
        SymmetryParameters_MakeMinimumImageVector3 ( self->symmetryParameters, translation, NULL ) ;
        Coordinates3_Translate ( coordinates3, translation, cSelection ) ;

        /* . Calculate the energy of the isolate at the new configuration. */
        eafter = NBModelMonteCarlo_MMMMEnergySingle ( self->nbModel, chosen, self->nbState ) ;

        /* . Calculate the total energy of the new configuration. */
        self->ecurrent = oldenergy + eafter - ebefore ;

        /* . Check to see if the move is rejected. */
        if ( IsMoveRejected ( self->beta * ( self->ecurrent - oldenergy ), self->random[rindex] ) )
        {
            /* . Increment nreject and nrejectm. */
            self->nreject  += 1 ;
            self->nrejectm += 1 ;

            /* . Reactivate the old configuration. */
            self->ecurrent = oldenergy ;
            Coordinates3_Scatter ( self->oldcoordinates3, coordinates3, cSelection ) ;
        }

        /* . Finish up. */
        status = Status_Success ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Do a volume move.
! . The volume is changed isotropically.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status MonteCarloState_MoveVolume ( MonteCarloState *self )
{
    Status status = Status_Null ;
    if ( self != NULL )
    {
        auto Real delta, gamma, oldenergy, newvolume, oldvolume ;
        auto Integer    i, rindex = 0                        ;
        auto Coordinates3 *coordinates3             ;
        auto Selection    *iselection               ;
        auto SymmetryParameters *symmetryParameters ;
        auto Vector3      *translation              ;

        /* . Set some aliases. */
        coordinates3       = self->coordinates3       ;
        translation        = self->translation        ;
        symmetryParameters = self->symmetryParameters ;

        /* . Increment the number of tries. */
        self->ntryv += 1 ;

        /* . Save the old configuration. */
        oldenergy = self->ecurrent ;
        oldvolume = self->volume   ;
        Coordinates3_CopyTo       ( coordinates3      , self->oldcoordinates3, NULL ) ;
        SymmetryParameters_CopyTo ( symmetryParameters, self->oldsymmetryParameters ) ;

        /* . Calculate the new volume and scale the symmetry parameters accordingly. */
        newvolume = oldvolume + 2.0e+00 * self->vmax * ( self->random[rindex] - 0.5e+00 ) ; rindex += 1 ;
        gamma     = pow ( ( newvolume / oldvolume ), 1.0e+00 / 3.0e+00 )   ;
        SymmetryParameters_IsotropicScale ( symmetryParameters, gamma )    ;
        self->volume = SymmetryParameters_Volume ( symmetryParameters ) ;

        /* . Check the minimum image convention. */
        if ( SymmetryParameters_IsMinimumImageConventionSatisfied ( symmetryParameters, self->nbModel->cutoff ) )
        {
            /* . Translate the coordinates of the particles in each isolate. */
            gamma -= 1.0e+00 ;
            for ( i = 0 ; i < self->isolates->nitems ; i++ )
            {
                iselection = self->isolates->items[i] ;
                Coordinates3_Center    ( coordinates3, iselection, NULL, &translation ) ;
                Vector3_Scale          ( translation, gamma ) ;
                Coordinates3_Translate ( coordinates3, translation, iselection ) ;
            }

            /* . Calculate the total energy for the configuration. */
            self->ecurrent = NBModelMonteCarlo_MMMMEnergyFull ( self->nbModel, self->nbState ) ;

            /* . Check to see if the move is rejected. */
            delta = self->beta * ( self->ecurrent - oldenergy + self->pressure * ( self->volume - oldvolume ) - self->tfact * log ( self->volume / oldvolume ) ) ;
            if ( IsMoveRejected ( delta, self->random[rindex] ) )
            {
                /* . Increment nreject and nrejectv. */
                self->nreject  += 1 ;
                self->nrejectv += 1 ;

                /* . Reactivate the old configuration. */
                self->ecurrent = oldenergy ;
                self->volume   = oldvolume ;
                Coordinates3_CopyTo       ( self->oldcoordinates3      , coordinates3, NULL ) ;
                SymmetryParameters_CopyTo ( self->oldsymmetryParameters, symmetryParameters ) ;
            }

            /* . Finish up. */
            status = Status_Success ;
        }
        else status = Status_ValueError ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Accumulate block statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloState_StatisticsBlockAccumulate ( MonteCarloState *self )
{
    if ( self != NULL )
    {
        auto Real e, h, v ;
        e = self->ecurrent ;
        v = self->volume   ;
        h = e + self->pressure * v ;
        self->eav += e ; self->eav2 += e * e ;
        self->hav += h ; self->hav2 += h * h ;
        self->vav += v ; self->vav2 += v * v ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Start block statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloState_StatisticsBlockStart ( MonteCarloState *self )
{
    if ( self != NULL )
    {
        self->nreject = 0 ;
        self->eav = 0.0e+00 ; self->eav2 = 0.0e+00 ;
        self->hav = 0.0e+00 ; self->hav2 = 0.0e+00 ;
        self->vav = 0.0e+00 ; self->vav2 = 0.0e+00 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Stop block statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloState_StatisticsBlockStop ( MonteCarloState *self )
{
    if ( self != NULL )
    {
        auto Real e2, h2, n, v2 ;

        /* . Accumulate run statistics. */
        self->nrejectt += self->nreject ;
        self->etot += self->eav ; self->etot2 += self->eav2 ;
        self->htot += self->hav ; self->htot2 += self->hav2 ;
        self->vtot += self->vav ; self->vtot2 += self->vav2 ;

        /* . Calculate block statistics. */
        n    = ( Real ) self->moves ;
        self->eav /= n ; e2 = self->eav2 / n - self->eav * self->eav ; self->eav2 = Maximum ( e2, 0.0e+00 ) ;
        self->hav /= n ; h2 = self->hav2 / n - self->hav * self->hav ; self->hav2 = Maximum ( h2, 0.0e+00 ) ;
        self->vav /= n ; v2 = self->vav2 / n - self->vav * self->vav ; self->vav2 = Maximum ( v2, 0.0e+00 ) ;

        /* . Accumulate run block statistics. */
        self->etotb += self->eav ; self->etotb2 += self->eav * self->eav ;
        self->htotb += self->hav ; self->htotb2 += self->hav * self->hav ;
        self->vtotb += self->vav ; self->vtotb2 += self->vav * self->vav ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Start statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloState_StatisticsStart ( MonteCarloState *self )
{
    if ( self != NULL )
    {
        self->nrejectt = 0 ;
        self->nrejectm = 0 ; self->ntrym = 0 ;
        self->nrejectv = 0 ; self->ntryv = 0 ;
        /* . Run. */
        self->etot  = 0.0e+00 ; self->etot2  = 0.0e+00 ;
        self->htot  = 0.0e+00 ; self->htot2  = 0.0e+00 ;
        self->vtot  = 0.0e+00 ; self->vtot2  = 0.0e+00 ;
        /* . Block. */
        self->etotb = 0.0e+00 ; self->etotb2 = 0.0e+00 ;
        self->htotb = 0.0e+00 ; self->htotb2 = 0.0e+00 ;
        self->vtotb = 0.0e+00 ; self->vtotb2 = 0.0e+00 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Start statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloState_StatisticsStop ( MonteCarloState *self )
{
    if ( ( self != NULL ) && ( self->blocks > 1 ) )
    {
        auto Real e2, h2, n, v2 ;

        /* . Run. */
        n = ( Real ) ( self->blocks * self->moves ) ;
        self->etot /= n ; e2 = self->etot2 / n - self->etot * self->etot ; self->etot2 = Maximum ( e2, 0.0e+00 ) ;
        self->htot /= n ; h2 = self->htot2 / n - self->htot * self->htot ; self->htot2 = Maximum ( h2, 0.0e+00 ) ;
        self->vtot /= n ; v2 = self->vtot2 / n - self->vtot * self->vtot ; self->vtot2 = Maximum ( v2, 0.0e+00 ) ;

        /* . Block - see A&T page 192. */
        n = ( Real ) self->blocks ;
        self->etotb /= n ; e2 = self->etotb2 / n + self->etot * ( self->etot - 2.0e+00 * self->etotb ) ; self->etotb2 = Maximum ( e2, 0.0e+00 ) ;
        self->htotb /= n ; h2 = self->htotb2 / n + self->htot * ( self->htot - 2.0e+00 * self->htotb ) ; self->htotb2 = Maximum ( h2, 0.0e+00 ) ;
        self->vtotb /= n ; v2 = self->vtotb2 / n + self->vtot * ( self->vtot - 2.0e+00 * self->vtotb ) ; self->vtotb2 = Maximum ( v2, 0.0e+00 ) ;
    }
}
