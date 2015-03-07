/*------------------------------------------------------------------------------
! . File      : SSBPModel.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Procedures for the SSBP boundary model.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "Coordinates3.h"
# include "Integer1DArray.h"
# include "Integer2DArray.h"
# include "Memory.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Selection.h"
# include "SSBPModel.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Angular potential. */
# define CAngular1  (  2.4085919002   * UNITS_ENERGY_CALORIES_TO_JOULES )
# define CAngular2  (  1.766839115555 * UNITS_ENERGY_CALORIES_TO_JOULES )
# define CAngular3  ( -3.067139576698 * UNITS_ENERGY_CALORIES_TO_JOULES )
# define CAngular4  ( -1.20064327666  * UNITS_ENERGY_CALORIES_TO_JOULES )
# define CAngular5  (  0.840661       * UNITS_ENERGY_CALORIES_TO_JOULES )
# define DRHA      -1.0

/* . Cavity potential. */
# define ACavity1  (  0.56198800      * UNITS_ENERGY_CALORIES_TO_JOULES )
# define ACavity2  ( -0.072798148     * UNITS_ENERGY_CALORIES_TO_JOULES )
# define ACavity3  (  0.00426122036   * UNITS_ENERGY_CALORIES_TO_JOULES )
# define ACavity4  ( -0.0000925233817 * UNITS_ENERGY_CALORIES_TO_JOULES )
# define ACavity5  ( -1.6649500       * UNITS_ENERGY_CALORIES_TO_JOULES )
# define ACavity6  (  0.0840          * UNITS_ENERGY_CALORIES_TO_JOULES )
# define ACavity7 15.39333
# define ACavity8  (  8.5             * UNITS_ENERGY_CALORIES_TO_JOULES )

# define BCavity1  1.319978287
# define BCavity2  ( -0.840953501               * UNITS_ENERGY_CALORIES_TO_JOULES )
# define BCavity3  ( -0.001602388122            * UNITS_ENERGY_CALORIES_TO_JOULES )
# define BCavity4  ( -8.392886499               * UNITS_ENERGY_CALORIES_TO_JOULES )
# define BCavity5  ( (-0.840953501-8.392886499) * UNITS_ENERGY_CALORIES_TO_JOULES )
# define BCavity6  ( 1.6                        * UNITS_ENERGY_CALORIES_TO_JOULES )
# define BCavity7  ( -8.4751210228              * UNITS_ENERGY_CALORIES_TO_JOULES )

# define BLower   -5.0

# define CavityDecrement 2.6

/* . Empirical correction. */
# define EightPt88 8.88

/* . Kirkwood parameters. */
# define OnePt6 1.6e+00
# define RSmall 1.0e-10

/* . Default parameters. */
# define DEFAULT_DOANGULARPOTENTIAL      True
# define DEFAULT_DOCAVITYPOTENTIAL       True
# define DEFAULT_DOEMPIRICALCORRECTION   True
# define DEFAULT_DOHARDSPHERERESTRICTION True
# define DEFAULT_DOKIRKWOOD              True
# define DEFAULT_FIXCAVITYRADIUS         False
# define DEFAULT_MAXIMUML                15
# define DEFAULT_CAVITYRADIUS             0.0
# define DEFAULT_CAVITYRADIUSINCREMENT    2.6
# define DEFAULT_DIELECTRICINSIDE         1.0
# define DEFAULT_DIELECTRICOUTSIDE       80.0
# define DEFAULT_EMPIRICAL1             ( 1.1     * UNITS_ENERGY_CALORIES_TO_JOULES )
# define DEFAULT_EMPIRICAL2               8.0e-03
# define DEFAULT_KIRKWOODRADIUSINCREMENT  2.8
# define DEFAULT_PRESSURE                 1.0
# define DEFAULT_SURFACETENSION         ( 3.3e-02 * UNITS_ENERGY_CALORIES_TO_JOULES )

/* . Need to handle fixed atoms explicitly? */

/*
! . There is some ambiguity in the original papers.
! . Here the angular and empirical terms can only be calculated for waters whereas the cavity term may be determined for any atom.
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Polynomial_CopyTo               ( const Integer l, const Real1DArray *a, Integer *dL, Real1DArray *dA ) ;
static void Polynomial_Differentiate        ( Integer *dL, Real1DArray *dA ) ;
static void Polynomial_GenerateCoefficients ( const Integer l, const Real1DArray *factorial, Real1DArray *a ) ;

static Real SSBPModel_AngularPotential      ( const SSBPModel *self, const Integer particleIndex, const Real radius, const Integer2DArray *waterAtomIndices,
                                                                      const Real1DArray *radii, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
static Real SSBPModel_CavityPotential       ( const SSBPModel *self, const Integer particleIndex, const Real radius, const Selection *cavitySelection,
                                                                const Real1DArray *radii, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
static void SSBPModel_DetermineRadii        ( const SSBPModel *self, SSBPModelState *state ) ;

static Real SSBPModel_EmpiricalCorrection   ( const SSBPModel *self, const Integer particleIndex, const Real radius, const Integer2DArray *waterAtomIndices,
                                                                       const Real1DArray *radii, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
static Real SSBPModel_HardSphereRestriction ( const SSBPModel *self, const Integer particleIndex, const Real radius, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;

static void SSBPModel_KirkwoodQLMs  ( const Integer l, const Integer m, const Integer lA, const Integer dL, const Boolean doGradients,
                                    			    const Real1DArray *factorial, const Real1DArray *a, const Real1DArray *dA,
				    			   const Real1DArray *rlr3, const Real1DArray *tvar12, const Real1DArray *xr2,
				    			       const Real1DArray *xrl, const Real1DArray *yr2, const Real1DArray *yrl,
				    			     const Real1DArray *zr2, const Real1DArray *zxr3, const Real1DArray *zyr3,
                                    			     const Real2DArray *comr, const Real2DArray *comi, const Real2DArray *rQs,
				    					       const Real2DArray *tvar, Real *bOut, Real1DArray *qLMr,
				    						 Real1DArray *qLMi, Coordinates3 *dQLMi, Coordinates3 *dQLMr ) ;
static void SSBPModel_KirkwoodSetup ( const SSBPModel *self, const Boolean doGradients, SSBPModelState *state ) ;

/*==================================================================================================================================
! .  Public functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
SSBPModel *SSBPModel_Allocate ( void )
{
    SSBPModel *self = NULL ;
    self = ( SSBPModel * ) Memory_Allocate ( sizeof ( SSBPModel ) ) ;
    if ( self != NULL )
    {
        self->doAngularPotential      = DEFAULT_DOANGULARPOTENTIAL      ;
        self->doCavityPotential       = DEFAULT_DOCAVITYPOTENTIAL       ;
        self->doEmpiricalCorrection   = DEFAULT_DOEMPIRICALCORRECTION   ;
        self->doHardSphereRestriction = DEFAULT_DOHARDSPHERERESTRICTION ;
        self->doKirkwood              = DEFAULT_DOKIRKWOOD              ;
        self->fixCavityRadius         = DEFAULT_FIXCAVITYRADIUS         ;
        self->maximumL                = DEFAULT_MAXIMUML                ;
        self->cavityRadius            = DEFAULT_CAVITYRADIUS            ;
        self->cavityRadiusIncrement   = DEFAULT_CAVITYRADIUSINCREMENT   ;
        self->dielectricInside        = DEFAULT_DIELECTRICINSIDE        ;
        self->dielectricOutside       = DEFAULT_DIELECTRICOUTSIDE       ;
        self->empirical1              = DEFAULT_EMPIRICAL1              ;
        self->empirical2              = DEFAULT_EMPIRICAL2              ;
        self->kirkwoodRadiusIncrement = DEFAULT_KIRKWOODRADIUSINCREMENT ;
        self->pressure                = DEFAULT_PRESSURE                ;
        self->surfaceTension          = DEFAULT_SURFACETENSION          ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
SSBPModel *SSBPModel_Clone ( const SSBPModel *self )
{
    SSBPModel *new = NULL ;
    if ( self != NULL )
    {
        new = SSBPModel_Allocate ( ) ;
        new->doAngularPotential      = self->doAngularPotential      ;
        new->doCavityPotential       = self->doCavityPotential       ;
        new->doEmpiricalCorrection   = self->doEmpiricalCorrection   ;
        new->doHardSphereRestriction = self->doHardSphereRestriction ;
        new->doKirkwood              = self->doKirkwood              ;
        new->fixCavityRadius         = self->fixCavityRadius         ;
        new->maximumL                = self->maximumL                ;
        new->cavityRadius            = self->cavityRadius            ;
        new->cavityRadiusIncrement   = self->cavityRadiusIncrement   ;
        new->dielectricInside        = self->dielectricInside        ;
        new->dielectricOutside       = self->dielectricOutside       ;
        new->empirical1              = self->empirical1              ;
        new->empirical2              = self->empirical2              ;
        new->kirkwoodRadiusIncrement = self->kirkwoodRadiusIncrement ;
        new->pressure                = self->pressure                ;
        new->surfaceTension          = self->surfaceTension          ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SSBPModel_Deallocate ( SSBPModel **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! .  Energy terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SSBPModel_Energy ( const SSBPModel *self, const Boolean doGradients, SSBPModelState *state )
{
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        /* . Determine radii. */
        SSBPModel_DetermineRadii ( self, state ) ;

        /* . Various energy terms. */
        if ( self->doAngularPotential      ) state->eAngular             = SSBPModel_AngularPotential      ( self, state->particleIndex, state->radius, state->waterAtomIndices,
                                                                                                                        state->radii, state->coordinates3, state->gradients3 ) ;
        if ( self->doCavityPotential       ) state->eCavity              = SSBPModel_CavityPotential       ( self, state->particleIndex, state->radius, state->cavitySelection,
                                                                                                                       state->radii, state->coordinates3, state->gradients3 ) ;
        if ( self->doEmpiricalCorrection   ) state->eEmpiricalCorrection = SSBPModel_EmpiricalCorrection   ( self, state->particleIndex, state->radius, state->waterAtomIndices,
                                                                                                                        state->radii, state->coordinates3, state->gradients3 ) ;
        if ( self->doHardSphereRestriction ) state->eHardSphere          = SSBPModel_HardSphereRestriction ( self, state->particleIndex, state->radius, state->coordinates3, state->gradients3 ) ;
        if ( self->doKirkwood              ) state->eKirkwood            = SSBPModel_Kirkwood              ( self, doGradients, state ) ;

        /* . Total energy. */
        state->eTotal = state->eAngular + state->eCavity + state->eEmpiricalCorrection + state->eHardSphere + state->eKirkwood ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Kirkwood electrostatic terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SSBPModel_Kirkwood ( const SSBPModel *self, const Boolean doGradients, SSBPModelState *state )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( state != NULL ) && ( self->doKirkwood ) )
    {
        auto Boolean doMM, doQC ;
        auto Integer dL, i, l, lA, m, s, t ;
        auto Real    b, c, coefficient, coefficientQ, coefficient0, coefficient0Q, cX = 0.0e+00, cY = 0.0e+00, cZ = 0.0e+00,
	                                       dCFactor = 0.0e+00, deltaR, e, f, gX, gY, gZ, iFactor, newRadius, newRadius2,
				      	                    qLMi = 0.0e+00, qLMiQ, qLMr = 0.0e+00, qLMrQ, rFactor, x, y, z ;

# ifdef SSBPTESTQC
{
if ( state->qcAtoms != NULL )
{
    auto Integer i, index ;
    for ( i = 0 ; i < state->qcAtoms->natoms ; i++ )
    {
        index = state->qcAtoms->data[i].index ;
        state->mmAtoms->data[index].QACTIVE = False ;
    }
}
# endif

        /* . Initialization. */
	doMM          = ( MMAtomContainer_NumberOfActiveAtoms ( state->mmAtoms ) > 0 ) ;
	doQC          = ( QCAtomContainer_Size                ( state->qcAtoms ) > 0 ) ;
	c             = - 0.5e+00 * ( self->dielectricOutside - self->dielectricInside ) ;
        coefficient0  = c * UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE ;
	coefficient0Q = c / UNITS_LENGTH_ANGSTROMS_TO_BOHRS ;
        deltaR        = self->kirkwoodRadiusIncrement - state->qTotal * OnePt6 * exp ( -0.5e+00 * state->radius ) ;
        newRadius     = state->radius + deltaR ;
        newRadius2    = newRadius * newRadius ;
        rFactor       = newRadius ;
        if ( ! self->fixCavityRadius ) dCFactor = ( ( self->kirkwoodRadiusIncrement - deltaR ) * 0.5e+00 + 1.0e+00 ) /
	                                                                               ( newRadius * state->radius ) ;

        /* . Determine some useful arrays for the calculation of multipoles and gradients. */
        SSBPModel_KirkwoodSetup ( self, doGradients, state ) ;

        /* . Loop over l. */
        for ( l = 0 ; l <= self->maximumL ; l++ )
        {
            /* . Initialization. */
            lA = l ;
            Polynomial_GenerateCoefficients ( lA, state->factorial, state->a ) ;
            Polynomial_CopyTo               ( lA, state->a, &dL, state->dA ) ;
            Polynomial_Differentiate        ( &dL, state->dA ) ;
            rFactor    /= newRadius2 ;
            /* . The CHARMM implementation has an error here as it uses an integer division (l/(l+1)) instead of a real division. */
            c = rFactor / ( ( self->dielectricOutside + ( ( Real ) l / ( ( Real ) ( l+1 ) ) ) *
	                                  self->dielectricInside ) * ( ( Real ) ( 2*l+1 ) ) ) ;
            coefficient  = coefficient0  * c ;
	    coefficientQ = coefficient0Q * c ;
            if ( ! self->fixCavityRadius )
            {
                Coordinates3_GetRow ( state->coordinates3, state->particleIndex, x, y, z ) ;
                c  = ( Real ) ( 2*l+1 ) * dCFactor ;
                cX = c * x ;
                cY = c * y ;
                cZ = c * z ;
            }
            /* . Loop over m. */
            for ( m = 0 ; m <= l ; m++ )
            {
	        e       = 0.0e+00 ;
                iFactor = 1.0e+00 ;
                if ( m > 0 )
                {
                    iFactor = 2.0e+00 ;
                    Polynomial_CopyTo        ( dL, state->dA, &lA, state->a ) ;
                    Polynomial_Differentiate ( &dL, state->dA ) ;
                }

                /* . MM atoms. */
		if ( doMM )
		{
                    SSBPModel_KirkwoodQLMs ( l, m, lA, dL, doGradients, state->factorial, state->a, state->dA, state->rlr3, state->tvar12,
                                                        	  state->xr2, state->xrl, state->yr2, state->yrl, state->zr2, state->zxr3,
							        	   state->zyr3, state->comr, state->comi, state->rQs, state->tvar,
							                	  &b, state->uS, state->vS, state->dQLMi, state->dQLMr ) ;

                    /* . Accumulation. */
		    qLMr = Real1DArray_Sum ( state->uS ) ;
		    qLMi = Real1DArray_Sum ( state->vS ) ;
		    e   += ( b * coefficient * iFactor * ( qLMr * qLMr + qLMi * qLMi ) ) ;
		    if ( doGradients )
		    {
        		Coordinates3_Set            ( state->dQLMs, 0.0e+00 ) ;
        		Coordinates3_AddScaledArray ( state->dQLMs, qLMr, state->dQLMr, NULL ) ;
        		Coordinates3_AddScaledArray ( state->dQLMs, qLMi, state->dQLMi, NULL ) ;
                    }
		}

                /* . QC atoms. */
		/*
		! . For potentials have:
		! . qLMr = Mr + Qr, etc. so:
		! . E = b ( Mr^2 + Mi^2 ) + 2b ( MrQr + MiQi ) + b ( Qr^2 + Qi^2 ) ;
		! . Qr = Sum_s q_s u_s and Qi = Sum_s q_s v_s ;
		! . So get:
		! . qcmm term [ q_s	  ] is : 2b ( Mr u_s  + Mi v_s  )
		! . qcqc term [ q_s * q_t ] is :  b ( u_s u_t + v_s v_t )
                !
		! . For gradients have:
		! . dE/dXm = 2b ( dMr/dXm Qr + dMi/dXm Qi )
		! . dE/dXq = 2b ( dQr/dXq ( Qr + Mr ) + dQi/dXq ( Qi + Mi ) )
		! . Note that the dMrs and dQrs are already scaled by 2b.
		*/
		if ( doQC )
		{
                    SSBPModel_KirkwoodQLMs ( l, m, lA, dL, doGradients, state->factorial, state->a, state->dA, state->rlr3Q, state->tvar12Q,
                                                              state->xr2Q, state->xrlQ, state->yr2Q, state->yrlQ, state->zr2Q, state->zxr3Q,
							     		state->zyr3Q, state->comrQ, state->comiQ, state->rQsQ, state->tvarQ,
							     			&b, state->uSQ, state->vSQ, state->dQLMiQ, state->dQLMrQ ) ;

		    /* . Accumulation. */
		    /* . Gradients (with energy as check). */
		    if ( doGradients )
		    {
			qLMrQ = Real1DArray_Dot ( state->qcCharges, state->uSQ, NULL ) ;
			qLMiQ = Real1DArray_Dot ( state->qcCharges, state->vSQ, NULL ) ;
			/* . Energy. */
			e += ( b * coefficientQ * iFactor * UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE *
                	                                   ( 2.0e+00 * ( qLMr  * qLMrQ + qLMi  * qLMiQ ) +
							               ( qLMrQ * qLMrQ + qLMiQ * qLMiQ ) ) ) ;
                        /* . MM. */
			if ( doMM )
			{
            	            Coordinates3_AddScaledArray ( state->dQLMs, qLMrQ, state->dQLMr, NULL ) ;
        	            Coordinates3_AddScaledArray ( state->dQLMs, qLMiQ, state->dQLMi, NULL ) ;
			    qLMrQ += qLMr ;
			    qLMiQ += qLMi ;
			}
                        /* . QC. */
        	        Coordinates3_Set            ( state->dQLMsQ, 0.0e+00 ) ;
        	        Coordinates3_AddScaledArray ( state->dQLMsQ, qLMrQ, state->dQLMrQ, NULL ) ;
        	        Coordinates3_AddScaledArray ( state->dQLMsQ, qLMiQ, state->dQLMiQ, NULL ) ;
			Coordinates3_ScaleRows      ( state->dQLMsQ, state->qcCharges, NULL ) ;
		    }
		    /* . QC potentials. */
		    else
		    {
		        f = b * coefficientQ * iFactor ;
			Real1DArray_AddScaledArray ( state->qcmmPotentials, 2.0e+00 * f * qLMr, state->uSQ, NULL ) ;
			Real1DArray_AddScaledArray ( state->qcmmPotentials, 2.0e+00 * f * qLMi, state->vSQ, NULL ) ;
			for ( s = 0 ; s < state->uSQ->length ; s++ )
			{
        		    for ( t = 0 ; t <= s ; t++ )
			    {
				SymmetricMatrix_Item ( state->qcqcPotentials, s, t ) += f * ( Real1DArray_Item ( state->uSQ, s ) *
				                                                              Real1DArray_Item ( state->uSQ, t ) +
                                                                        		      Real1DArray_Item ( state->vSQ, s ) *
											      Real1DArray_Item ( state->vSQ, t ) ) ;
        		    }
			}
		    }
		}

                /* . Energy and gradients. */
		energy += e ;
                if ( doGradients )
		{
		    if ( doMM )
		    {
                	Coordinates3_Scale ( state->dQLMs, coefficient * iFactor ) ;
                	for ( i = s = 0 ; i < state->mmAtoms->natoms ; i++ )
                	{
		            if ( state->mmAtoms->data[i].QACTIVE )
			    {
                        	Coordinates3_GetRow   ( state->dQLMs     , s, gX, gY, gZ ) ;
                        	Coordinates3_IncrementRow ( state->gradients3, i, gX, gY, gZ ) ;
				s += 1 ;
			    }
                	}
		    }
		    if ( doQC )
		    {
                	Coordinates3_AddScaledArray ( state->qcGradients3, coefficientQ * iFactor, state->dQLMsQ, NULL ) ;
		    }
                    if ( ! self->fixCavityRadius )
		    {
		        Coordinates3_IncrementRow ( state->gradients3, state->particleIndex, -e*cX, -e*cY, -e*cZ ) ;
		    }
		}
            }
        }

        /* . Finish up the QC terms. */
	if ( doQC )
	{
	    /* . Gradients. */
	    if ( doGradients )
	    {
                Coordinates3_Scale ( state->qcGradients3, UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE ) ;
    	        for ( i = 0 ; i < state->qcAtoms->natoms ; i++ )
	        {
	            Coordinates3_GetRow ( state->qcGradients3, i, gX, gY, gZ ) ;
           	    QCAtomContainer_SetAtomGradients3 ( state->qcAtoms, i, state->coordinates3, gX, gY, gZ, state->gradients3 ) ;
		}
	    }
	}

# ifdef SSBPTESTQC
{
/* . Finish the energy calculation. */
if ( doQC && ( ! doGradients ) )
{
    auto Real1DArray *work ;
    work = Real1DArray_Allocate    ( state->qcCharges->length, NULL ) ;
    SymmetricMatrix_VectorMultiply ( state->qcqcPotentials, state->qcCharges, work, NULL ) ;
    Real1DArray_Add                ( work, state->qcmmPotentials, NULL ) ;
    energy += ( Real1DArray_Dot    ( state->qcCharges, work, NULL ) * UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE ) ;
    Real1DArray_Deallocate         ( &work ) ;
}
if ( doQC )
{
    auto Integer i, index ;
    for ( i = 0 ; i < state->qcAtoms->natoms ; i++ )
    {
        index = state->qcAtoms->data[i].index ;
        state->mmAtoms->data[index].QACTIVE = False ;
    }
}
# endif

        /* . Finish up. */
	state->eKirkwoodCheck = energy ;
    }
    return energy ;
}

/*==================================================================================================================================
! .  Private functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy a Legendre polynomial.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Polynomial_CopyTo ( const Integer l, const Real1DArray *a, Integer *dL, Real1DArray *dA )
{
    (*dL) = l ;
    Real1DArray_CopyTo ( a, dA, NULL ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Differentiate a Legendre polynomial of degree l.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Polynomial_Differentiate ( Integer *dL, Real1DArray *dA )
{
    auto Integer l ;
    l = (*dL) ;
    if ( l > 0 )
    {
        auto Integer i ;
        (*dL) = l - 1 ;
        for ( i = 0 ; i <= l ; i++ ) Real1DArray_Item ( dA, i ) *= ( ( Real ) ( l - i ) ) ;
    }
    else
    {
        (*dL) = 0 ;
        Real1DArray_Set ( dA, 0.0e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate the coefficients of a Legendre polynomial of degree l.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . This uses Rodigues representation in which a[2*k] is the coefficient for x^(l-2*k), k = [0,floor (l/2)]. */
/* . a must be of at least length l+1. */

static void Polynomial_GenerateCoefficients ( const Integer l, const Real1DArray *factorial, Real1DArray *a )
{
    Integer k, kUpper ;
    Real1DArray_Set ( a, 0.0e+00 ) ;
    if ( l % 2 == 0 ) kUpper =   l       / 2 ;
    else              kUpper = ( l - 1 ) / 2 ;
    for ( k = 0 ; k <= kUpper ; k++ ) Real1DArray_Item ( a, 2*k ) = ( Real1DArray_Item ( factorial, 2*(l-k) ) * pow ( -1.0e+00, ( Real ) k ) ) /
                                                                    ( Real1DArray_Item ( factorial, k ) * Real1DArray_Item ( factorial, l-k ) * Real1DArray_Item ( factorial, l-2*k ) ) ;
    Real1DArray_Scale ( a, 1.0e+00 / pow ( 2.0e+00, ( Real ) l ) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! .  Angular potential term.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . This contribution was developed empirically to make more isotropic the orientational distribution function
! . of the waters near the outer shell. Warning: it works only for a 3 site water models such as tip3p.
*/

static Real SSBPModel_AngularPotential ( const SSBPModel *self, const Integer particleIndex, const Real radius, const Integer2DArray *waterAtomIndices,
                                                                   const Real1DArray *radii, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( waterAtomIndices != NULL )
    {
        Boolean doGradients ;
        Integer c, H1, H2, j, O, s ;
        Real    c1, c2, c3, c4, dE, dot, dR, e, f, fC, g, gTotal, r, rHO2, rMaximum, t1, t2, t3, x, y, z ;
        Real    dHO[3], gRH[3][3], gRO[3], rH[3][3], rO[3] ;

        /* . Initialization. */
        doGradients = ( gradients3 != NULL );
        gTotal      = 0.0e+00 ;
        if ( self->fixCavityRadius ) rMaximum = radius ;
        else                         rMaximum = Real1DArray_Item ( radii, particleIndex ) ;

        /* . Loop over atoms. */
        for ( s = 0 ; s < waterAtomIndices->length0 ; s++ )
        {
            O  = Integer2DArray_Item ( waterAtomIndices, s, 0 ) ;
            r  = Real1DArray_Item ( radii, O ) ;
            dR = r - rMaximum - DRHA ;
            if ( dR > 0.0e+00 )
            {
                /* . Coordinates. */
                H1 = Integer2DArray_Item ( waterAtomIndices, s, 1 ) ;
                H2 = Integer2DArray_Item ( waterAtomIndices, s, 2 ) ;
                Coordinates3_GetRow ( coordinates3, O , rO[0]   , rO[1]   , rO[2]    ) ;
                Coordinates3_GetRow ( coordinates3, H1, rH[0][0], rH[0][1], rH[0][2] ) ;
                Coordinates3_GetRow ( coordinates3, H2, rH[1][0], rH[1][1], rH[1][2] ) ;

                /* . Energy terms - loop over hydrogens. */
                e = 0.0e+00 ;
                for ( c = 0 ; c < 3 ; c++ ) gRO[c] = 0.0e+00 ;

                for ( j = 0 ; j < 2 ; j++ )
                {
                    dot = 0.0e+00 ;
                    for ( c = 0 ; c < 3 ; c++ )
                    {
                        dHO[c] = rH[j][c] - rO[c] ;
                        dot += dHO[c] * rO[c] ;
                    }
                    rHO2 = dHO[0]*dHO[0] + dHO[1]*dHO[1] + dHO[2]*dHO[2] ;
                    t1 = 1.0e+00 / ( r * sqrt ( rHO2 ) ) ;
                    t2 = dot / rHO2 ;
                    t3 = dot / ( r * r ) ;
                    c1 = dot * t1 ;
                    c2 = c1 * c1 ;
                    c3 = c2 * c1 ;
                    c4 = c3 * c1 ;
                    dE = 2.0e+00 * ( 4.0e+00 * CAngular1 * c3 + 3.0e+00 * CAngular2 * c2 + 2.0e+00 * CAngular3 * c1 + CAngular4 ) ;
                    e += 2.0e+00 * (           CAngular1 * c4 +           CAngular2 * c3 +           CAngular3 * c2 + CAngular4 * c1 + CAngular5 ) ;
                    for ( c = 0 ; c < 3 ; c++ )
                    {
                        gRH[j][c] = dE * t1 * ( -t2 * dHO[c] + rO[c] ) ;
                        gRO[c]   += dE * t1 * ( -rO[c] * ( t3 + 1.0e+00 ) + dHO[c] * ( 1.0e+00 + t2 ) ) ;
                    }
                }
                g  = dR * e ;
                f  = 0.5e+00 *  g * dR ;
                fC = 0.5e+00 * dR * dR ;

                /* . Energy and gradients. */
                energy += f ;
                if ( doGradients )
                {
                    gTotal += g ;
                    g      /= r ;
                    Coordinates3_IncrementRow ( gradients3, O, g*rO[0]+fC*gRO[0]   , g*rO[1]+fC*gRO[1]   , g*rO[2]+fC*gRO[2]    ) ;
                    Coordinates3_IncrementRow ( gradients3, H1,        fC*gRH[0][0],         fC*gRH[0][1],         fC*gRH[0][2] ) ;
                    Coordinates3_IncrementRow ( gradients3, H2,        fC*gRH[1][0],         fC*gRH[1][1],         fC*gRH[1][2] ) ;
                }
            }
        }
        if ( ( ! self->fixCavityRadius ) && doGradients )
        {
            g = - gTotal / radius ;
            Coordinates3_GetRow     ( coordinates3, particleIndex,   x,   y,   z ) ;
            Coordinates3_IncrementRow ( gradients3  , particleIndex, g*x, g*y, g*z ) ;
        }
    }
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! .  Cavity potential term.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . This constribution was calculated from the cavity potential of mean force between a LJ TIP3P-like oxygen atom and a large
! . hard-sphere of radius RMAX. The polynomial approximation has been fitted and yields a boundary potential that is very similar
! . to the SBOUND method of C. L. Brooks.
*/

static Real SSBPModel_CavityPotential ( const SSBPModel *self, const Integer particleIndex, const Real radius, const Selection *cavitySelection,
                                                            const Real1DArray *radii, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( cavitySelection != NULL )
    {
        Boolean doGradients ;
        Integer i, s ;
        Real    a, a2, a3, a4, d, f, g, gTotal, r, rFactor, x, y, z ;

        /* . Initialization. */
        doGradients = ( gradients3 != NULL ) ;
        gTotal      = 0.0e+00 ;
        rFactor     = radius + self->cavityRadiusIncrement - CavityDecrement ;

        /* . Loop over atoms. */
        for ( s = 0 ; s < cavitySelection->nindices ; s++ )
        {
            i = cavitySelection->indices[s] ;
            r = Real1DArray_Item ( radii, i ) ;
            a = r - rFactor ;
            if ( a < BLower ) { f = BCavity7 ; g = 0.0e+00 ; }
            else
            {
                a2 = a * a ;
                if ( a > 0.0e+00 )
                {
                    f = BCavity5 + BCavity6 * a2 ;
                    g = 2.0e+00 * a * BCavity6   ;
                }
                else
                {
                    d = ( 1.0e+00 + a2 / BCavity1 ) ;
                    f = BCavity2 / d + BCavity3 * a2 + BCavity4 ;
                    g = 2.0e+00 * a * ( BCavity3 - BCavity2 / ( BCavity1 * d * d ) ) ;
                }
            }
            energy += f ;
            if ( doGradients )
            {
                gTotal += g ;
                g      /= Real1DArray_Item ( radii, i ) ;
                Coordinates3_GetRow     ( coordinates3, i,   x,   y,   z ) ;
                Coordinates3_IncrementRow ( gradients3  , i, g*x, g*y, g*z ) ;
            }
        }

        /* . Function-axis shift of the approximated cavity potential. */
        if ( ! self->fixCavityRadius )
        {
            a = rFactor ;
            if ( a > ACavity7 ) { f = ACavity6 ; g = 0.0e+00 ; }
            else
            {
                a2 = a  * a ;
                a3 = a2 * a ;
                a4 = a3 * a ;
                f  = ACavity5 + ACavity1 * a +           ACavity2 * a2 +           ACavity3 * a3 +           ACavity4 * a4 ;
                g  =            ACavity1     + 2.0e+00 * ACavity2 * a  + 3.0e+00 * ACavity3 * a2 + 4.0e+00 * ACavity4 * a3 ;
            }
            f += ACavity8 ;
            f *= ( Real ) ( cavitySelection->nindices ) ;

            /* . Energy and gradients. */
            energy += f ;
            if ( doGradients )
            {
                g *= ( Real ) ( cavitySelection->nindices ) ;
                g -= gTotal ;
                g /= radius ;
                Coordinates3_GetRow     ( coordinates3, particleIndex,   x,   y,   z ) ;
                Coordinates3_IncrementRow ( gradients3  , particleIndex, g*x, g*y, g*z ) ;
            }
        }
    }
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Define radii.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void SSBPModel_DetermineRadii ( const SSBPModel *self, SSBPModelState *state )
{
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        auto Integer i, p, s ;
        auto Real    m, r, x, y, z ;

        /* . Calculate radii. */
        for ( i = 0 ; i < state->radii->length ; i++ )
        {
            Coordinates3_GetRow ( state->coordinates3, i, x, y, z ) ;
            Real1DArray_Item ( state->radii, i ) = sqrt ( x*x + y*y + z*z ) ;
        }

        /* . Find maximum radius. */
        if ( state->radiusSelection != NULL )
        {
            p = -1 ;
            if ( self->fixCavityRadius ) m = self->cavityRadius ;
            else                         m = 0.0e+00 ;
            for ( s = 0 ; s < state->radiusSelection->nindices ; s++ )
            {
                i = state->radiusSelection->indices[s] ;
                r = Real1DArray_Item ( state->radii, i ) ;
                if ( r > m )
                {
                    if ( self->fixCavityRadius ) state->atomOutsideCavity = True ;
                    else { p = i ; m = r ; }
                }
            }

            /* . Finish up. */
            state->particleIndex = p ;
            state->radius        = m ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! .  Empirical correction term.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real SSBPModel_EmpiricalCorrection ( const SSBPModel *self, const Integer particleIndex, const Real radius, const Integer2DArray *waterAtomIndices,
                                                                      const Real1DArray *radii, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( waterAtomIndices != NULL )
    {
        Boolean doGradients ;
        Integer i, s ;
        Real    f, g, h, p1, p2, r, radiusSquared, x, y, z ;

        /* . Initialization. */
        doGradients   = ( gradients3 != NULL ) ;
        radiusSquared = radius * radius ;
        p1 = self->empirical1 * EightPt88             / radius ;
        p2 = self->empirical2 * EightPt88 * EightPt88 / radiusSquared ;

        /* . Loop over atoms. */
        h = 0.0e+00 ;
        for ( s = 0 ; s < waterAtomIndices->length0 ; s++ )
        {
            i = Integer2DArray_Item ( waterAtomIndices, s, 0 ) ;
            r = Real1DArray_Item ( radii, i ) ;
            f = p1 * exp ( - p2 * r * r ) ;
            energy += f ;
            if ( doGradients )
            {
                g  = - 2.0e+00 * p2 * f  ;
                h += ( f * r * r ) ;
                Coordinates3_GetRow     ( coordinates3, i,   x,   y,   z ) ;
                Coordinates3_IncrementRow ( gradients3  , i, g*x, g*y, g*z ) ;
            }
        }
        if ( ( ! self->fixCavityRadius ) && doGradients )
        {
            g = ( 2.0e+00 * p2 * h - energy ) / radiusSquared ;
            Coordinates3_GetRow     ( coordinates3, particleIndex,   x,   y,   z ) ;
            Coordinates3_IncrementRow ( gradients3  , particleIndex, g*x, g*y, g*z ) ;
        }
    }
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! .  Hard sphere restriction term ( P * V + Surface Tension * Surface Area ).
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real SSBPModel_HardSphereRestriction ( const SSBPModel *self, const Integer particleIndex, const Real radius, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    Real g, pvTerm, radiusSquared, stTerm, x, y, z ;

    /* . Energy. */
    radiusSquared = radius * radius ;
    pvTerm        = self->pressure       * ( 4.0e+00 * M_PI * radiusSquared * radius / 3.0e+00 ) * UNITS_PRESSURE_ATMOSPHERES_TO_KILOJOULES_PER_MOLE ;
    stTerm        = self->surfaceTension * ( 4.0e+00 * M_PI * radiusSquared ) ;
    energy        = pvTerm + stTerm ;

    /* . Gradients. */
    if ( ( ! self->fixCavityRadius ) && ( gradients3 != NULL ) )
    {
        g = ( 3.0e+00 * pvTerm + 2.0e+00 * stTerm ) / ( radiusSquared ) ;
        Coordinates3_GetRow     ( coordinates3, particleIndex,   x,   y,   z ) ;
        Coordinates3_IncrementRow ( gradients3  , particleIndex, g*x, g*y, g*z ) ;
    }

    /* . Finish up. */
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate quantities needed for the |Qlm|^2.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void SSBPModel_KirkwoodQLMs ( const Integer l, const Integer m, const Integer lA, const Integer dL, const Boolean doGradients,
                                                           const Real1DArray *factorial, const Real1DArray *a, const Real1DArray *dA,
							  const Real1DArray *rlr3, const Real1DArray *tvar12, const Real1DArray *xr2,
							      const Real1DArray *xrl, const Real1DArray *yr2, const Real1DArray *yrl,
							    const Real1DArray *zr2, const Real1DArray *zxr3, const Real1DArray *zyr3,
                                                            const Real2DArray *comr, const Real2DArray *comi, const Real2DArray *rQs,
						    	                      const Real2DArray *tvar, Real *bOut, Real1DArray *qLMr,
								   	        Real1DArray *qLMi, Coordinates3 *dQLMi, Coordinates3 *dQLMr )
{
    if ( qLMr != NULL )
    {
	auto Integer i, n ;
	auto Real    b, b1, bi1, br1, d, d1, d2, qLMiX, qLMiY, qLMiZ,  qLMrX, qLMrY, qLMrZ, rQ,
                	 yLMi1, yLMi2, yLMi3, yLMic, yLMit, yLMr1, yLMr2, yLMr3, yLMrc, yLMrt ;

	/* . Energy and derivative terms. */
	b = ( ( Real ) ( 2*l+1 ) ) * Real1DArray_Item ( factorial, l-m ) / Real1DArray_Item ( factorial, l+m ) ;
	for ( i = 0 ; i < qLMr->length; i++ )
	{
            d = 0.0e+00 ;
            for ( n = 0 ; n <= lA ; n++ ) d += ( Real1DArray_Item ( a, n ) * Real2DArray_Item ( tvar, lA-n, i ) ) ;
            yLMit = Real2DArray_Item ( comi, m, i ) ;
            yLMrt = Real2DArray_Item ( comr, m, i ) ;
            yLMic = d * yLMit ;
            yLMrc = d * yLMrt ;
            rQ    = Real2DArray_Item ( rQs, l, i ) ;
            Real1DArray_Item ( qLMi, i ) = rQ * yLMic ;
            Real1DArray_Item ( qLMr, i ) = rQ * yLMrc ;
            if ( doGradients )
            {
        	d1 = 0.0e+00 ;
        	for ( n = 0 ; n <= dL ; n++ ) d1 += ( Real1DArray_Item ( dA, n ) * Real2DArray_Item ( tvar, dL-n, i ) ) ;
        	d2 = d1 - d * Real1DArray_Item ( tvar12, i ) * m ;
        	b1 = d2 * Real1DArray_Item ( rlr3, i ) ;
        	yLMr3 = yLMrt * b1 ;
        	yLMi3 = yLMit * b1 ;
        	br1   = -d2 * Real1DArray_Item ( zyr3, i ) ;
        	bi1   =   m * Real1DArray_Item ( xrl , i ) ;
        	yLMr2 = yLMrt * br1 - yLMic * bi1 ;
        	yLMi2 = yLMit * br1 + yLMrc * bi1 ;
        	br1   = -d2 * Real1DArray_Item ( zxr3, i ) ;
        	bi1   = - m * Real1DArray_Item ( yrl , i ) ;
        	yLMr1 = yLMrt * br1 - yLMic * bi1 ;
        	yLMi1 = yLMit * br1 + yLMrc * bi1 ;
        	rQ   *= 2.0e+00 * b ;
        	b1    = l * Real1DArray_Item ( xr2, i ) ;
        	qLMrX = rQ * ( b1*yLMrc + yLMr1 ) ;
        	qLMiX = rQ * ( b1*yLMic + yLMi1 ) ;
        	b1    = l * Real1DArray_Item ( yr2, i ) ;
        	qLMrY = rQ * ( b1*yLMrc + yLMr2 ) ;
        	qLMiY = rQ * ( b1*yLMic + yLMi2 ) ;
        	b1    = l * Real1DArray_Item ( zr2, i ) ;
        	qLMrZ = rQ * ( b1*yLMrc + yLMr3 ) ;
        	qLMiZ = rQ * ( b1*yLMic + yLMi3 ) ;
        	Coordinates3_SetRow ( dQLMr, i, qLMrX, qLMrY, qLMrZ ) ;
        	Coordinates3_SetRow ( dQLMi, i, qLMiX, qLMiY, qLMiZ ) ;
            }
	}
	(*bOut)= b ;
    }
    else (*bOut) = 0.0e+00 ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Kirkwood setup of temporary arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void SSBPModel_KirkwoodSetup ( const SSBPModel *self, const Boolean doGradients, SSBPModelState *state )
{
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        auto Integer i, l, s ;
        auto Real    ar, ai, bi, br, r, r2, r2XY, r3, t, u, x, y, z ;

        /* . MM atoms. */
	if ( MMAtomContainer_NumberOfActiveAtoms ( state->mmAtoms ) > 0 )
	{
            for ( i = s = 0 ; i < state->mmAtoms->natoms ; i++ )
            {
		if ( state->mmAtoms->data[i].QACTIVE )
		{
        	    r = Real1DArray_Item ( state->radii, i ) ;
        	    Coordinates3_GetRow  ( state->coordinates3, i, x, y, z ) ;
        	    u = x ;
        	    if ( ( u < RSmall ) && ( y < RSmall ) ) { u += RSmall ; r = sqrt ( u*u + y*y + z*z ) ; }
        	    t = z / r ;
        	    Real2DArray_Item ( state->comr, 0, s ) = 1.0e+00 ;
        	    Real2DArray_Item ( state->comi, 0, s ) = 0.0e+00 ;
        	    Real2DArray_Item ( state->rQs , 0, s ) = state->mmAtoms->data[i].charge ;
        	    Real2DArray_Item ( state->tvar, 0, s ) = 1.0e+00 ;
        	    if ( self->maximumL > 0 )
        	    {
                	br = u / r ; bi = y / r ;
                	for ( l = 1 ; l <= self->maximumL ; l++ )
                	{
                	    ar = Real2DArray_Item ( state->comr, l-1, s ) * br - Real2DArray_Item ( state->comi, l-1, s ) * bi ;
                	    ai = Real2DArray_Item ( state->comr, l-1, s ) * bi + Real2DArray_Item ( state->comi, l-1, s ) * br ;
                	    Real2DArray_Item ( state->comr, l, s ) = ar ;
                	    Real2DArray_Item ( state->comi, l, s ) = ai ;
                	    Real2DArray_Item ( state->rQs , l, s ) = Real2DArray_Item ( state->rQs , l-1, s ) * r ;
                	    Real2DArray_Item ( state->tvar, l, s ) = Real2DArray_Item ( state->tvar, l-1, s ) * t  ;
                	}
        	    }
		    if ( doGradients )
		    {
			r2XY = u  * u + y * y ;
        		r2   = r  * r ;
        		r3   = r2 * r ;
        		Real1DArray_Item ( state->rlr3  , s ) = r2XY  / r3 ;
        		Real1DArray_Item ( state->zyr3  , s ) = z * y / r3 ;
        		Real1DArray_Item ( state->zxr3  , s ) = z * u / r3 ;
        		Real1DArray_Item ( state->xr2   , s ) = u / r2 ;
        		Real1DArray_Item ( state->yr2   , s ) = y / r2 ;
        		Real1DArray_Item ( state->zr2   , s ) = z / r2 ;
        		Real1DArray_Item ( state->xrl   , s ) = u / r2XY ;
        		Real1DArray_Item ( state->yrl   , s ) = y / r2XY ;
        		Real1DArray_Item ( state->tvar12, s ) = t / ( 1.0e+00 - t*t ) ;
		    }
		    s += 1 ;
		}
	    }
	}

        /* . QC atoms. */
	if ( QCAtomContainer_Size ( state->qcAtoms ) > 0 )
	{
            for ( s = 0 ; s < state->qcAtoms->natoms ; s++ )
            {
                QCAtomContainer_GetAtomCoordinates3 ( state->qcAtoms, s, state->coordinates3, &x, &y, &z ) ; /* . In Angstroms. */
		/* . Always recalculate radius as boundary atoms may be present. */
        	u = x ;
        	if ( ( u < RSmall ) && ( y < RSmall ) ) u += RSmall ;
		r = sqrt ( u*u + y*y + z*z ) ;
        	t = z / r ;
        	Real2DArray_Item ( state->comrQ, 0, s ) = 1.0e+00 ;
        	Real2DArray_Item ( state->comiQ, 0, s ) = 0.0e+00 ;
        	Real2DArray_Item ( state->rQsQ , 0, s ) = 1.0e+00 ;
        	Real2DArray_Item ( state->tvarQ, 0, s ) = 1.0e+00 ;
        	if ( self->maximumL > 0 )
        	{
                    br = u / r ; bi = y / r ;
                    for ( l = 1 ; l <= self->maximumL ; l++ )
                    {
                	ar = Real2DArray_Item ( state->comrQ, l-1, s ) * br - Real2DArray_Item ( state->comiQ, l-1, s ) * bi ;
                	ai = Real2DArray_Item ( state->comrQ, l-1, s ) * bi + Real2DArray_Item ( state->comiQ, l-1, s ) * br ;
                	Real2DArray_Item ( state->comrQ, l, s ) = ar ;
                	Real2DArray_Item ( state->comiQ, l, s ) = ai ;
                	Real2DArray_Item ( state->rQsQ , l, s ) = Real2DArray_Item ( state->rQsQ , l-1, s ) * r ;
                	Real2DArray_Item ( state->tvarQ, l, s ) = Real2DArray_Item ( state->tvarQ, l-1, s ) * t  ;
                    }
        	}
		if ( doGradients )
		{
        	    r2XY = u  * u + y * y ;
        	    r2   = r  * r ;
        	    r3   = r2 * r ;
        	    Real1DArray_Item ( state->rlr3Q  , s ) = r2XY  / r3 ;
        	    Real1DArray_Item ( state->zyr3Q  , s ) = z * y / r3 ;
        	    Real1DArray_Item ( state->zxr3Q  , s ) = z * u / r3 ;
        	    Real1DArray_Item ( state->xr2Q   , s ) = u / r2 ;
        	    Real1DArray_Item ( state->yr2Q   , s ) = y / r2 ;
        	    Real1DArray_Item ( state->zr2Q   , s ) = z / r2 ;
        	    Real1DArray_Item ( state->xrlQ   , s ) = u / r2XY ;
        	    Real1DArray_Item ( state->yrlQ   , s ) = y / r2XY ;
        	    Real1DArray_Item ( state->tvar12Q, s ) = t / ( 1.0e+00 - t*t ) ;
                }
            }
        }
    }
}
