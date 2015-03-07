/*------------------------------------------------------------------------------
! . File      : SSBPModel.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SSBPMODEL
# define _SSBPMODEL

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "SSBPModelState.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The model type. */
typedef struct {
    Boolean doAngularPotential      ;
    Boolean doCavityPotential       ;
    Boolean doEmpiricalCorrection   ;
    Boolean doHardSphereRestriction ;
    Boolean doKirkwood              ;
    Boolean fixCavityRadius         ;
    Integer maximumL                ;
    Real    cavityRadius            ;
    Real    cavityRadiusIncrement   ;
    Real    dielectricInside        ;
    Real    dielectricOutside       ;
    Real    empirical1              ;
    Real    empirical2              ;
    Real    kirkwoodRadiusIncrement ;
    Real    pressure                ;
    Real    surfaceTension          ;
} SSBPModel ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern SSBPModel *SSBPModel_Allocate   ( void ) ;
extern SSBPModel *SSBPModel_Clone      ( const SSBPModel  *self ) ;
extern void       SSBPModel_Deallocate (       SSBPModel **self ) ;
extern void       SSBPModel_Energy     ( const SSBPModel  *self, const Boolean doGradients, SSBPModelState *state ) ;
extern Real       SSBPModel_Kirkwood   ( const SSBPModel  *self, const Boolean doGradients, SSBPModelState *state ) ;

# endif
