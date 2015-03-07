/*------------------------------------------------------------------------------
! . File      : SSBPModelState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SSBPMODELSTATE
# define _SSBPMODELSTATE

# include "Boolean.h"
# include "Coordinates3.h"
# include "Integer.h"
# include "Integer2DArray.h"
# include "MMAtomContainer.h"
# include "QCAtomContainer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Selection.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The model type. */
typedef struct {
    Boolean	     atomOutsideCavity    ;
    Integer	     numberOfAtoms	  ;
    Integer	     particleIndex	  ;
    Real  	     eAngular		  ;
    Real  	     eCavity		  ;
    Real  	     eEmpiricalCorrection ;
    Real  	     eHardSphere	  ;
    Real  	     eKirkwood  	  ;
    Real             eKirkwoodCheck       ;
    Real  	     eTotal		  ;
    Real  	     qTotal		  ;
    Real  	     radius		  ;
    /* . Aliases. */
    Coordinates3    *coordinates3	  ;
    Coordinates3    *gradients3 	  ;
    MMAtomContainer *mmAtoms  	          ;
    QCAtomContainer *qcAtoms              ;
    Real1DArray     *qcCharges            ;
    Real1DArray     *qcmmPotentials       ;
    SymmetricMatrix *qcqcPotentials       ;
    /* . Arrays. */
    Integer2DArray  *waterAtomIndices	  ;
    Real1DArray     *radii		  ;
    Selection	    *cavitySelection	  ;
    Selection	    *radiusSelection	  ;
    /* . Kirkwood - MM. */
    Coordinates3    *dQLMi		  ;
    Coordinates3    *dQLMr		  ;
    Coordinates3    *dQLMs		  ;
    Real1DArray     *a  		  ;
    Real1DArray     *dA 		  ;
    Real1DArray     *factorial  	  ;
    Real1DArray     *tvar12		  ;
    Real1DArray     *rlr3		  ;
    Real1DArray     *uS   		  ;
    Real1DArray     *vS 		  ;
    Real1DArray     *xr2		  ;
    Real1DArray     *xrl		  ;
    Real1DArray     *yr2		  ;
    Real1DArray     *yrl		  ;
    Real1DArray     *zr2		  ;
    Real1DArray     *zxr3		  ;
    Real1DArray     *zyr3		  ;
    Real2DArray     *comr		  ;
    Real2DArray     *comi		  ;
    Real2DArray     *rQs		  ;
    Real2DArray     *tvar		  ;
    /* . Kirkwood - QC. */
    Coordinates3    *dQLMiQ		  ;
    Coordinates3    *dQLMrQ		  ;
    Coordinates3    *dQLMsQ		  ;
    Coordinates3    *qcGradients3         ;
    Real1DArray     *tvar12Q		  ;
    Real1DArray     *rlr3Q		  ;
    Real1DArray     *uSQ   		  ;
    Real1DArray     *vSQ 		  ;
    Real1DArray     *xr2Q		  ;
    Real1DArray     *xrlQ		  ;
    Real1DArray     *yr2Q		  ;
    Real1DArray     *yrlQ		  ;
    Real1DArray     *zr2Q		  ;
    Real1DArray     *zxr3Q		  ;
    Real1DArray     *zyr3Q		  ;
    Real2DArray     *comrQ		  ;
    Real2DArray     *comiQ		  ;
    Real2DArray     *rQsQ		  ;
    Real2DArray     *tvarQ		  ;
} SSBPModelState ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern SSBPModelState *SSBPModelState_Allocate   ( void ) ;
extern void            SSBPModelState_Deallocate ( SSBPModelState **self ) ;
extern void            SSBPModelState_Initialize ( SSBPModelState  *self, Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern SSBPModelState *SSBPModelState_Setup      ( const Boolean doKirkwood, const Integer maximumL, MMAtomContainer *mmAtoms,
                                                                         QCAtomContainer *qcAtoms, Selection *cavitySelection,
                                                                 Selection *radiusSelection, Integer2DArray *waterAtomIndices,
                                                                          Real1DArray *qcCharges, Real1DArray *qcmmPotentials,
                                                                             SymmetricMatrix *qcqcPotentials, Status *status ) ;

# endif
