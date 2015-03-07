/*------------------------------------------------------------------------------
! . File      : ChargeConstraintContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _CHARGECONSTRAINTCONTAINER
# define _CHARGECONSTRAINTCONTAINER

# include "Boolean.h"
# include "Integer.h"
# include "Integer1DArray.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Selection.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The charge constraint type. */
typedef struct {
    Integer         numberOfCharges ;
    Integer         numberOfSpins   ;
    Real            target          ;
    Integer1DArray *chargeIndices   ;
    Integer1DArray *spinIndices     ;
    Real1DArray    *chargeWeights   ;
    Real1DArray    *spinWeights     ;
} ChargeConstraint ;

/* . The grid type. */
typedef struct {
    Boolean            hasCharges          ;
    Boolean            hasSpins            ;
    Integer            highestIndex        ;
    Integer            numberOfConstraints ;
    ChargeConstraint **constraints         ;
} ChargeConstraintContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern ChargeConstraintContainer *ChargeConstraintContainer_Allocate     ( const Integer numberOfConstraints, Status *status ) ;
extern ChargeConstraintContainer *ChargeConstraintContainer_Clone        ( const ChargeConstraintContainer  *self, Status *status ) ;
extern void                       ChargeConstraintContainer_Deallocate   (       ChargeConstraintContainer **self ) ;
extern void                       ChargeConstraintContainer_Deviations   ( const ChargeConstraintContainer  *self       ,
                                                                           const Real1DArray                *charges    ,
                                                                           const Real1DArray                *spins      ,
                                                                                 Real1DArray                *deviations ,
                                                                                 Status                     *status     ) ;
extern Integer                    ChargeConstraintContainer_HighestIndex ( const ChargeConstraintContainer  *self ) ;
extern ChargeConstraintContainer *ChargeConstraintContainer_Prune        ( const ChargeConstraintContainer  *self, Selection *selection, Status *status ) ;
extern void                       ChargeConstraintContainer_SetItem      (       ChargeConstraintContainer  *self, const Integer index, ChargeConstraint **constraint, Status *status ) ;
extern Integer                    ChargeConstraintContainer_Size         ( const ChargeConstraintContainer  *self ) ;

extern ChargeConstraint *ChargeConstraint_Allocate   ( const Integer numberOfCharges, const Integer numberOfSpins, Status *status ) ;
extern ChargeConstraint *ChargeConstraint_Clone      ( const ChargeConstraint  *self, Status *status ) ;
extern void              ChargeConstraint_Deallocate (       ChargeConstraint **self ) ;
extern ChargeConstraint *ChargeConstraint_Prune      ( const ChargeConstraint  *self, const Selection *selection, Status *status ) ;

# endif
