/*------------------------------------------------------------------------------
! . File      : QCMMInteractionState.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCMMINTERACTIONSTATE
# define _QCMMINTERACTIONSTATE

# include "Boolean.h"
# include "Integer.h"
# include "Real1DArray.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The QC/MM interaction state type. */
typedef struct {
    Real1DArray     *qcCharges      ;
    Real1DArray     *qcmmPotentials ;
    SymmetricMatrix *qcqcPotentials ;
} QCMMInteractionState ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern QCMMInteractionState *QCMMInteractionState_Allocate               ( const Integer n, const Boolean includeQCQC, Status *status ) ;
extern void                  QCMMInteractionState_AllocateQCQCPotentials ( QCMMInteractionState  *self, Status *status ) ;
extern void                  QCMMInteractionState_Deallocate             ( QCMMInteractionState **self ) ;
extern void                  QCMMInteractionState_Initialize             ( QCMMInteractionState  *self ) ;

# endif
