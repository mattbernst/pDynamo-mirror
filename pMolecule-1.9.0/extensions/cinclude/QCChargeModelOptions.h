/*------------------------------------------------------------------------------
! . File      : QCChargeModelOptions.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCCHARGEMODELOPTIONS
# define _QCCHARGEMODELOPTIONS

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . QC charge model options. */
typedef enum {
    QCChargeModel_CoulombFitting = 0,
    QCChargeModel_Lowdin         = 1,
    QCChargeModel_Mulliken       = 2
} QCChargeModel ;

# endif
