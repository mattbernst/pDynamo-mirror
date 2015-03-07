/*------------------------------------------------------------------------------
! . File      : QCMMLinkAtomCouplingOptions.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCMMLINKATOMCOUPLINGOPTIONS
# define _QCMMLINKATOMCOUPLINGOPTIONS

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . QC/MM coupling options for the electrostatic interactions involving link atoms. */
typedef enum {
    QCMMLinkAtomCoupling_MM = 0,
    QCMMLinkAtomCoupling_RC = 1,
    QCMMLinkAtomCoupling_RD = 2
} QCMMLinkAtomCoupling ;

# endif
