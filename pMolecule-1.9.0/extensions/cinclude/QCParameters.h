/*------------------------------------------------------------------------------
! . File      : QCParameters.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCPARAMETERS
# define _QCPARAMETERS

# include "Definitions.h"
# include "GaussianBasis.h"
# include "MNDOParameters.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The qccenter type. */
typedef struct {
    Integer         atomicNumber   ;
    GaussianBasis  *densitybasis   ;
    GaussianBasis  *orbitalbasis   ;
    GaussianBasis  *poissonbasis   ;
    MNDOParameters *mndoparameters ;
} QCCenter ;

/* . The qcparameter type. */
typedef struct {
    Integer   ncenters ;
    QCCenter *centers  ;
} QCParameter ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern QCParameter *QCParameters_Allocate   ( const Integer ncenters ) ;
extern QCParameter *QCParameters_Clone      ( const QCParameter  *self ) ;
extern void         QCParameters_Deallocate (       QCParameter **self ) ;

# endif
