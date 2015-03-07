/*------------------------------------------------------------------------------
! . File      : StatusHandler.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _STATUSHANDLER
# define _STATUSHANDLER

# include "Boolean.h"
# include "Integer.h"
# include "Status.h"

/* . Need (fixed) status stack of given length - can use circular queue. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The status handler type. */
typedef struct {
    Integer events ;
    Integer flags  ;
} StatusHandler ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Status  StatusHandler_Allocate   ( StatusHandler **self ) ;
extern Boolean StatusHandler_Continue   ( StatusHandler  *self ) ;
extern void    StatusHandler_Deallocate ( StatusHandler **self ) ;
extern void    StatusHandler_Initialize ( StatusHandler  *self ) ;
extern void    StatusHandler_Record     ( StatusHandler  *self, Status status ) ;

# endif
