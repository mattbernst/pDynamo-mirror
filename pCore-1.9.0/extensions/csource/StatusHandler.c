/*------------------------------------------------------------------------------
! . File      : StatusHandler.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Very basic status handling.
!=================================================================================================================================*/

# define STATUSHANDLERPRINTING

# include <stdio.h>

# include "Memory.h"
# include "StatusHandler.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status StatusHandler_Allocate ( StatusHandler **self )
{
    StatusHandler *new = NULL ;
    if ( ( self != NULL ) && ( (*self) == NULL ) )
    {
        MEMORY_ALLOCATE  ( new, StatusHandler ) ;
        StatusHandler_Initialize ( new ) ;
        (*self) = new ;
    }
    return ( new == NULL ) ? Status_MemoryAllocationFailure : Status_Success ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Continuation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean StatusHandler_Continue ( StatusHandler *self )
{
    return ( ( self == NULL ) || ( ( self != NULL ) && ( self->flags == 0 ) ) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void StatusHandler_Deallocate ( StatusHandler **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) ) MEMORY_DEALLOCATE ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void StatusHandler_Initialize ( StatusHandler *self )
{
    if ( self != NULL )
    {
        self->events = 0 ;
        self->flags  = 0 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Record an event.
!---------------------------------------------------------------------------------------------------------------------------------*/
void StatusHandler_Record ( StatusHandler *self, Status status )
{
    if ( self != NULL )
    {
        self->events += 1 ;
        if ( status != Status_Success )
        {
            self->flags += 1 ;
# ifdef STATUSHANDLERPRINTING
            printf ( "\nStatus Record Error> events: %d ; flags: %d ; status: %d.", self->events, self->flags, status ) ;
# endif
        }
    }
}
