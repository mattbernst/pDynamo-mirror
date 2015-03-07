/*------------------------------------------------------------------------------
! . File      : Timing.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _TIMER
# define _TIMER

/*define USETIMER*/

#ifdef USETIMER

# include "Real.h"

/* . A better clock often than default timer. */
# include <sys/time.h>

# define Timer_CurrentTime \
    static Real Timer_Current( void ) \
    { \
        struct timeval current ; \
        gettimeofday ( &current, NULL ) ; \
        Real time = current.tv_sec + ( current.tv_usec / 1000000.0 ) ; \
        return time ; \
    }

# endif
# endif


