/*------------------------------------------------------------------------------
! . File      : Stack.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _STACK
# define _STACK

# include "Definitions.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The stack element type. */
struct StackElementStructure
{
   struct StackElementStructure *previous ;
   void                         *object   ;
} ;
typedef struct StackElementStructure StackElement ;

/* . The stack type. */
typedef struct
{
   Integer       nelements ;
   StackElement *last      ;
} Stack ;

/*------------------------------------------------------------------------------
! . Procedures.
!-----------------------------------------------------------------------------*/
extern Stack  *Stack_Allocate   ( void ) ;
extern void    Stack_Deallocate ( Stack *stack ) ;
extern void   *Stack_Pop        ( Stack *stack ) ;
extern void    Stack_Push       ( Stack *stack, void *object ) ;
extern Integer Stack_Size       ( Stack *stack ) ;

# endif
