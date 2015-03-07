/*------------------------------------------------------------------------------
! . File      : Stack.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . This module provides procedures for handling a FIFO stack of objects.
!=============================================================================*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "Definitions.h"
# include "Stack.h"

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
Stack *Stack_Allocate ( void )
{
   Stack *stack ;
   stack = ( Stack * ) malloc ( sizeof ( Stack ) ) ;
   stack->last      = NULL ;
   stack->nelements = 0    ;
   return stack ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void Stack_Deallocate ( Stack *stack )
{
   free ( stack ) ;
   stack = NULL ;
}

/*------------------------------------------------------------------------------
! . Pop an object from the stack.
!-----------------------------------------------------------------------------*/
void *Stack_Pop ( Stack *stack )
{
   StackElement *old ;
   void               *object ;
   if ( stack->last == NULL ) return NULL ;
   else
   {
      old         = stack->last   ;
      stack->last = old->previous ;
      stack->nelements --         ;
      object      = old->object   ;
      old->object = NULL          ;
      free ( old ) ;
      return object ;
   }
}

/*------------------------------------------------------------------------------
! . Push an object onto the stack.
!-----------------------------------------------------------------------------*/
void Stack_Push ( Stack *stack, void *object )
{
   StackElement *new ;
   new           = ( StackElement * ) malloc ( sizeof ( StackElement ) ) ;
   new->previous = NULL ;
   new->object   = object ;
   if ( stack->last != NULL ) new->previous = stack->last ;
   stack->last   = new ;
   stack->nelements ++ ;
}

/*------------------------------------------------------------------------------
! . Return the size of a stack.
!-----------------------------------------------------------------------------*/
int Stack_Size ( Stack *stack )
{
   return stack->nelements ;
}
