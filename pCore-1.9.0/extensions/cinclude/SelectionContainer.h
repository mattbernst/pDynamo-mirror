/*------------------------------------------------------------------------------
! . File      : SelectionContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SELECTIONCONTAINER
# define _SELECTIONCONTAINER

# include "Definitions.h"
# include "Selection.h"
# include "Status.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The container type. */
typedef struct {
    Boolean     QOWNER     ;
    Integer     nitems     ;
    Integer     upperBound ;
    Selection **items      ;
} SelectionContainer ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern SelectionContainer *SelectionContainer_Allocate       ( const Integer nitems ) ;
extern SelectionContainer *SelectionContainer_Clone          ( const SelectionContainer  *self ) ;
extern void                SelectionContainer_Deallocate     (       SelectionContainer **self ) ;
extern Selection          *SelectionContainer_GetAllIndices  ( const SelectionContainer  *self, Selection *indices  ) ;
extern Status              SelectionContainer_MergeIsolates  (       SelectionContainer  *self, Selection *toMerge  ) ;
extern Status              SelectionContainer_RemoveIsolates (       SelectionContainer  *self, Selection *toRemove ) ;
extern Integer             SelectionContainer_UpperBound     (       SelectionContainer  *self ) ;

# endif
