/*------------------------------------------------------------------------------
! . File      : Selection.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SELECTION
# define _SELECTION

# include "Boolean1DArray.h"
# include "Definitions.h"
# include "Integer.h"
# include "Status.h"

/* . To be redone.

Is ImmutableOrderedIntegerSet. Leave as ordered.

Have I1DA for indices or have plain array or pointer for easy access?

positions is better implemented as a dictionary with default -1.

flags is just look up in position keys.

*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The atom selection type. */
typedef struct {
    Boolean  QSORTED    ;
    Integer  nflags     ;
    Integer  nindices   ;
    Integer  npositions ;
    Boolean *flags      ;
    Integer *indices    ;
    Integer *positions  ;
} Selection ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Selection *Selection_Allocate             ( const Integer nindices ) ;
extern void       Selection_ClearFlags           (       Selection  *self ) ;
extern void       Selection_ClearPositions       (       Selection  *self ) ;
extern void       Selection_ClearRepresentations (       Selection  *self ) ;
extern Selection *Selection_Clone                ( const Selection  *self ) ;
extern Integer    Selection_Compare              (       Selection  *self, Selection *other ) ;
extern Selection *Selection_Complement           (       Selection  *self, const Integer upperBound ) ;
extern void       Selection_Deallocate           (       Selection **self ) ;
extern Status     Selection_FromBoolean1DArray   (       Selection **self, const Boolean1DArray *flags ) ;
Selection        *Selection_FromFlags            ( const Integer nflags, const Boolean *flags ) ;
extern Selection *Selection_FromIntegerArray     ( const Integer nindices, const Integer *indices ) ;
extern Selection *Selection_Intersection         ( const Integer nselections, ... ) ;
extern void       Selection_MakeFlags            (       Selection  *self, const Integer upperBound ) ;
extern void       Selection_MakePositions        (       Selection  *self, const Integer upperBound ) ;
extern Selection *Selection_Merge                ( const Selection  *self, const Selection *other, const Integer indexincrement ) ;
extern Selection *Selection_Prune                (       Selection  *self, Selection *selection ) ;
extern Integer    Selection_Size                 ( const Selection  *self ) ;
extern void       Selection_Sort                 (       Selection  *self ) ;
extern Selection *Selection_Union                ( const Integer nselections, Selection **selections ) ;
extern Integer    Selection_UpperBound           (       Selection  *self ) ;

# endif
