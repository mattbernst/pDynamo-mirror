/*------------------------------------------------------------------------------
! . File      : IndexedSelection.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _INDEXEDSELECTION
# define _INDEXEDSELECTION

# include "Integer.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The indexed selection type. */
typedef struct {
   Integer  index    ;
   Integer  nindices ;
   Integer *indices  ;
} IndexedSelection ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
/* . IndexedSelection. */
extern IndexedSelection *IndexedSelection_Allocate       ( const Integer index, const Integer nindices, const Integer *indices ) ;
extern void              IndexedSelection_ListDeallocate ( void *vself ) ;
extern void              IndexedSelection_Sort           ( IndexedSelection *self ) ;

# endif
