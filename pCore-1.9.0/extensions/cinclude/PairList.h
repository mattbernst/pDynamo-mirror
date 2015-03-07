/*------------------------------------------------------------------------------
! . File      : PairList.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _PAIRLIST
# define _PAIRLIST

# include "Coordinates3.h"
# include "Definitions.h"
# include "IndexedSelection.h"
# include "List.h"
# include "Real1DArray.h"
# include "RegularGrid.h"
# include "RegularGridOccupancy.h"
# include "Selection.h"
# include "SelectionContainer.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The pairlist type. */
typedef struct {
    Boolean            QSELF           ; /* . A flag indicating whether the pairlist is cross or self. */
    Integer            nconnectionsi   ; /* . The upperbounds for the connection representation. */
    Integer            npairs          ; /* . The number of pairs in the list. */
    Integer            numberOfRecords ; /* . The number of records (indexed selections) in the list. */
    Integer            upperBoundi     ; /* . The sizes of the sets from which the pairlist is drawn (equal to the greatest index in the list + 1). */
    Integer            upperBoundj     ;
    Integer           *connectionsi    ; /* . The indices into connectionsj for each index. */
    Integer           *connectionsj    ; /* . The number of pairs * 2. */
    IndexedSelection **records         ; /* . An array containing the indexed selections. */
    List              *pairs           ;
} PairList ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . General procedures. */
extern PairList           *PairList_Allocate             ( const Boolean QSELF ) ;
extern void                PairList_ClearRepresentations (       PairList  *self ) ;
extern void                PairList_Deallocate           (       PairList **self ) ;
extern PairList           *PairList_FromIntegerPairArray ( const Boolean QSELF, const Integer npairs, Integer *indices ) ;
extern IndexedSelection   *PairList_Iterate              (       PairList  *self ) ;
extern Integer             PairList_Length               ( const PairList  *self ) ;
extern void                PairList_MakeRecords          (       PairList  *self ) ;
extern Integer             PairList_NumberOfRecords      (       PairList  *self ) ;
extern Integer            *PairList_ToIntegerPairArray   (       PairList  *self ) ;

/* . Cross-pairlist procedures. */
extern Status              CrossPairList_MakeConnections            ( PairList *self, const Integer upperBound ) ;

/* . Self-pairlist procedures. */
extern PairList           *SelfPairList_FromSelfPairList            ( PairList  *self, Selection *andSelection, Selection *orSelection, const Boolean QRENUMBER ) ;
extern Status              SelfPairList_MakeConnections             ( PairList  *self, const Integer upperBound ) ;
extern PairList           *SelfPairList_ToCrossPairList             ( PairList  *self, Selection *andSelection1, Selection *andSelection2, Selection *orSelection,
                                                                                                           const Boolean QEXCLUDESELF, const Boolean QINCLUDESELF,
                                                                                                              const Boolean QRENUMBER1, const Boolean QRENUMBER2 ) ;
extern SelectionContainer *SelfPairList_ToIsolateSelectionContainer ( PairList  *self, const Integer upperBound ) ;
extern Integer             SelfPairList_UpperBound                  ( PairList  *self ) ;

# endif
