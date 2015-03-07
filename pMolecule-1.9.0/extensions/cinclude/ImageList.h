/*------------------------------------------------------------------------------
! . File      : ImageList.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _IMAGELIST
# define _IMAGELIST

# include "Definitions.h"
# include "List.h"
# include "PairList.h"
# include "Transformation3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The image type. */
typedef struct {
    Boolean          QOWNER          ;
    Integer          a               ;
    Integer          b               ;
    Integer          c               ;
    Real             scale           ;
    PairList        *pairlist        ;
    Transformation3 *transformation3 ;
} Image ;

/* . The image list type. */
typedef struct {
    Integer   nimages         ;
    Integer   npairs          ;
    Integer   numberOfRecords ;
    Image   **records         ;
    List     *images          ;
} ImageList ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Image. */
extern Image *Image_Allocate   ( void ) ;
extern void   Image_Deallocate ( void *vimage ) ;

/* . Image list. */
extern ImageList *ImageList_Allocate        ( void ) ;
extern Boolean    ImageList_CreateImage     (       ImageList  *self, const Integer a, const Integer b, const Integer c, const Real scale, Transformation3 *transformation3, PairList **pairlist ) ;
extern void       ImageList_Deallocate      (       ImageList **self ) ;
extern Image     *ImageList_Iterate         (       ImageList  *self ) ;
extern void       ImageList_MakeRecords     (       ImageList  *self ) ;
extern Integer    ImageList_NumberOfImages  ( const ImageList  *self ) ;
extern Integer    ImageList_NumberOfPairs   ( const ImageList  *self ) ;
extern Integer    ImageList_NumberOfRecords (       ImageList  *self ) ;

# endif
