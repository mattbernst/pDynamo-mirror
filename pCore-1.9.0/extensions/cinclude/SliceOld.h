/*------------------------------------------------------------------------------
! . File      : SliceOld.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SLICEOLD
# define _SLICEOLD

# include "Boolean.h"
# include "Integer.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The slice type. */
typedef struct {
    Boolean isScalar ;
    Integer extent   ;
    Integer start    ;
    Integer stop     ;
    Integer stride   ;
} SliceX ;

/* . The multislice type. */
typedef struct {
    Integer  rank  ;
    Integer  size  ;
    SliceX  *items ;
} MultiSliceX ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Slice. */
extern Status  SliceX_Allocate      ( SliceX **self ) ;
extern void    SliceX_Deallocate    ( SliceX **self ) ;
extern void    SliceX_Initialize    ( SliceX  *self ) ;
extern Status  SliceX_SetFromScalar ( SliceX  *self, Integer index, const Integer rawExtent ) ;
extern Status  SliceX_SetFromSlice  ( SliceX  *self, Integer start, Integer stop, const Integer stride, const Integer rawExtent ) ;

/* . Multislice. */
extern Status  MultiSliceX_Allocate    ( MultiSliceX **self, const Integer rank ) ;
extern Status  MultiSliceX_AllocateRaw ( MultiSliceX **self ) ;
extern void    MultiSliceX_Deallocate  ( MultiSliceX **self ) ;
extern void    MultiSliceX_SetRank     ( MultiSliceX  *self ) ;

extern Integer SliceIndices_GetExtent ( const Integer start, const Integer stop, const Integer stride, const Integer extent, Status *status ) ;
extern Integer SliceIndices_GetLength ( const Integer start, const Integer stop, const Integer stride, const Integer length, Status *status ) ;

# endif
