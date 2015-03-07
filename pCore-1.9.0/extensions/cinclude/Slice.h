/*------------------------------------------------------------------------------
! . File      : Slice.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SLICE
# define _SLICE

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
} Slice ;

/* . The multislice type. */
typedef struct {
    Integer  rank  ;
    Integer  size  ;
    Slice   *items ;
} MultiSlice ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Slice. */
extern Slice      *Slice_Allocate            ( Status *status ) ;
extern Boolean     Slice_ApplyExtent         ( Slice  *self, const Integer extent, Status *status ) ;
extern void        Slice_Deallocate          ( Slice **self ) ;
extern void        Slice_Initialize          ( Slice  *self ) ;
extern void        Slice_SetFromScalar       ( Slice  *self, const Integer index ) ;
extern void        Slice_SetFromSliceIndices ( Slice  *self, Integer *start, Integer *stop, Integer *stride, Status *status ) ;

/* . Multislice. */
extern MultiSlice *MultiSlice_Allocate       ( const Integer size, Status *status ) ;
extern MultiSlice *MultiSlice_AllocateRaw    ( Status *status ) ;
extern void        MultiSlice_Deallocate     ( MultiSlice **self ) ;
extern Boolean     MultiSlice_IsConformable  ( MultiSlice  *self, const Integer rank, const Integer *extents, Status *status ) ;
extern void        MultiSlice_SetRank        ( MultiSlice  *self ) ;

# endif
