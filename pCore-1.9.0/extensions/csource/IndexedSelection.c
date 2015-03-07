/*------------------------------------------------------------------------------
! . File      : IndexedSelection.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
! . Indexed selections are used in pairlists.
!=============================================================================*/

# include "Definitions.h"
# include "IndexedSelection.h"
# include "Memory.h"

/*------------------------------------------------------------------------------
! . Local procedures.
!-----------------------------------------------------------------------------*/
static int Index_Compare ( const void *vterm1, const void *vterm2 ) ;

/*------------------------------------------------------------------------------
! . Allocation.
! . This is always done even if the number of indices is zero.
!-----------------------------------------------------------------------------*/
IndexedSelection *IndexedSelection_Allocate ( const int index, const int nindices, const int *indices )
{
    IndexedSelection *self = NULL ;
    self = ( IndexedSelection * ) malloc ( sizeof ( IndexedSelection ) ) ;
    if ( self != NULL )
    {
        self->index    = index                   ;
        self->indices  = NULL                    ;
        self->nindices = Maximum ( nindices, 0 ) ;
        if ( nindices > 0 )
        {
            self->indices  = Memory_Allocate_Array_Integer_Initialize ( nindices, -1 ) ;
            if ( self->indices == NULL )
            {
                Memory_Deallocate ( self ) ;
            }
            else if ( indices != NULL )
            {
                auto int i ;
                for ( i = 0 ; i < self->nindices ; i++ ) self->indices[i] = indices[i] ;
            }
        }
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Deallocation (as part of a list).
!-----------------------------------------------------------------------------*/
void IndexedSelection_ListDeallocate ( void *vself )
{
    if ( vself != NULL )
    {
        IndexedSelection *self ;
        self = ( IndexedSelection * ) vself ;
        free ( self->indices ) ;
        free ( self ) ;
    }
}

/*------------------------------------------------------------------------------
! . Sorting.
!-----------------------------------------------------------------------------*/
void IndexedSelection_Sort ( IndexedSelection *self )
{
    if ( ( self != NULL ) && ( self->nindices > 1 ) && ( self->indices != NULL ) )
    {
        auto int i, n ;
        /* . Sort. */
        qsort ( ( void * ) self->indices, ( size_t ) self->nindices, sizeof ( int ), ( void * ) Index_Compare ) ;
        /* . Remove redundant elements. */
        for ( i = 1, n = 1 ; i < self->nindices ; i++ )
        {
            if ( self->indices[i-1] != self->indices[i] )
            {
                if ( n < i ) self->indices[n] = self->indices[i] ;
                n++ ;
            }
        }
        self->nindices = n ;
    }
}

/*==============================================================================
! . Private procedures.
!============================================================================*/
static int Index_Compare ( const void *vterm1, const void *vterm2 )
{
    int  i ;
    int *term1, *term2 ;
    term1 = ( int * ) vterm1 ;
    term2 = ( int * ) vterm2 ;
         if ( (*term1) < (*term2) ) i = -1 ;
    else if ( (*term1) > (*term2) ) i =  1 ;
    else i = 0 ;
    return i ;
}
