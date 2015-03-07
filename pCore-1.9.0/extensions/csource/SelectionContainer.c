/*------------------------------------------------------------------------------
! . File      : SelectionContainer.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include <math.h>
# include <stdio.h>

# include "SelectionContainer.h"
# include "Memory.h"

/* . Upperbound must be set explicitly by whatever creates the container. */

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
SelectionContainer *SelectionContainer_Allocate ( const Integer nitems )
{
    SelectionContainer *self = NULL ;
    if ( nitems >= 0 )
    {
        auto Integer i ;
        self             = ( SelectionContainer * ) Memory_Allocate ( sizeof ( SelectionContainer ) ) ;
	self->items      = ( Selection ** ) Memory_Allocate_Array ( nitems, sizeof ( Selection * ) ) ;
	self->nitems     = nitems ;
        self->upperBound = -1     ;
        self->QOWNER     = True   ; /* . May need to make this an array? */
	for ( i = 0 ; i < nitems ; i++ ) self->items[i] = NULL ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
SelectionContainer *SelectionContainer_Clone ( const SelectionContainer  *self )
{
    SelectionContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Boolean QOK ;
        new = SelectionContainer_Allocate ( self->nitems ) ;
        QOK = ( new != NULL ) && ( new->items != NULL ) ;
        if ( QOK )
        {
            auto Integer i ;
            for ( i = 0 ; i < new->nitems ; i++ )
            {
                new->items[i] = Selection_Clone ( self->items[i] ) ;
                if ( new->items[i] == NULL )
                {
                    QOK = False ;
                    break ;
                }
            }
        }
        if ( ! QOK ) SelectionContainer_Deallocate ( &new ) ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void SelectionContainer_Deallocate ( SelectionContainer **self )
{
    if ( (*self) != NULL )
    {
        if ( (*self)->QOWNER )
        {
            auto Integer i ;
	    for ( i = 0 ; i < (*self)->nitems ; i++ ) Selection_Deallocate ( &((*self)->items[i]) ) ;
        }
        free ( (*self)->items ) ;
        free ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*------------------------------------------------------------------------------
! . Get the union of all container selections that contain any of the indices
! . in |indices|.
! . This algorithm would be easier if the container was a partition. However,
! . this is not necessarily the case (or though it could be in the future).
!-----------------------------------------------------------------------------*/
Selection *SelectionContainer_GetAllIndices ( const SelectionContainer *self, Selection *indices )
{
    Selection *total = NULL ;
    if ( ( self != NULL ) && ( indices != NULL ) )
    {
        auto Boolean *QFLAG ;

        /* . Allocate space. */
        QFLAG = Memory_Allocate_Array_Boolean_Initialize ( self->nitems, False ) ;
        if ( QFLAG != NULL )
        {
            auto Integer             i, n, s, t ;
            auto Selection *selection  ;

            /* . Flag all selections containing an index. */
            for ( i = n = 0 ; i < self->nitems ; i++ )
            {
                selection = self->items[i] ;
                for ( s = 0 ; s < selection->nindices ; s++ )
                {
                    for ( t = 0 ; t < indices->nindices ; t++ )
                    {
                        if ( selection->indices[s] == indices->indices[t] )
                        {
                            n += selection->nindices ;
                            QFLAG[i] = True ;
                            break ;
                        }
                    }
                    if ( QFLAG[i] ) break ;
                }
            }

            /* . Allocate the selection. */
            total = Selection_Allocate ( n ) ;
            if ( total != NULL )
            {
                /* . Create the selection. */
                for ( i = n = 0 ; i < self->nitems ; i++ )
                {
                    if ( QFLAG[i] )
                    {
                        selection = self->items[i] ;
                        for ( s = 0 ; s < selection->nindices ; n++, s++ ) total->indices[n] = selection->indices[s] ;
                    }
                }
                Selection_Sort ( total ) ;
            }

            /* . Finish up. */
            Memory_Deallocate_Boolean ( &QFLAG ) ;
        }
    }
    return total ;
}

/*------------------------------------------------------------------------------
! . Merge isolates.
!-----------------------------------------------------------------------------*/
Status SelectionContainer_MergeIsolates ( SelectionContainer *self, Selection *toMerge )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( toMerge != NULL ) )
    {
        auto Boolean            *QMERGE = NULL, QOK = False      ;
        auto Integer              i, nitems, nmerges, s           ;
        auto Selection **merges = NULL, *new, *selection ;
        /* . Allocate space.*/
	merges = ( Selection ** ) Memory_Allocate_Array ( self->nitems, sizeof ( Selection * ) ) ;
        QMERGE = Memory_Allocate_Array_Boolean_Initialize ( self->nitems, False ) ;
        /* . Make the flags representation of toMerge. */
        Selection_MakeFlags ( toMerge, SelectionContainer_UpperBound ( self ) ) ;
        /* . Check for memory.*/
        if ( ( merges != NULL ) && ( QMERGE != NULL ) )
        {
            /* . Loop over selections. */
            for ( i = nmerges = 0 ; i < self->nitems ; i++ )
            {
                selection = self->items[i] ;
                /* . Loop over indices. */
                for ( s = 0 ; s < selection->nindices ; s++ )
                {
                    if ( toMerge->flags[selection->indices[s]] ) { merges[nmerges] = selection ; QMERGE[i] = True ; nmerges++ ; break ; }
                }
            }
            /* . Only modify if there are more than one selections to merge. */
            if ( nmerges > 1 )
            {
                /* . Get the union of the selections. */
                new = Selection_Union ( nmerges, &(merges[0]) ) ;
                if ( new != NULL )
                {
                    /* . Compress the existing array and remove old selections. */
                    for ( i = nitems = 0 ; i < self->nitems ; i++ )
                    {
                        if ( QMERGE[i] )
                        {
                            if ( self->QOWNER ) Selection_Deallocate ( &(self->items[i]) ) ;
                        }
                        else
                        {
                            self->items[nitems] = self->items[i] ;
                            nitems++ ;
                        }
                    }
                    /* . Append the merged selection. */
                    self->items[nitems] = new ;
                    self->nitems = nitems + 1 ;
                    /* . Finish up. */
                    QOK = True ;
                }
            }
            else QOK = True ;
        }
        /* . Finish up.*/
        free ( merges ) ;
        free ( QMERGE ) ;
        if ( QOK ) status = Status_Success     ;
        else       status = Status_OutOfMemory ;
    }
    return status ;
}

/*------------------------------------------------------------------------------
! . Remove isolates.
!-----------------------------------------------------------------------------*/
Status SelectionContainer_RemoveIsolates ( SelectionContainer *self, Selection *toRemove )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( toRemove != NULL ) )
    {
        auto Boolean    isFlagged ;
        auto Integer    i, n = 0, s ;
        auto Selection *selection ;
        /* . Make the flags representation of toRemove. */
        Selection_MakeFlags ( toRemove, SelectionContainer_UpperBound ( self ) ) ;
        /* . Loop over the indices of each selection. */
        for ( i = 0 ; i < self->nitems ; i++ )
        {
            isFlagged = False ;
            selection = self->items[i] ;
            for ( s = 0 ; s < selection->nindices ; s++ )
            {
                if ( toRemove->flags[selection->indices[s]] ) { isFlagged = True ; break ; }
            }
            /* . Remove the selection. */
            if ( isFlagged )
            {
                if ( self->QOWNER ) Selection_Deallocate ( &(self->items[i]) ) ;
            }
            /* . Save the selection. */
            else
            {
                self->items[n] = self->items[i] ;
                n++ ;
            }
        }
        self->nitems     =  n ;
        self->upperBound = -1 ;
        SelectionContainer_UpperBound ( self ) ;
    }
    return status ;
}

/*------------------------------------------------------------------------------
! . Determine the upper bound of the container.
!-----------------------------------------------------------------------------*/
Integer SelectionContainer_UpperBound ( SelectionContainer *self )
{
    Integer upperBound = -1 ;
    if ( self != NULL )
    {
        if ( self->upperBound >= 0 ) { upperBound = self->upperBound ; }
        else
        {
            auto Integer i ;
            for ( i = 0 ; i < self->nitems ; i++ ) { upperBound = Maximum ( Selection_UpperBound ( self->items[i] ), upperBound ) ; }
            self->upperBound = upperBound ;
        }
    }
    return upperBound ;
}
