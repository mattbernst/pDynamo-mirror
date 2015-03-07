/*------------------------------------------------------------------------------
! . File      : Selection.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include <stdarg.h>
# include <stdio.h>
# include <stdlib.h>

# include "Memory.h"
# include "Selection.h"

/* . This needs to be rationalized. */

/*------------------------------------------------------------------------------
! . Local procedures.
!-----------------------------------------------------------------------------*/
static Integer SelectionIndex_Compare ( const void *vterm1, const void *vterm2 ) ;

/*==============================================================================
! . Standard procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
Selection *Selection_Allocate ( const Integer nindices )
{
    Selection *self = NULL ;
    self             = ( Selection * ) Memory_Allocate ( sizeof ( Selection ) ) ;
    self->QSORTED    = False ;
    self->nflags     = -1   ;
    self->nindices   =  0   ;
    self->npositions = -1   ;
    self->flags      = NULL ;
    self->indices    = NULL ;
    self->positions  = NULL ;
    if ( nindices > 0 )
    {
        self->nindices = nindices ;
        self->indices  = Memory_Allocate_Array_Integer_Initialize ( nindices, -1 ) ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Clear the flags representation.
!-----------------------------------------------------------------------------*/
void Selection_ClearFlags ( Selection *self )
{
    if ( self != NULL )
    {
        free ( self->flags ) ;
        self->nflags = -1    ;
        self->flags  = NULL  ;
    }
}

/*------------------------------------------------------------------------------
! . Clear the position representation.
!-----------------------------------------------------------------------------*/
void Selection_ClearPositions ( Selection *self )
{
    if ( self != NULL )
    {
        free ( self->positions ) ;
        self->npositions = -1    ;
        self->positions  = NULL  ;
    }
}

/*------------------------------------------------------------------------------
! . Clear representations.
!-----------------------------------------------------------------------------*/
void Selection_ClearRepresentations ( Selection *self )
{
    Selection_ClearFlags     ( self ) ;
    Selection_ClearPositions ( self ) ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
Selection *Selection_Clone ( const Selection *self )
{
    Selection *new = NULL ;
    if ( self != NULL )
    {
        new = Selection_Allocate ( self->nindices ) ;
        if ( new->indices != NULL )
        {
            auto Integer i ;
            for ( i = 0 ; i < new->nindices ; i++ ) new->indices[i] = self->indices[i] ;
        }
        new->QSORTED = self->QSORTED ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Comparison.
!-----------------------------------------------------------------------------*/
int Selection_Compare ( Selection *self, Selection *other )
{
    Integer comparison = 0 ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i, nindices ;
        /* . Both selections should be sorted. */
        Selection_Sort ( self  ) ;
        Selection_Sort ( other ) ;
        /* . Loop over the indices in both selections. */
        nindices = self->nindices ;
        if ( other->nindices < nindices ) nindices = other->nindices ;
        for ( i = 0 ; i < nindices ; i++ )
        {
            if ( self->indices[i] == other->indices[i] ) continue ;
            else
            {
                     if ( self->indices[i] > other->indices[i] ) comparison =  1 ;
                else if ( self->indices[i] < other->indices[i] ) comparison = -1 ;
                break ;
            }
        }
        /* . If the indices checked are the same, the longer selection wins.*/
        if ( comparison == 0 )
        {
                 if ( self->nindices > other->nindices ) comparison =  1 ;
            else if ( self->nindices < other->nindices ) comparison = -1 ;
        }
    }
    return comparison ;
}

/*------------------------------------------------------------------------------
! . Return the complement of the selection.
!-----------------------------------------------------------------------------*/
Selection *Selection_Complement ( Selection *self, const Integer upperBound )
{
    Selection *new = NULL ;
    if ( self != NULL )
    {
        auto Integer n, range ;
        range = Selection_UpperBound ( self ) ;
        if ( range < upperBound ) range = upperBound ;
        n   = range - self->nindices ;
        new = Selection_Allocate ( n ) ;
        if ( ( new->nindices > 0 ) && ( new->indices != NULL ) )
        {
            auto Integer i ;
            Selection_MakeFlags ( self, range ) ;
            for ( i = 0, n = 0 ; i < range ; i++ )
            {
                if ( ! self->flags[i] ) { new->indices[n] = i ; n++ ; }
            }
        }
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void Selection_Deallocate ( Selection **self )
{
    if ( (*self) != NULL )
    {
        Selection_ClearRepresentations ( (*self) ) ;
        Memory_Deallocate ( (*self)->indices   ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Constructor from boolean array.
! . The selection is automatically sorted.
!-----------------------------------------------------------------------------*/
Status Selection_FromBoolean1DArray ( Selection **self, const Boolean1DArray *flags )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( flags != NULL ) )
    {
        auto Integer i, n ;
        /* . Find the number of True flags. */
        for ( i = n = 0 ; i < flags->length ; i++ ) { if ( Boolean1DArray_Item ( flags, i ) ) n++ ; }
        /* . Allocate space. */
        (*self) = Selection_Allocate ( n ) ;
        /* . Fill the selection. */
        if ( (*self) == NULL ) status = Status_OutOfMemory ;
        else
        {
            for ( i = n = 0 ; i < flags->length ; i++ ) { if ( Boolean1DArray_Item ( flags, i ) ) { (*self)->indices[n] = i ; n++ ; } }
            status = Status_Success ;
        }
        (*self)->QSORTED = True ;
    }
    return status ;
}

/*------------------------------------------------------------------------------
! . Constructor from flags.
! . The selection is automatically sorted.
!-----------------------------------------------------------------------------*/
Selection *Selection_FromFlags ( const Integer nflags, const Boolean *flags )
{
    Selection *self = NULL ;
    if ( ( nflags > 0 ) && ( flags != NULL ) )
    {
        auto Integer i, n ;
        for ( i = 0, n = 0 ; i < nflags ; i++ )
        {
            if ( flags[i] ) n++ ;
        }
        self = Selection_Allocate ( n ) ;
        if ( self->indices != NULL )
        {
            for ( i = 0, n = 0 ; i < nflags ; i++ )
            {
                if ( flags[i] ) { self->indices[n] = i ; n++ ; }
            }
        }
        self->QSORTED = True ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Constructor from an integer array.
! . The selection is sorted.
!-----------------------------------------------------------------------------*/
Selection *Selection_FromIntegerArray ( const Integer nindices, const Integer *indices )
{
    Selection *self = NULL ;
    if ( ( nindices == 0 ) || ( ( nindices > 0 ) && ( indices != NULL ) ) )
    {
        self = Selection_Allocate ( nindices ) ;
        if ( self->indices != NULL )
        {
            auto Integer i ;
            for ( i = 0 ; i < nindices ; i++ ) self->indices[i] = indices[i] ;
        }
        Selection_Sort ( self ) ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Return the intersection of a set of selections.
!-----------------------------------------------------------------------------*/
Selection *Selection_Intersection ( const Integer nselections, ... )
{
    Selection *self = NULL ;
    if ( nselections > 0 )
    {
        Integer             i, n, nflags = 0, upperBound ;
        Selection *other = NULL ;
        va_list         argp ;
        /* . Determine the value of the smallest upperbound. */
        va_start ( argp, nselections ) ;
        for ( i = 0 ; i < nselections ; i++ )
        {
            other      = va_arg ( argp, Selection * ) ;
            upperBound = Selection_UpperBound ( other ) ;
            if ( i == 0 ) nflags = upperBound ;
            else if ( upperBound < nflags ) nflags = upperBound ;
        }
        va_end ( argp ) ;
        /* . Non-empty selection. */
        if ( nflags > 0 )
        {
            auto Boolean *flags = NULL ;
            flags = Memory_Allocate_Array_Boolean_Initialize ( nflags, True ) ;
            va_start ( argp, nselections ) ;
            for ( i = 0 ; i < nselections ; i++ )
            {
                other = va_arg ( argp, Selection * ) ;
                Selection_MakeFlags ( other, nflags ) ;
                for ( n = 0 ; n < nflags ; n++ )
                {
                    if ( flags[n] && other->flags[n] ) flags[n] = True ;
                    else flags[n] = False ;
                }
            }
            va_end ( argp ) ;
            /* . Create the selection. */
            self = Selection_FromFlags ( nflags, flags ) ;
            Memory_Deallocate ( flags ) ;
        }
        /* . Empty selection. */
        else self = Selection_Allocate ( 0 ) ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Make the flags representation of the selection.
! . The representation will have a minimum size of self's upperbound but it can
! . be made bigger by using the argument upperbound.
!-----------------------------------------------------------------------------*/
void Selection_MakeFlags ( Selection *self, const Integer upperBound )
{
    if ( self != NULL )
    {
        auto Integer range ;
        /* . Get the size of the representation. */
        range = Selection_UpperBound ( self ) ;
        if ( ( upperBound > 0 ) && ( range < upperBound ) ) range = upperBound ;
        /* . Remove representations that are too small. */
        if ( range > self->nflags ) Selection_ClearFlags ( self ) ;
        /* . Create the representation. */
        if ( ( range > 0 ) && ( ( self->nflags < 0 ) || ( self->flags == NULL ) ) )
        {
            auto Integer i ;
            self->nflags = range ;
            self->flags  = Memory_Allocate_Array_Boolean_Initialize ( range, False ) ;
            for ( i = 0 ; i < self->nindices ; i++ ) self->flags[self->indices[i]] = True ;
        }
    }
}

/*------------------------------------------------------------------------------
! . Make the positions representation of the selection.
!-----------------------------------------------------------------------------*/
void Selection_MakePositions ( Selection *self, const Integer upperBound )
{
    if ( self != NULL )
    {
        auto Integer range ;
        /* . Get the size of the representation. */
        range = Selection_UpperBound ( self ) ;
        if ( ( upperBound > 0 ) && ( range < upperBound ) ) range = upperBound ;
        /* . Remove representations that are too small. */
        if ( range > self->npositions ) Selection_ClearPositions ( self ) ;
        /* . Create the representation - the selection must be sorted. */
        if ( ( range > 0 ) && ( ( self->npositions < 0 ) || ( self->positions == NULL ) ) )
        {
            auto Integer i, n ;
            self->npositions = range ;
            self->positions  = Memory_Allocate_Array_Integer_Initialize ( range, -1 ) ;
            for ( i = 0, n = 0 ; i < self->nindices ; i++ ) { self->positions[self->indices[i]] = n ; n++ ; }
        }
    }
}

/*------------------------------------------------------------------------------
! . Merging.
!-----------------------------------------------------------------------------*/
Selection *Selection_Merge ( const Selection *self, const Selection *other, const Integer indexincrement )
{
    Selection *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer   i ;
        new = Selection_Allocate ( self->nindices + other->nindices ) ;
        for ( i = 0 ; i < self->nindices  ; i++ ) new->indices[i]                = self->indices[i] ;
        for ( i = 0 ; i < other->nindices ; i++ ) new->indices[i+self->nindices] = other->indices[i] + indexincrement ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Pruning.
!-----------------------------------------------------------------------------*/
Selection *Selection_Prune ( Selection *self, Selection *selection )
{
    Selection *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Integer i, n  ;
        n = Selection_UpperBound ( self ) ;
        Selection_MakeFlags     ( selection, n ) ;
        Selection_MakePositions ( selection, n ) ;
        /* . Find the number of indices. */
        for ( i = 0, n = 0 ; i < self->nindices ; i++ )
        {
            if ( selection->flags[self->indices[i]] ) n++ ;
        }
        /* . Create the new selection. */
        if ( n >= 0 )
        {
            new = Selection_Allocate ( n ) ;
            if ( n > 0 )
            {
                for ( i = 0, n = 0 ; i < self->nindices ; i++ )
                {
                    if ( selection->flags[self->indices[i]] )
                    {
                        new->indices[n] = selection->positions[self->indices[i]] ;
                        n++ ;
                    }
                }
                new->QSORTED = self->QSORTED ;
            }
        }
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Size.
!-----------------------------------------------------------------------------*/
int Selection_Size ( const Selection *self )
{
    Integer size = 0 ;
    if ( self != NULL ) size = self->nindices ;
    return size ;
}

/*------------------------------------------------------------------------------
! . Sorting.
! . Indices are sorted in ascending order with no duplicates.
!-----------------------------------------------------------------------------*/
void Selection_Sort ( Selection *self )
{
    if ( ( self != NULL ) && ( ! self->QSORTED ) )
    {
        if ( self->nindices > 1 )
        {
            auto Integer i, n ;
            /* . Order the terms within the container. */
            qsort ( ( void * ) self->indices, ( size_t ) self->nindices, sizeof ( Integer ), ( void * ) SelectionIndex_Compare ) ;
            /* . Remove redundant terms. */
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
        self->QSORTED = True ;
    }
}

/*------------------------------------------------------------------------------
! . Return the union of a set of selections.
!-----------------------------------------------------------------------------*/
Selection *Selection_Union ( const Integer nselections, Selection **selections )
{
    Selection *self = NULL ;
    if ( nselections > 0 )
    {
        Integer             i, n, nflags = 0, upperBound ;
        Selection *other = NULL ;
        /* . Determine the value of the largest upperbound. */
        for ( i = 0 ; i < nselections ; i++ )
        {
            other      = selections[i] ;
            upperBound = Selection_UpperBound ( other ) ;
            if ( i == 0 ) nflags = upperBound ;
            else if ( upperBound > nflags ) nflags = upperBound ;
        }
        /* . Non-empty selection. */
        if ( nflags > 0 )
        {
            auto Boolean *flags = NULL ;
            flags = Memory_Allocate_Array_Boolean_Initialize ( nflags, False ) ;
            for ( i = 0 ; i < nselections ; i++ )
            {
                other = selections[i] ;
                Selection_MakeFlags ( other, nflags ) ;
                for ( n = 0 ; n < nflags ; n++ )
                {
                    if ( other->flags[n] ) flags[n] = True ;
                }
            }
            /* . Create the selection. */
            self = Selection_FromFlags ( nflags, flags ) ;
            Memory_Deallocate ( flags ) ;
        }
        /* . Empty selection. */
        else self = Selection_Allocate ( 0 ) ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one and will be the minimal
! . size of flag and position representations.
!-----------------------------------------------------------------------------*/
int Selection_UpperBound ( Selection *self )
{
    Integer upperBound = 0 ;
    if ( ( self != NULL ) && ( self->nindices > 0 ) )
    {
        Selection_Sort ( self ) ;
        upperBound = self->indices[self->nindices-1] + 1 ;
    }
    return upperBound ;
}

/*==============================================================================
! . Private procedures.
!============================================================================*/
static Integer SelectionIndex_Compare ( const void *vterm1, const void *vterm2 )
{
    Integer  i ;
    Integer *term1, *term2 ;
    term1 = ( Integer * ) vterm1 ;
    term2 = ( Integer * ) vterm2 ;
         if ( (*term1) < (*term2) ) i = -1 ;
    else if ( (*term1) > (*term2) ) i =  1 ;
    else i = 0 ;
    return i ;
}
