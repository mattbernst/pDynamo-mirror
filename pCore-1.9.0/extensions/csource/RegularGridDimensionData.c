/*------------------------------------------------------------------------------
! . File      : RegularGridDimensionData.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . 1-D indexed value arrays.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "Macros.h"
# include "Memory.h"
# include "RegularGridDimensionData.h"
# include "SliceOld.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer Value_CompareAscending  ( const void *vterm1, const void *vterm2 ) ;
static Integer Value_CompareDescending ( const void *vterm1, const void *vterm2 ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGridDimensionData *RegularGridDimensionData_Allocate ( const Integer length, Status *status )
{
    RegularGridDimensionData *self = NULL ;
    MEMORY_ALLOCATE ( self, RegularGridDimensionData ) ;
    if ( self != NULL )
    {
        self->items   = NULL   ;
        self->isOwner = True   ;
        self->isSlice = False  ;
        self->length  = Maximum ( length, 0 ) ;
        self->offset  = 0      ;
        self->size    = Maximum ( length, 0 ) ;
        self->stride  = 1      ;
        if ( length > 0 )
        {
            MEMORY_ALLOCATEARRAY ( self->items, length, RegularGridDimensionDatum ) ;
            if ( self->items == NULL ) RegularGridDimensionData_Deallocate ( &self ) ;
        }
        else if ( length < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the data indices to ensure that they fall within a specific range.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGridDimensionData_CheckIndices ( RegularGridDimensionData *self, const Integer range, const Boolean isPeriodic )
{
    Integer numberOutsideRange = 0 ;
    if ( self != NULL )
    {
        auto Integer i, index ;
        if ( isPeriodic )
        {
            for ( i = 0 ; i < self->length ; i++ )
            {
                index = self->items[i].index ;
                if ( index < 0 )
                {
                    while ( index <  0     ) index += range ;
                    self->items[i].index = index ;
                }
                else if ( index >= range )
                {
                    while ( index >= range ) index -= range ;
                    self->items[i].index = index ;
                }
            }
        }
        else
        {
            for ( i = 0 ; i < self->length ; i++ )
            {
                index = self->items[i].index ;
                if ( index < 0 )
                {
                    numberOutsideRange += 1 ;
                    self->items[i].index = 0 ;
                }
                else if ( index >= range )
                {
                    numberOutsideRange += 1 ;
                    self->items[i].index = range - 1 ;
                }
            }
        }
    }
    return numberOutsideRange ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridDimensionData_Deallocate ( RegularGridDimensionData **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        if ( (*self)->isOwner ) MEMORY_DEALLOCATE ( (*self)->items ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGridDimensionData_Length ( const RegularGridDimensionData *self )
{
    if ( self == NULL ) return 0 ;
    else                return self->length ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting by value - ascending in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridDimensionData_SortByValue ( RegularGridDimensionData *self, const Boolean doAscending )
{
    if ( ( self != NULL ) && ( self->length > 1 ) )
    {
        if ( doAscending ) qsort ( ( void * ) RegularGridDimensionData_Items ( self ), ( size_t ) self->length, self->stride * sizeof ( RegularGridDimensionDatum ), ( void * ) Value_CompareAscending  ) ;
        else               qsort ( ( void * ) RegularGridDimensionData_Items ( self ), ( size_t ) self->length, self->stride * sizeof ( RegularGridDimensionDatum ), ( void * ) Value_CompareDescending ) ;
    }
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
static Integer Value_CompareAscending ( const void *vterm1, const void *vterm2 )
{
    Integer i ;
    RegularGridDimensionDatum *term1, *term2 ;
    term1 = ( RegularGridDimensionDatum * ) vterm1 ;
    term2 = ( RegularGridDimensionDatum * ) vterm2 ;
         if ( term1->value < term2->value ) i = -1 ;
    else if ( term1->value > term2->value ) i =  1 ;
    else i = 0 ;
    return i ;
}

static Integer Value_CompareDescending ( const void *vterm1, const void *vterm2 )
{
    Integer i ;
    RegularGridDimensionDatum *term1, *term2 ;
    term1 = ( RegularGridDimensionDatum * ) vterm1 ;
    term2 = ( RegularGridDimensionDatum * ) vterm2 ;
         if ( term1->value < term2->value ) i =  1 ;
    else if ( term1->value > term2->value ) i = -1 ;
    else i = 0 ;
    return i ;
}
