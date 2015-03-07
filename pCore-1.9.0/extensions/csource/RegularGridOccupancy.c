/*------------------------------------------------------------------------------
! . File      : RegularGridOccupancy.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Regular grid occupancy procedures.
!=================================================================================================================================*/

# include <stdio.h>

# include "Memory.h"
# include "RegularGridOccupancy.h"

/* # define DEBUGPRINTING */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGridOccupancy *RegularGridOccupancy_Allocate ( const Integer numberOfCells, const Integer numberOfPoints, Status *status )
{
    RegularGridOccupancy *self = NULL ;
# ifdef DEBUGPRINTING
printf ( "\nRegular Grid Occupancy Allocate:\n" ) ;
printf ( "Number Of Cells      = %6d\n", numberOfCells  ) ;
printf ( "Number Of Points     = %6d\n", numberOfPoints ) ;
# endif
    if ( ( numberOfCells > 0 ) && ( numberOfPoints > 0 ) )
    {
        MEMORY_ALLOCATE ( self, RegularGridOccupancy ) ;
        if ( self != NULL )
        {
            /* . Initialization. */
            self->numberOfCells   = numberOfCells  ;
            self->numberOfPoints  = numberOfPoints ;
            self->cellFirstPoints = NULL ;
            self->cellPoints      = NULL ;
            self->cellTotalPoints = NULL ;
            self->pointCells      = NULL ;
            /* . Array allocation. */
            self->cellFirstPoints = Integer1DArray_Allocate ( numberOfCells , NULL ) ;
            self->cellPoints      = Integer1DArray_Allocate ( numberOfPoints, NULL ) ;
            self->cellTotalPoints = Integer1DArray_Allocate ( numberOfCells , NULL ) ;
            self->pointCells      = Integer1DArray_Allocate ( numberOfPoints, NULL ) ;
            /* . Check. */
            if ( ( self->cellFirstPoints == NULL ) || ( self->cellPoints == NULL ) ||
                 ( self->cellTotalPoints == NULL ) || ( self->pointCells == NULL ) ) RegularGridOccupancy_Deallocate ( &self ) ;
            /* . Array initialization. */
            RegularGridOccupancy_Initialize ( self ) ;
        }
        if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    else Status_Set ( status, Status_InvalidArgument ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridOccupancy_Deallocate ( RegularGridOccupancy **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Integer1DArray_Deallocate ( &((*self)->cellFirstPoints) ) ;
        Integer1DArray_Deallocate ( &((*self)->cellPoints     ) ) ;
        Integer1DArray_Deallocate ( &((*self)->cellTotalPoints) ) ;
        Integer1DArray_Deallocate ( &((*self)->pointCells     ) ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fill the structure given a grid and a set of points.
! . Note that points that are not on the grid are ignored.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGridOccupancy_Fill ( RegularGridOccupancy *self, const RegularGrid *grid, const Real2DArray *points, Status *status )
{
    Integer totalOffGrid = -1 ;
    if ( ( self != NULL ) && ( grid != NULL ) && ( points != NULL ) )
    {
# ifdef DEBUGPRINTING
printf ( "\nRegular Grid Occupancy Fill Start:\n" ) ;
printf ( "Number Of Cells      = %6d (%6d)\n", self->numberOfCells  , RegularGrid_NumberOfGridPoints ( grid ) ) ;
printf ( "Number Of Dimensions = %6d (%6d)\n", grid->ndimensions    , Real2DArray_Length ( points, 1 )        ) ;
printf ( "Number Of Points     = %6d (%6d)\n", self->numberOfPoints , Real2DArray_Length ( points, 0 )        ) ;
# endif
        if ( ( self->numberOfCells  == RegularGrid_NumberOfGridPoints ( grid ) ) &&
             ( self->numberOfPoints == Real2DArray_Length ( points, 0 )        ) &&
             ( grid->ndimensions    == Real2DArray_Length ( points, 1 )        ) )
        {
            auto Integer c, n, p ;

            /* . Determine the cell position of each point and the total number of points in each cell. */
            for ( p = totalOffGrid = 0 ; p < self->numberOfPoints ; p++ )
            {
                c = RegularGrid_FindCellIDOfPoint ( grid, Real2DArray_RowPointer ( points, p ) ) ;
                if ( c >= 0 )
                {
                    Integer1DArray_Item ( self->cellTotalPoints , c ) += 1 ;
                    Integer1DArray_Item ( self->pointCells      , p )  = c ;
                }
                else totalOffGrid++ ;
            }

            /* . Determine the index of the first point for each cell. */
            for ( c = n = 0 ; c < self->numberOfCells ; c++ )
            {
                Integer1DArray_Item ( self->cellFirstPoints, c ) = n ;
                n += Integer1DArray_Item ( self->cellTotalPoints, c ) ;
            }

            /* . Fill the cell points index array. */
            for ( p = 0 ; p < self->numberOfPoints ; p++ )
            {
                c = Integer1DArray_Item ( self->pointCells , p ) ;
                if ( c >= 0 )
                {
                    n = Integer1DArray_Item ( self->cellFirstPoints , c ) ;
                    Integer1DArray_Item ( self->cellPoints          , n )  = p ;
                    Integer1DArray_Item ( self->cellFirstPoints     , c ) += 1 ;
                }
            }

            /* . Reset cellFirstPoints. */
            Integer1DArray_AddScaledArray ( self->cellFirstPoints, -1, self->cellTotalPoints, NULL ) ;
# ifdef DEBUGPRINTING
printf ( "\nRegular Grid Occupancy Fill Stop:\n" ) ;
printf ( "Number Of Cells      = %6d\n", self->numberOfCells  ) ;
printf ( "Number Of Dimensions = %6d\n", grid->ndimensions    ) ;
printf ( "Number Of Points     = %6d\n", self->numberOfPoints ) ;
printf ( "Total Off Grid       = %6d\n", totalOffGrid         ) ;
printf ( "Total On  Grid       = %6d\n", Integer1DArray_Sum ( self->cellTotalPoints ) ) ;
printf ( "cellFirstPoints\n" ) ;
Integer1DArray_Print ( self->cellFirstPoints  ) ;
printf ( "\ncellPoints\n" ) ;
Integer1DArray_Print ( self->cellPoints       ) ;  
printf ( "\ncellTotalPoints\n" ) ;
Integer1DArray_Print ( self->cellTotalPoints  ) ;
printf ( "\npointCells\n" ) ;
Integer1DArray_Print ( self->pointCells       ) ;
# endif
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return totalOffGrid ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor given a grid and associated points.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGridOccupancy *RegularGridOccupancy_FromGridAndPoints ( const RegularGrid *grid, const Real2DArray *points, Status *status )
{
    RegularGridOccupancy *self = NULL ;
    if ( ( grid != NULL ) && ( points != NULL ) )
    {
        self = RegularGridOccupancy_Allocate ( RegularGrid_NumberOfGridPoints ( grid ), Real2DArray_Rows ( points ), status ) ;
        RegularGridOccupancy_Fill ( self, grid, points, status ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridOccupancy_Initialize ( RegularGridOccupancy *self )
{
    if ( self != NULL )
    {
        Integer1DArray_Set ( self->cellFirstPoints ,  0 ) ;
        Integer1DArray_Set ( self->cellPoints      , -1 ) ;
        Integer1DArray_Set ( self->cellTotalPoints ,  0 ) ;
        Integer1DArray_Set ( self->pointCells      , -1 ) ;
    }
}
