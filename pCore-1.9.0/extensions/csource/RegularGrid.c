/*------------------------------------------------------------------------------
! . File      : RegularGrid.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Defines a regular n-dimensional grid.
!===================================================================================================================================
!
! . Notes:
!
! . Each dimension, i, is independent and is divided into ni cells numbered from 0 to ni-1.
!
! . The origin, oi, is the extreme left-hand side of the box and the box width is ni * hi.
!
! . The position of the center of the first cell is oi + hi/2 and of the last oi + ( ni - 1/2 ) * hi. The general position is
!   oi + ( i + 1/2 ) * hi.
!
! . Cells are stored sequentially with the rightmost dimension changing quickest (i.e. the same as C's row-major ordering).
!
! . Periodic dimensions are allowed.
!
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>
# include <stdio.h>

# include "Macros.h"
# include "Memory.h"
# include "RegularGrid.h"

/*==================================================================================================================================
! . Regular Grid Functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGrid *RegularGrid_Allocate ( const Integer ndimensions, Status *status )
{
    RegularGrid *self = NULL ;
    MEMORY_ALLOCATE ( self, RegularGrid ) ;
    if ( self != NULL )
    {
        self->dimensions  = NULL ;
        self->ndimensions = Maximum ( ndimensions, 0 ) ;
        if ( ndimensions > 0 )
        {
            MEMORY_ALLOCATEARRAY ( self->dimensions, ndimensions, RegularGridDimension ) ;
            self->workI = Memory_Allocate_Array_Integer ( 2 * ndimensions ) ; /* . Temporary? */
            self->workR = Memory_Allocate_Array_Real    (     ndimensions ) ;
            if ( ( self->dimensions == NULL ) || ( self->workI == NULL ) || ( self->workR == NULL ) ) RegularGrid_Deallocate ( &self ) ;
            else
            {
                auto Integer i ;
                for ( i = 0 ; i < self->ndimensions ; i++ ) RegularGridDimension_Initialize ( &(self->dimensions[i]) ) ;
            }
        }
        else if ( ndimensions < 0 ) Status_Set ( status, Status_NegativeArrayLength ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert a cell ID to cell indices (non-compact).
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGrid_CellIDToIndices ( const RegularGrid *self, const Integer cellID, Integer *indices, Status *status )
{
    auto Boolean isOK = False ;
    if ( ( self != NULL ) && ( indices != NULL ) && ( cellID >= 0 ) ) /* . Check also maximum cellID? */
    {
        auto Integer d, i, n = cellID ;
        for ( d = 0, isOK = True ; d < self->ndimensions ; d++ )
        {
            i = n / self->dimensions[d].stride ;
            n = n % self->dimensions[d].stride ;
            RegularGridDimension_RegularizeIndexToMinusOne ( self->dimensions[d], i ) ;
            if ( i < 0 ) isOK = False ;
            indices[d] = i ;
        }
        if ( n != 0 ) isOK = False ;
    }
    if ( ! isOK ) Status_Set ( status, Status_InvalidArgument ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert cell indices to cell ID (non-compact).
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGrid_CellIndicesToID ( const RegularGrid *self, const Integer *indices )
{
    Integer cellID = -1 ;
    if ( ( self != NULL ) && ( indices != NULL ) )
    {
        auto Integer outsideGrid ;
        RegularGrid_CopyDimensionData ( self, indices, self->workI ) ;
        RegularGrid_RegularizeIndices ( self, outsideGrid ) ;
        if ( outsideGrid == 0 ) RegularGrid_MakeCellID ( self, cellID ) ;
    }
    return cellID ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGrid *RegularGrid_Clone ( const RegularGrid *self, Status *status )
{
    RegularGrid *clone = NULL ;
    if ( self != NULL )
    {
        clone = RegularGrid_Allocate ( self->ndimensions, status ) ;
        if ( clone != NULL )
        {
            auto Integer i ;
            for ( i = 0 ; i < self->ndimensions ; i++ ) RegularGridDimension_CopyTo ( &(self->dimensions[i]), &(clone->dimensions[i]) ) ;
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGrid_Deallocate ( RegularGrid **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Memory_Deallocate_Integer ( &((*self)->workI) ) ;
        Memory_Deallocate_Real    ( &((*self)->workR) ) ;
        MEMORY_DEALLOCATE ( (*self)->dimensions ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the maximum and minimum distances squared between a cell and a point.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGrid_DistancesSquaredBetweenCellAndPoint ( const RegularGrid *self, const Integer *indices, const Real *coordinates, 
                                                                     Real *maximumDistanceSquared, Real *minimumDistanceSquared,
                                                                                                                Status *status )
{
    if ( maximumDistanceSquared != NULL ) (*maximumDistanceSquared) = 0.0e+00 ;
    if ( minimumDistanceSquared != NULL ) (*minimumDistanceSquared) = 0.0e+00 ;
    if ( ( self != NULL ) && ( indices != NULL ) && ( coordinates != NULL ) && ( ( maximumDistanceSquared != NULL ) || ( minimumDistanceSquared != NULL ) ) )
    {
        auto Integer d, n ;
        auto Real    dL, dU, h, lowerN, maxR2, maxX, minR2, minX, upperN, x ;
        maxR2 = 0.0e+00 ;
        minR2 = 0.0e+00 ;
        for ( d = 0 ; d < self->ndimensions ; d++ )
        {
            h = self->dimensions[d].binSize ;
            n = indices[d] ;
            x = coordinates[d] ;
            if ( self->dimensions[d].isPeriodic )
            {
                dL = x - ( ( Real ) n * h + self->dimensions[d].lower ) ;
                dU = dL - h ;
                if ( ( dL * dU ) < 0.0e+00 ) { maxX = Maximum ( fabs ( dL ), fabs ( dU ) ) ; minX = 0.0e+00 ; }
                else
                {
                    RegularGridDimension_AdjustPeriodicDistance ( self->dimensions[d], dL ) ;
                    RegularGridDimension_AdjustPeriodicDistance ( self->dimensions[d], dU ) ;
                    minX = Minimum ( fabs ( dL ), fabs ( dU ) ) ;
                    maxX = minX + h ;
                }
            }
            else
            {
                lowerN = ( Real ) n * h + self->dimensions[d].lower ;
                upperN = lowerN + h ;
                     if ( x > upperN ) { maxX = x - lowerN ; minX = x - upperN ; }
                else if ( x < lowerN ) { maxX = upperN - x ; minX = lowerN - x ; }
                else { maxX = Maximum ( ( x - lowerN ), ( upperN - x ) ) ; minX = 0.0e+00 ; }
            }
            maxR2 += maxX * maxX ;
            minR2 += minX * minX ;
        }
        if ( maximumDistanceSquared != NULL ) (*maximumDistanceSquared) = maxR2 ;
        if ( minimumDistanceSquared != NULL ) (*minimumDistanceSquared) = minR2 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the maximum and minimum distances squared between cells.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGrid_DistancesSquaredBetweenCells ( const RegularGrid *self, const Integer *indices1, const Integer *indices2,
                                                               Real *maximumDistanceSquared, Real *minimumDistanceSquared,
                                                                                                          Status *status )
{
    if ( maximumDistanceSquared != NULL ) (*maximumDistanceSquared) = 0.0e+00 ;
    if ( minimumDistanceSquared != NULL ) (*minimumDistanceSquared) = 0.0e+00 ;
    if ( ( self != NULL ) && ( indices1 != NULL ) && ( indices2 != NULL ) && ( ( maximumDistanceSquared != NULL ) || ( minimumDistanceSquared != NULL ) ) )
    {
        auto Integer d, d12, n ;
        auto Real    delta, dL1U2, dU1L2, h, maxR2, maxX, minR2, minX ;
        auto RegularGridDimension axis ;
        maxR2 = 0.0e+00 ;
        minR2 = 0.0e+00 ;
        for ( d = 0 ; d < self->ndimensions ; d++ )
        {
            axis = self->dimensions[d] ;
            h    = axis.binSize ;
            if ( axis.isPeriodic )
            {
                n   = axis.bins ;
                d12 = abs ( Modulo ( indices1[d] - indices2[d], n ) ) ;
                     if   ( d12 == 0 )                       { maxX =           h ; minX = 0.0e+00 ; }
                else if ( ( d12 == 1 )  || ( d12 = n - 1 ) ) { maxX = 2.0e+00 * h ; minX = 0.0e+00 ; }
                else
                {
                    delta = ( Real ) d12 ;
                    dL1U2 = ( delta + 1.0e+00 ) * h ;
                    dU1L2 =   dL1U2 - 2.0e+00   * h ;
                    RegularGridDimension_AdjustPeriodicDistance ( axis, dL1U2 ) ;
                    RegularGridDimension_AdjustPeriodicDistance ( axis, dU1L2 ) ;
                    maxX = Minimum ( Maximum ( dL1U2, dU1L2 ) + 2.0e+00 * h, axis.period ) ;
                    minX =           Minimum ( dL1U2, dU1L2 ) ;
                }
            }
            else
            {
                delta = fabs ( ( Real ) ( indices1[d] - indices2[d] ) ) ;
                maxX  = ( delta + 1.0e+00 ) * h ;
                minX  = Maximum ( ( delta - 1.0e+00 ) * h, 0.0e+00 ) ;
            }
            maxR2 += maxX * maxX ;
            minR2 += minX * minX ;
        }
        if ( maximumDistanceSquared != NULL ) (*maximumDistanceSquared) = maxR2 ;
        if ( minimumDistanceSquared != NULL ) (*minimumDistanceSquared) = minR2 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the ID of the cell on the grid within which a point lies.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . No checking of dimensions is done. */
Integer RegularGrid_FindCellIDOfPoint ( const RegularGrid *self, const Real *coordinates )
{
    Integer cellID = -1 ;
    if ( ( self != NULL ) && ( coordinates != NULL ) )
    {
        auto Integer outsideGrid ;
        RegularGrid_DecomposePoint    ( self, coordinates ) ;
        RegularGrid_RegularizeIndices ( self, outsideGrid ) ;
        if ( outsideGrid == 0 ) RegularGrid_MakeCellID ( self, cellID ) ;
    }
    return cellID ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the indices of the cell on the grid within which a point lies.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGrid_FindCellIndicesOfPoint ( const RegularGrid *self, const Real *coordinates, const Boolean regularizeToBoundary, Integer *indices, Real *fractional )
{
    Integer outsideGrid = 0 ;
    if ( ( self != NULL ) && ( coordinates != NULL ) && ( indices != NULL ) )
    {
        auto Integer outsideGrid ;
        RegularGrid_DecomposePoint ( self, coordinates ) ;
        if ( regularizeToBoundary ) { RegularGrid_RegularizeIndicesToBoundary ( self, outsideGrid ) ; }
        else                        { RegularGrid_RegularizeIndices           ( self, outsideGrid ) ; }
        if ( fractional != NULL ) RegularGrid_CopyDimensionData ( self, self->workR, fractional ) ;
        if ( indices    != NULL ) RegularGrid_CopyDimensionData ( self, self->workI, indices    ) ;
    }
    return outsideGrid ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find cells within range of a given cell.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGrid_FindCellsWithinRangeOfCell ( const RegularGrid *self, const Integer cellID, RegularGridSearchRange *range, Status *status )
{
    Integer n = 0 ;
    if ( ( self != NULL ) && ( range != NULL ) )
    {
        auto Integer c, d, *indices, *newIndices, outsideGrid ;
        /* . Get cell indices. */
        indices = range->workI ;
        RegularGrid_CellIDToIndices ( self, cellID, indices, status ) ;
        /* . Loop over possible cells within range. */
        newIndices = self->workI ;
        for ( c = 0 ; c < range->numberOfCells0 ; c++ )
        {
            for ( d = 0 ; d < self->ndimensions ; d++ ) newIndices[d] = indices[d] + Integer2DArray_Item ( range->cellIndices0, c, d ) ;
            RegularGrid_RegularizeIndices ( self, outsideGrid ) ;
            if ( outsideGrid == 0 )
            {
                RegularGrid_MakeCellID ( self, Integer1DArray_Item ( range->cellIDs, n ) ) ;
                Boolean1DArray_Item ( range->isFullyWithinRange, n ) = Boolean1DArray_Item ( range->isFullyWithinRange0, c ) ;
                n++ ;
            }
        }
        /* . Finish up. */
        range->numberOfCells = n ;
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find cells within range of a given point.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGrid_FindCellsWithinRangeOfPoint ( const RegularGrid *self, const Integer cellID, const Real *coordinates, RegularGridSearchRange *range, Status *status )
{
    Integer n = 0 ;
    if ( ( self != NULL ) && ( coordinates != NULL ) && ( range != NULL ) )
    {
        auto Integer c, d, *indices, *newIndices, outsideGrid ;
        auto Real    maximumDistanceSquared, minimumDistanceSquared ;
        /* . Get cell indices. */
        indices = range->workI ;
        RegularGrid_CellIDToIndices ( self, cellID, indices, status ) ;
        /* . Loop over possible cells within range. */
        newIndices = self->workI ;
        for ( c = 0 ; c < range->numberOfCells0 ; c++ )
        {
            for ( d = 0 ; d < self->ndimensions ; d++ ) newIndices[d] = indices[d] + Integer2DArray_Item ( range->cellIndices0, c, d ) ;
            RegularGrid_RegularizeIndices ( self, outsideGrid ) ;
            if ( outsideGrid == 0 )
            {
                RegularGrid_DistancesSquaredBetweenCellAndPoint ( self, newIndices, coordinates, &maximumDistanceSquared, &minimumDistanceSquared, status ) ;
                if ( minimumDistanceSquared < range->cutoffSquared ) /* <= ? */
                {
                    RegularGrid_MakeCellID ( self, Integer1DArray_Item ( range->cellIDs, n ) ) ;
                    Boolean1DArray_Item ( range->isFullyWithinRange, n ) = ( maximumDistanceSquared <= range->cutoffSquared ) ;
                    n++ ;
                }
            }
        }
        /* . Finish up. */
        range->numberOfCells = n ;
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find cells on a conforming grid within range of a given cell.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGrid_FindConformingCellsWithinRangeOfCell ( const RegularGrid *self, const Integer cellID, const RegularGrid *other, const Integer1DArray *offset,
                                                                                                                  RegularGridSearchRange *range, Status *status )
{
    Integer n = 0 ;
    if ( ( self != NULL ) && ( other != NULL ) && ( offset != NULL ) && ( range != NULL ) )
    {
        auto Integer c, d, *indices, *newIndices, outsideGrid ;
        /* . Get cell indices. */
        indices = range->workI ;
        RegularGrid_CellIDToIndices ( self, cellID, indices, status ) ;
        /* . Loop over possible cells within range. */
        newIndices = other->workI ;
        for ( c = 0 ; c < range->numberOfCells0 ; c++ )
        {
            for ( d = 0 ; d < self->ndimensions ; d++ ) newIndices[d] = indices[d] + Integer2DArray_Item ( range->cellIndices0, c, d ) - Integer1DArray_Item ( offset, d ) ;
            RegularGrid_RegularizeIndices ( other, outsideGrid ) ;
            if ( outsideGrid == 0 )
            {
                RegularGrid_MakeCellID ( other, Integer1DArray_Item ( range->cellIDs, n ) ) ;
                Boolean1DArray_Item ( range->isFullyWithinRange, n ) = Boolean1DArray_Item ( range->isFullyWithinRange0, c ) ;
                n++ ;
            }
        }
        /* . Finish up. */
        range->numberOfCells = n ;
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find cells on a conforming grid within range of a given point.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGrid_FindConformingCellsWithinRangeOfPoint ( const RegularGrid *self, const Integer cellID, const Real *coordinates, const RegularGrid *other,
                                                                                const Integer1DArray *offset, RegularGridSearchRange *range, Status *status )
{
    Integer n = 0 ;
    if ( ( self != NULL ) && ( coordinates != NULL ) && ( other != NULL ) && ( offset != NULL ) && ( range != NULL ) )
    {
        auto Integer c, d, *indices, *newIndices, outsideGrid ;
        auto Real    maximumDistanceSquared, minimumDistanceSquared ;
        /* . Get cell indices. */
        indices = range->workI ;
        RegularGrid_CellIDToIndices ( self, cellID, indices, status ) ;
        /* . Loop over possible cells within range. */
        newIndices = other->workI ;
        for ( c = 0 ; c < range->numberOfCells0 ; c++ )
        {
            for ( d = 0 ; d < self->ndimensions ; d++ ) newIndices[d] = indices[d] + Integer2DArray_Item ( range->cellIndices0, c, d ) - Integer1DArray_Item ( offset, d ) ;
            RegularGrid_RegularizeIndices ( other, outsideGrid ) ;
            if ( outsideGrid == 0 )
            {
                RegularGrid_DistancesSquaredBetweenCellAndPoint ( other, newIndices, coordinates, &maximumDistanceSquared, &minimumDistanceSquared, status ) ;
                if ( minimumDistanceSquared < range->cutoffSquared ) /* <= ? */
                {
                    RegularGrid_MakeCellID ( other, Integer1DArray_Item ( range->cellIDs, n ) ) ;
                    Boolean1DArray_Item ( range->isFullyWithinRange, n ) = ( maximumDistanceSquared <= range->cutoffSquared ) ;
                    n++ ;
                }
            }
        }
        /* . Finish up. */
        range->numberOfCells = n ;
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the coordinates of a cell given by an array of indices.
! . The indices do not have to be on the grid.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . No checking of dimensions is done. */
void RegularGrid_GetGridPointCoordinates ( const RegularGrid  *self, const Integer *indices, Real *coordinates )
{
    if ( ( self != NULL ) && ( indices != NULL ) && ( coordinates != NULL ) )
    {
        auto Integer d ;
        for ( d = 0 ; d < self->ndimensions ; d++ )
        {
            coordinates[d] = self->dimensions[d].lower + self->dimensions[d].binSize * ( ( Real ) indices[d] + 0.5e+00 ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Does the array have the same shape as the grid?
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean RegularGrid_IsConformingRealNDArray ( const RegularGrid *self,const RealNDArray *data )
{
    Boolean isSameShape = False ;
    if ( ( self != NULL ) && ( data != NULL ) )
    {
        if ( self->ndimensions == data->view->rank )
        {
            auto Integer d ;
            for ( d = 0 ; d < self->ndimensions ; d++ )
            {
                isSameShape = ( self->dimensions[d].bins == data->view->extents[d] ) ;
                if ( ! isSameShape ) break ;
            }
        }
    }
    return isSameShape ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Is the grid cubic?
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean RegularGrid_IsCubic ( const RegularGrid *self )
{
    Boolean isCubic = False ;
    if ( self != NULL )
    {
        auto Integer d ;
        for ( d = 1 ; d < self->ndimensions ; d++ )
        {
            isCubic = ( self->dimensions[d].binSize == self->dimensions[0].binSize ) ;
            if ( ! isCubic ) break ;
        }
    }
    return isCubic ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Is the grid periodic in at least one dimension?
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean RegularGrid_IsPeriodic ( const RegularGrid *self )
{
    Boolean isPeriodic = False ;
    if ( self != NULL )
    {
        auto Integer d ;
        for ( d = 0 ; d < self->ndimensions ; d++ )
        {
            isPeriodic = self->dimensions[d].isPeriodic ;
            if ( isPeriodic ) break ;
        }
    }
    return isPeriodic ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a periodic grid.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGrid *RegularGrid_MakePeriodicGrid3 ( const Vector3 *boxSize, const Real approximateGridSize, Status *status )
{
    RegularGrid *self = NULL ;
    if ( boxSize != NULL )
    {
        /* . Allocation. */
        self = RegularGrid_Allocate ( 3, status ) ;
        if ( self != NULL )
        {
            auto Integer d, stride ;
            auto Real    cells, cellsL, cellsU, extent, gridSize, lower, upper ;

            /* . Basic data. */
            for ( d = 0 ; d < 3 ; d++ )
            {
                lower    = 0.0e+00 ;
                upper    = Vector3_Item ( boxSize, d ) ;
                extent   = upper - lower ;
                cellsL   = floor ( extent / approximateGridSize ) ;
                cellsU   = ceil  ( extent / approximateGridSize ) ;
                cells    = ( ( fabs ( approximateGridSize - extent / cellsL ) < fabs ( approximateGridSize - extent / cellsU ) ) ? cellsL : cellsU ) ;
                gridSize = extent / cells ;
                self->dimensions[d].bins          = ( Integer ) cells ;
                self->dimensions[d].binSize       = gridSize          ;
                self->dimensions[d].isPeriodic    = True              ;
                self->dimensions[d].lower         = lower             ;
                self->dimensions[d].midPointLower = lower + 0.5e+00 * gridSize ;
                self->dimensions[d].period        = extent            ;
                self->dimensions[d].upper         = upper             ;
            }

            /* . Remaining data (do in reverse order for the strides). */
            stride = 1 ;
            for ( d = 2 ; d >= 0 ; d-- ) { self->dimensions[d].stride = stride ; stride *= self->dimensions[d].bins ; }
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a search range given a cutoff.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGridSearchRange *RegularGrid_MakeSearchRange ( const RegularGrid *self, const Real cutoff, Status *status )
{
    RegularGridSearchRange *range = NULL ;
    if ( ( self != NULL ) && ( cutoff > 0.0 ) )
    {
        auto Integer d, *lower, numberOfCells0, *upper ;
        /* . Determine the maximum number of possible cells within the cutoff. */
        lower = &(self->workI[0]) ;
        upper = &(self->workI[3]) ;
        for ( d = 0, numberOfCells0 = 1 ; d < self->ndimensions ; d++ )
        {
            RegularGridDimension_BinsWithinCutoff ( self->dimensions[d], cutoff, lower[d], upper[d] ) ;
            numberOfCells0 *= ( upper[d] - lower[d] + 1 ) ;
        }
        /* . Create the range and fill it with zeroth-order data. */
        range = RegularGridSearchRange_Allocate ( numberOfCells0, self->ndimensions, status ) ;
        if ( range != NULL )
        {
            auto Integer c, index, *indices, n, *origin ;
            auto Real    maximumDistanceSquared, minimumDistanceSquared ;
            /* . Set the cutoffs. */
            range->cutoff        = cutoff ;
            range->cutoffSquared = cutoff * cutoff ;
            /* . Set the origin - use last available slot in cellIndices0. */
            origin = Integer2DArray_RowPointer ( range->cellIndices0, numberOfCells0 - 1 ) ;
            for ( d = 0 ; d < self->ndimensions ; d++ ) origin[d] = 0 ;
            /* . Set the data for cells within range. */
            indices = range->workI ;
            for ( d = 0 ; d < self->ndimensions ; d++ ) indices[d] = lower[d] ;
            for ( c = n = 0 ; c < numberOfCells0 ; c++ )
            {
                RegularGrid_DistancesSquaredBetweenCells ( self, origin, indices, &maximumDistanceSquared, &minimumDistanceSquared, status ) ;
                if ( minimumDistanceSquared < range->cutoffSquared ) /* <= ? */
                {
                    for ( d = 0 ; d < self->ndimensions ; d++ ) Integer2DArray_Item ( range->cellIndices0, n, d ) = indices[d] ;
                    Boolean1DArray_Item ( range->isFullyWithinRange0, n ) = ( maximumDistanceSquared <= range->cutoffSquared ) ;
                    n++ ;
                }
                for ( d = self->ndimensions-1 ; d >= 0 ; d-- )
                {
                   index = indices[d] + 1 ;
                   if ( index > upper[d] ) { indices[d] = lower[d] ; }
                   else                    { indices[d] = index ; break ;   }
                }
            }
            range->numberOfCells0 = n ;
        }
    }
    return range ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Number of grid points.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGrid_NumberOfGridPoints ( const RegularGrid *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0, n = 1 ; i < self->ndimensions ; i++ ) n *= self->dimensions[i].bins ;
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGrid_Print ( const RegularGrid *self )
{
    if ( self != NULL )
    {
        auto Integer c, n ;
        auto RegularGridDimension *d ;
        n = self->ndimensions ;
        printf ( "\nRegular Grid:\n" ) ;
        for ( c = 0 ; c < n ; c++ )
        {
            d = &(self->dimensions[c]) ;
            printf ( "%10d %10d %10d %10.3f %10.3f %10.3f %10.3f %10.3f", c, d->bins, d->stride, d->binSize, d->lower, d->midPointLower, d->period, d->upper ) ;
            if ( d->isPeriodic ) printf ( "    Yes\n" ) ; else printf ( "    No\n" ) ;
        }
        printf ( "\n" ) ;
    }
}

/*==================================================================================================================================
! . Regular Grid Dimension Functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridDimension_CopyTo ( const RegularGridDimension *self, RegularGridDimension *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        other->isPeriodic    = self->isPeriodic    ;
        other->bins          = self->bins          ;
        other->stride        = self->stride        ;
        other->binSize       = self->binSize       ;
        other->lower         = self->lower         ;
        other->midPointLower = self->midPointLower ;
        other->period        = self->period        ;
        other->upper         = self->upper         ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridDimension_Initialize ( RegularGridDimension *self )
{
    if ( self != NULL )
    {
        self->isPeriodic    = False ;
        self->bins          =     0 ;
        self->stride        =     0 ;
        self->binSize       =   0.0 ;
        self->lower         =   0.0 ;
        self->midPointLower =   0.0 ;
        self->period        =   0.0 ;
        self->upper         =   0.0 ;
    }
}

/*==================================================================================================================================
! . Regular Grid Search Range Functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGridSearchRange *RegularGridSearchRange_Allocate ( const Integer cells, const Integer dimensions, Status *status )
{
    RegularGridSearchRange *self = NULL ;
    MEMORY_ALLOCATE ( self, RegularGridSearchRange ) ;
    if ( self != NULL )
    {
        /* . Initialization. */
        self->numberOfCells  = 0 ;
        self->numberOfCells0 = 0 ;
        self->cutoff         = 0.0e+00 ;
        self->cutoffSquared  = 0.0e+00 ;
        /* . Allocate space. */
        self->isFullyWithinRange  = Boolean1DArray_Allocate ( cells, status ) ;
        self->isFullyWithinRange0 = Boolean1DArray_Allocate ( cells, status ) ;
        self->cellIDs             = Integer1DArray_Allocate ( cells, status ) ;
        self->cellIndices0        = Integer2DArray_Allocate ( cells, dimensions, status ) ;
        self->workI               = Memory_Allocate_Array_Integer ( dimensions ) ;
        if ( ( self->isFullyWithinRange == NULL ) || ( self->isFullyWithinRange0 == NULL ) || ( self->cellIDs == NULL ) || 
             ( self->cellIndices0 == NULL ) || ( self->workI == NULL ) ) RegularGridSearchRange_Deallocate ( &self ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridSearchRange_Deallocate ( RegularGridSearchRange **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Boolean1DArray_Deallocate ( &((*self)->isFullyWithinRange ) ) ;
        Boolean1DArray_Deallocate ( &((*self)->isFullyWithinRange0) ) ;
        Integer1DArray_Deallocate ( &((*self)->cellIDs            ) ) ;
        Integer2DArray_Deallocate ( &((*self)->cellIndices0       ) ) ;
        Memory_Deallocate_Integer ( &((*self)->workI              ) ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridSearchRange_Print ( const RegularGridSearchRange *self )
{
    if ( self != NULL )
    {
        auto Integer c, i, n ;
        n = self->cellIndices0->length1 ;
        printf ( "\nRegular Grid Search Range:\n" ) ;
        printf ( "Number of Cells = %20d\n"  , self->numberOfCells0 ) ;
        printf ( "Cutoff          = %20.5f\n", self->cutoff         ) ;
        printf ( "Cutoff Squared  = %20.5f\n", self->cutoffSquared  ) ;
        for ( c = 0 ; c < self->numberOfCells0 ; c++ )
        {
            printf ( "%10d", c ) ;
            for ( i = 0 ; i < n ; i++ ) printf ( "%10d", Integer2DArray_Item ( self->cellIndices0, c, i ) ) ;
            if ( Boolean1DArray_Item ( self->isFullyWithinRange0, c ) ) printf ( "    Yes\n" ) ; else printf ( "    No\n" ) ;
        }
        printf ( "\n" ) ;
    }
}
