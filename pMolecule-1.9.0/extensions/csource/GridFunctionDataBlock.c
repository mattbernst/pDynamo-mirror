/*------------------------------------------------------------------------------
! . File      : GridFunctionDataBlock.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module defines a data structure for storing function and derivative values on a grid.
!=================================================================================================================================*/

# include "GridFunctionDataBlock.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
GridFunctionDataBlock *GridFunctionDataBlock_Allocate ( const Integer numberOfFunctions, const Integer numberOfPoints, const Integer order, Status *status )
{
    GridFunctionDataBlock *self = NULL ;
    MEMORY_ALLOCATE ( self, GridFunctionDataBlock ) ;
    if ( self != NULL )
    {
        auto Integer f, o, p ;
        auto Status  localStatus = Status_Continue ;
        f = Maximum ( numberOfFunctions   , 0 ) ; self->numberOfFunctions = f ;
        p = Maximum ( numberOfPoints      , 0 ) ; self->numberOfPoints    = p ;
        o = Minimum ( Maximum ( order, 0 ), 3 ) ; self->order             = o ;
        self->indices = Integer1DArray_Allocate ( f,    &localStatus ) ;
        self->f       = Real2DArray_Allocate    ( f, p, &localStatus ) ;
        if ( o > 0 )
        {
           self->fX   = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fY   = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fZ   = Real2DArray_Allocate ( f, p, &localStatus ) ;
        }
        else { self->fX = NULL ; self->fY = NULL ; self->fZ = NULL ; }
        if ( o > 1 )
        {
           self->fXX  = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fXY  = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fXZ  = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fYY  = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fYZ  = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fZZ  = Real2DArray_Allocate ( f, p, &localStatus ) ;
        }
        else
        {
            self->fXX = NULL ; self->fXY = NULL ; self->fXZ = NULL ;
            self->fYY = NULL ; self->fYZ = NULL ; self->fZZ = NULL ;
        }
        if ( o > 2 )
        {
           self->fXXX = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fXXY = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fXXZ = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fXYY = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fXYZ = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fXZZ = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fYYY = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fYYZ = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fYZZ = Real2DArray_Allocate ( f, p, &localStatus ) ;
           self->fZZZ = Real2DArray_Allocate ( f, p, &localStatus ) ;
        }
        else
        {
            self->fXXX = NULL ; self->fXXY = NULL ; self->fXXZ = NULL ; self->fXYY = NULL ; self->fXYZ = NULL ;
            self->fXZZ = NULL ; self->fYYY = NULL ; self->fYYZ = NULL ; self->fYZZ = NULL ; self->fZZZ = NULL ;
        }
        if ( localStatus != Status_Continue ) GridFunctionDataBlock_Deallocate ( &self ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The size in bytes of the data block.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real GridFunctionDataBlock_ByteSize ( const GridFunctionDataBlock *self )
{
    Real size = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer f, n, p ;
        f = self->f->length0 ;
        p = self->f->length1 ;
        n = 1 ;
        if ( self->order > 0 ) n +=  3 ;
        if ( self->order > 1 ) n +=  6 ;
        if ( self->order > 2 ) n += 10 ;
        size  = sizeof ( GridFunctionDataBlock ) + ( sizeof ( Integer1DArray ) + sizeof ( Integer ) * f ) + n * ( sizeof ( Real2DArray ) + sizeof ( Real ) * f * p ) ;
    }
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GridFunctionDataBlock_Deallocate ( GridFunctionDataBlock **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Integer1DArray_Deallocate ( &((*self)->indices) ) ;
        Real2DArray_Deallocate    ( &((*self)->f      ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fX     ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fY     ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fZ     ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fXX    ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fXY    ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fXZ    ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fYY    ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fYZ    ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fZZ    ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fXXX   ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fXXY   ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fXXZ   ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fXYY   ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fXYZ   ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fXZZ   ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fYYY   ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fYYZ   ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fYZZ   ) ) ;
        Real2DArray_Deallocate    ( &((*self)->fZZZ   ) ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Filtering.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GridFunctionDataBlock_FilterValues ( GridFunctionDataBlock *self, const Integer fStart, const Real *tolerance )
{
    if ( ( self != NULL ) && ( tolerance != NULL ) && ( (*tolerance) > 0.0e+00 ) ) ;
    {
        auto Integer     f, f0, f1, n, o ;
        auto Real        t   ;
        auto Real1DArray row ;
        f0 = Maximum ( fStart, 0 ) ;
        f1 = self->numberOfFunctions ;
        o  = self->order ;
        t  = (*tolerance) ;
        for ( f = n = f0 ; f < f1 ; f++ )
        {
            Real2DArray_RowSlice ( self->f, f, &row, NULL ) ;
            if ( Real1DArray_AbsoluteMaximum ( &row ) > t )
            {
                if ( f != n )
                {
                    Integer1DArray_Item ( self->indices, n ) = Integer1DArray_Item ( self->indices, f ) ;
                    Real2DArray_RowCopyTo ( self->f, f, self->f, n, NULL ) ;
                    if ( o > 0 )
                    {
                        Real2DArray_RowCopyTo ( self->fX, f, self->fX, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fY, f, self->fY, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fZ, f, self->fZ, n, NULL ) ;
                    }
                    if ( o > 1 )
                    {
                        Real2DArray_RowCopyTo ( self->fXX, f, self->fXX, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fXY, f, self->fXY, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fXZ, f, self->fXZ, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fYY, f, self->fYY, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fYZ, f, self->fYZ, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fZZ, f, self->fZZ, n, NULL ) ;
                    }
                    if ( o > 2 )
                    {
                        Real2DArray_RowCopyTo ( self->fXXX, f, self->fXXX, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fXXY, f, self->fXXY, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fXXZ, f, self->fXXZ, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fXYY, f, self->fXYY, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fXYZ, f, self->fXYZ, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fXZZ, f, self->fXZZ, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fYYY, f, self->fYYY, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fYYZ, f, self->fYYZ, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fYZZ, f, self->fYZZ, n, NULL ) ;
                        Real2DArray_RowCopyTo ( self->fZZZ, f, self->fZZZ, n, NULL ) ;
                    }
                }
                n++ ;
            }
        }
        self->numberOfFunctions = n ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GridFunctionDataBlock_Initialize ( GridFunctionDataBlock *self )
{
    if ( self != NULL )
    {
        auto Integer f ;
        self->numberOfFunctions = 0 ;
        for ( f = 0 ; f < self->indices->length ; f++ ) Integer1DArray_Item ( self->indices, f ) = f ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GridFunctionDataBlock_Resize ( GridFunctionDataBlock *self, const Integer numberOfFunctions, Status *status )
{
    if ( ( self != NULL ) && ( self->f->length0 != numberOfFunctions ) )
    {
        self->numberOfFunctions = Minimum ( self->numberOfFunctions, numberOfFunctions ) ;
        Integer1DArray_Resize ( self->indices, numberOfFunctions, NULL, status ) ;
        Real2DArray_Resize    ( self->f      , numberOfFunctions, NULL, status ) ;
        if ( self->order > 0 )
        {
            Real2DArray_Resize ( self->fX  , numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fY  , numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fZ  , numberOfFunctions, NULL, status ) ;
        }
        if ( self->order > 1 )
        {
            Real2DArray_Resize ( self->fXX , numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fXY , numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fXZ , numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fYY , numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fYZ , numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fZZ , numberOfFunctions, NULL, status ) ;
        }
        if ( self->order > 2 )
        {
            Real2DArray_Resize ( self->fXXX, numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fXXY, numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fXXZ, numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fXYY, numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fXYZ, numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fXZZ, numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fYYY, numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fYYZ, numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fYZZ, numberOfFunctions, NULL, status ) ;
            Real2DArray_Resize ( self->fZZZ, numberOfFunctions, NULL, status ) ;
        }
    }
}
