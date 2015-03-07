/*------------------------------------------------------------------------------
! . File      : PolygonalSurface.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Simple surface type.
!=================================================================================================================================*/

# include <stdio.h>

# include "Macros.h"
# include "Memory.h"
# include "PolygonalSurface.h"
# include "Real1DArray.h"

# define INDEXINITIALIZER -1
# define VALUEINITIALIZER 0.0e+00

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PolygonalSurface *PolygonalSurface_Allocate ( const Integer ndimensions, const Integer nvertices, const Integer npolygons, Status *status )
{
    PolygonalSurface *self = NULL ;
    /* . Check the dimension (must be > 1). */
    if ( ndimensions > 1 )
    {
        MEMORY_ALLOCATE ( self, PolygonalSurface ) ;
        if ( self != NULL )
        {
            auto Integer n ;
            /* . Initialization. */
            self->polygons = NULL ;
            self->normals  = NULL ;
            self->vertices = NULL ;
            /* . Allocation. */
            n = Maximum ( npolygons, 0 ) ;
            self->polygons = Integer2DArray_Allocate ( n, ndimensions, status ) ;
            n = Maximum ( nvertices, 0 ) ;
            self->normals  = Real2DArray_Allocate    ( n, ndimensions, status ) ;
            self->vertices = Real2DArray_Allocate    ( n, ndimensions, status ) ;
            if ( ( self->polygons == NULL ) || ( self->normals == NULL ) || ( self->vertices == NULL ) ) PolygonalSurface_Deallocate ( &self ) ;
        }
        if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    else Status_Set ( status, Status_InvalidDimension ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
PolygonalSurface *PolygonalSurface_Clone ( const PolygonalSurface *self, Status *status )
{
    PolygonalSurface *clone = NULL ;
    if ( self != NULL )
    {
        clone = PolygonalSurface_Allocate ( PolygonalSurface_Rank             ( self ),
                                            PolygonalSurface_NumberOfVertices ( self ),
                                            PolygonalSurface_NumberOfPolygons ( self ), status ) ;
        PolygonalSurface_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Construct the normals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PolygonalSurface_ConstructNormals ( PolygonalSurface *self )
{
    if ( self != NULL )
    {
        ; /* . Do nothing for now. */
/* . Normals can be constructed from grid data (gradients). Alternatively from averages of the polygon normals with (optional)
weights proportional to their areas or their inverse areas. For triangle, if f = axb + bxc + cxa, the normal is f(normalized)
and the area is |f|/2.
*/
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PolygonalSurface_CopyTo ( const PolygonalSurface *self, PolygonalSurface *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        Integer2DArray_CopyTo ( self->polygons, other->polygons, status ) ;
        Real2DArray_CopyTo    ( self->normals,  other->normals,  status ) ;
        Real2DArray_CopyTo    ( self->vertices, other->vertices, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PolygonalSurface_Deallocate ( PolygonalSurface **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Integer2DArray_Deallocate ( &((*self)->polygons) ) ;
        Real2DArray_Deallocate    ( &((*self)->normals ) ) ;
        Real2DArray_Deallocate    ( &((*self)->vertices) ) ;
        MEMORY_DEALLOCATE ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialize the arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PolygonalSurface_InitializeArrays ( PolygonalSurface *self )
{
    if ( self != NULL )
    {
        Integer2DArray_Set ( self->polygons, INDEXINITIALIZER ) ;
        Real2DArray_Set    ( self->normals,  VALUEINITIALIZER ) ;
        Real2DArray_Set    ( self->vertices, VALUEINITIALIZER ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalize the normals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PolygonalSurface_NormalizeNormals ( PolygonalSurface *self )
{
    if ( self != NULL )
    {
        auto Integer i ;
        auto Real1DArray row ;
        for ( i = 0 ; i < PolygonalSurface_NumberOfVertices ( self ) ; i++ )
        {
            Real2DArray_RowSlice ( self->normals, i, &row, NULL ) ;
            Real1DArray_Normalize ( &row, NULL, NULL ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Number of polygons.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer PolygonalSurface_NumberOfPolygons ( const PolygonalSurface *self )
{
    if ( self == NULL ) return 0 ;
    else                return Integer2DArray_Length ( self->polygons, 0 ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Number of vertices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer PolygonalSurface_NumberOfVertices ( const PolygonalSurface *self )
{
    if ( self == NULL ) return 0 ;
    else                return Real2DArray_Length ( self->vertices, 0 ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Origin and extents.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PolygonalSurface_OriginAndExtents ( const PolygonalSurface *self, Real1DArray *origin, Real1DArray *extents, Status *status )
{
    if ( ( self != NULL ) && ( ( extents != NULL ) || ( origin != NULL ) ) )
    {
        auto Integer n ;
        n = PolygonalSurface_Rank ( self ) ;
        if ( ( ( extents != NULL ) && ( n != Real1DArray_Length ( extents ) ) ) || ( ( origin != NULL ) && ( n != Real1DArray_Length ( origin ) ) ) ) Status_Set ( status, Status_ArrayDimensionMismatch ) ;
        else
        {
            auto Integer d ;
            auto Real lower, upper ;
            auto Real1DArray column ;
            for ( d = 0 ; d < n ; d++ )
            {
                Real2DArray_ColumnSlice ( self->vertices, d, &column, NULL ) ;
                Real1DArray_Range ( &column, &lower, &upper ) ;
                Real1DArray_Item ( extents, d ) = upper - lower ;
                Real1DArray_Item ( origin , d ) =         lower ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Rank.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer PolygonalSurface_Rank ( const PolygonalSurface *self )
{
    if ( self == NULL ) return 0 ;
    else                return Real2DArray_Length ( self->vertices, 1 ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PolygonalSurface_Resize ( PolygonalSurface *self, const Integer nvertices, const Integer npolygons, const Boolean initialize, Status *status )
{
    if ( self != NULL )
    {
        auto Integer i = INDEXINITIALIZER, *ivalue = NULL ;
        auto Real    r = VALUEINITIALIZER, *rvalue = NULL ;
        if ( initialize ) { ivalue = &i ; rvalue = &r ; }
        Integer2DArray_Resize ( self->polygons, npolygons, ivalue, status ) ;
        Real2DArray_Resize    ( self->normals,  nvertices, rvalue, status ) ;
        Real2DArray_Resize    ( self->vertices, nvertices, rvalue, status ) ;
    }
}

# undef INDEXINITIALIZER
# undef VALUEINITIALIZER
