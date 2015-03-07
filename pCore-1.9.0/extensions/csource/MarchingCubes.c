/*------------------------------------------------------------------------------
! . File      : MarchingCubes.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . Marching cubes algorithm.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
!
! . Code for marching cubes adapted from the following:
!
! * @file    MarchingCubes.cpp
! * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
! * @author  Math Dept, PUC-Rio
! * @version 0.2
! * @date    12/08/2002
! *
! * @brief   MarchingCubes Algorithm
!
! . Counterclockwise vertex order for the triangles.
!
!---------------------------------------------------------------------------------------------------------------------------------*/

# include <math.h>
# include <stdlib.h>
# include <stdio.h>

# include "Boolean.h"
# include "Integer.h"
# include "Integer2DArrayN3.h"
# include "IntegerNDArray.h"
# include "Macros.h"
# include "MarchingCubes.h"
# include "Real2DArrayN3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A small value for face and interior testing. */
# define EPSILON 1.0e-10

/* . Factor for estimating the starting number of polygons. */
# define POLYGONFACTOR0 4

/* . Polygon number increments. */
# define POLYGONINCREMENT 5000

/* . A small value to reduce numerical problems in linear interpolation. */
# define SAFEMINIMUM 1.0e-10

/* . Starting vertex count. */
# define VERTEXCOUNT0 10000

/* . Vertex number increments. */
# define VERTEXINCREMENT1 10000
# define VERTEXINCREMENT2  5000

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void    AddTriangle         ( Integer *ntriangles, Integer2DArray *triangles, const IntegerNDArray *intersections, const Integer i, const Integer j, const Integer k, const Integer8 *trig, const Integer8 n, const Integer v12  ) ;
static Integer AddInteriorVertex   ( const Integer i, const Integer j, const Integer k, const IntegerNDArray *intersections, Integer *nvertices, Real2DArray *vertices, Real2DArray *normals ) ;
static void    GetGradient         ( const RealNDArray *data, const Integer i, const Integer j, const Integer k, const Integer ni, const Integer nj, const Integer nk, Real *gx, Real *gy, Real *gz ) ;
static Real    LinearlyInterpolate ( const Real f0, const Real f1, const Real u0 ) ;
static void    PrintCube           ( const Real *cube ) ;
static void    ProcessCube         ( const Integer i, const Integer j, const Integer k, const Real *cube, const Integer8 cubecase, const Integer8 configuration, const IntegerNDArray *intersections,
                                                                                  Integer *nvertices, Real2DArray *vertices, Real2DArray *normals, Integer *ntriangles, Integer2DArray *triangles ) ;
static Boolean TestFace            ( const Real *cube, const Integer8 face ) ;
static Boolean TestInterior        ( const Real *cube, const Integer8 cubecase, const Integer8 configuration, const Integer8 subconfiguration, const Integer8 s ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate an isosurface given data on a regular grid.
!---------------------------------------------------------------------------------------------------------------------------------*/
PolygonalSurface *MarchingCubes_GenerateIsosurface ( const RegularGrid *grid, const RealNDArray *data, const Real isovalue, Status *status )
{
    PolygonalSurface *surface = NULL ;
    if ( ( grid != NULL ) && ( data != NULL ) )
    {
        /* . Check that the grid has three dimensions and the data is compatible with the grid. */
        if ( ( grid->ndimensions == 3 ) && RegularGrid_IsConformingRealNDArray ( grid, data ) )
        {
            auto Boolean         doResize ;
            auto Cardinal8       tableentry ;
            auto Integer         d, i, index, indices[4], j, k, ncubes, ni, nj, nk, nold, npolygons, npolygons0, nvertices, nvertices0, p ;
            auto Integer8        cubecase ;
            auto Integer2DArray *polygons = NULL ;
            auto IntegerNDArray *intersections = NULL ;
            auto Real            cube[8], f0, faxis[3], gx, gy, gz, hx, hy, hz, u ;
            auto Real1DArray     column ;
            auto Real2DArray    *normals = NULL, *vertices = NULL ;
            auto Status          localstatus = Status_Continue ;

            /* . Get the grid extents. */
            ni = data->view->extents[0] ;
            nj = data->view->extents[1] ;
            nk = data->view->extents[2] ;

            /* . Estimate the number of vertices for the surface - be conservative. */
            ncubes     = RegularGrid_NumberOfGridPoints ( grid ) ;
            nvertices0 = Minimum ( 3 * ncubes, VERTEXCOUNT0 ) ;

            /* . Allocate an initial surface with an approximate size. */
            surface = PolygonalSurface_Allocate ( 3, nvertices0, 0, &localstatus ) ;
            PolygonalSurface_InitializeArrays ( surface ) ;

            /* . Allocate the intersections. */
            indices[0] = ni ;
            indices[1] = nj ;
            indices[2] = nk ;
            indices[3] =  3 ;
            intersections = IntegerNDArray_Allocate ( 4, indices, &localstatus ) ;
            IntegerNDArray_Set ( intersections, -1 ) ;

            /* . Ensure that everything is OK. */
            if ( ! Status_OK ( &localstatus ) ) goto FinishUp ;

            /* . Aliases. */
            normals  = surface->normals  ;
            polygons = surface->polygons ;
            vertices = surface->vertices ;

            /* . Compute intersections for each cube along the cube edges (almost all of them). */
            for ( i = nvertices = 0 ; i < ni ; i++ )
            {
                for ( j = 0 ; j < nj ; j++ )
                {
                    for ( k = 0 ; k < nk ; k++ )
                    {
                        /* . Make sure that there is enough space for the maximum three vertices that can be added here. */
                        if ( nvertices + 3 > nvertices0 )
                        {
                            nvertices0 += VERTEXINCREMENT1 ;
                            PolygonalSurface_Resize ( surface, nvertices0, 0, True, &localstatus ) ;
                            if ( ! Status_OK ( &localstatus ) ) goto FinishUp ;
                        }
                        /* . Get the function values along the lower corner of the cube. */
                        f0 = RealNDArray_Item3D ( data, i, j, k ) - isovalue ;
                        if ( i < ni - 1 ) faxis[0] = RealNDArray_Item3D ( data, i+1, j, k ) - isovalue ;
                        else              faxis[0] = f0 ;
                        if ( j < nj - 1 ) faxis[1] = RealNDArray_Item3D ( data, i, j+1, k ) - isovalue ;
                        else              faxis[1] = f0 ;
                        if ( k < nk - 1 ) faxis[2] = RealNDArray_Item3D ( data, i, j ,k+1 ) - isovalue ;
                        else              faxis[2] = f0 ;
                        /* . Determine if there are any vertices. */
                        nold = nvertices ;
                        if ( f0 < 0 )
                        {
                            for ( d = 0 ; d < 3 ; d++ )
                            {
                                if ( faxis[d] >= 0 )
                                {
                                    Real2DArray_Item ( vertices, nvertices, d ) = LinearlyInterpolate ( f0, faxis[d], 1.0e+00 ) ;
                                    IntegerNDArray_Item4D ( intersections, i, j, k, d ) = nvertices ;
                                    nvertices ++ ;
                                }
                            }
                        }
                        else
                        {
                            for ( d = 0 ; d < 3 ; d++ )
                            {
                                if ( faxis[d] < 0 )
                                {
                                    Real2DArray_Item ( vertices, nvertices, d ) = LinearlyInterpolate ( f0, faxis[d], 0.0e+00 ) ;
                                    IntegerNDArray_Item4D ( intersections, i, j, k, d ) = nvertices ;
                                    nvertices ++ ;
                                }
                            }
                        }
                        /* . Process the vertex information. */
                        if ( nvertices > nold )
                        {
                            /* . Get the gradient at i, j, k. */
                            GetGradient ( data, i, j, k, ni, nj, nk, &gx, &gy, &gz ) ;
                            /* . Loop over the dimensions. */
                            for ( d = 0 ; d < 3 ; d++ )
                            {
                                index = IntegerNDArray_Item4D ( intersections, i, j, k, d ) ;
                                if ( index >= 0 )
                                {
                                    u = Real2DArray_Item ( vertices, index, d ) ;
                                    Real2DArray_IncrementRowN3 ( vertices, index, ( Real ) i, ( Real ) j, ( Real ) k ) ;
                                    /* . Get the gradient along the appropriate axis. */
                                    GetGradient ( data, i + IJKTERMS[d][0], j + IJKTERMS[d][1], k + IJKTERMS[d][2], ni, nj, nk, &hx, &hy, &hz ) ;
                                    /* . Determine the normal. */
                                    Real2DArray_SetRowN3 ( normals, index, ( 1.0e+00 - u ) * gx + u * hx, ( 1.0e+00 - u ) * gy + u * hy, ( 1.0e+00 - u ) * gz + u * hz ) ;
                                }
                            }
                        }
                    }
                }
            }

            /* . Allocate space for the polygons. */
            npolygons0 = POLYGONFACTOR0 * nvertices ;
            PolygonalSurface_Resize ( surface, nvertices0, npolygons0, True, &localstatus ) ;
            if ( ! Status_OK ( &localstatus ) ) goto FinishUp ;

            /* . Process each cube. */
            for ( i = npolygons = 0 ; i < ni - 1 ; i++ )
            {
                for ( j = 0 ; j < nj - 1 ; j++ )
                {
                    for ( k = 0 ; k < nk - 1 ; k++ )
                    {
                        /* . Make sure that there is enough space for the maximum one vertex and twelve polygons that can be added here. */
                        doResize = False ;
                        if ( npolygons + 12 > npolygons0 ) { npolygons0 += POLYGONINCREMENT ; doResize = True ; }
                        if ( nvertices +  1 > nvertices0 ) { nvertices0 += VERTEXINCREMENT2 ; doResize = True ; }
                        if ( doResize )
                        {
                            PolygonalSurface_Resize ( surface, nvertices0, npolygons0, True, &localstatus ) ;
                            if ( ! Status_OK ( &localstatus ) ) goto FinishUp ;
                        }
                        /* . Determine to which case the cube belongs (values between 0 and 255). */
                        tableentry = 0 ;
                        for ( p = 0 ; p < 8 ; ++p )
                        {
                            cube[p] = RealNDArray_Item3D ( data,  i+((p^(p>>1))&1), j+((p>>1)&1), k+((p>>2)&1) ) - isovalue ;
                            if ( cube[p] >= 0 ) tableentry += 1 << p ;
                        }
                        /* . Process the cube if necessary. */
                        cubecase = CUBECASES[tableentry][0] ;
                        if ( cubecase > 0 ) ProcessCube ( i, j, k, cube, cubecase, CUBECASES[tableentry][1], intersections, &nvertices, vertices, normals, &npolygons, polygons ) ;
                    }
                }
            }

            /* . Remove unused space from the surface data structure. */
            PolygonalSurface_Resize ( surface, nvertices, npolygons, False, &localstatus ) ;
            if ( ! Status_OK ( &localstatus ) ) goto FinishUp ;

            /* . Scale the vertices and normals to the grid coordinates. */
            for ( d = 0 ; d < 3 ; d++ )
            {
                Real2DArray_ColumnSlice ( normals,  d, &column, NULL ) ;
                Real1DArray_Scale ( &column, grid->dimensions[d].binSize ) ;
                Real2DArray_ColumnSlice ( vertices, d, &column, NULL ) ;
                Real1DArray_Scale ( &column, grid->dimensions[d].binSize ) ;
            }

            /* . Translate the vertex coordinates to the correct origin. */
            for ( d = 0 ; d < 3 ; d++ )
            {
                Real2DArray_ColumnSlice ( vertices, d, &column, NULL ) ;
                Real1DArray_AddScalar  ( &column, grid->dimensions[d].midPointLower ) ;
            }

            /* . Normalize the normals. */
            PolygonalSurface_NormalizeNormals ( surface ) ;

            /* . Finish up. */
        FinishUp:
            IntegerNDArray_Deallocate ( &intersections ) ;
            if ( ! Status_OK ( &localstatus ) ) PolygonalSurface_Deallocate ( &surface ) ;
            Status_Set ( status, localstatus ) ;
        }
        /* . Argument error. */
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return surface ;
}

/*==================================================================================================================================
! . Marching cubes auxiliary procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Add triangles.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void AddTriangle ( Integer *ntriangles, Integer2DArray *triangles, const IntegerNDArray *intersections, const Integer i, const Integer j, const Integer k, const Integer8 *trig, const Integer8 n, const Integer v12  )
{
    Integer t, tv[3] ;
    for( t = 0 ; t < 3*n ; t++ )
    {
        switch ( trig[t] )
        {
            case  0 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i  , j  , k  , 0 ) ; break ;
            case  1 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i+1, j  , k  , 1 ) ; break ;
            case  2 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i  , j+1, k  , 0 ) ; break ;
            case  3 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i  , j  , k  , 1 ) ; break ;
            case  4 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i  , j  , k+1, 0 ) ; break ;
            case  5 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i+1, j  , k+1, 1 ) ; break ;
            case  6 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i  , j+1, k+1, 0 ) ; break ;
            case  7 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i  , j  , k+1, 1 ) ; break ;
            case  8 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i  , j  , k  , 2 ) ; break ;
            case  9 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i+1, j  , k  , 2 ) ; break ;
            case 10 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i+1, j+1, k  , 2 ) ; break ;
            case 11 : tv[ t % 3 ] = IntegerNDArray_Item4D ( intersections, i  , j+1, k  , 2 ) ; break ;
            case 12 : tv[ t % 3 ] = v12 ; break ;
            default : break ;
        }
        if ( tv[t%3] == -1 ) { printf ( "Marching Cubes: invalid triangle %d %d %d %d\n", (*ntriangles) + 1, i, j, k ) ; }
        if ( t%3 == 2 )
        {
            Integer2DArray_SetRowN3 ( triangles, (*ntriangles), tv[0], tv[1], tv[2] ) ;
            (*ntriangles)++ ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add an interior vertex using the average of the intersection points of a cube.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer AddInteriorVertex ( const Integer i, const Integer j, const Integer k, const IntegerNDArray *intersections, Integer *nvertices, Real2DArray *vertices, Real2DArray *normals )
{
    Integer c, current, di, dj, dk, n, t, vid ;
    Real    scale, x, y, z ;
    current = (*nvertices) ;
    for ( c = n = 0 ; c < 3 ; c++ )
    {
        for ( t = 0 ; t < 4 ; t++ )
        {
            di  = i + IVERTEXTERMS[c][t][0] ;
            dj  = j + IVERTEXTERMS[c][t][1] ;
            dk  = k + IVERTEXTERMS[c][t][2] ;
            vid = IntegerNDArray_Item4D ( intersections, di, dj, dk, c ) ;
            if ( vid != -1 )
            {
                Real2DArray_GetRowN3       ( vertices, vid,     x, y, z ) ;
                Real2DArray_IncrementRowN3 ( vertices, current, x, y, z ) ;
                Real2DArray_GetRowN3       ( normals,  vid,     x, y, z ) ;
                Real2DArray_IncrementRowN3 ( normals,  current, x, y, z ) ;
                n++ ;
            }
        }
    }
    scale = 1.0e+00 / ( Real ) n ;
    Real2DArray_ScaleRowN3 ( vertices, current, scale ) ;
    Real2DArray_ScaleRowN3 ( vertices, current, scale ) ;
    (*nvertices)++ ;
    return current ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gradient calculation at a grid point using finite differences.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GetGradient ( const RealNDArray *data, const Integer i, const Integer j, const Integer k, const Integer ni, const Integer nj, const Integer nk, Real *gx, Real *gy, Real *gz )
{
    if ( i > 0 )
    {
        if ( i < ni - 1 ) (*gx) = 0.5e+00 * ( RealNDArray_Item3D ( data, i+1, j, k ) - RealNDArray_Item3D ( data, i-1, j, k ) ) ;
        else              (*gx) =             RealNDArray_Item3D ( data, i  , j, k ) - RealNDArray_Item3D ( data, i-1, j, k )   ;
    }
    else                  (*gx) =             RealNDArray_Item3D ( data, i+1, j, k ) - RealNDArray_Item3D ( data, i  , j, k )   ;
    if ( j > 0 )
    {
        if ( j < nj - 1 ) (*gy) = 0.5e+00 * ( RealNDArray_Item3D ( data, i, j+1, k ) - RealNDArray_Item3D ( data, i, j-1, k ) ) ;
        else              (*gy) =             RealNDArray_Item3D ( data, i, j  , k ) - RealNDArray_Item3D ( data, i, j-1, k )   ;
    }
    else                  (*gy) =             RealNDArray_Item3D ( data, i, j+1, k ) - RealNDArray_Item3D ( data, i, j  , k )   ;
    if ( k > 0 )
    {
        if ( k < nk - 1 ) (*gz) = 0.5e+00 * ( RealNDArray_Item3D ( data, i, j, k+1 ) - RealNDArray_Item3D ( data, i, j, k-1 ) ) ;
        else              (*gz) =             RealNDArray_Item3D ( data, i, j, k   ) - RealNDArray_Item3D ( data, i, j, k-1 )   ;
    }
    else                  (*gz) =             RealNDArray_Item3D ( data, i, j, k+1 ) - RealNDArray_Item3D ( data, i, j, k   )   ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Linear interpolation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real LinearlyInterpolate ( const Real f0, const Real f1, const Real u0 )
{
    Real delta, u = u0 ;
    delta = f0 - f1 ;
    if ( fabs ( delta ) > SAFEMINIMUM ) u = f0 / delta ;
    return u ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Print a cube for debugging.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PrintCube ( const Real *cube ) { printf ( "\t%f %f %f %f %f %f %f %f\n", cube[0], cube[1], cube[2], cube[3], cube[4], cube[5], cube[6], cube[7] ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Tesselate a cube.
! . This procedure adds at most one vertex and at most 12 triangles. As no checks are made on storage here, it should be ensured
! . that the vertex and triangle arrays have at least this amount of space available on entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void ProcessCube ( const Integer i, const Integer j, const Integer k, const Real *cube, const Integer8 cubecase, const Integer8 configuration, const IntegerNDArray *intersections,
                                                                         Integer *nvertices, Real2DArray *vertices, Real2DArray *normals, Integer *ntriangles, Integer2DArray *triangles )
{
    Integer8 subconfiguration = 0 ;
    Integer   v12 = -1 ;
    switch ( cubecase )
    {
/* . Done earlier.
    case  0 :
        break ;
*/
    case  1 :
        AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING1[configuration], 1, -1 ) ; break ;

    case  2 :
        AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING2[configuration], 2, -1 ) ; break ;

    case  3 :
        if ( TestFace ( cube, TEST3[configuration] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING3_2[configuration], 4, -1 ) ; /* .  3.2 */
        else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING3_1[configuration], 2, -1 ) ; /* .  3.1 */
        break ;

    case  4 :
        if ( TestInterior ( cube, cubecase, configuration, subconfiguration, TEST4[configuration] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING4_1[configuration], 2, -1 ) ; /* . 4.1.1 */
        else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING4_2[configuration], 6, -1 ) ; /* .  4.1.2 */
        break ;

    case  5 :
        AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING5[configuration], 3, -1 ) ;
        break ;

    case  6 :
        if ( TestFace ( cube, TEST6[configuration][0] ) )
            AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING6_2[configuration], 5, -1 ) ; /* .  6.2 */
        else
        {
            if ( TestInterior ( cube, cubecase, configuration, subconfiguration, TEST6[configuration][1] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING6_1_1[configuration], 3, -1 ) ; /* .  6.1.1 */
            else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING6_1_2[configuration], 7, -1 ) ; /* .  6.1.2 */
        }
        break ;

    case  7 :
        if ( TestFace ( cube, TEST7[configuration][0] ) ) subconfiguration +=  1 ;
        if ( TestFace ( cube, TEST7[configuration][1] ) ) subconfiguration +=  2 ;
        if ( TestFace ( cube, TEST7[configuration][2] ) ) subconfiguration +=  4 ;
        switch ( subconfiguration )
        {
            case 0 :
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING7_1[configuration], 3, -1 ) ; break ;
            case 1 :
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING7_2[configuration][0], 5, -1 ) ; break ;
            case 2 :
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING7_2[configuration][1], 5, -1 ) ; break ;
            case 3 :
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING7_3[configuration][0], 9, v12 ) ; break ;
            case 4 :
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING7_2[configuration][2], 5, -1 ) ; break ;
            case 5 :
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING7_3[configuration][1], 9, v12 ) ; break ;
            case 6 :
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING7_3[configuration][2], 9, v12 ) ; break ;
            case 7 :
                if ( TestInterior ( cube, cubecase, configuration, subconfiguration, TEST7[configuration][3] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING7_4_2[configuration], 9, -1 ) ;
                else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING7_4_1[configuration], 5, -1 ) ;
              break ;
        } ;
        break ;

    case  8 :
        AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING8[configuration], 2, -1 ) ;
        break ;

    case  9 :
        AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING9[configuration], 4, -1 ) ;
        break ;

    case 10 :
        if ( TestFace ( cube, TEST10[configuration][0] ) )
        {
            if ( TestFace ( cube, TEST10[configuration][1] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING10_1_1_[configuration], 4, -1 ) ; /* .  10.1.1 */
            else
            {
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING10_2[configuration], 8, v12 ) ; /* .  10.2 */
            }
        }
        else
        {
            if ( TestFace ( cube, TEST10[configuration][1] ) )
            {
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING10_2_[configuration], 8, v12 ) ; /* .  10.2 */
            }
            else
            {
                if ( TestInterior ( cube, cubecase, configuration, subconfiguration, TEST10[configuration][2] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING10_1_1[configuration], 4, -1 ) ; /* .  10.1.1 */
                else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING10_1_2[configuration], 8, -1 ) ; /* .  10.1.2 */
            }
        }
        break ;

    case 11 :
        AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING11[configuration], 4, -1 ) ;
        break ;

    case 12 :
        if ( TestFace ( cube, TEST12[configuration][0] ) )
        {
            if ( TestFace ( cube, TEST12[configuration][1] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING12_1_1_[configuration], 4, -1 ) ; /* .  12.1.1 */
            else
            {
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING12_2[configuration], 8, v12 ) ; /* .  12.2 */
            }
        }
        else
        {
            if ( TestFace ( cube, TEST12[configuration][1] ) )
            {
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING12_2_[configuration], 8, v12 ) ; /* .  12.2 */
            }
            else
            {
                if ( TestInterior ( cube, cubecase, configuration, subconfiguration, TEST12[configuration][2] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING12_1_1[configuration], 4, -1 ) ; /* .  12.1.1 */
                else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING12_1_2[configuration], 8, -1 ) ; /* .  12.1.2 */
            }
        }
        break ;

    case 13 :
        if ( TestFace ( cube, TEST13[configuration][0] ) ) subconfiguration +=  1 ;
        if ( TestFace ( cube, TEST13[configuration][1] ) ) subconfiguration +=  2 ;
        if ( TestFace ( cube, TEST13[configuration][2] ) ) subconfiguration +=  4 ;
        if ( TestFace ( cube, TEST13[configuration][3] ) ) subconfiguration +=  8 ;
        if ( TestFace ( cube, TEST13[configuration][4] ) ) subconfiguration += 16 ;
        if ( TestFace ( cube, TEST13[configuration][5] ) ) subconfiguration += 32 ;
        switch ( SUBCONFIGURATION13[subconfiguration] )
        {
            case 0 :/* 13.1 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_1[configuration], 4, -1 ) ; break ;

            case 1 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2[configuration][0], 6, -1 ) ; break ;
            case 2 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2[configuration][1], 6, -1 ) ; break ;
            case 3 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2[configuration][2], 6, -1 ) ; break ;
            case 4 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2[configuration][3], 6, -1 ) ; break ;
            case 5 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2[configuration][4], 6, -1 ) ; break ;
            case 6 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2[configuration][5], 6, -1 ) ; break ;

            case 7 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][0], 10, v12 ) ; break ;
            case 8 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][1], 10, v12 ) ; break ;
            case 9 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][2], 10, v12 ) ; break ;
            case 10 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][3], 10, v12 ) ; break ;
            case 11 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][4], 10, v12 ) ; break ;
            case 12 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][5], 10, v12 ) ; break ;
            case 13 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][6], 10, v12 ) ; break ;
            case 14 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][7], 10, v12 ) ; break ;
            case 15 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][8], 10, v12 ) ; break ;
            case 16 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][9], 10, v12 ) ; break ;
            case 17 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][10], 10, v12 ) ; break ;
            case 18 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3[configuration][11], 10, v12 ) ; break ;

            case 19 :/* 13.4 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_4[configuration][0], 12, v12 ) ; break ;
            case 20 :/* 13.4 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_4[configuration][1], 12, v12 ) ; break ;
            case 21 :/* 13.4 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_4[configuration][2], 12, v12 ) ; break ;
            case 22 :/* 13.4 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_4[configuration][3], 12, v12 ) ; break ;

            case 23 :/* 13.5 */
                subconfiguration = 0 ;
                if ( TestInterior ( cube, cubecase, configuration, subconfiguration, TEST13[configuration][6] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_5_1[configuration][0], 6, -1 ) ;
                else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_5_2[configuration][0], 10, -1 ) ;
                break ;
            case 24 :/* 13.5 */
                subconfiguration = 1 ;
                if ( TestInterior ( cube, cubecase, configuration, subconfiguration, TEST13[configuration][6] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_5_1[configuration][1], 6, -1 ) ;
                else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_5_2[configuration][1], 10, -1 ) ;
                break ;
            case 25 :/* 13.5 */
                subconfiguration = 2 ;
                if ( TestInterior ( cube, cubecase, configuration, subconfiguration, TEST13[configuration][6] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_5_1[configuration][2], 6, -1 ) ;
                else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_5_2[configuration][2], 10, -1 ) ;
                break ;
            case 26 :/* 13.5 */
                subconfiguration = 3 ;
                if ( TestInterior ( cube, cubecase, configuration, subconfiguration, TEST13[configuration][6] ) ) AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_5_1[configuration][3], 6, -1 ) ;
                else AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_5_2[configuration][3], 10, -1 ) ;
                break ;

            case 27 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][0], 10, v12 ) ; break ;
            case 28 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][1], 10, v12 ) ; break ;
            case 29 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][2], 10, v12 ) ; break ;
            case 30 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][3], 10, v12 ) ; break ;
            case 31 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][4], 10, v12 ) ; break ;
            case 32 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][5], 10, v12 ) ; break ;
            case 33 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][6], 10, v12 ) ; break ;
            case 34 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][7], 10, v12 ) ; break ;
            case 35 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][8], 10, v12 ) ; break ;
            case 36 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][9], 10, v12 ) ; break ;
            case 37 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][10], 10, v12 ) ; break ;
            case 38 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nvertices, vertices, normals ) ;
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][11], 10, v12 ) ; break ;

            case 39 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][0], 6, -1 ) ; break ;
            case 40 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][1], 6, -1 ) ; break ;
            case 41 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][2], 6, -1 ) ; break ;
            case 42 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][3], 6, -1 ) ; break ;
            case 43 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][4], 6, -1 ) ; break ;
            case 44 :/* 13.2 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][5], 6, -1 ) ; break ;

            case 45 :/* 13.1 */
                AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING13_1_[configuration], 4, -1 ) ; break ;

            default :
                printf ( "Marching Cubes: Impossible case 13?\n" ) ;  PrintCube ( cube ) ;
        }
        break ;
    case 14 :
        AddTriangle ( ntriangles, triangles, intersections, i, j, k, TILING14[configuration], 4, -1 ) ;
        break ;
    } ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Tests if the s of the tesselation of the cube should be connected by the interior of an ambiguous face.
! . Returns True if the face contains a part of the surface.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean TestFace ( const Real *cube, const Integer8 face )
{
    Real A,B,C,D ;
    switch ( face )
    {
        case -1 : case 1 :  A = cube[0] ;  B = cube[4] ;  C = cube[5] ;  D = cube[1] ;  break ;
        case -2 : case 2 :  A = cube[1] ;  B = cube[5] ;  C = cube[6] ;  D = cube[2] ;  break ;
        case -3 : case 3 :  A = cube[2] ;  B = cube[6] ;  C = cube[7] ;  D = cube[3] ;  break ;
        case -4 : case 4 :  A = cube[3] ;  B = cube[7] ;  C = cube[4] ;  D = cube[0] ;  break ;
        case -5 : case 5 :  A = cube[0] ;  B = cube[3] ;  C = cube[2] ;  D = cube[1] ;  break ;
        case -6 : case 6 :  A = cube[4] ;  B = cube[7] ;  C = cube[6] ;  D = cube[5] ;  break ;
        default : printf ( "Invalid face code %d\n", face ) ;  PrintCube ( cube ) ;  A = B = C = D = 0 ;
    } ;
    if ( fabs ( A*C - B*D ) < EPSILON ) return face >= 0 ;
    return face * A * ( A*C - B*D ) >= 0  ;  /* . face and A invert signs. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Tests if the components of the tesselation of the cube should be connected through the interior of the cube.
! . If the interior is empty returns True for s = 7 and False for s = -7.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean TestInterior ( const Real *cube, const Integer8 cubecase, const Integer8 configuration, const Integer8 subconfiguration, const Integer8 s )
{
    Integer8 edge = -1, test = 0 ; /* . edge is the reference edge of the triangulation. */
    Real     At = 0, Bt = 0, Ct = 0, Dt = 0 ;
    Real     a, b, t ;
    switch ( cubecase )
    {
        case  4 :
        case 10 :
            a = ( cube[4] - cube[0] ) * ( cube[6] - cube[2] ) - ( cube[7] - cube[3] ) * ( cube[5] - cube[1] ) ;
            b =  cube[2] * ( cube[4] - cube[0] ) + cube[0] * ( cube[6] - cube[2] ) - cube[1] * ( cube[7] - cube[3] ) - cube[3] * ( cube[5] - cube[1] ) ;
            t = - b / ( 2 * a ) ;
            if ( ( t < 0 ) || ( t > 1 ) ) return s > 0 ;
            At = cube[0] + ( cube[4] - cube[0] ) * t ;
            Bt = cube[3] + ( cube[7] - cube[3] ) * t ;
            Ct = cube[2] + ( cube[6] - cube[2] ) * t ;
            Dt = cube[1] + ( cube[5] - cube[1] ) * t ;
            break ;
        case  6 :
        case  7 :
        case 12 :
        case 13 :
            switch ( cubecase )
            {
                case  6 : edge = TEST6 [configuration][2] ; break ;
                case  7 : edge = TEST7 [configuration][4] ; break ;
                case 12 : edge = TEST12[configuration][3] ; break ;
                case 13 : edge = TILING13_5_1[configuration][subconfiguration][0] ; break ;
            }
            switch ( edge )
            {
            case  0 :
                t  = cube[0] / ( cube[0] - cube[1] ) ;
                At = 0 ;
                Bt = cube[3] + ( cube[2] - cube[3] ) * t ;
                Ct = cube[7] + ( cube[6] - cube[7] ) * t ;
                Dt = cube[4] + ( cube[5] - cube[4] ) * t ;
                break ;
            case  1 :
                t  = cube[1] / ( cube[1] - cube[2] ) ;
                At = 0 ;
                Bt = cube[0] + ( cube[3] - cube[0] ) * t ;
                Ct = cube[4] + ( cube[7] - cube[4] ) * t ;
                Dt = cube[5] + ( cube[6] - cube[5] ) * t ;
                break ;
            case  2 :
                t  = cube[2] / ( cube[2] - cube[3] ) ;
                At = 0 ;
                Bt = cube[1] + ( cube[0] - cube[1] ) * t ;
                Ct = cube[5] + ( cube[4] - cube[5] ) * t ;
                Dt = cube[6] + ( cube[7] - cube[6] ) * t ;
                break ;
            case  3 :
                t  = cube[3] / ( cube[3] - cube[0] ) ;
                At = 0 ;
                Bt = cube[2] + ( cube[1] - cube[2] ) * t ;
                Ct = cube[6] + ( cube[5] - cube[6] ) * t ;
                Dt = cube[7] + ( cube[4] - cube[7] ) * t ;
                break ;
            case  4 :
                t  = cube[4] / ( cube[4] - cube[5] ) ;
                At = 0 ;
                Bt = cube[7] + ( cube[6] - cube[7] ) * t ;
                Ct = cube[3] + ( cube[2] - cube[3] ) * t ;
                Dt = cube[0] + ( cube[1] - cube[0] ) * t ;
                break ;
            case  5 :
                t  = cube[5] / ( cube[5] - cube[6] ) ;
                At = 0 ;
                Bt = cube[4] + ( cube[7] - cube[4] ) * t ;
                Ct = cube[0] + ( cube[3] - cube[0] ) * t ;
                Dt = cube[1] + ( cube[2] - cube[1] ) * t ;
                break ;
            case  6 :
                t  = cube[6] / ( cube[6] - cube[7] ) ;
                At = 0 ;
                Bt = cube[5] + ( cube[4] - cube[5] ) * t ;
                Ct = cube[1] + ( cube[0] - cube[1] ) * t ;
                Dt = cube[2] + ( cube[3] - cube[2] ) * t ;
                break ;
            case  7 :
                t  = cube[7] / ( cube[7] - cube[4] ) ;
                At = 0 ;
                Bt = cube[6] + ( cube[5] - cube[6] ) * t ;
                Ct = cube[2] + ( cube[1] - cube[2] ) * t ;
                Dt = cube[3] + ( cube[0] - cube[3] ) * t ;
                break ;
            case  8 :
                t  = cube[0] / ( cube[0] - cube[4] ) ;
                At = 0 ;
                Bt = cube[3] + ( cube[7] - cube[3] ) * t ;
                Ct = cube[2] + ( cube[6] - cube[2] ) * t ;
                Dt = cube[1] + ( cube[5] - cube[1] ) * t ;
                break ;
            case  9 :
                t  = cube[1] / ( cube[1] - cube[5] ) ;
                At = 0 ;
                Bt = cube[0] + ( cube[4] - cube[0] ) * t ;
                Ct = cube[3] + ( cube[7] - cube[3] ) * t ;
                Dt = cube[2] + ( cube[6] - cube[2] ) * t ;
                break ;
            case 10 :
                t  = cube[2] / ( cube[2] - cube[6] ) ;
                At = 0 ;
                Bt = cube[1] + ( cube[5] - cube[1] ) * t ;
                Ct = cube[0] + ( cube[4] - cube[0] ) * t ;
                Dt = cube[3] + ( cube[7] - cube[3] ) * t ;
                break ;
            case 11 :
                t  = cube[3] / ( cube[3] - cube[7] ) ;
                At = 0 ;
                Bt = cube[2] + ( cube[6] - cube[2] ) * t ;
                Ct = cube[1] + ( cube[5] - cube[1] ) * t ;
                Dt = cube[0] + ( cube[4] - cube[0] ) * t ;
                break ;
            default : printf ( "Invalid edge %d\n", edge ) ;  PrintCube ( cube ) ;  break ;
        }
        break ;
        default : printf ( "Invalid ambiguous case %d\n", cubecase ) ;  PrintCube ( cube ) ;  break ;
    }
    if ( At >= 0 ) test ++   ;
    if ( Bt >= 0 ) test += 2 ;
    if ( Ct >= 0 ) test += 4 ;
    if ( Dt >= 0 ) test += 8 ;
    switch ( test )
    {
        case  0 : return s > 0 ;
        case  1 : return s > 0 ;
        case  2 : return s > 0 ;
        case  3 : return s > 0 ;
        case  4 : return s > 0 ;
        case  5 : if ( ( At * Ct - Bt * Dt ) <  EPSILON ) return s > 0 ; break ;
        case  6 : return s > 0 ;
        case  7 : return s < 0 ;
        case  8 : return s > 0 ;
        case  9 : return s > 0 ;
        case 10 : if ( ( At * Ct - Bt * Dt ) >= EPSILON ) return s > 0 ; break ;
        case 11 : return s < 0 ;
        case 12 : return s > 0 ;
        case 13 : return s < 0 ;
        case 14 : return s < 0 ;
        case 15 : return s < 0 ;
    }
    return ( s < 0 ) ;
}
