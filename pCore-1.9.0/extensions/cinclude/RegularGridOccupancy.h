/*------------------------------------------------------------------------------
! . File      : RegularGridOccupancy.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _REGULARGRIDOCCUPANCY
# define _REGULARGRIDOCCUPANCY

# include "Integer.h"
# include "Integer1DArray.h"
# include "Real2DArray.h"
# include "RegularGrid.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The occupancy type. */
typedef struct {
    Integer         numberOfCells   ;
    Integer         numberOfPoints  ;
    Integer1DArray *cellFirstPoints ; /* . The indices of the first points in each cell in the cellPoints array. */
    Integer1DArray *cellPoints      ; /* . The indices of the points within each cell. */
    Integer1DArray *cellTotalPoints ; /* . The total number of points per cell. */
    Integer1DArray *pointCells      ; /* . The grid cell within which a point is found. */
} RegularGridOccupancy ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern RegularGridOccupancy *RegularGridOccupancy_Allocate          ( const Integer numberOfCells, const Integer numberOfPoints, Status *status ) ;
extern void                  RegularGridOccupancy_Deallocate        (       RegularGridOccupancy **self ) ;
extern Integer               RegularGridOccupancy_Fill              (       RegularGridOccupancy  *self, const RegularGrid *grid, const Real2DArray *points, Status *status ) ;
extern RegularGridOccupancy *RegularGridOccupancy_FromGridAndPoints ( const RegularGrid *grid, const Real2DArray *points, Status *status ) ;
extern void                  RegularGridOccupancy_Initialize        (       RegularGridOccupancy  *self ) ;

# endif
