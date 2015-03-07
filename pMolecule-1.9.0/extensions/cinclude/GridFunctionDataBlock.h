/*------------------------------------------------------------------------------
! . File      : GridFunctionDataBlock.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _GRIDFUNCTIONDATABLOCK
# define _GRIDFUNCTIONDATABLOCK

# include "Boolean.h"
# include "Integer.h"
# include "Integer1DArray.h"
# include "Real2DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The grid function data block type. */
typedef struct {
    Integer numberOfFunctions ;
    Integer numberOfPoints    ;
    Integer order             ;
    Integer1DArray *indices   ;
    Real2DArray    *f         ;
    Real2DArray    *fX        ;
    Real2DArray    *fY        ;
    Real2DArray    *fZ        ;
    Real2DArray    *fXX       ;
    Real2DArray    *fXY       ;
    Real2DArray    *fXZ       ;
    Real2DArray    *fYY       ;
    Real2DArray    *fYZ       ;
    Real2DArray    *fZZ       ;
    Real2DArray    *fXXX      ;
    Real2DArray    *fXXY      ;
    Real2DArray    *fXXZ      ;
    Real2DArray    *fXYY      ;
    Real2DArray    *fXYZ      ;
    Real2DArray    *fXZZ      ;
    Real2DArray    *fYYY      ;
    Real2DArray    *fYYZ      ;
    Real2DArray    *fYZZ      ;
    Real2DArray    *fZZZ      ;
} GridFunctionDataBlock ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern GridFunctionDataBlock *GridFunctionDataBlock_Allocate     ( const Integer numberOfFunctions, const Integer numberOfPoints, const Integer order, Status *status ) ;
extern Real                   GridFunctionDataBlock_ByteSize     ( const GridFunctionDataBlock  *self ) ;
extern void                   GridFunctionDataBlock_Deallocate   (       GridFunctionDataBlock **self ) ;
extern void                   GridFunctionDataBlock_FilterValues (       GridFunctionDataBlock  *self, const Integer fStart, const Real *tolerance ) ;
extern void                   GridFunctionDataBlock_Initialize   (       GridFunctionDataBlock  *self ) ;
extern void                   GridFunctionDataBlock_Resize       (       GridFunctionDataBlock  *self, const Integer numberOfFunctions, Status *status ) ;

# endif
