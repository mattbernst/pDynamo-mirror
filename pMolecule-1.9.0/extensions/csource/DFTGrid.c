/*------------------------------------------------------------------------------
! . File      : DFTGrid.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module handles the DFT grid.
!
! . This module uses a modified Mura-Knowles method for radial integration and Lebedev grids for angular integration.
!
! . All units are atomic.
!===================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "DFTGrid.h"
# include "Lebedev.h"
# include "Memory.h"
# include "Units.h"

/* . Nothing complicated done here until decide what to do about the grid. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The elemental bragg radii. */
/* . Units are Angstroms and must be converted to Bohrs before use. */
# define NELEMENTS 110
const Real BraggRadii[NELEMENTS] = { 0.75, 0.35, 0.35, 1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50 ,
                                     0.50, 1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.00, 2.20 ,
                                     1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35 ,
                                     1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.15, 2.35, 2.00, 1.80 ,
                                     1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55 ,
                                     1.45, 1.45, 1.40, 1.40, 1.40, 2.60, 2.15, 1.95, 1.85, 1.85 ,
                                     1.85, 1.85, 1.85, 1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75 ,
                                     1.75, 1.75, 1.55, 1.45, 1.35, 1.35, 1.30, 1.35, 1.35, 1.35 ,
                                     1.50, 1.90, 1.80, 1.60, 1.90, 1.90, 1.90, 2.60, 2.15, 1.95 ,
                                     1.80, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75 ,
                                     1.75, 1.75, 1.75, 1.75, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55 } ;

/* . Tolerances for basis function calculation. */
const Real BFTolerances[NDFTGRID_ACCURACY] = { 1.0e-8, 1.0e-9, 1.0e-10, 1.0e-12, 1.0e-15 } ;

/* . Tolerances for density calculation. */
const Real RhoTolerances[NDFTGRID_ACCURACY] = { 1.0e-13, 1.0e-14, 1.0e-15, 1.0e-17, 1.0e-20 } ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void               DFTGrid_Atom_Parameters ( const DFTGrid *self, const Integer ni, Integer *nr, Integer *lvalue, Real *maximumRadius ) ;
static Real               DFTGrid_Bragg_Radius    ( const Integer atomicNumber ) ;
static void               DFTGrid_Radial_Points   ( const Integer nr, const Real range, Real *r, Real *w ) ;

static DFTGridPointBlock *DFTGridPointBlock_Allocate   ( const Integer numberOfPoints, const Integer atom, Coordinates3 **rG, Real1DArray **wG ) ;
static void               DFTGridPointBlock_Deallocate ( void *vSelf ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTGrid *DFTGrid_Allocate ( const DFTGridAccuracy accuracy )
{
    DFTGrid *self = NULL ;
    self = ( DFTGrid * ) Memory_Allocate ( sizeof ( DFTGrid ) ) ;
    self->accuracy          = accuracy ;
    self->blockSize         = 128 ;
    self->numberOfPoints    =   0 ;
    self->numberOfRecords   =  -1 ;
    self->bfTolerance       = BFTolerances [accuracy] ;
    self->rhoTolerance      = RhoTolerances[accuracy] ;
    self->points            = List_Allocate ( ) ;
    self->points->Element_Deallocate = DFTGridPointBlock_Deallocate ;
    self->records           = NULL ;
    self->weights           = NULL ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Construct a grid.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define MINIMUM_LVALUE       9
# define RADIAL_CUTOFF_FACTOR 0.2e+00
# define WEIGHT_TOLERANCE     1.0e-30

DFTGrid *DFTGrid_Construct ( const DFTGridAccuracy accuracy, const QCAtomContainer *qcAtoms, const Coordinates3 *qcCoordinates3 )
{
    DFTGrid *self = NULL ;
    if ( ( qcAtoms != NULL ) && ( qcCoordinates3 != NULL ) && ( qcAtoms->natoms > 0 ) )
    {
        auto Integer            ia, ipt, iqm, ir, lmax = -1, lold, lStart, lStop, lval, n, nangmax, nangmin, napts = 0, nLocal, nr = 0 ;
        auto Real               maximumRadius = 0.0e+00, range, rCutoff, w, wfac, xqm, yqm, zqm ;
        auto Real               pg[3], *rr, *wa, *work1, *work2, *wr, *xa, *ya, *za ;
        auto Coordinates3      *rG, *rLocal, rSlice ;
        auto DFTGridPointBlock *block ;
        auto Real1DArray       *wG, *wLocal, wSlice ;

        /* . Allocate temporary space. */
        work1 = Memory_Allocate_Array_Real ( qcAtoms->natoms ) ;
        work2 = Memory_Allocate_Array_Real ( qcAtoms->natoms ) ;

        /* . Temporarily fill work1 with the Bragg radii. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ ) work1[iqm] = DFTGrid_Bragg_Radius ( qcAtoms->data[iqm].atomicNumber ) ;

        /* . Allocation. */
        self = DFTGrid_Allocate ( accuracy ) ;
        self->weights = DFTGridWeights_Allocate ( qcCoordinates3, work1 ) ;

        /* . Loop over the atoms. */
        for ( iqm = 0 ; iqm < qcAtoms->natoms ; iqm++ )
        {
            /* . Get the data for the atom. */
            Coordinates3_GetRow ( qcCoordinates3, iqm, xqm, yqm, zqm ) ;
            DFTGrid_Atom_Parameters ( self, qcAtoms->data[iqm].atomicNumber, &nr, &lmax, &maximumRadius ) ;
            nangmax = LebedevLaikov_Number_Of_Points ( lmax           ) ;
            nangmin = LebedevLaikov_Number_Of_Points ( MINIMUM_LVALUE ) ;
            range   = DFTGrid_Bragg_Radius ( qcAtoms->data[iqm].atomicNumber ) ;

            /* . Allocate space. */
            rr = Memory_Allocate_Array_Real ( nr      ) ;
            wr = Memory_Allocate_Array_Real ( nr      ) ;
            wa = Memory_Allocate_Array_Real ( nangmax ) ;
            xa = Memory_Allocate_Array_Real ( nangmax ) ;
            ya = Memory_Allocate_Array_Real ( nangmax ) ;
            za = Memory_Allocate_Array_Real ( nangmax ) ;
            rG = Coordinates3_Allocate ( nangmax * nr ) ;
            wG = Real1DArray_Allocate  ( nangmax * nr, NULL ) ;

            /* . Get the radial grid points. */
            rCutoff = RADIAL_CUTOFF_FACTOR * range ;
            DFTGrid_Radial_Points ( nr, range, rr, wr ) ;

            /* . Loop  over the radial grid points. */
            for ( ir = 0, ipt = 0, lold = -1 ; ir < nr ; ir++ )
            {
                /* . Check for the maximum value of r. */
                if ( rr[ir] > maximumRadius ) break ;

                /* . Get the angular points. */
                if ( rr[ir] > rCutoff )
                {
                    lval  = lmax ;
                }
                else
                {
                    n = ( Integer ) ceil ( ( ( Real ) nangmax ) * rr[ir] / rCutoff ) ;
                    if ( n < nangmin ) lval = MINIMUM_LVALUE ;
                    else               lval = LebedevLaikov_Angular_Momentum_Value ( n ) ;
                }

                /* . Get the angular points. */
                if ( lval != lold )
                {
                    napts = LebedevLaikov_Number_Of_Points ( lval ) ;
                    napts = LebedevLaikov_Points ( napts, xa, ya, za, wa ) ;
                    lold  = lval ;
                }

                /* . Construct the integration points. */
                {
                    auto Real sum = 0.0e+00, xdev ;
                    wfac = 4.0e+00 * M_PI * wr[ir] ;
                    for ( ia = 0 ; ia < napts ; ia++ )
                    {
                        /* . Get the point. */
                        pg[0] = xa[ia] * rr[ir] + xqm ;
                        pg[1] = ya[ia] * rr[ir] + yqm ;
                        pg[2] = za[ia] * rr[ir] + zqm ;
                        w = wfac * wa[ia] * DFTGridWeights_Weight ( self->weights, iqm, pg, work1, work2 ) ;
                        if ( fabs ( w ) > WEIGHT_TOLERANCE )
                        {
                            Coordinates3_SetRow ( rG, ipt, pg[0], pg[1], pg[2] ) ;
                            Real1DArray_Item ( wG, ipt ) = w ;
                            ipt++ ;
                        }
                        xdev = fabs ( 1.0e+00 - xa[ia]*xa[ia] - ya[ia]*ya[ia] - za[ia]*za[ia] ) ;
                        if ( xdev > 1.0e-8 ) printf ( "Node Inaccuracy = %5d %5d %5d %25.15f\n", ir, napts, lval, xdev ) ;
                        sum += wa[ia] ;
                    }
                    xdev = fabs ( sum - 1.0e+00 ) ;
                    if ( xdev > 1.0e-9 ) printf ( "Weight Inaccuracy = %5d %5d %5d %25.15f\n", ir, napts, lval, xdev ) ;
                }
            }

            /* . Save the grid points in blocks of the appropriate size. */
            lStart = 0 ;
            do
            {
                lStop  = Minimum ( lStart + self->blockSize, ipt ) ;
                nLocal = lStop - lStart ;
                self->numberOfPoints += nLocal ;
                rLocal = Coordinates3_Allocate ( nLocal ) ;
                wLocal = Real1DArray_Allocate  ( nLocal, NULL ) ;
                Coordinates3_2DSlice ( rG, lStart, lStop, 1, 0, 3, 1, &rSlice, NULL ) ;
                Coordinates3_CopyTo  ( &rSlice, rLocal, NULL ) ;
                Real1DArray_Slice    ( wG, lStart, lStop, 1, &wSlice, NULL ) ;
                Real1DArray_CopyTo   ( &wSlice, wLocal, NULL ) ;
                block  = DFTGridPointBlock_Allocate ( nLocal, iqm, &rLocal, &wLocal ) ;
                List_Element_Append ( self->points, ( void * ) block ) ;
                lStart += self->blockSize ;
            }
            while ( lStart < ipt ) ;

            /* . Deallocate all space. */
            Coordinates3_Deallocate ( &rG ) ;
            Memory_Deallocate_Real  ( &rr ) ;
            Memory_Deallocate_Real  ( &wr ) ;
            Memory_Deallocate_Real  ( &wa ) ;
            Memory_Deallocate_Real  ( &xa ) ;
            Memory_Deallocate_Real  ( &ya ) ;
            Memory_Deallocate_Real  ( &za ) ;
            Real1DArray_Deallocate  ( &wG ) ;
        }

        /* . Deallocate space. */
        Memory_Deallocate_Real ( &work1 ) ;
        Memory_Deallocate_Real ( &work2 ) ;
    }
    return self ;
}

# undef MINIMUM_LVALUE
# undef RADIAL_CUTOFF_FACTOR
# undef WEIGHT_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTGrid_Deallocate ( DFTGrid **self )
{
    if ( (*self) != NULL )
    {
        DFTGrid_DeallocateFunctionData ( (*self) ) ;
        DFTGridWeights_Deallocate ( &((*self)->weights) ) ;
        List_Deallocate           ( &((*self)->points)  ) ;
        Memory_Deallocate ( (*self)->records ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation of function data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTGrid_DeallocateFunctionData ( DFTGrid *self )
{
    if ( ( self != NULL ) && DFTGrid_HasFunctionData ( self ) )
    {
        auto Integer r ;
        for ( r = 0 ; r < self->numberOfRecords ; r++ ) GridFunctionDataBlock_Deallocate ( &(self->records[r]->functionData) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Size of the block storage in (decimal) GB.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define TOGIGABYTES 1.0e+09 /* . Assumes sizeof returns bytes. */
Real DFTGrid_FunctionDataSize ( DFTGrid *self )
{
    Real size = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->numberOfPoints > 0 ) && DFTGrid_HasFunctionData ( self ) )
    {
        auto Integer r ;
        for ( r = 0 ; r < self->numberOfRecords ; r++ ) size += GridFunctionDataBlock_ByteSize ( self->records[r]->functionData ) ;
        size /= TOGIGABYTES ;
    }
    return size ;
}
# undef TOGIGABYTES

/*----------------------------------------------------------------------------------------------------------------------------------
! . Has the grid function data?
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean DFTGrid_HasFunctionData ( DFTGrid *self )
{
    Boolean hasData = False ;
    if ( ( self != NULL ) && ( self->numberOfPoints > 0 ) )
    {
        DFTGrid_MakeRecords ( self ) ;
        if ( self->numberOfRecords > 0 ) hasData = ( self->records[0]->functionData != NULL ) ;
    }
    return hasData ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Iteration.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTGridPointBlock *DFTGrid_Iterate ( DFTGrid *self )
{
   if ( self == NULL ) return NULL ;
   else                return ( DFTGridPointBlock * ) List_Iterate ( self->points ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the records representation of the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTGrid_MakeRecords ( DFTGrid *self )
{
    if ( ( self != NULL ) && ( self->records == NULL ) )
    {
        auto Integer n ;
        n = DFTGrid_NumberOfRecords ( self ) ;
        MEMORY_ALLOCATEARRAY ( self->records, n, DFTGridPointBlock * ) ;
        if ( self->records != NULL )
        {
            auto DFTGridPointBlock *record ;
            n = 0 ;
            List_Iterate_Initialize ( self->points ) ;
            while ( ( record = DFTGrid_Iterate ( self ) ) != NULL ) { self->records[n] = record ; n += 1 ; }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the number of stored function values.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer DFTGrid_NumberOfFunctionValues ( DFTGrid *self )
{
    Integer n = 0 ;
    if ( ( self != NULL ) && ( self->numberOfPoints > 0 ) )
    {
        auto Integer r ;
        DFTGrid_MakeRecords ( self ) ;
        for ( r = 0 ; r < self->numberOfRecords ; r++ )
        {
            if ( self->records[r]->functionData != NULL ) n += ( self->records[r]->functionData->numberOfFunctions * self->records[r]->numberOfPoints ) ;
        }
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the number of records in the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer DFTGrid_NumberOfRecords ( DFTGrid *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        if ( self->numberOfRecords >= 0 ) n = self->numberOfRecords ;
        else
        {
            auto DFTGridPointBlock *record ;
            List_Iterate_Initialize ( self->points ) ;
            while ( ( record = DFTGrid_Iterate ( self ) ) != NULL ) n += 1 ;
            self->numberOfRecords = n ;
        }
    }
    return n ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the grid parameters for an element.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTGrid_Atom_Parameters ( const DFTGrid *self, const Integer ni, Integer *nr, Integer *lvalue, Real *maximumRadius )
{
    switch ( self->accuracy )
    {
        case DFTGridAccuracy_VeryLow:
                 if ( ni < 11 ) (*nr) = 21 ;
            else if ( ni < 19 ) (*nr) = 42 ;
            else if ( ni < 37 ) (*nr) = 75 ;
            else                (*nr) = 84 ;
            (*maximumRadius) = 20.0e+00 ;
            (*lvalue)        = 23 ;
            break ;
        case DFTGridAccuracy_Low:
                 if ( ni < 11 ) (*nr) =  35 ;
            else if ( ni < 19 ) (*nr) =  70 ;
            else if ( ni < 37 ) (*nr) =  95 ;
            else                (*nr) = 104 ;
            (*maximumRadius) = 25.0e+00 ;
            (*lvalue)        = 35 ;
            break ;
        case DFTGridAccuracy_Medium:
                 if ( ni < 11 ) (*nr) =  49 ;
            else if ( ni < 19 ) (*nr) =  88 ;
            else if ( ni < 37 ) (*nr) = 112 ;
            else                (*nr) = 123 ;
            (*maximumRadius) = 30.0e+00 ;
            if ( ni < 11 ) (*lvalue) = 35 ;
            else           (*lvalue) = 41 ;
            break ;
        case DFTGridAccuracy_High:
                 if ( ni < 11 ) (*nr) =  70 ;
            else if ( ni < 19 ) (*nr) = 123 ;
            else if ( ni < 37 ) (*nr) = 130 ;
            else                (*nr) = 155 ;
            (*maximumRadius) = 35.0e+00 ;
                 if ( ni < 11 ) (*lvalue) = 41 ;
            else if ( ni < 19 ) (*lvalue) = 47 ;
            else if ( ni < 89 ) (*lvalue) = 53 ;
            else                (*lvalue) = 59 ;
            break ;
        case DFTGridAccuracy_VeryHigh:
                 if ( ni < 11 ) (*nr) = 100 ;
            else if ( ni < 19 ) (*nr) = 125 ;
            else if ( ni < 37 ) (*nr) = 160 ;
            else                (*nr) = 205 ;
            (*maximumRadius) = 35.0e+00 ;
            (*lvalue)        = 65 ;
            break ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the bragg radius for an atom.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real DFTGrid_Bragg_Radius ( const Integer atomicNumber )
{
    Integer ni ;
    ni = atomicNumber ;
    if      ( ni <  0         ) ni = 0 ;
    else if ( ni >= NELEMENTS ) ni = NELEMENTS - 1 ;
    return UNITS_LENGTH_ANGSTROMS_TO_BOHRS * BraggRadii[ni] ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the radial integration points.
! . Follows the Mura and Knowles scheme (JCP 104, 9848, 1996).
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifdef DFTGRID_MURAKNOWLES
# define EXPONENT 3.0e+00
# define SCALING  3.3e+00
static void DFTGrid_Radial_Points ( const Integer nr, const Real range, Real *r, Real *w )
{
    Integer i ;
    Real    alpha, fmn, qi, ri, wi ;
    alpha = SCALING * range ;
    fmn   = EXPONENT / ( 1.0e+00 + ( Real ) nr ) ;
    for ( i = 0 ; i < nr ; i++ )
    {
       qi = ( Real ) i / ( ( Real ) nr + 1.0e+00 ) ;
       ri = - alpha * log ( 1.0e+00 - pow ( qi, EXPONENT ) ) ;
       wi = fmn * alpha * ( ri * ri ) / ( 1.0e+00 - pow ( qi, EXPONENT ) ) * pow ( qi, EXPONENT - 1.0e+00 ) ;
       r[i] = ri ;
       w[i] = wi ;
    }
}
# undef EXPONENT
# undef SCALING
# else
static void DFTGrid_Radial_Points ( const Integer nr, const Real range, Real *r, Real *w )
{
    Integer i ;
    Real    dfac, ri, wi, xnode ;
    dfac = ( 1.0e+00 + ( Real ) nr ) ;
    for ( i = 0 ; i < nr ; i++ )
    {
       xnode = ( Real ) ( i + 1 ) / dfac ;
       ri    = range * pow ( ( xnode / ( 1.0e+00 - xnode ) ), 2.0e+00 ) ;
       wi    = 2.0e+00 * pow ( range, 3.0e+00 ) * pow ( xnode, 5.0e+00 ) / pow ( ( 1.0e+00 - xnode ), 7.0e+00 ) ;
       r[i] = ri ;
       w[i] = wi / dfac ;
    }
}
# endif

/*==================================================================================================================================
! . Private grid point block procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate a block.
!---------------------------------------------------------------------------------------------------------------------------------*/
static DFTGridPointBlock *DFTGridPointBlock_Allocate ( const Integer numberOfPoints, const Integer atom, Coordinates3 **rG, Real1DArray **wG )
{
   DFTGridPointBlock *self = NULL ;
   if ( numberOfPoints > 0 )
   {
       self = ( DFTGridPointBlock * ) Memory_Allocate ( sizeof ( DFTGridPointBlock ) ) ;
       self->atom           = atom  ;
       self->coordinates3   = (*rG) ;
       self->functionData   = NULL  ;
       self->numberOfPoints = numberOfPoints ;
       self->weights        = (*wG) ;
       /* . Take ownership of the arrays. */
       (*rG) = NULL ;
       (*wG) = NULL ;
   }
   return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocate a block.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTGridPointBlock_Deallocate ( void *vSelf )
{
   DFTGridPointBlock *self ;
   self = ( DFTGridPointBlock * ) vSelf ;
   Coordinates3_Deallocate          ( &(self->coordinates3) ) ;
   GridFunctionDataBlock_Deallocate ( &(self->functionData) ) ;
   Real1DArray_Deallocate           ( &(self->weights     ) ) ;
   free ( self ) ;
}
