/*------------------------------------------------------------------------------
! . File      : QCAtomContainer.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "Definitions.h"
# include "GaussianBasisCore.h"
# include "GaussianBasisGrid.h"
# include "Memory.h"
# include "QCAtomContainer.h"
# include "Units.h"

/*==============================================================================
! . Parameters.
!=============================================================================*/
/* . Nuclear charge distribution parameters. */
# define NUCLEAR_WIDTH  1.0e-4
# define NUCLEAR_WIDTHE ( 4.0e+00 / ( NUCLEAR_WIDTH * NUCLEAR_WIDTH * M_PI ) )
# define NUCLEAR_WIDTHN sqrt ( pow ( ( NUCLEAR_WIDTHE / M_PI ), 3 ) )

/*==============================================================================
! . Standard procedures.
!=============================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCAtomContainer *QCAtomContainer_Allocate ( const Integer natoms )
{
    QCAtomContainer *self = NULL ;
    if ( natoms > 0 )
    {
        auto Integer i ;
        self                 = ( QCAtomContainer * ) Memory_Allocate ( sizeof ( QCAtomContainer ) ) ;
        self->data           = ( QCAtom * ) Memory_Allocate_Array ( natoms, sizeof ( QCAtom ) ) ;
        self->natoms         = natoms ;
        self->nboundary      = 0 ;
        self->ndbasis        = 0 ;
        self->nfbasis        = 0 ;
        self->nobasis        = 0 ;
        self->npbasis        = 0 ;
        self->ndbasisw       = 0 ;
        self->nfbasisw       = 0 ;
        self->nobasisw       = 0 ;
        self->npbasisw       = 0 ;
        self->nuclearCharge  = 0 ;
        self->energybaseline = 0.00e+00 ;
        self->QLINKRATIO     = False    ;
        self->QTOSPHERICAL   = False    ;
        for ( i = 0 ; i < natoms ; i++ )
        {
            self->data[i].QBOUNDARY    = False ;
            self->data[i].atomicNumber = -1 ;
            self->data[i].center       = -1 ;
            self->data[i].index        = -1 ;
            self->data[i].qcpartner    = -1 ;
            self->data[i].dstart       =  0 ;
            self->data[i].fstart       =  0 ;
            self->data[i].ostart       =  0 ;
            self->data[i].pstart       =  0 ;
            self->data[i].ndbasis      =  0 ;
            self->data[i].nfbasis      =  0 ;
            self->data[i].nobasis      =  0 ;
            self->data[i].npbasis      =  0 ;
            self->data[i].dstartw      =  0 ;
            self->data[i].fstartw      =  0 ;
            self->data[i].ostartw      =  0 ;
            self->data[i].pstartw      =  0 ;
            self->data[i].ndbasisw     =  0 ;
            self->data[i].nfbasisw     =  0 ;
            self->data[i].nobasisw     =  0 ;
            self->data[i].npbasisw     =  0 ;
            self->data[i].linkfactor   = 0.0e+00 ;
            self->data[i].widthe       = NUCLEAR_WIDTHE ;
            self->data[i].widthn       = NUCLEAR_WIDTHN ;
            self->data[i].nmmpartners  =  0 ;
            self->data[i].mmpartners   = NULL ;
        }
        /* . Representations. */
        self->baselection  = NULL ;
        self->mmpselection = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCAtomContainer *QCAtomContainer_Clone ( const QCAtomContainer *self )
{
    QCAtomContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i, m ;
        new = QCAtomContainer_Allocate ( self->natoms ) ;
        new->natoms         = self->natoms          ;
        new->nboundary      = self->nboundary       ;
        new->ndbasis        = self->ndbasis         ;
        new->nfbasis        = self->nfbasis         ;
        new->nobasis        = self->nobasis         ;
        new->npbasis        = self->npbasis         ;
        new->ndbasisw       = self->ndbasisw        ;
        new->nfbasisw       = self->nfbasisw        ;
        new->nobasisw       = self->nobasisw        ;
        new->npbasisw       = self->npbasisw        ;
        new->nuclearCharge  = self->nuclearCharge ;
        new->energybaseline = self->energybaseline  ;
        new->QLINKRATIO     = self->QLINKRATIO      ;
        new->QTOSPHERICAL   = self->QTOSPHERICAL    ;
        for ( i = 0 ; i < self->natoms ; i++ )
        {
            new->data[i].QBOUNDARY    = self->data[i].QBOUNDARY    ;
            new->data[i].atomicNumber = self->data[i].atomicNumber ;
            new->data[i].center       = self->data[i].center       ;
            new->data[i].index        = self->data[i].index        ;
            new->data[i].qcpartner    = self->data[i].qcpartner    ;
            new->data[i].dstart       = self->data[i].dstart       ;
            new->data[i].fstart       = self->data[i].fstart       ;
            new->data[i].ostart       = self->data[i].ostart       ;
            new->data[i].pstart       = self->data[i].pstart       ;
            new->data[i].ndbasis      = self->data[i].ndbasis      ;
            new->data[i].nfbasis      = self->data[i].nfbasis      ;
            new->data[i].nobasis      = self->data[i].nobasis      ;
            new->data[i].npbasis      = self->data[i].npbasis      ;
            new->data[i].dstartw      = self->data[i].dstartw      ;
            new->data[i].fstartw      = self->data[i].fstartw      ;
            new->data[i].ostartw      = self->data[i].ostartw      ;
            new->data[i].pstartw      = self->data[i].pstartw      ;
            new->data[i].ndbasisw     = self->data[i].ndbasisw     ;
            new->data[i].nfbasisw     = self->data[i].nfbasisw     ;
            new->data[i].nobasisw     = self->data[i].nobasisw     ;
            new->data[i].npbasisw     = self->data[i].npbasisw     ;
            new->data[i].linkfactor   = self->data[i].linkfactor   ;
            new->data[i].widthe       = self->data[i].widthe       ;
            new->data[i].widthn       = self->data[i].widthn       ;
            new->data[i].nmmpartners  = self->data[i].nmmpartners  ;
            if ( self->data[i].nmmpartners > 0 )
            {
                new->data[i].mmpartners = Memory_Allocate_Array_Integer ( self->data[i].nmmpartners ) ;
                for ( m = 0 ; m < self->data[i].nmmpartners ; m++ ) new->data[i].mmpartners[m] = self->data[i].mmpartners[m] ;
            }
        }
        /* . Representations. */
        new->baselection  = Selection_Clone ( self->baselection  ) ;
        new->mmpselection = Selection_Clone ( self->mmpselection ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCAtomContainer_Deallocate ( QCAtomContainer **self )
{
    if ( (*self) != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < (*self)->natoms ; i++ ) Memory_Deallocate_Integer ( &((*self)->data[i].mmpartners) ) ;
        Memory_Deallocate ( (*self)->data ) ;
        Selection_Deallocate ( &((*self)->baselection)  ) ;
        Selection_Deallocate ( &((*self)->mmpselection) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the boundary atom representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCAtomContainer_MakeBoundaryAtomSelection ( QCAtomContainer *self )
{
    if ( ( self != NULL ) && ( self->baselection == NULL ) )
    {
        auto Integer i, n ;
        self->baselection = Selection_Allocate ( self->nboundary ) ;
        for ( i = n = 0 ; i < self->natoms ; i++ )
        {
            if ( self->data[i].QBOUNDARY )
            {
                self->baselection->indices[n] = self->data[i].index ;
                n++ ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a full selection (including boundary atoms).
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *QCAtomContainer_MakeFullSelection ( const QCAtomContainer *self )
{
    Selection *selection = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        selection = Selection_Allocate ( self->natoms ) ;
        for ( i = 0 ; i < self->natoms ; i++ ) selection->indices[i] = self->data[i].index ;
    }
    return selection ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert a pairlist into one good for MM atom indices.
! . Always return a valid pairlist even if there are no interactions.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *QCAtomContainer_MakeMMPairList ( QCAtomContainer  *self, const MMAtomContainer *mmAtoms, PairList *pairlist )
{
    PairList *new = NULL ;
    if ( ( self != NULL ) && ( pairlist != NULL ) )
    {
        auto Boolean QEMPTY = True ;

        /* . Check for an empty list. */
        if ( PairList_Length ( pairlist ) > 0 )
        {
            auto Integer  i, j, n, nnew, nold, p, s, t ;
            auto Integer  *inew = NULL, *iold, *mapping, *partners ;

            /* . Convert the old list to an array of indices. */
            iold = PairList_ToIntegerPairArray ( pairlist ) ;
            nold = pairlist->npairs ;

            /* . Allocate space. */
            mapping  = Memory_Allocate_Array_Integer_Initialize ( mmAtoms->natoms, -1 ) ;
            partners = Memory_Allocate_Array_Integer_Initialize ( mmAtoms->natoms,  0 ) ;

            if ( ( iold != NULL ) && ( mapping != NULL ) && ( partners != NULL ) )
            {
                /* . Make the pure MM selection. */
                QCAtomContainer_MakePureMMSelection ( self, mmAtoms->natoms ) ;

                /* . Create the mapping - pure MM atoms. */
                for ( n = s = 0 ; s < self->mmpselection->nindices ; n++, s++ )
                {
                    mapping[self->mmpselection->indices[s]]  = n ;
                    partners[self->mmpselection->indices[s]] = 1 ;
                }
                /* . Boundary atoms - without partners. */
                for ( i = 0 ; i < self->natoms ; i++ )
                {
                    if ( self->data[i].QBOUNDARY )
                    {
                        mapping[self->data[i].index]  = n ;
                        partners[self->data[i].index] = self->data[i].nmmpartners ;
                        n += self->data[i].nmmpartners ;
                    }
                }

                /* . Find the number of new pairs. */
                for ( n = nnew = p = 0 ; p < nold ; n += 2, p++ ) nnew += partners[iold[n]] * partners[iold[n+1]] ;

                /* . There are pairs. */
                if ( nnew > 0 )
                {
                    /* . Allocate space. */
                    inew = Memory_Allocate_Array_Integer_Initialize ( 2 * nnew, -1 ) ;

                    /* . Fill the array. */
                    for ( n = nnew = p = 0 ; p < nold ; n += 2, p++ )
                    {
                        i = iold[n]   ;
                        j = iold[n+1] ;
                        if ( ( mapping[i] >= 0 ) && ( mapping[j] >= 0 ) )
                        {
                            for ( s = 0 ; s < partners[i] ; s++ )
                            {
                                for ( t = 0 ; t < partners[j] ; t++ )
                                {
                                    inew[2*nnew]   = mapping[i] + s;
                                    inew[2*nnew+1] = mapping[j] + t;
                                    nnew++ ;
                                }
                            }
                        }
                    }

                    /* . Make the new list.*/
                    new    = PairList_FromIntegerPairArray ( True, nnew, inew ) ;
                    QEMPTY = False ;
                }
            }

            /* . Deallocate space. */
            Memory_Deallocate_Integer ( &inew     ) ;
            Memory_Deallocate_Integer ( &iold     ) ;
            Memory_Deallocate_Integer ( &mapping  ) ;
            Memory_Deallocate_Integer ( &partners ) ;

        }
        /* . An empty list. */
        if ( QEMPTY ) new = PairList_Allocate ( True ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the pure MM representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCAtomContainer_MakePureMMSelection ( QCAtomContainer *self, const Integer upperBound )
{
    if ( self != NULL )
    {
        /* . Existing representation. */
        if ( self->mmpselection != NULL )
        {
            if ( Selection_UpperBound ( self->mmpselection ) < upperBound ) Selection_Deallocate ( &(self->mmpselection) ) ;
        }
        /* . Create the representation. */
        if ( self->mmpselection == NULL )
        {
            auto Integer n ;
            auto Selection *qcselection = NULL ;
            n = Maximum ( upperBound, self->data[self->natoms-1].index + 1 ) ;
            qcselection        = QCAtomContainer_MakeFullSelection ( self ) ;
            self->mmpselection = Selection_Complement ( qcselection, n ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a pure selection (excluding boundary atoms).
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *QCAtomContainer_MakePureSelection ( const QCAtomContainer *self )
{
    Selection *selection = NULL ;
    if ( self != NULL )
    {
        auto Integer i, n ;
        selection = Selection_Allocate ( self->natoms - self->nboundary ) ;
        for ( i = n = 0 ; i < self->natoms ; i++ )
        {
            if ( ! self->data[i].QBOUNDARY )
            {
                selection->indices[n] = self->data[i].index ;
                n++ ;
            }
        }
    }
    return selection ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return an array of atom indices for the orbital basis functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer1DArray *QCAtomContainer_OrbitalBasisAtomIndices ( const QCAtomContainer *self, Status *status )
{
    Integer1DArray *indices = NULL ;
    if ( self != NULL )
    {
        indices = Integer1DArray_Allocate ( self->nobasisw, status ) ;
        if ( indices != NULL )
        {
            auto Integer a, b ;
            for ( a = 0 ; a < self->natoms ; a++ )
            {
                for ( b = 0 ; b < self->data[a].nobasisw ; b++ ) Integer1DArray_Item ( indices, b+self->data[a].ostartw ) = a ;
            }
        }
    }
    return indices ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate orbital basis function values and their derivatives for all QC atoms at a set of grid points.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCAtomContainer_OrbitalBasisGridPointValues ( const QCAtomContainer       *self             ,
                                                   const QCParameter           *qcParameters     ,
                                                   const Coordinates3          *qcCoordinates3   ,
                                                   const Coordinates3          *gridCoordinates3 ,
                                                   const Boolean                resize           ,
                                                   const Real                  *tolerance        ,
                                                         GridFunctionDataBlock *data             ,
                                                         Status                *status           )
{ 
    if ( ( self != NULL ) && ( qcParameters != NULL ) && ( qcCoordinates3 != NULL ) && ( data != NULL ) )
    {
        auto Integer a, c  ;
        GridFunctionDataBlock_Initialize ( data ) ;
        for ( a = 0 ; a < self->natoms ; a++ )
        {
            c  = self->data[a].center ;
            GaussianBasis_GridPointValues ( qcParameters->centers[c].orbitalbasis, Coordinates3_RowPointer ( qcCoordinates3, a ), gridCoordinates3, data ) ;
        }
        if ( ( resize ) && ( tolerance != NULL ) && ( (*tolerance) > 0.0e+00 ) )
        {
            GridFunctionDataBlock_FilterValues ( data, 0, tolerance ) ;
            GridFunctionDataBlock_Resize ( data, data->numberOfFunctions, status ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Size.
!---------------------------------------------------------------------------------------------------------------------------------*/
int QCAtomContainer_Size ( const QCAtomContainer *self )
{
    Integer size = 0 ;
    if ( self != NULL ) size = self->natoms ;
    return size ;
}

/*==============================================================================
! . Non-standard procedures.
!=============================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the coordinates of a QC atom without unit conversion.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCAtomContainer_GetAtomCoordinates3 ( const QCAtomContainer *self, const Integer q, const Coordinates3 *coordinates3, Real *x, Real *y, Real *z )
{
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
        if ( self->data[q].QBOUNDARY )
        {
            auto Real d, xij, yij, zij ;
            Coordinates3_GetRow ( coordinates3, self->data[q].qcpartner, (*x), (*y), (*z) ) ;
            Coordinates3_DifferenceRow ( coordinates3, self->data[q].index, self->data[q].qcpartner, xij, yij, zij ) ;
            d = self->data[q].linkfactor ;
            if ( ! self->QLINKRATIO ) d /= sqrt ( xij * xij + yij * yij + zij * zij ) ;
            (*x) += d * xij ;
            (*y) += d * yij ;
            (*z) += d * zij ;
        }
        else Coordinates3_GetRow ( coordinates3, self->data[q].index, (*x), (*y), (*z) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert atom coordinates into QC atom coordinates with optional unit
! . conversion.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status QCAtomContainer_GetCoordinates3 ( const QCAtomContainer *self, const Coordinates3 *coordinates3, const Boolean QBOHR, Coordinates3 **qccoordinates3 )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
        auto Real x, y, z ;
        auto Integer    iqm ;
        /* . Allocate space if necessary. */
        if ( (*qccoordinates3) == NULL ) (* qccoordinates3) = Coordinates3_Allocate ( self->natoms ) ;
        if ( (*qccoordinates3) == NULL ) status = Status_OutOfMemory ;
        else
        {
            for ( iqm = 0 ; iqm < self->natoms ; iqm++ )
            {
                QCAtomContainer_GetAtomCoordinates3 ( self, iqm, coordinates3, &x, &y, &z ) ;
                Coordinates3_SetRow ( (*qccoordinates3), iqm, x, y, z ) ;
            }
            if ( QBOHR ) Coordinates3_Scale ( (*qccoordinates3), UNITS_LENGTH_ANGSTROMS_TO_BOHRS ) ;
            status = Status_Success ;
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the MM charges for the QC/MM interactions (RC/RD coupling only).
!---------------------------------------------------------------------------------------------------------------------------------*/
Real1DArray *QCAtomContainer_GetMMCharges ( QCAtomContainer *self, const MMAtomContainer *mmAtoms, const Boolean QRCCOUPLING )
{
    Real1DArray *mmCharges = NULL ;
    if ( ( self != NULL ) && ( mmAtoms != NULL ) )
    {
        auto Real q0 ;
        auto Integer    b, i, n, np, p, s, t ;
        /* . Make the boundary atom and pure MM selections. */
        QCAtomContainer_MakeBoundaryAtomSelection ( self                  ) ;
        QCAtomContainer_MakePureMMSelection       ( self, mmAtoms->natoms ) ;
        /* . Find the number of centers. */
        /* . Pure MM atoms. */
        n = self->mmpselection->nindices ;
        /* . Boundary MM atoms and their partners. */
        for ( i = 0 ; i < self->natoms ; i++ )
        {
            if ( self->data[i].QBOUNDARY ) n += self->data[i].nmmpartners ;
        }
        /* . Allocate space. */
        mmCharges = Real1DArray_Allocate ( n, NULL ) ;
        if ( mmCharges != NULL )
        {
            /* . Make the correct representations of the selections. */
            Selection_MakeFlags     ( self->baselection,  mmAtoms->natoms ) ;
            Selection_MakePositions ( self->mmpselection, mmAtoms->natoms ) ;
            /* . Initialize charges. */
            Real1DArray_Set ( mmCharges, 0.0e+00 ) ;
            /* . Pure MM atoms. */
            for ( n = s = 0 ; s < self->mmpselection->nindices ; n++, s++ )
            {
                i = self->mmpselection->indices[s] ;
                Real1DArray_Item ( mmCharges, n ) = mmAtoms->data[i].charge ;
            }
            /* . Boundary MM atoms and their partners. */
            for ( i = 0 ; i < self->natoms ; i++ )
            {
                if ( self->data[i].QBOUNDARY )
                {
                    b = self->data[i].index ;
                    /* . Loop over the MM partners - np should never be zero. */
                    np = self->data[i].nmmpartners ;
                    q0 = mmAtoms->data[b].charge / ( Real ) np ;
                    for ( t = 0 ; t < self->data[i].nmmpartners ; n++, t++ )
                    {
                        p = self->data[i].mmpartners[t] ;
                        if ( self->baselection->flags[p] || QRCCOUPLING )
                        {
                            Real1DArray_Item ( mmCharges, n ) = q0 ;
                        }
                        else
                        {
                            Real1DArray_Item ( mmCharges, n )                                 = 2.0e+00 * q0 ;
                            Real1DArray_Item ( mmCharges, self->mmpselection->positions[p] ) -=           q0 ;
                        }
                    }
                }
            }
        }
    }
    return mmCharges ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the MM coordinates for the QC/MM interactions (RC/RD coupling only).
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCAtomContainer_GetMMCoordinates3 ( const QCAtomContainer *self, const Coordinates3 *coordinates3, Coordinates3 *mmcoordinates3 )
{
    if ( ( self != NULL ) && ( coordinates3 != NULL ) && ( mmcoordinates3 != NULL ) )
    {
        auto Real xm, xp, ym, yp, zm, zp ;
        auto Integer    b, i, m, n, p, s, t       ;
        /* . Pure MM atoms. */
        for ( n = s = 0 ; s < self->mmpselection->nindices ; n++, s++ )
        {
            m = self->mmpselection->indices[s] ;
            Coordinates3_GetRow ( coordinates3,   m, xm, ym, zm ) ;
            Coordinates3_SetRow ( mmcoordinates3, n, xm, ym, zm ) ;
        }
        /* . Boundary MM atoms and their partners. */
        for ( i = 0 ; i < self->natoms ; i++ )
        {
            if ( self->data[i].QBOUNDARY )
            {
                b = self->data[i].index ;
                Coordinates3_GetRow ( coordinates3, b, xm, ym, zm ) ;
                for ( t = 0 ; t < self->data[i].nmmpartners ; n++, t++ )
                {
                    p = self->data[i].mmpartners[t] ;
                    Coordinates3_GetRow ( coordinates3,   p, xp, yp, zp ) ;
                    Coordinates3_SetRow ( mmcoordinates3, n, 0.5e+00 * ( xm + xp ), 0.5e+00 * ( ym + yp ), 0.5e+00 * ( zm + zp ) ) ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the nuclear charges.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCAtomContainer_GetNuclearCharges ( const QCAtomContainer *self, const QCParameter *qcParameters, Real1DArray *qcCharges )
{
    if ( ( self != NULL ) && ( qcParameters != NULL ) && ( qcCharges != NULL ) && ( self->natoms <= qcCharges->length ) )
    {
        auto Integer iatom, ic ;
        auto Real    q ;
        for ( iatom = 0 ; iatom < self->natoms ; iatom++ )
        {
            ic = self->data[iatom].center ;
            if ( qcParameters->centers[ic].mndoparameters == NULL ) q = ( Real ) qcParameters->centers[ic].atomicNumber ;
            else                                                    q = ( Real ) qcParameters->centers[ic].mndoparameters->zcore ;
            Real1DArray_Item ( qcCharges, iatom ) = q ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the gradients for an atom without unit conversion.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCAtomContainer_SetAtomGradients3 ( const QCAtomContainer *self, const Integer q, const Coordinates3 *coordinates3,
                                              const Real gx, const Real gy, const Real gz, Coordinates3 *gradients3 )
{
    if ( ( self != NULL ) && ( coordinates3 != NULL ) && ( gradients3 != NULL ) )
    {
        if ( self->data[q].QBOUNDARY )
        {
            auto Real d, dg, tx, ty, tz, xij, yij, zij ;
            /* . First contribution. */
            Coordinates3_IncrementRow ( gradients3, self->data[q].qcpartner, gx, gy, gz ) ;
            /* . Ratio of distances. */
            if ( self->QLINKRATIO )
            {
                /* . Scale gradient. */
                tx = gx * self->data[q].linkfactor ;
                ty = gy * self->data[q].linkfactor ;
                tz = gz * self->data[q].linkfactor ;
                /* . Second contribution. */
                Coordinates3_IncrementRow ( gradients3, self->data[q].index,     tx, ty, tz ) ;
                Coordinates3_DecrementRow ( gradients3, self->data[q].qcpartner, tx, ty, tz ) ;
            }
            /* . Absolute distance. */
            else
            {
                /* . Normalized interatomic vector. */
                Coordinates3_DifferenceRow ( coordinates3, self->data[q].index, self->data[q].qcpartner, xij, yij, zij ) ;
                d = sqrt ( xij * xij + yij * yij + zij * zij ) ;
                xij /= d ;
                yij /= d ;
                zij /= d ;
                /* . Scale gradient. */
                tx = gx * ( self->data[q].linkfactor / d ) ;
                ty = gy * ( self->data[q].linkfactor / d ) ;
                tz = gz * ( self->data[q].linkfactor / d ) ;
                /* . Dot product. */
                dg = xij * tx + yij * ty + zij * tz ;
                /* . Modified gradient. */
                tx -= dg * xij ;
                ty -= dg * yij ;
                tz -= dg * zij ;
                /* . Second contribution. */
                Coordinates3_IncrementRow ( gradients3, self->data[q].index,     tx, ty, tz ) ;
                Coordinates3_DecrementRow ( gradients3, self->data[q].qcpartner, tx, ty, tz ) ;
            }
        }
        else Coordinates3_IncrementRow ( gradients3, self->data[q].index, gx, gy, gz ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert QC atom gradients into atom gradients with optional unit conversion
! . from atomic units to regular units.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define CONVERSION_FACTOR ( UNITS_LENGTH_ANGSTROMS_TO_BOHRS * UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE )
Status QCAtomContainer_SetGradients3 ( const QCAtomContainer *self, const Coordinates3 *coordinates3, const Coordinates3 *qcgradients3, const Boolean QCONVERT, Coordinates3 **gradients3 )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) && ( qcgradients3 != NULL ) )
    {
        auto Real gx, gy, gz ;
        auto Integer    iqm ;
        /* . Allocate space if necessary. */
        if ( (*gradients3) == NULL )
        {
            (*gradients3) = Coordinates3_Allocate ( self->natoms ) ;
            Coordinates3_Set ( (*gradients3), 0.0e+00 ) ;
        }
        if ( (*gradients3) == NULL ) status = Status_OutOfMemory ;
        else
        {
            for ( iqm = 0 ; iqm < self->natoms ; iqm++ )
            {
                Coordinates3_GetRow ( qcgradients3, iqm, gx, gy, gz ) ;
                if ( QCONVERT )
                {
                    gx *= CONVERSION_FACTOR ;
                    gy *= CONVERSION_FACTOR ;
                    gz *= CONVERSION_FACTOR ;
                }
                QCAtomContainer_SetAtomGradients3 ( self, iqm, coordinates3, gx, gy, gz, (*gradients3) ) ;
            }
            status = Status_Success ;
        }
    }
    return status ;
}
# undef CONVERSION_FACTOR

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the MM gradients for the QC/MM interactions (RC/RD coupling only).
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCAtomContainer_SetMMGradients3 ( const QCAtomContainer *self, const Coordinates3 *mmgradients3, Coordinates3 *gradients3 )
{
    if ( ( self != NULL ) && ( mmgradients3 != NULL ) && ( gradients3 != NULL ) )
    {
        auto Real gx, gy, gz ;
        auto Integer    b, i, m, n, p, s, t ;
        /* . Pure MM atoms. */
        for ( n = s = 0 ; s < self->mmpselection->nindices ; n++, s++ )
        {
            m = self->mmpselection->indices[s] ;
            Coordinates3_GetRow       ( mmgradients3, n, gx, gy, gz ) ;
            Coordinates3_IncrementRow ( gradients3,   m, gx, gy, gz ) ;
        }
        /* . Boundary MM atoms and their partners. */
        for ( i = 0 ; i < self->natoms ; i++ )
        {
            if ( self->data[i].QBOUNDARY )
            {
                b = self->data[i].index ;
                for ( t = 0 ; t < self->data[i].nmmpartners ; n++, t++ )
                {
                    p = self->data[i].mmpartners[t] ;
                    Coordinates3_GetRow ( mmgradients3, n, gx, gy, gz ) ;
                    gx *= 0.5e+00 ;
                    gy *= 0.5e+00 ;
                    gz *= 0.5e+00 ;
                    Coordinates3_IncrementRow ( gradients3, b, gx, gy, gz ) ;
                    Coordinates3_IncrementRow ( gradients3, p, gx, gy, gz ) ;
                }
            }
        }
    }
}
