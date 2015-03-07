/*------------------------------------------------------------------------------
! . File      : QCOnePDM.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . One-particle density matrix procedures.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "Memory.h"
# include "QCOnePDM.h"

# include "LAPACKEigenvalueSolver.h"
# include "NMREigenvalueSolver.h"

# define NMREIGENVALUESOLVERNOCHECK

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Defaults. */
# define DefaultDensityType         QCOnePDMDensityType_Total
# define DefaultOccupancyType       QCOnePDMOccupancyType_Cardinal
# define DefaultFermiBroadening     1.0e+03
# define DefaultOccupancyFactor     2.0e+00

/* . Tolerance for cardinal occupancy check. */
# define CardinalOccupancyTolerance 1.0e-10

/* . Parameters for determination of variable occupancies. */
# define FermiBisectionTolerance    1.0e-12
# define FermiChargeTolerance       1.0e-10
# define FermiEnergyLower          -1000.0e+00
# define FermiEnergyUpper           1000.0e+00
# define FermiMaximumExponent       500.0e+00
# define FermiMaximumIterations     500

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real QCOnePDM_FermiOccupancies ( QCOnePDM *self ) ;
static void QCOnePDM_Occupancies      ( QCOnePDM *self ) ;
static void QCOnePDM_OccupancyEnergy  ( QCOnePDM *self ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCOnePDM *QCOnePDM_Allocate ( const Integer length, Status *status )
{
    QCOnePDM *self = NULL ;
    if ( length >= 0 )
    {
        MEMORY_ALLOCATE ( self, QCOnePDM ) ;
        if ( self != NULL )
        {
            /* . Basic initialization. */
            self->isValid         = False ;
            self->numberOccupied  = 0     ;
            self->numberOrbitals  = 0     ;
            self->fermiBroadening = DefaultFermiBroadening ;
            self->fermiEnergy     = 0.0e+00 ;
            self->occupancyEnergy = 0.0e+00 ;
            self->occupancyFactor = DefaultOccupancyFactor ;
            self->totalCharge     = 0.0e+00 ;
            self->densityType     = DefaultDensityType     ;
            self->occupancyType   = DefaultOccupancyType   ;
            self->energies        = NULL ;
            self->occupancies     = NULL ;
            self->orbitals        = NULL ;
            self->density         = NULL ;
            self->fock            = NULL ;
            /* . Allocate arrays . */
            self->occupancies     = Real1DArray_Allocate      ( length, status ) ;
            self->density         = SymmetricMatrix_AllocateN ( length, status ) ;
            if ( ( self->density == NULL ) || ( self->occupancies == NULL ) ) QCOnePDM_Deallocate ( &self ) ;
        }
        if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
    }
    else Status_Set ( status, Status_InvalidDimension ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate the orbital data for a density.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCOnePDM_AllocateOrbitalData ( QCOnePDM *self, Real2DArray *orthogonalizer, Status *status )
{
    if ( self != NULL )
    {
        if ( orthogonalizer == NULL ) self->numberOrbitals = self->density->dimension ;
        else                          self->numberOrbitals = orthogonalizer->length1  ;
        self->fermiEnergy     = 0.0e+00 ;
        self->occupancyEnergy = 0.0e+00 ;
        self->energies        = Real1DArray_Allocate      ( self->numberOrbitals, status ) ;
        self->orbitals        = Real2DArray_Allocate      ( self->density->dimension, self->numberOrbitals, status ) ;
        self->fock            = SymmetricMatrix_AllocateN ( self->density->dimension, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert total and spin densities to alpha and beta densities.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCOnePDM_AlphaBetaFromTotalSpin ( QCOnePDM *self, QCOnePDM *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        auto Real    qa, qb ;
        qa = 0.5e+00 * ( self->totalCharge + other->totalCharge ) ;
        qb = 0.5e+00 * ( self->totalCharge - other->totalCharge ) ;
        self ->totalCharge = qa ;
        other->totalCharge = qb ;
        for ( i = 0 ; i < self->density->size ; i++ )
        {
            qa = 0.5e+00 * ( self->density->data[i] + other->density->data[i] ) ;
            qb = 0.5e+00 * ( self->density->data[i] - other->density->data[i] ) ;
            self ->density->data[i] = qa ;
            other->density->data[i] = qb ;
        }
        self ->densityType     = QCOnePDMDensityType_Alpha ;
        other->densityType     = QCOnePDMDensityType_Beta  ;
        self ->occupancyFactor = 1.0e+00 ;
        other->occupancyFactor = 1.0e+00 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert alpha and beta densities to total and spin densities.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCOnePDM_AlphaBetaToTotalSpin ( QCOnePDM *self, QCOnePDM *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        auto Real    qs, qt ;
        qt = self->totalCharge + other->totalCharge ;
        qs = self->totalCharge - other->totalCharge ;
        self ->totalCharge = qt ;
        other->totalCharge = qs ;
        for ( i = 0 ; i < self->density->size ; i++ )
        {
            qt = self->density->data[i] + other->density->data[i] ;
            qs = self->density->data[i] - other->density->data[i] ;
            self ->density->data[i] = qt ;
            other->density->data[i] = qs ;
        }
        self ->densityType = QCOnePDMDensityType_Total ;
        other->densityType = QCOnePDMDensityType_Spin  ;
        self ->occupancyFactor = 2.0e+00 ;
        other->occupancyFactor = 1.0e+00 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCOnePDM *QCOnePDM_Clone ( const QCOnePDM *self, Status *status )
{
    QCOnePDM *other = NULL ;
    if ( self != NULL )
    {
        other = QCOnePDM_Allocate ( self->density->dimension, status ) ;
        if ( other != NULL )
        {
            /* . Scalar data. */
            other->isValid         = self->isValid         ;
            other->numberOccupied  = self->numberOccupied  ;
            other->numberOrbitals  = self->numberOrbitals  ;
            other->fermiBroadening = self->fermiBroadening ;
            other->fermiEnergy     = self->fermiEnergy     ;
            other->occupancyEnergy = self->occupancyEnergy ;
            other->occupancyFactor = self->occupancyFactor ;
            other->totalCharge     = self->totalCharge     ;
            other->densityType     = self->densityType     ;
            other->occupancyType   = self->occupancyType   ;
            /* . Array data. */
            Real1DArray_CopyTo     ( self->occupancies, other->occupancies, status ) ;
            SymmetricMatrix_CopyTo ( self->density    , other->density     ) ;
            other->energies        = Real1DArray_Clone     ( self->energies, status ) ;
            other->orbitals        = Real2DArray_Clone     ( self->orbitals, status ) ;
            other->fock            = SymmetricMatrix_Clone ( self->fock             ) ;
        }
    }
    return other ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCOnePDM_Deallocate ( QCOnePDM **self )
{
    if ( (*self) != NULL )
    {
        QCOnePDM_DeallocateOrbitalData (   (*self) ) ;
        Real1DArray_Deallocate         ( &((*self)->occupancies) ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->density    ) ) ;
        Memory_Deallocate              (   (*self) ) ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocate or reset orbital-related data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCOnePDM_DeallocateOrbitalData ( QCOnePDM *self )
{
    if ( self != NULL )
    {
        if ( self->occupancyType == QCOnePDMOccupancyType_FractionalVariable ) self->numberOccupied = 0 ;
        self->numberOrbitals  = 0 ;
        self->fermiEnergy     = 0.0e+00 ;
        self->occupancyEnergy = 0.0e+00 ;
        Real1DArray_Deallocate     ( &(self->energies) ) ;
        Real2DArray_Deallocate     ( &(self->orbitals) ) ;
        SymmetricMatrix_Deallocate ( &(self->fock    ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the density value at a grid point given a set of function values.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCOnePDM_DensityGridValues ( const QCOnePDM *self, const GridFunctionDataBlock *basisData, Real1DArray *rho, Status *status )
{
    Real1DArray_Set ( rho, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( basisData != NULL ) && ( rho != NULL ) )
    {
        auto Integer          f, g ;
        auto Real2DArray     *b = basisData->f, *new, *work ;
        auto SymmetricMatrix *old ;
        f   = b->length0 ;
        g   = b->length1 ;
        old = self->density ;
        if ( ( f <= old->dimension ) && ( g == rho->length ) )
        {
            new  = Real2DArray_Allocate ( f, f, status ) ;
            work = Real2DArray_Allocate ( f, g, status ) ;
            if ( ( new != NULL ) && ( work != NULL ) )
            {
                SymmetricMatrix_IndexedCopyToReal2DArray ( old, basisData->indices, new, NULL ) ;
                Real2DArray_MatrixMultiply               ( False, False, 1.0e+00, new, b, 0.0e+00, work, NULL ) ;
                Real2DArray_ColumnDotProducts            ( True, b, work, rho ) ;
            }
            Real2DArray_Deallocate ( &new  ) ;
            Real2DArray_Deallocate ( &work ) ;
        }
        else Status_Set ( status, Status_ArrayLengthMismatch ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the occupancies and total charge given a Fermi energy.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real QCOnePDM_FermiOccupancies ( QCOnePDM *self )
{
    Real q = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer i ;
        auto Real    e ;
        for ( i = 0 ; i < self->numberOrbitals ; i++ )
        {
            e = - self->fermiBroadening * ( self->fermiEnergy - Real1DArray_Item ( self->energies, i ) ) ;
                 if ( e >   FermiMaximumExponent ) Real1DArray_Item ( self->occupancies, i ) = 0.0e+00 ;
            else if ( e < - FermiMaximumExponent ) Real1DArray_Item ( self->occupancies, i ) = 1.0e+00 ;
            else Real1DArray_Item ( self->occupancies, i ) = 1.0e+00 / ( 1.0e+00 + exp ( e ) ) ;
            q += Real1DArray_Item ( self->occupancies, i ) ;
        }
    }
    return q ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Obtain a density with a diagonal guess.
! . The occupancies are assigned for the fixed occupancy options.
! . For fixed fractional occupancies, the charge in the X HOOs is spread evenly over these and the Y LUOs.
!---------------------------------------------------------------------------------------------------------------------------------*/
QCOnePDM *QCOnePDM_FromDiagonalGuess ( const QCOnePDMDensityType densityType, const QCOnePDMOccupancyType occupancyType, const Integer length, const Real totalCharge, const Integer numberFractionalHOOs, const Integer numberFractionalLUOs, Status *status )
{
    Boolean   isOK = True ;
    Integer   numberCardinal = 0, numberFractional = 0, numberOccupied = 0 ;
    Real      fractionalCharge = 0.0e+00, occupancyFactor ;
    QCOnePDM *self = NULL ;

    /* . Basic checks. */
    /* . Get the occupancy factor and estimate the number of occupied orbitals given the charge. */
    if ( densityType == QCOnePDMDensityType_Total ) occupancyFactor = 2.0e+00 ;
    else                                            occupancyFactor = 1.0e+00 ;
    numberOccupied = ( Integer ) Round ( ceil ( totalCharge / occupancyFactor ) ) ;

    /* . Cardinal occupancies. */
    if ( occupancyType == QCOnePDMOccupancyType_Cardinal )
    {
        isOK           = ( fabs (  ( Real ) numberOccupied * occupancyFactor - totalCharge ) < CardinalOccupancyTolerance ) ;
        numberCardinal = numberOccupied ;
    }
    /* . Fixed fractional occupancies. */
    else if ( occupancyType == QCOnePDMOccupancyType_FractionalFixed )
    {
        numberCardinal   = numberOccupied - numberFractionalHOOs       ;
        numberFractional = numberFractionalHOOs + numberFractionalLUOs ;
        numberOccupied   = numberCardinal + numberFractional ;
        fractionalCharge = ( ( totalCharge / occupancyFactor ) - ( Real ) numberCardinal ) / ( Real ) ( Maximum ( numberFractional, 1 ) ) ;
        isOK             = ( numberFractionalHOOs > 0 ) && ( numberFractionalLUOs >= 0 ) && ( numberCardinal >= 0 ) && ( ( numberCardinal + numberFractional ) <= length ) && ( fractionalCharge >= 0.0e+00 ) && ( fractionalCharge <= 1.0e+00 ) ;
    }

    /* . Allocation. */
    if ( isOK && ( length >= 0 ) && ( totalCharge >= 0.0e+00 ) )
    {
        self = QCOnePDM_Allocate ( length, status ) ;
        if ( self != NULL )
        {
            auto Integer i ;

            /* . Basic initialization. */
            self->densityType     = densityType     ;
            self->occupancyType   = occupancyType   ;
            self->occupancyFactor = occupancyFactor ;
            self->totalCharge     = totalCharge     ;

            /* . Density and occupancies. */
            if ( length > 0 )
            {
                /* . Density. */
                SymmetricMatrix_Set          ( self->density, 0.0e+00 ) ;
                SymmetricMatrix_Set_Diagonal ( self->density, 0, self->density->dimension, totalCharge / ( ( Real ) length ) ) ;

                /* . Occupancies. */
                Real1DArray_Set ( self->occupancies, 0.0e+00 ) ;
                if ( occupancyType != QCOnePDMOccupancyType_FractionalVariable )
                {
                    self->numberOccupied = numberOccupied ;
                    for ( i = 0 ; i < numberCardinal   ; i++ ) Real1DArray_Item ( self->occupancies, i                  ) =                    occupancyFactor ;
                    for ( i = 0 ; i < numberFractional ; i++ ) Real1DArray_Item ( self->occupancies, i + numberCardinal ) = fractionalCharge * occupancyFactor ;
                }
            }
        }
    }
    else Status_Set ( status, Status_InvalidArgument ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the index of the HOMO.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer QCOnePDM_HOMOIndex ( const QCOnePDM *self )
{
    Integer index = -1 ;
    if ( ( self != NULL ) && ( self->orbitals != NULL ) ) index = self->numberOccupied - 1 ;
    return index ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the index of the LUMO.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer QCOnePDM_LUMOIndex ( const QCOnePDM *self )
{
    Integer index = -1 ;
    if ( ( self != NULL ) && ( self->orbitals != NULL ) ) index = Minimum ( self->numberOccupied, self->orbitals->length1 ) ;
    return index ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a density matrix from orbitals.
! . The RMS difference between the old and new densities is returned.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real QCOnePDM_Make ( QCOnePDM *self, SymmetricMatrix *scratch, Status *status )
{
    Real rmsDifference = 0.0e+00 ;
    if ( self != NULL )
    {
        if ( self->numberOccupied > 0 )
        {
            auto SymmetricMatrix *new = self->density, *old ;
            if ( scratch == NULL ) old = SymmetricMatrix_AllocateN ( new->dimension, status ) ;
            else                   old = scratch ;
            SymmetricMatrix_CopyTo              ( new, old ) ;
            SymmetricMatrix_MakeFromEigensystem ( new, True, self->numberOccupied, self->occupancies, self->orbitals, status ) ;
            SymmetricMatrix_AddScaledMatrix     ( old, -1.0e+00, new ) ;
            rmsDifference = sqrt ( SymmetricMatrix_Trace2 ( old, old, NULL ) ) / ( Real ) ( old->dimension ) ;
            if ( scratch == NULL ) SymmetricMatrix_Deallocate ( &old ) ;
        }
        else SymmetricMatrix_Set ( self->density, 0.0e+00 ) ;
        /* . The density matrix is valid even if it has no electrons. */
        self->isValid = True ;
    }
    return rmsDifference ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a density from the individual components.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCOnePDM_MakeElementary ( SymmetricMatrix *self, const Integer numberOccupied, const Real1DArray *occupancies, const Real2DArray *orbitals, Status *status )
{
    SymmetricMatrix_MakeFromEigensystem ( self, True, numberOccupied, occupancies, orbitals, status ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Compute a new density given a Fock matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real QCOnePDM_MakeFromFock ( QCOnePDM *self, Real2DArray *orthogonalizer, Status *status )
{
    Real rmsDifference = 0.0e+00 ;
    if ( self != NULL )
    {
        if ( self->numberOrbitals > 0 )
        {
            auto Real2DArray     *eigenVectors = NULL ;
            auto SymmetricMatrix *fock         = NULL ;

            /* . Transform the Fock matrix to the orthogonal basis if necessary. */
            if ( orthogonalizer == NULL )
            {
                eigenVectors = self->orbitals ;
                fock         = self->fock ;
            }
            else
            {
                eigenVectors = Real2DArray_Allocate      ( orthogonalizer->length1, orthogonalizer->length1, status ) ;
                fock         = SymmetricMatrix_AllocateN ( orthogonalizer->length1, status ) ;
                SymmetricMatrix_Transform ( self->fock, orthogonalizer, False, fock ) ;
            }

            /* . Diagonalize the Fock matrix. */
# ifdef USENMREIGENVALUESOLVER
# ifdef NMREIGENVALUESOLVERNOCHECK
            NMREigenvalueSolver_Solve ( fock, self->energies, eigenVectors, status ) ;
# else
            NMREigenvalueSolver_FullCheck ( fock, self->energies, eigenVectors, status ) ;
# endif
# else
            LAPACKEigenvalueSolver_Diagonalize ( fock, self->energies, eigenVectors, status ) ;
# endif

            /* . Back transform the orbitals - this could be done in place later. */
            if ( orthogonalizer != NULL )
            {
                Real2DArray_MatrixMultiply ( False, False, 1.0e+00, orthogonalizer, eigenVectors, 0.0e+00, self->orbitals, status ) ;
                Real2DArray_Deallocate     ( &eigenVectors ) ;
                SymmetricMatrix_Deallocate ( &fock ) ;
            }

            /* . Get the orbital occupancies and occupancy energy. */
            QCOnePDM_Occupancies ( self ) ;

            /* . Form the density matrix. */
            rmsDifference = QCOnePDM_Make ( self, NULL, NULL ) ;
        }
    }
    return rmsDifference ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a weighted density matrix from orbitals and their energies.
! . A symmetric matrix is returned.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Care needed here as the length of energies can be less than that of occupancies (e.g. for spherical basis sets). */
static void QCOnePDM_MakeWeightedOne ( const QCOnePDM        *self    , 
                                             Real1DArray     *scratch , 
                                             SymmetricMatrix *wDM     )
{
    if ( self != NULL )
    {
        Real1DArray_CopyTo ( self->energies , scratch , NULL ) ;
        if ( self->energies->length < self->occupancies->length )
        {
            auto Real1DArray slice ;
            Real1DArray_Slice ( self->occupancies, 0, self->energies->length, 1, &slice, NULL ) ;
            Real1DArray_Multiply ( scratch, &slice, NULL ) ;
        }
        else
        {
            Real1DArray_Multiply ( scratch, self->occupancies, NULL ) ;
        }
        SymmetricMatrix_MakeFromEigensystem ( wDM, False, self->numberOccupied, scratch, self->orbitals, NULL ) ;
    }
}

SymmetricMatrix *QCOnePDM_MakeWeighted ( const QCOnePDM *self, const QCOnePDM *other, Status *status )
{
    SymmetricMatrix *wDM = NULL ;
    if ( self != NULL )
    {
        auto Real1DArray *scratch = NULL ;
        scratch = Real1DArray_Allocate      ( self->energies->length  , status ) ;
        wDM     = SymmetricMatrix_AllocateN ( self->density->dimension, status ) ;
        if ( ( scratch != NULL ) && ( wDM != NULL ) )
        {
            SymmetricMatrix_Set ( wDM, 0.0e+00 ) ;
            QCOnePDM_MakeWeightedOne ( self  , scratch , wDM ) ; 
            QCOnePDM_MakeWeightedOne ( other , scratch , wDM ) ; 
            SymmetricMatrix_Scale ( wDM, -2.0e+00 ) ;
        }
        Real1DArray_Deallocate ( &scratch ) ;
     }
    return wDM ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Compute orbital occupancies and the occupancy energy when there are variable fractional occupancies.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define ProcedureName "QCOnePDM_Occupancies"
static void QCOnePDM_Occupancies ( QCOnePDM *self )
{
    if ( self != NULL )
    {
        /* . Initialization. */
        self->fermiEnergy     = 0.0e+00 ;
        self->occupancyEnergy = 0.0e+00 ;

        /* . Variable fractional orbitals. */
        if ( self->occupancyType == QCOnePDMOccupancyType_FractionalVariable )
        {
            auto Real targetCharge ;

            /* . Initialization. */
            self->numberOccupied = 0 ;
            targetCharge         = self->totalCharge / self->occupancyFactor ;
            Real1DArray_Set ( self->occupancies, 0.0e+00 ) ;
            if ( targetCharge > 0.0e+00 )
            {
                auto Integer i ;
                auto Real    dx, eflower, efupper, flower, fmiddle, fupper ;

                /* . Calculate the function at lower and upper bounds to the Fermi energy. */
                eflower = FermiEnergyLower ; self->fermiEnergy = eflower ; flower = QCOnePDM_FermiOccupancies ( self ) - targetCharge ;
                efupper = FermiEnergyUpper ; self->fermiEnergy = efupper ; fupper = QCOnePDM_FermiOccupancies ( self ) - targetCharge ;

                /* . Check to see if a root is bracketed. */
                if ( ( fupper < 0.0e+00 ) || ( flower > 0.0e+00 ) ) printf ( "%s%s\n", ProcedureName, "Extreme orbital energies do not satisfy Fermi function." ) ;

                /* . Iterate. */
                for ( i = 0, dx = ( efupper - eflower ) ; i < FermiMaximumIterations ; i++ )
                {
	           dx                = 0.5e+00 * dx ;
	           self->fermiEnergy = eflower + dx ;
	           fmiddle           = QCOnePDM_FermiOccupancies ( self ) - targetCharge ;
	           if ( fmiddle <= 0.0e+00 ) eflower = self->fermiEnergy ;
	           if ( fabs ( fmiddle ) <= FermiBisectionTolerance ) break ;
                }

                /* . There have been too many bisections. */
                if ( fabs ( fmiddle ) > FermiBisectionTolerance ) printf ( "%s%s\n", ProcedureName, "Unable to find Fermi energy to required accuracy." ) ;

                /* . Get the number of occupied orbitals and the occupancy energy. */
                QCOnePDM_OccupancyEnergy ( self ) ;

                /* . Scale the occupancies and occupancy energy if necessary. */
                if ( self->occupancyFactor != 1.0e+00 )
                {
                    self->occupancyEnergy *= self->occupancyFactor ;
                    Real1DArray_Scale ( self->occupancies, self->occupancyFactor ) ;
                }
            }
        }
        /* . Other cases. */
        else if ( self->numberOccupied > 0 )
        {
            /* . Set the Fermi level as the average of the HOMO/LUMO energy if possible. */
            if ( self->numberOccupied < self->numberOrbitals ) self->fermiEnergy = 0.5e+00 * ( Real1DArray_Item ( self->energies, self->numberOccupied-1 ) + Real1DArray_Item ( self->energies, self->numberOccupied ) ) ;
            else                                               self->fermiEnergy =             Real1DArray_Item ( self->energies, self->numberOccupied-1 ) ;
        }
    }
}
# undef ProcedureName

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the number of occupied orbitals and the occupancy energy (variable occupancies only).
!---------------------------------------------------------------------------------------------------------------------------------*/
static void QCOnePDM_OccupancyEnergy ( QCOnePDM *self )
{
    if ( self != NULL )
    {
        self->occupancyEnergy = 0.0e+00 ;
        if ( self->occupancyType == QCOnePDMOccupancyType_FractionalVariable )
        {
            auto Integer i ;
            auto Real    e, p, q ;
            self->numberOccupied = 0 ;
            for ( i = 0, e = 0.0e+00 ; i < self->numberOrbitals ; i++ )
            {
                p = Real1DArray_Item ( self->occupancies, i ) ;
                q = 1.0e+00 - p ;
                if ( p > FermiChargeTolerance )
                {
                    self->numberOccupied ++ ;
                    e += p * log ( p ) ;
                    if ( q > FermiChargeTolerance ) e += q * log ( q ) ;
                }
            }
            self->occupancyEnergy = e / self->fermiBroadening ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Compute a new density given an old density and an orthogonalizer.
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCOnePDM_Reorthogonalize ( QCOnePDM *self, Real2DArray *orthogonalizer, Status *status )
{
    if ( ( self != NULL ) && ( orthogonalizer != NULL ) )
    {
        auto Real2DArray     *eigenVectors  ;
        auto SymmetricMatrix *density       ;

        /* . Transform the density matrix to the orthogonal basis. */
        density = SymmetricMatrix_AllocateN ( orthogonalizer->length1, status ) ;
        SymmetricMatrix_Transform ( self->density, orthogonalizer, False, density ) ;

        /* . Diagonalize the transformed density. */
        eigenVectors = Real2DArray_Allocate ( orthogonalizer->length1, orthogonalizer->length1, status ) ;
        SymmetricMatrix_Diagonalize ( density, self->energies, eigenVectors, status ) ;
        Real1DArray_Scale ( self->energies, -1.0e+00 ) ;

        /* . Back transform the orbitals. */
        Real2DArray_MatrixMultiply ( False, False, 1.0e+00, orthogonalizer, eigenVectors, 0.0e+00, self->orbitals, status ) ;
        Real2DArray_Deallocate     ( &eigenVectors ) ;
        SymmetricMatrix_Deallocate ( &density      ) ;

        /* . Get the number of occupied orbitals and the occupancy energy. */
        QCOnePDM_Occupancies ( self ) ;

        /* . Form the density matrix. */
        QCOnePDM_MakeElementary ( self->density, self->numberOccupied, self->occupancies, self->orbitals, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Reset the density from another density.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Could also test for idempotency here or, more generally, 0<=Tr(PSP)<=Tr(PS). */
/* . The assumption is made that the density is valid. */

# define ChargeTolerance 1.0e-6
Boolean QCOnePDM_ResetDensityFromDensity ( QCOnePDM *self, SymmetricMatrix *density, SymmetricMatrix *overlap, Status *status )
{
    Boolean isOK = False ;
    if ( ( self != NULL ) && ( density != NULL ) )
    {
        auto Real totalCharge ;

        /* . Trace. */
        if ( overlap == NULL ) totalCharge = SymmetricMatrix_Trace ( density ) ;
        else                   totalCharge = SymmetricMatrix_Multiply2_Trace ( density, overlap ) ;

        /* . Basic checks only. */
        isOK = ( self->density->dimension == density->dimension ) && ( fabs ( self->totalCharge - totalCharge ) <= ChargeTolerance ) ;

        /* . Copy over if OK. */
        if ( isOK ) { SymmetricMatrix_CopyTo ( density, self->density ) ; self->isValid = True ; }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return isOK ;
}
# undef ChargeTolerance

/*----------------------------------------------------------------------------------------------------------------------------------
! . Reset the density from a set of orbitals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The assumption is made that the density is valid. */

Boolean QCOnePDM_ResetDensityFromOrbitals ( QCOnePDM *self, Real2DArray *orbitals, Status *status )
{
    Boolean isOK = False ;
    if ( ( self != NULL ) && ( orbitals != NULL ) )
    {
        auto Status localStatus ;
        Status_Set ( &localStatus, Status_Continue ) ;

        /* . Get occupancies for variable fractional. */
        if ( self->occupancyType == QCOnePDMOccupancyType_FractionalVariable )
        {
            auto Integer i, numberCardinal ;
            auto Real    fractionalCharge ;
            numberCardinal = ( Integer ) Round ( floor ( self->totalCharge / self->occupancyFactor ) ) ;
            self->numberOccupied = numberCardinal ;
            for ( i = 0 ; i < numberCardinal ; i++ ) Real1DArray_Item ( self->occupancies, i ) = self->occupancyFactor ;
            fractionalCharge = ( self->totalCharge - self->occupancyFactor * ( ( Real ) numberCardinal ) ) ;
            if ( fractionalCharge > 0.0e+00 )
            {
                self->numberOccupied += 1 ;
                Real1DArray_Item ( self->occupancies, numberCardinal ) = fractionalCharge ;
            }
        }

        /* . Make the density. */
        QCOnePDM_MakeElementary ( self->density, self->numberOccupied, self->occupancies, orbitals, &localStatus ) ;

        /* . Reset occupancies for variable fractional. */
        if ( self->occupancyType == QCOnePDMOccupancyType_FractionalVariable )
        {
            self->numberOccupied = 0 ;
            Real1DArray_Set ( self->occupancies, 0.0e+00 ) ;
        }

        /* . Finish up. */
        isOK = Status_OK ( &localStatus ) ;
        if ( isOK ) self->isValid = True ;
        else Status_Set ( status, localStatus ) ;
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Compute spin expectation values for a density.
! . This procedure requires alpha and beta densities.
! . The beta density should be present even if there is no beta charge
! . (this applies to few cases).
!---------------------------------------------------------------------------------------------------------------------------------*/
void QCOnePDM_SpinExpectationValues ( QCOnePDM *self, QCOnePDM *other, SymmetricMatrix *overlap, Real *Sz, Real *S2 )
{
    /* . Initialization. */
    (*Sz) = 0.0e+00 ;
    (*S2) = 0.0e+00 ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Real trace ;
        if ( overlap == NULL ) trace = SymmetricMatrix_Multiply2_Trace    ( self->density, other->density ) ;
        else                   trace = SymmetricMatrix_MultiplyASBS_Trace ( self->density, other->density, overlap, NULL ) ;

        /* . Calculate Sz. */
        (*Sz) = 0.5e+00 * ( self->totalCharge - other->totalCharge ) ;

        /* . Calculate S2. */
        (*S2) = ( (*Sz) * (*Sz) ) + 0.5e+00 * ( self->totalCharge + other->totalCharge ) - trace ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Undefine local parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# undef DefaultDensityType
# undef DefaultFermiBroadening
# undef DefaultOccupancyType
# undef DefaultOccupancyFactor

# undef CardinalOccupancyTolerance

# undef FermiBisectionTolerance
# undef FermiChargeTolerance
# undef FermiEnergyLower
# undef FermiEnergyUpper
# undef FermiMaximumExponent
# undef FermiMaximumIterations
