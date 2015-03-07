/*------------------------------------------------------------------------------
! . File      : DIISSCFConvergerState.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . The state for DIIS SCF convergers.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "DIISSCFConvergerState.h"
# include "Memory.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DIISSCFConvergerState *DIISSCFConvergerState_Allocate ( void )
{
    DIISSCFConvergerState *self = NULL ;
    self = ( DIISSCFConvergerState * ) Memory_Allocate ( sizeof ( DIISSCFConvergerState ) ) ;
    if ( self != NULL )
    {
        /* . Scalars. */
        self->isConverged           = False ;
        self->QDAMP                 = False ;
        self->QDIIS                 = False ;
        self->QRCA                  = False ;
        self->QRCAPREDICTED         = False ;
        self->iteration             = 0 ;
        self->matnum                = 0 ;
        self->ndiis                 = 0 ;
        self->damp                  = 0.0e+00 ;
        self->deavg                 = 0.0e+00 ;
        self->deltaeold             = 0.0e+00 ;
        self->diiserror             = 0.0e+00 ;
        self->dmpsav                = 0.0e+00 ;
        self->eold                  = 0.0e+00 ;
        self->rcamu                 = 1.0e+00 ; /* . Default first value. */
        self->rmsdifference         = 0.0e+00 ;
        /* . Pointers. */
        self->densityp       = NULL ;
        self->densityq       = NULL ;
        self->overlap        = NULL ;
        self->orthogonalizer = NULL ;
        /* . The DIIS matrix indices. */
        self->matind         = NULL ;
        /* . The DIIS error vector overlaps. */
        self->bcoeff         = NULL ;
        /* . Error vectors and Fock matrices. */
        /* . Initialization. */
        /* . Damping. */
        self->dfockp1   = NULL ;
        self->dfockp2   = NULL ;
        self->dfockq1   = NULL ;
        self->dfockq2   = NULL ;
        /* . DIIS. */
/*        self->asmwork   = NULL ;*/
        self->errorp    = NULL ;
        self->errorq    = NULL ;
        self->work      = NULL ;
        self->xdensityp = NULL ;
        self->xdensityq = NULL ;
        self->xfockp    = NULL ;
        self->xfockq    = NULL ;
        /* . RCA. */
        self->rdensityp = NULL ;
        self->rfockp    = NULL ;
        self->rdensityq = NULL ;
        self->rfockq    = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for convergence.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean DIISSCFConvergerState_Converged ( const DIISSCFConvergerState *self )
{
    Boolean isConverged = False ;
    if ( self != NULL )
    {
        isConverged = self->isConverged ;
    }
    return isConverged ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DIISSCFConvergerState_Deallocate ( DIISSCFConvergerState **self )
{
    if ( (*self) != NULL )
    {
        auto Integer i ;
        Real2DArray_Deallocate ( &((*self)->bcoeff)  ) ;
        Memory_Deallocate      (   (*self)->matind   ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->dfockp1  ) ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->dfockp2  ) ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->dfockq1  ) ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->dfockq2  ) ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->rdensityp) ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->rfockp   ) ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->rdensityq) ) ;
        SymmetricMatrix_Deallocate     ( &((*self)->rfockq   ) ) ;
/*        AntisymmetricMatrix_Deallocate ( &((*self)->asmwork  ) ) ;*/
        if ( (*self)->errorp != NULL )
        {
            for ( i = 0 ; i < (*self)->ndiis ; i++ ) { AntisymmetricMatrix_Deallocate ( &((*self)->errorp[i]   ) ) ;
                                                       SymmetricMatrix_Deallocate     ( &((*self)->xdensityp[i]) ) ;
                                                       SymmetricMatrix_Deallocate     ( &((*self)->xfockp[i]   ) ) ; }
            Memory_Deallocate ( (*self)->errorp    ) ;
            Memory_Deallocate ( (*self)->xdensityp ) ;
            Memory_Deallocate ( (*self)->xfockp    ) ;
        }
        if ( (*self)->errorq != NULL )
        {
            for ( i = 0 ; i < (*self)->ndiis ; i++ ) { AntisymmetricMatrix_Deallocate ( &((*self)->errorq[i]   ) ) ;
                                                       SymmetricMatrix_Deallocate     ( &((*self)->xdensityq[i]) ) ;
                                                       SymmetricMatrix_Deallocate     ( &((*self)->xfockq[i]   ) ) ; }
            Memory_Deallocate ( (*self)->errorq    ) ;
            Memory_Deallocate ( (*self)->xdensityq ) ;
            Memory_Deallocate ( (*self)->xfockq    ) ;
        }
        if ( (*self)->work != NULL )
        {
            for ( i = 0 ; i < 3 ; i++ ) Real2DArray_Deallocate ( &((*self)->work[i]) ) ;
            Memory_Deallocate ( (*self)->work ) ;
        }
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get table data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DIISSCFConvergerState_GetTableData ( const DIISSCFConvergerState *self, Integer *iteration, Real *deltaeold, Real *rmsdifference )
{
    if ( self != NULL )
    {
        (*iteration)     = self->iteration     ;
        (*deltaeold)     = self->deltaeold     ;
        (*rmsdifference) = self->rmsdifference ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set up a state.
!---------------------------------------------------------------------------------------------------------------------------------*/
DIISSCFConvergerState *DIISSCFConvergerState_SetUp ( const Boolean QUSERCA, const Integer ndiis, QCOnePDM *densityp, QCOnePDM *densityq, SymmetricMatrix *overlap, Real2DArray *orthogonalizer )
{
    DIISSCFConvergerState *self = NULL ;
    if ( densityp != NULL )
    {
        auto Integer  i, m, n ;
        /* . Allocate the structure. */
        self = DIISSCFConvergerState_Allocate ( ) ;
        /* . Options. */
        if ( QUSERCA ) self->QRCA  = True ;
        else           self->QDAMP = True ;
        self->ndiis          = ndiis          ;
        /* . Pointers. */
        self->densityp       = densityp       ;
        self->overlap        = overlap        ;
        self->orthogonalizer = orthogonalizer ;
        /* . Orthogonalizer specific items. */
        if ( ( overlap != NULL ) && ( orthogonalizer != NULL ) ) { m = orthogonalizer->length0      ; n = orthogonalizer->length1 ; }
        else                                                     { m = densityp->density->dimension ; n = m ; }
        /* . Is a second density required? */
        if ( densityq != NULL ) self->densityq = densityq ;
        /* . The DIIS matrix indices. */
        self->matind = Memory_Allocate_Array_Integer ( self->ndiis ) ;
        for ( i = 0 ; i < self->ndiis ; i++ ) self->matind[i] = i ;
        /* . The DIIS error vector overlaps. */
        self->bcoeff = Real2DArray_Allocate ( self->ndiis, self->ndiis, NULL ) ;
        Real2DArray_Set ( self->bcoeff, 0.0e+00 ) ;
        /* . Allocate space - first density. */
        /* . RCA. */
        if ( QUSERCA )
        {
            self->rdensityp = SymmetricMatrix_Allocate ( densityp->density->dimension ) ;
            self->rfockp    = SymmetricMatrix_Allocate ( densityp->density->dimension ) ;
        }
        /* . Damping. */
        else
        {
            self->dfockp1   = SymmetricMatrix_Allocate ( densityp->density->dimension ) ;
            self->dfockp2   = SymmetricMatrix_Allocate ( densityp->density->dimension ) ;
        }
        /* . DIIS. */
        if ( densityp->density->dimension > 1 )
        {
            self->errorp    = ( AntisymmetricMatrix ** ) Memory_Allocate_Array ( self->ndiis, sizeof ( AntisymmetricMatrix * ) ) ;
            self->xdensityp = ( SymmetricMatrix ** ) Memory_Allocate_Array ( self->ndiis, sizeof ( SymmetricMatrix * ) ) ;
            self->xfockp    = ( SymmetricMatrix ** ) Memory_Allocate_Array ( self->ndiis, sizeof ( SymmetricMatrix * ) ) ;
            for ( i = 0 ; i < self->ndiis ; i++ ) { self->errorp[i]    = AntisymmetricMatrix_Allocate ( n, NULL ) ;
                                                    self->xdensityp[i] = SymmetricMatrix_Allocate     ( densityp->density->dimension ) ;
                                                    self->xfockp[i]    = SymmetricMatrix_Allocate     ( densityp->density->dimension ) ; }
            /* . Work space. */
            self->work = ( Real2DArray ** ) Memory_Allocate_Array ( 3, sizeof ( Real2DArray * ) ) ;
            for ( i = 0 ; i < 3 ; i++ ) self->work[i] = Real2DArray_Allocate ( m, m, NULL ) ;
        }
        /* . Allocate space - second density. */
        if ( densityq != NULL )
        {
            if ( densityq->density->dimension > 1 )
            {
               /* . RCA. */
               if ( QUSERCA )
               {
                  self->rdensityq = SymmetricMatrix_Allocate ( densityq->density->dimension ) ;
                  self->rfockq    = SymmetricMatrix_Allocate ( densityq->density->dimension ) ;
               }
               /* . Damping. */
               {
                  self->dfockq1   = SymmetricMatrix_Allocate ( densityp->density->dimension ) ;
                  self->dfockq2   = SymmetricMatrix_Allocate ( densityp->density->dimension ) ;
               }
               /* . DIIS. */
               self->errorq    = (  AntisymmetricMatrix ** ) Memory_Allocate_Array ( self->ndiis, sizeof ( AntisymmetricMatrix * ) ) ;
               self->xdensityq = ( SymmetricMatrix ** ) Memory_Allocate_Array ( self->ndiis, sizeof ( SymmetricMatrix * ) ) ;
               self->xfockq    = ( SymmetricMatrix ** ) Memory_Allocate_Array ( self->ndiis, sizeof ( SymmetricMatrix * ) ) ;
               for ( i = 0 ; i < self->ndiis ; i++ ) { self->errorq[i]    = AntisymmetricMatrix_Allocate ( n, NULL ) ;
                                                       self->xdensityq[i] = SymmetricMatrix_Allocate     ( densityq->density->dimension ) ;
                                                       self->xfockq[i]    = SymmetricMatrix_Allocate     ( densityq->density->dimension ) ; }
            }
        }
    }
    return self ;
}

