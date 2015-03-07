/*------------------------------------------------------------------------------
! . File      : QCParameters.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "Memory.h"
# include "QCParameters.h"

/*==============================================================================
! . Standard procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
QCParameter *QCParameters_Allocate ( const Integer ncenters )
{
    QCParameter *self = NULL ;
    if ( ncenters > 0 )
    {
        auto Integer i ;
        self           = ( QCParameter * ) Memory_Allocate ( sizeof ( QCParameter ) ) ;
        self->centers  = ( QCCenter    * ) Memory_Allocate_Array ( ncenters, sizeof ( QCCenter ) ) ;
        self->ncenters = ncenters ;
        for ( i = 0 ; i < ncenters ; i++ )
        {
            self->centers[i].atomicNumber   = -1 ;
            self->centers[i].densitybasis   = NULL ;
            self->centers[i].orbitalbasis   = NULL ;
            self->centers[i].poissonbasis   = NULL ;
            self->centers[i].mndoparameters = NULL ;
        }
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
QCParameter *QCParameters_Clone ( const QCParameter *self )
{
    QCParameter *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = QCParameters_Allocate ( self->ncenters ) ;
        for ( i = 0 ; i < new->ncenters ; i++ )
        {
            new->centers[i].atomicNumber   = self->centers[i].atomicNumber ;
            new->centers[i].densitybasis   = GaussianBasis_Clone  ( self->centers[i].densitybasis   ) ;
            new->centers[i].orbitalbasis   = GaussianBasis_Clone  ( self->centers[i].orbitalbasis   ) ;
            new->centers[i].poissonbasis   = GaussianBasis_Clone  ( self->centers[i].poissonbasis   ) ;
            new->centers[i].mndoparameters = MNDOParameters_Clone ( self->centers[i].mndoparameters ) ;
        }
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void QCParameters_Deallocate ( QCParameter **self )
{
    if ( (*self) != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < (*self)->ncenters ; i++ )
        {
            GaussianBasis_Deallocate  ( &((*self)->centers[i].densitybasis)   ) ;
            GaussianBasis_Deallocate  ( &((*self)->centers[i].orbitalbasis)   ) ;
            GaussianBasis_Deallocate  ( &((*self)->centers[i].poissonbasis)   ) ;
            MNDOParameters_Deallocate ( &((*self)->centers[i].mndoparameters) ) ;
        }
        Memory_Deallocate ( (*self)->centers ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}
