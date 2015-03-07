/*------------------------------------------------------------------------------
! . File      : Transformation3.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include "Memory.h"
# include "Transformation3.h"

/*------------------------------------------------------------------------------
! . Allocation.
! . The matrix is not initialized.
!-----------------------------------------------------------------------------*/
Transformation3 *Transformation3_Allocate ( void )
{
    Transformation3 *self = NULL ;
    self              = ( Transformation3 * ) Memory_Allocate ( sizeof ( Transformation3 ) ) ;
    self->rotation    = Matrix33_Allocate ( ) ;
    self->translation = Vector3_Allocate  ( ) ;
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
Transformation3 *Transformation3_Clone ( const Transformation3 *self )
{
    Transformation3 *new = NULL ;
    if ( self != NULL )
    {
        new = Transformation3_Allocate ( ) ;
	Matrix33_CopyTo ( self->rotation,    new->rotation   , NULL ) ;
	Vector3_CopyTo  ( self->translation, new->translation, NULL ) ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Copying.
!-----------------------------------------------------------------------------*/
void Transformation3_Copy ( Transformation3 *self, const Transformation3 *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        Matrix33_CopyTo ( other->rotation,    self->rotation   , NULL ) ;
        Vector3_CopyTo  ( other->translation, self->translation, NULL ) ;
    }
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void Transformation3_Deallocate ( Transformation3 **self )
{
    if ( (*self) != NULL )
    {
        Matrix33_Deallocate ( &((*self)->rotation)    ) ;
        Vector3_Deallocate  ( &((*self)->translation) ) ;
        Memory_Deallocate   ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Test for equality of contents.
!-----------------------------------------------------------------------------*/
Boolean Transformation3_IsEqual ( const Transformation3 *self, const Transformation3 *other )
{
    Boolean QEQUAL = False ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        QEQUAL = ( Matrix33_IsEqual ( self->rotation, other->rotation ) && Vector3_IsEqual ( self->translation, other->translation ) ) ;
    }
    return QEQUAL ;
}

/*------------------------------------------------------------------------------
! . Test for identity.
!-----------------------------------------------------------------------------*/
Boolean Transformation3_IsIdentity ( const Transformation3 *self )
{
    Boolean QIDENTITY = False ;
    if ( self != NULL )
    {
        QIDENTITY = ( Matrix33_IsIdentity ( self->rotation ) && Vector3_IsNull ( self->translation ) ) ;
    }
    return QIDENTITY ;
}

/*------------------------------------------------------------------------------
! . Orthogonalization (transform to M * R * M^-1 + M t).
!-----------------------------------------------------------------------------*/
void Transformation3_Orthogonalize ( Transformation3 *self, const Matrix33 *A, const Matrix33 *B )
{
    if ( ( self != NULL ) && ( A != NULL ) && ( B != NULL ) )
    {
        Matrix33_PreMultiplyBy  ( self->rotation, A ) ;
        Matrix33_PostMultiplyBy ( self->rotation, B ) ;
        Matrix33_ApplyToVector3 ( A, self->translation ) ;
    }
}

