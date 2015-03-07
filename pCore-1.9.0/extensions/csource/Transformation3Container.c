/*------------------------------------------------------------------------------
! . File      : Transformation3Container.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include <math.h>
# include <stdio.h>

# include "Transformation3Container.h"
# include "Memory.h"

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
Transformation3Container *Transformation3Container_Allocate ( const int nitems )
{
    Transformation3Container *self = NULL ;
    if ( nitems >= 0 )
    {
        auto int i ;
        self           = ( Transformation3Container * ) Memory_Allocate ( sizeof ( Transformation3Container ) ) ;
        self->identity =   -1 ;
        self->inverses = NULL ;
	self->items    = ( Transformation3 ** ) Memory_Allocate_Array ( nitems, sizeof ( Transformation3 * ) ) ;
	self->nitems   = nitems ;
        self->QOWNER   = False  ;
	for ( i = 0 ; i < nitems ; i++ ) self->items[i] = NULL ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void Transformation3Container_Deallocate ( Transformation3Container **self )
{
    if ( (*self) != NULL )
    {
        if ( (*self)->QOWNER )
        {
            auto int i ;
	    for ( i = 0 ; i < (*self)->nitems ; i++ ) Transformation3_Deallocate ( &((*self)->items[i]) ) ;
        }
        free ( (*self)->inverses ) ;
        free ( (*self)->items    ) ;
        free ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*------------------------------------------------------------------------------
! . Find the identity.
!-----------------------------------------------------------------------------*/
void Transformation3Container_FindIdentity ( Transformation3Container *self )
{
    if ( ( self != NULL ) && ( self->nitems > 0 ) )
    {
        auto int i ;
        for ( i = 0 ; i < self->nitems ; i++ )
        {
            if ( Transformation3_IsIdentity ( self->items[i] ) )
            {
                self->identity = i ;
                break ;
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Find the inverse translation for a transformation (in terms of a, b, c).
!-----------------------------------------------------------------------------*/
# define BIGNUM    -999999
# define TOLERANCE  1.0e-4

void Transformation3Container_FindInverseIntegerTranslation ( const Transformation3Container *self, const int t, const int a, const int b, const int c,
                                                                                    Vector3 *translation, int *ainverse, int *binverse, int *cinverse )
{
    (*ainverse) = BIGNUM ;
    (*binverse) = BIGNUM ;
    (*cinverse) = BIGNUM ;
    if ( ( self != NULL ) && ( translation != NULL ) && ( t >= 0 ) && ( self->inverses[t] >= 0 ) )
    {
        auto int i, n[3], tinverse, v ;
        tinverse = self->inverses[t] ;
        /* . Find the inverse translation minus the inverse's translation. */
        Vector3_CopyTo ( self->items[t]->translation, translation, NULL ) ;
        translation->data[0] += ( double ) a ;
        translation->data[1] += ( double ) b ;
        translation->data[2] += ( double ) c ;
        Matrix33_ApplyToVector3 ( self->items[tinverse]->rotation, translation ) ;
        Vector3_Scale ( translation, -1.0e+00 ) ;
        Vector3_AddScaledVector ( translation, -1.0e+00, self->items[tinverse]->translation, NULL ) ;
        /* . Convert the result to integers. */
        for ( i = 0 ; i < 3 ; i++ )
        {
            v = Round ( translation->data[i] ) ;
            if ( fabs ( translation->data[i] - ( double ) v ) < TOLERANCE ) n[i] = v      ;
            else                                                            n[i] = BIGNUM ;
if ( v != n[i] )
printf ( "i, t, tinverse, v, n[i], translation->data[i] = %d %d %d %d %d %f\n", i, t, tinverse, v, n[i], translation->data[i] ) ;
        }
        (*ainverse) = n[0] ;
        (*binverse) = n[1] ;
        (*cinverse) = n[2] ;
    }
}

# undef BIGNUM
# undef TOLERANCE

/*------------------------------------------------------------------------------
! . Find the inverses - only the rotational part!
! . This should probably be changed (have both rotational and full inverses).
!-----------------------------------------------------------------------------*/
void Transformation3Container_FindInverses ( Transformation3Container *self )
{
    /* . Initialization. */
    free ( self->inverses ) ;
    self->inverses = NULL   ;
    /* . Find the inverses. */
    if ( ( self != NULL ) && ( self->nitems > 0 ) )
    {
        auto int            i, j ;
        auto Matrix33 *m    ;
        m = Matrix33_Allocate   ( ) ;
        Matrix33_Set ( m, 0.0e+00 ) ;
        self->inverses = Memory_Allocate_Array_Integer_Initialize ( self->nitems, -1 ) ;
        for ( i = 0 ; i < self->nitems ; i++ )
        {
            if ( self->inverses[i] < 0 )
            {
                Matrix33_Invert ( m, self->items[i]->rotation ) ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    if ( self->inverses[j] < 0 )
                    {
                        if ( Matrix33_IsEqual ( m, self->items[j]->rotation ) )
                        {
                            self->inverses[i] = j ;
                            self->inverses[j] = i ;
                            break ;
                        }
                    }
                }
            }
        }
        Matrix33_Deallocate ( &m ) ;
    }
}
