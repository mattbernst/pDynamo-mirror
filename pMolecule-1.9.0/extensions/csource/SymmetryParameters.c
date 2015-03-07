/*------------------------------------------------------------------------------
! . File      : SymmetryParameters.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/

# include <math.h>

# include "Memory.h"
# include "SymmetryParameters.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetryParameters *SymmetryParameters_Allocate ( void )
{
    SymmetryParameters *self = NULL ;
    self           = ( SymmetryParameters * ) Memory_Allocate ( sizeof ( SymmetryParameters ) ) ;
    self->QM       = False   ;
    self->a        = 0.0e+00 ;
    self->b        = 0.0e+00 ;
    self->c        = 0.0e+00 ;
    self->alpha    = 0.0e+00 ;
    self->beta     = 0.0e+00 ;
    self->gamma    = 0.0e+00 ;
    self->inverseM = NULL ;
    self->M        = NULL ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Center coordinates within the primary image by index.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SymmetryParameters_CenterCoordinates3ByIndex ( const SymmetryParameters *self, const Selection *selection, Coordinates3 *coordinates3 )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
        auto Coordinates3 *fractional ;
        auto Integer       n          ;
        /* . Allocate space. */
        if ( selection == NULL ) n = coordinates3->length0 ;
        else                     n = selection->nindices   ;
        fractional = Coordinates3_Allocate ( n ) ;
        if ( fractional != NULL )
        {
            auto Real x, y, z ;
            auto Integer    i ;
            /* . Get the fractional coordinates. */
            Coordinates3_Gather ( fractional, coordinates3  , selection ) ;
            Coordinates3_Rotate ( fractional, self->inverseM, NULL      ) ;
            /* . Shift the coordinates. */
            for ( i = 0 ; i < fractional->length0 ; i++ )
            {
                Coordinates3_GetRow ( fractional, i, x, y, z ) ;
                x -= floor ( x ) ;
                y -= floor ( y ) ;
                z -= floor ( z ) ;
                Coordinates3_SetRow ( fractional, i, x, y, z ) ;
            }
            /* . Back transform and copy back. */
            Coordinates3_Rotate  ( fractional, self->M     , NULL      ) ;
            Coordinates3_Scatter ( fractional, coordinates3, selection ) ;
            /* . Success. */
            status = Status_Success ;
        }
        else status = Status_OutOfMemory ;
        /* . Finish up. */
        Coordinates3_Deallocate ( &fractional ) ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Center coordinates within the primary image by isolate.
! . Only isolates with selected members are centered.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SymmetryParameters_CenterCoordinates3ByIsolate ( const SymmetryParameters *self, const SelectionContainer *isolates, Selection *selection, Coordinates3 *coordinates3 )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( isolates != NULL ) && ( coordinates3 != NULL ) )
    {
        auto Boolean QOK = True, QSELECTION ;
        auto Integer  i, s ;
        auto Selection *isolate ;
        auto Real1DArray    *center, *translation ;
        /* . Allocate space. */
        center      = Vector3_Allocate ( ) ;
        translation = Vector3_Allocate ( ) ;
        /* . Check for a selection. */
        QSELECTION = ( selection != NULL ) ;
        if ( QSELECTION ) Selection_MakeFlags ( selection, isolates->upperBound ) ;
        /* . Check memory. */
        if ( ( center != NULL ) && ( translation != NULL ) )
        {
            /* . Loop over isolates. */
            for ( i = 0 ; i < isolates->nitems ; i++, QOK = True )
            {
                isolate = isolates->items[i] ;
                /* . Check whether to include the isolate. */
                if ( QSELECTION )
                {
                    for ( s = 0 ; s < isolate->nindices ; s++ )
                    {
                        if ( ! selection->flags[isolate->indices[s]] ) { QOK = False ; break ; }
                    }
                }
                /* . Do the centering. */
                if ( QOK )
                {
                    Coordinates3_Center ( coordinates3, isolate, NULL, &center ) ;
                    SymmetryParameters_FindCenteringTranslation ( self, center, translation ) ;
                    Coordinates3_Translate ( coordinates3, translation, isolate ) ;
                }
            }
            /* . Success. */
            status = Status_Success ;
        }
        else status = Status_OutOfMemory ;
        /* . Finish up. */
        Vector3_Deallocate ( &center      ) ;
        Vector3_Deallocate ( &translation ) ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Clear the M/inverseM representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_ClearM ( SymmetryParameters *self )
{
    if ( self != NULL )
    {
        Matrix33_Deallocate ( &(self->inverseM) ) ;
        Matrix33_Deallocate ( &(self->M)        ) ;
        self->QM = False ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transfer data from one structure to another.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_CopyTo ( const SymmetryParameters  *self, SymmetryParameters *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        other->a     = self->a     ;
        other->b     = self->b     ;
        other->c     = self->c     ;
        other->alpha = self->alpha ;
        other->beta  = self->beta  ;
        other->gamma = self->gamma ;
        other->QM    = False       ;
        SymmetryParameters_MakeM ( other ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_Deallocate ( SymmetryParameters **self )
{
    if ( (*self) != NULL )
    {
        SymmetryParameters_ClearM ( (*self) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a displacement (given in terms of a, b, c which are the columns
! . of the matrix M).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_Displacement ( const SymmetryParameters *self, const Integer a, const Integer b, const Integer c, Vector3 *displacement )
{
    if ( ( self != NULL ) && ( displacement != NULL ) )
    {
        auto Real d, na, nb, nc ;
        auto Integer    i ;
        na = ( Real ) a ; nb = ( Real ) b ; nc = ( Real ) c ;
        for ( i = 0 ; i < 3 ; i++ )
        {
            d = na * Matrix33_Item ( self->M, i, 0 ) + nb * Matrix33_Item ( self->M, i, 1 ) + nc * Matrix33_Item ( self->M, i, 2 ) ;
            Vector3_Item ( displacement, i ) = d ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the range of a, b and c values for which an image box overlaps with
! . a central box. This procedure makes use of the special properties of M
! . and the boxes are assumed to be orthorhombic.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Limits method. */
static Boolean GetLimits ( const Real bl, const Real bu, Real il, Real iu, const Real t, Integer *low, Integer *high )
{
    Boolean QOK = False ;
    Integer  n   = 0     ;
    /* . Move left until iu < bl. */
    while ( iu >= bl ) { il -= t ; iu -= t ; n-- ; }
    /* . Move right until iu >= bl. */
    while ( iu <  bl ) { il += t ; iu += t ; n++ ; }
    /* . Check that there is some overlap. */
    if ( il <= bu ) /* . && iu >= bl. */
    {
        (*low)  = n ;
        /* . Move right until il > bu. */
        while ( il <= bu ) { il += t ; iu += t ; n++ ; }
        (*high) = n - 1 ;
        /* . Everything is OK. */
        QOK     = True ;
    }
    return QOK ;
}

/* . Box method. */
void SymmetryParameters_FindBoxSearchLimits ( const SymmetryParameters *self   ,
                                              const Vector3            *lower  ,
                                              const Vector3            *upper  ,
                                              const Vector3            *ilower ,
                                              const Vector3            *iupper ,
                                                    Integer            *alow   ,
                                                    Integer            *ahigh  ,
                                                    Integer            *blow   ,
                                                    Integer            *bhigh  ,
                                                    Integer            *clow   ,
                                                    Integer            *chigh  )
{
    /* . Initialization. */
    (*alow)  = (*blow)  = (*clow)  =  0 ;
    (*ahigh) = (*bhigh) = (*chigh) = -1 ;
    if ( ( self != NULL ) && ( lower != NULL ) && ( upper != NULL ) && ( ilower != NULL ) && ( iupper != NULL ) )
    {
        auto Boolean   QOK ;
        auto Real bl, bu, d1, d2 ;
        /* . Do c first as this is the only lattice vector that contributes to z. */
        QOK = GetLimits ( lower->data[2], upper->data[2], ilower->data[2], iupper->data[2], Matrix33_Item ( self->M, 2, 2 ), clow, chigh ) ;
        /* . Now do b which contributes to y (along with c). */
        if ( QOK )
        {
            d1 = (*clow)  * Matrix33_Item ( self->M, 1, 2 ) ;
            d2 = (*chigh) * Matrix33_Item ( self->M, 1, 2 ) ;
            bl = lower->data[1] - Maximum ( d1, d2 ) ;
            bu = upper->data[1] - Minimum ( d1, d2 ) ;
            QOK = GetLimits ( bl, bu, ilower->data[1], iupper->data[1], Matrix33_Item ( self->M, 1, 1 ), blow, bhigh ) ;
            /* . Now do a which contributes to x (along with b and c). */
            if ( QOK )
            {
                d1  = (*clow)  * Matrix33_Item ( self->M, 0, 2 ) ;
                d2  = (*chigh) * Matrix33_Item ( self->M, 0, 2 ) ;
                bl  = lower->data[0] - Maximum ( d1, d2 ) ;
                bu  = upper->data[0] - Minimum ( d1, d2 ) ;
                d1  = (*blow)  * Matrix33_Item ( self->M, 0, 1 ) ;
                d2  = (*bhigh) * Matrix33_Item ( self->M, 0, 1 ) ;
                bl -= Maximum ( d1, d2 ) ;
                bu -= Minimum ( d1, d2 ) ;
                QOK = GetLimits ( bl, bu, ilower->data[0], iupper->data[0], Matrix33_Item ( self->M, 0, 0 ), alow, ahigh ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the translation that puts the point (in real space) inside the
! . primary image (defined as between [0,1] in fractional coordinates).
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SymmetryParameters_FindCenteringTranslation ( const SymmetryParameters *self, const Vector3 *point, Vector3 *translation )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( point != NULL ) && ( translation != NULL ) )
    {
        auto Integer a, b, c ;
        /* . Find the point in fractional coordinates using translation as temporary storage. */
        Vector3_CopyTo ( point, translation, NULL ) ;
        Matrix33_ApplyToVector3 ( self->inverseM, translation ) ;
        /* . Find the integral parts of the translation. */
        a = ( Integer ) floor ( translation->data[0] ) ;
        b = ( Integer ) floor ( translation->data[1] ) ;
        c = ( Integer ) floor ( translation->data[2] ) ;
        /* . Find the displacement. */
        SymmetryParameters_Displacement ( self, -a, -b, -c, translation ) ;
        /* . Success. */
        status = Status_Success ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check to see if the minimum image convention is satisfied given a length.
! . THIS NEEDS TO BE GENERALIZED TO NON-ORTHOGONAL PARAMETERS.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SymmetryParameters_IsMinimumImageConventionSatisfied ( const SymmetryParameters *self, const Real length )
{
    Boolean QOK = False ;
    if ( self != NULL )
    {
        auto Real d ;
        d   = 2.0 * length ;
        QOK = ( d < self->a ) && ( d < self->b ) && ( d < self->c ) ;
    }
    return QOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check to see if the symmetry parameters are orthorhombic.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SymmetryParameters_IsOrthorhombic ( const SymmetryParameters *self )
{
    Boolean QOK = False ;
    if ( self != NULL ) QOK = ( self->alpha == 90.0e+00 ) && ( self->beta == 90.0e+00 ) && ( self->gamma == 90.0e+00 ) ;
    return QOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the symmetry parameters isotropically.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_IsotropicScale ( SymmetryParameters *self, const Real scale )
{
    if ( self != NULL )
    {
        self->a *= scale ;
        self->b *= scale ;
        self->c *= scale ;
        Matrix33_Scale ( self->M, scale ) ;
        if ( scale == 0.0e+00 ) Matrix33_Set   ( self->inverseM, 0.0e+00         ) ;
        else                    Matrix33_Scale ( self->inverseM, 1.0e+00 / scale ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply the minimum image convention to a vector (in the form of a vector3).
! . THIS NEEDS TO BE GENERALIZED TO NON-ORTHOGONAL PARAMETERS.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_MakeMinimumImageVector3 ( const SymmetryParameters *self, Vector3 *r, Vector3 *dr )
{
    if ( ( self != NULL ) && ( r != NULL ) )
    {
        auto Boolean   QDR ;
        auto Real dc  ;
        QDR = ( dr != NULL ) ;
        dc = - self->a * Round ( r->data[0] / self->a ) ; r->data[0] += dc ; if ( QDR ) dr->data[0] = dc ;
        dc = - self->b * Round ( r->data[1] / self->b ) ; r->data[1] += dc ; if ( QDR ) dr->data[1] = dc ;
        dc = - self->c * Round ( r->data[2] / self->c ) ; r->data[2] += dc ; if ( QDR ) dr->data[2] = dc ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply the minimum image convention to a vector (in the form of xyz).
! . THIS NEEDS TO BE GENERALIZED TO NON-ORTHOGONAL PARAMETERS.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_MakeMinimumImageXYZ ( const SymmetryParameters *self, Real *x, Real *y, Real *z, Real *dx, Real *dy, Real *dz )
{
    if ( ( self != NULL ) && ( x != NULL ) && ( y != NULL ) && ( z != NULL ) )
    {
        auto Real dc ;
        dc = - self->a * Round ( (*x) / self->a ) ; (*x) += dc ; if ( dx != NULL ) (*dx) = dc ;
        dc = - self->b * Round ( (*y) / self->b ) ; (*y) += dc ; if ( dy != NULL ) (*dy) = dc ;
        dc = - self->c * Round ( (*z) / self->c ) ; (*z) += dc ; if ( dz != NULL ) (*dz) = dc ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the M/inverseM representation.
! . The columns of the matrix M correspond to the lattice vectors a, b, c.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_MakeM ( SymmetryParameters *self )
{
    if ( ( self != NULL ) && ( ( ! self->QM ) || ( self->inverseM == NULL ) || ( self->M == NULL ) ) )
    {
        auto Real alpha, beta, cosalpha, cosbeta, cosgamma, gamma, singamma ;

        /* . Clear and reallocate only if necessary. */
        if ( ( self->inverseM == NULL ) || ( self->M == NULL ) )
        {
            SymmetryParameters_ClearM ( self ) ;
            self->inverseM = Matrix33_Allocate ( ) ;
            self->M        = Matrix33_Allocate ( ) ;
        }

        /* . Some factors. */
        alpha    = self->alpha * UNITS_ANGLE_DEGREES_TO_RADIANS ;
        beta     = self->beta  * UNITS_ANGLE_DEGREES_TO_RADIANS ;
        gamma    = self->gamma * UNITS_ANGLE_DEGREES_TO_RADIANS ;
        cosalpha = cos ( alpha ) ;
        cosbeta  = cos ( beta  ) ;
        cosgamma = cos ( gamma ) ;
        singamma = sin ( gamma ) ;

        /* . Create the M matrix - standard orientation. */
        Matrix33_Set  ( self->M, 0.0e+00 ) ;
        Matrix33_Item ( self->M, 0, 0 ) = self->a ;
        Matrix33_Item ( self->M, 0, 1 ) = self->b * cosgamma ;
        Matrix33_Item ( self->M, 1, 1 ) = self->b * singamma ;
        Matrix33_Item ( self->M, 0, 2 ) = self->c * cosbeta  ;
        Matrix33_Item ( self->M, 1, 2 ) = self->c * ( cosalpha - cosbeta * cosgamma ) / singamma ;
        Matrix33_Item ( self->M, 2, 2 ) = self->c * sqrt ( 1.0e+00 - cosalpha * cosalpha - cosbeta * cosbeta - cosgamma * cosgamma +
                                                           2.0e+00 * cosalpha * cosbeta * cosgamma ) / singamma ;

        /* . Invert it. */
        Matrix33_Invert ( self->inverseM, self->M ) ;

        /* . Finish up. */
        self->QM = True ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the orthorhombic widths corresponding to the symmetry parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_MakeOrthorhombicWidths ( const SymmetryParameters *self, Vector3 *widths )
{
    if ( ( self != NULL ) && ( widths != NULL ) )
    {
        if ( SymmetryParameters_IsOrthorhombic ( self ) )
        {
            Vector3_Item ( widths, 0 ) = self->a ;
            Vector3_Item ( widths, 1 ) = self->b ;
            Vector3_Item ( widths, 2 ) = self->c ;
        }
        else
        {
            auto Integer  i, j, x, y ;
            auto Vector3 *u, *v ;
            u = Vector3_Allocate ( ) ;
            v = Vector3_Allocate ( ) ;
            for ( i = 0 ; i < 3 ; i++ )
            {
                x = ( i + 1 ) % 3 ;
                y = ( i + 2 ) % 3 ;
                for ( j = 0 ; j < 3 ; j++ )
                {
                    Vector3_Item ( u, j ) = Matrix33_Item ( self->M, j, x ) ;
                    Vector3_Item ( v, j ) = Matrix33_Item ( self->M, j, y ) ;
                }
                Vector3_CrossProduct ( u, v ) ;
                for ( j = 0 ; j < 3 ; j++ ) { Vector3_Item ( v, j ) = Matrix33_Item ( self->M, j, i ) ; }
                Vector3_Item ( widths, i ) = fabs ( Vector3_Dot ( u, v, NULL ) ) / Vector3_Norm2 ( u ) ;
            }
            Vector3_Deallocate ( &u ) ;
            Vector3_Deallocate ( &v ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the parameters appropriate for a crystal.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_SetCrystalParameters (       SymmetryParameters *self  ,
                                               const Real                a     ,
                                               const Real                b     ,
                                               const Real                c     ,
                                               const Real                alpha ,
                                               const Real                beta  ,
                                               const Real                gamma )
{
    if ( self != NULL )
    {
        self->a     = a     ;
        self->b     = b     ;
        self->c     = c     ;
        self->alpha = alpha ;
        self->beta  = beta  ;
        self->gamma = gamma ;
        self->QM    = False ;
        SymmetryParameters_MakeM ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the volume which is the product of the diagonal elements of M.
!---------------------------------------------------------------------------------------------------------------------------------*/
double SymmetryParameters_Volume ( const SymmetryParameters *self )
{
    Real volume = 0.0e+00 ;
    if ( self != NULL )
    {
        volume = Matrix33_Item ( self->M, 0, 0 ) * Matrix33_Item ( self->M, 1, 1 ) * Matrix33_Item ( self->M, 2, 2 ) ;
    }
    return volume ;
}
