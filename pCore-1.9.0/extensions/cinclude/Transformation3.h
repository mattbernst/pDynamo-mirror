/*------------------------------------------------------------------------------
! . File      : Transformation3.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _TRANSFORMATION3
# define _TRANSFORMATION3

# include "Definitions.h"
# include "Matrix33.h"
# include "Vector3.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
typedef struct {
   Matrix33 *rotation    ;
   Vector3  *translation ;
} Transformation3 ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern Transformation3 *Transformation3_Allocate      ( void ) ;
extern Transformation3 *Transformation3_Clone         ( const Transformation3  *self ) ;
extern void             Transformation3_Copy          (       Transformation3  *self, const Transformation3 *other ) ;
extern void             Transformation3_Deallocate    (       Transformation3 **self ) ;
extern Boolean          Transformation3_IsEqual       ( const Transformation3  *self, const Transformation3 *other ) ;
extern Boolean          Transformation3_IsIdentity    ( const Transformation3  *self ) ;
extern void             Transformation3_Orthogonalize (       Transformation3  *self, const Matrix33 *A, const Matrix33 *B ) ;

# endif
