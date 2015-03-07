/*------------------------------------------------------------------------------
! . File      : AntisymmetricMatrix.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _ANTISYMMETRICMATRIX
# define _ANTISYMMETRICMATRIX

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The array type. */
typedef struct {
    Boolean  isOwner ;
    Boolean  isView  ;
    Integer  length  ;
    Integer  offset  ;
    Integer  size    ;
    Integer  stride  ;
    Real    *data    ;
    Real1DArray view ;
} AntisymmetricMatrix ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define AntisymmetricMatrix_Data( self ) ( &((self)->data[(self)->offset]) )

/* . An item (i > j). */
# define AntisymmetricMatrix_Item( self, i, j ) ( (self)->data[(self)->offset+(self)->stride*(((i)*(i-1))/2+j)] )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real                 AntisymmetricMatrix_AbsoluteMaximum        ( const AntisymmetricMatrix  *self ) ;
extern void                 AntisymmetricMatrix_AddScaledMatrix        (       AntisymmetricMatrix  *self, const Real value, const AntisymmetricMatrix *other, Status *status ) ;
extern AntisymmetricMatrix *AntisymmetricMatrix_Allocate               ( const Integer length, Status *status ) ;
extern void                 AntisymmetricMatrix_AnticommutatorAS       ( const AntisymmetricMatrix  *self, const SymmetricMatrix *a, AntisymmetricMatrix *result, Status *status ) ;
extern void                 AntisymmetricMatrix_CommutatorAS           ( const AntisymmetricMatrix  *self, const SymmetricMatrix *a,     SymmetricMatrix *result, Status *status ) ;
extern void                 AntisymmetricMatrix_CommutatorSS_Fast      (       AntisymmetricMatrix  *self   ,
                                                                         const SymmetricMatrix      *a      ,
                                                                         const SymmetricMatrix      *b      ,
                                                                               Real2DArray          *mA     ,
                                                                               Real2DArray          *mB     ,
                                                                               Real2DArray          *mC     ,
                                                                               Status               *status ) ;
extern void                 AntisymmetricMatrix_CommutatorSS_Reference (       AntisymmetricMatrix  *self, const SymmetricMatrix *a, const SymmetricMatrix *b, Status *status ) ;
extern void                 AntisymmetricMatrix_CommutatorSSS          (       AntisymmetricMatrix  *self, const SymmetricMatrix *a, const SymmetricMatrix *b, const SymmetricMatrix *c, Status *status ) ;
extern void                 AntisymmetricMatrix_CommutatorTSSST        (       AntisymmetricMatrix  *self       ,
                                                                         const SymmetricMatrix      *a          ,
                                                                         const SymmetricMatrix      *b          ,
                                                                         const SymmetricMatrix      *c          ,
                                                                         const Real2DArray          *m          ,
                                                                         const Boolean               mTranspose ,
                                                                               Real2DArray          *u          ,
                                                                               Real2DArray          *v          ,
                                                                               Real2DArray          *w          ,
                                                                               Status               *status     ) ;
extern void                 AntisymmetricMatrix_CopyFromReal2DArray    (       AntisymmetricMatrix  *self, const Real2DArray *other, const Boolean scale, Status *status ) ;
extern void                 AntisymmetricMatrix_CopyTo                 ( const AntisymmetricMatrix  *self, AntisymmetricMatrix *other, Status *status ) ;
extern void                 AntisymmetricMatrix_CopyToReal2DArray      ( const AntisymmetricMatrix  *self, Real2DArray *other, Status *status ) ;
extern void                 AntisymmetricMatrix_Deallocate             (       AntisymmetricMatrix **self ) ;
extern void                 AntisymmetricMatrix_GetColumn              ( const AntisymmetricMatrix  *self, const Integer n, Real1DArray *column, Status *status ) ;
extern Real                 AntisymmetricMatrix_GetItem                ( const AntisymmetricMatrix  *self, const Integer i, const Integer j, Status *status ) ;
extern void                 AntisymmetricMatrix_GetItemIndexAndSign    ( const AntisymmetricMatrix  *self, const Integer i, const Integer j, Integer *index, Real *sign, Status *status ) ;
extern void                 AntisymmetricMatrix_Print                  ( const AntisymmetricMatrix  *self ) ;
extern void                 AntisymmetricMatrix_Scale                  (       AntisymmetricMatrix  *self, const Real value ) ;
extern void                 AntisymmetricMatrix_Set                    (       AntisymmetricMatrix  *self, const Real value ) ;
extern void                 AntisymmetricMatrix_SetItem                (       AntisymmetricMatrix  *self, const Integer i, const Integer j, const Real value, Status *status ) ;
extern void                 AntisymmetricMatrix_SymmetricTransform     ( const AntisymmetricMatrix  *self, const SymmetricMatrix *matrix, AntisymmetricMatrix *result, Status *status ) ;
extern Real                 AntisymmetricMatrix_Trace2                 ( const AntisymmetricMatrix  *self, const AntisymmetricMatrix *other, Status *status ) ;
extern void                 AntisymmetricMatrix_Transform              ( const AntisymmetricMatrix  *self, const Real2DArray *matrix, const Boolean useTranspose, AntisymmetricMatrix *result, Status *status ) ;
extern void                 AntisymmetricMatrix_Transpose              (       AntisymmetricMatrix  *self ) ;
extern void                 AntisymmetricMatrix_ViewOfRaw              (       AntisymmetricMatrix  *self, const Integer length, const Integer stride, Real *data, const Integer rawstride, Status *status ) ;

# endif
