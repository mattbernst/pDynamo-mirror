/*------------------------------------------------------------------------------
! . File      : Transformation3Container.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _TRANSFORMATION3CONTAINER
# define _TRANSFORMATION3CONTAINER

# include "Definitions.h"
# include "Transformation3.h"
# include "Vector3.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The container type. */
typedef struct {
    Boolean           QOWNER   ;
    Integer           identity ;
    Integer           nitems   ;
    Integer          *inverses ;
    Transformation3 **items    ;
} Transformation3Container ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern Transformation3Container *Transformation3Container_Allocate                      ( const Integer nitems ) ;
extern void                      Transformation3Container_Deallocate                    (       Transformation3Container **self ) ;
extern void                      Transformation3Container_FindIdentity                  (       Transformation3Container  *self ) ;
extern void                      Transformation3Container_FindInverseIntegerTranslation ( const Transformation3Container  *self, const Integer t,
                                                                                               const Integer a, const Integer b, const Integer c,
                                                                                                                            Vector3 *translation,
                                                                                          Integer *ainverse, Integer *binverse, Integer *cinverse ) ;
extern void                      Transformation3Container_FindInverses                  (       Transformation3Container  *self ) ;

# endif
