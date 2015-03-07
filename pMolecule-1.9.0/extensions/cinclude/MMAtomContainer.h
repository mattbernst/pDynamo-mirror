/*------------------------------------------------------------------------------
! . File      : MMAtomContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _MMATOMCONTAINER
# define _MMATOMCONTAINER

# include "Boolean.h"
# include "Coordinates3.h"
# include "Definitions.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Selection.h"
# include "Status.h"
# include "Vector3.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The atom type. */
typedef struct {
    Boolean QACTIVE  ;
    Integer atomtype ;
    Integer ljtype   ;
    Real    charge   ;
} MMAtom ;

/* . The atom container type. */
typedef struct {
   Integer natoms ;
   MMAtom *data   ;
} MMAtomContainer ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern MMAtomContainer *MMAtomContainer_Allocate            ( const Integer natoms ) ;
extern MMAtomContainer *MMAtomContainer_Clone               ( const MMAtomContainer  *self ) ;
extern void             MMAtomContainer_Deallocate          (	    MMAtomContainer **self ) ;
extern Vector3         *MMAtomContainer_DipoleMoment        ( const MMAtomContainer  *self, const Coordinates3 *coordinates3,
                                                                                                    const Vector3 *center ) ;
extern Real1DArray     *MMAtomContainer_GetCharges          ( const MMAtomContainer  *self, const Boolean activeOnly, Status *status ) ;
extern MMAtomContainer *MMAtomContainer_Merge               ( const MMAtomContainer  *self, const MMAtomContainer *other,
                                                                                         const Integer atomtypeincrement,
						            		                const Integer ljtypeincrement ) ;
extern Integer          MMAtomContainer_NumberOfActiveAtoms ( const MMAtomContainer  *self ) ;
extern MMAtomContainer *MMAtomContainer_Prune               ( const MMAtomContainer  *self, Selection *selection ) ;
extern Integer          MMAtomContainer_Size                ( const MMAtomContainer  *self ) ;
extern Real             MMAtomContainer_TotalCharge         ( const MMAtomContainer  *self ) ;

# endif
