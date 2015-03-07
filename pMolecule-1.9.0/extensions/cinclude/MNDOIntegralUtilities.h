/*------------------------------------------------------------------------------
! . File      : MNDOIntegralUtilities.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _MNDOINTEGRALUTILITIES
# define _MNDOINTEGRALUTILITIES

# include "BlockStorage.h"
# include "Boolean.h"
# include "Integer.h"
# include "MNDOParameters.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Type definition for function pointer. */
typedef Real ( * ChargeInteractionFunction ) ( const Real r, const Integer l1, const Integer l2, const Integer m, const Real da, const Real db, const Real add ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real MNDOIntegralUtilities_CoreCore                  ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real rij ) ;
extern Real MNDOIntegralUtilities_CoreCoreD                 ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real rij ) ;
extern void MNDOIntegralUtilities_CoreCoreFD                ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real rij, Real *fCore, Real *dCore ) ;
extern Real MNDOIntegralUtilities_GetTransformationMatrices ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real *coordi, const Real *coordj, Real2DArray **itransformation, Real2DArray **jtransformation,
                                                                                                                Real2DArray **itransformationX, Real2DArray **itransformationY, Real2DArray **itransformationZ,
                                                                                                                Real2DArray **jtransformationX, Real2DArray **jtransformationY, Real2DArray **jtransformationZ,
                                                                                                                                                                                Real *x, Real *y, Real *z ) ;
extern void MNDOIntegralUtilities_LocalFrame2CTEIs          ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real r,
                                                                          Real2DArray  *lfteis, Real1DArray  *core1b, Real1DArray  *core2a,
                                                                          Real2DArray *dlfteis, Real1DArray *dcore1b, Real1DArray *dcore2a ) ;
extern void MNDOIntegralUtilities_LocalFrame2COEIsSP        ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real r, const Boolean swapped, Real1DArray *core, Real1DArray *dcore ) ;
extern void MNDOIntegralUtilities_LocalFrame2CTEIsSP        ( const MNDOParameters *idata, const MNDOParameters *jdata, const Real r, Real2DArray *lfteis, Real2DArray *dlfteis ) ;
# ifdef MNDODORBITALS
extern Real MNDOIntegralUtilities_LocalFrame2CTEI           ( const ChargeInteractionFunction Evaluate, const MNDOParameters *idata, const MNDOParameters *jdata, const Integer ij, const Integer kl,
                                                                                                                  const Integer i, const Integer j, const Integer k, const Integer l, const Integer c, const Real r ) ;
# endif

# endif
