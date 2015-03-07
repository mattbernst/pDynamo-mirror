/*------------------------------------------------------------------------------
! . File      : PairwiseInteraction.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _PAIRWISEINTERACTION
# define _PAIRWISEINTERACTION

# include "Boolean.h"
# include "Coordinates3.h"
# include "CubicSpline.h"
# include "Integer.h"
# include "LJParameterContainer.h"
# include "Macros.h"
# include "MMAtomContainer.h"
# include "PairList.h"
# include "QCAtomContainer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The ABFS interaction type. */
typedef struct {
    Boolean      useAnalyticForm     ;
    Integer      splinePointDensity  ;
    Real         dampingCutoff       ;
    Real         innerCutoff         ;
    Real         outerCutoff         ;
    CubicSpline *electrostaticSpline ;
    CubicSpline *lennardJonesASpline ;
    CubicSpline *lennardJonesBSpline ;
} PairwiseInteractionABFS ;

/* . The ABFS interaction factors type for analytic evaluation. */
typedef struct {
    /* . Cutoff factors. */
    Real r2Damp   ;
    Real r2Off    ;
    Real r2On     ;
    /* . Electrostatic factors. */
    Real a        ;
    Real b        ;
    Real c        ;
    Real d        ;
    Real qShift1  ;
    Real qShift2  ;
    Real qF0      ;
    Real qAlpha   ;
    /* . Lennard-Jones A factors. */
    Real aF6      ;
    Real aK12     ;
    Real aShift12 ;
    Real aF0      ;
    Real aAlpha   ;
    /* . Lennard-Jones B factors. */
    Real bF3      ;
    Real bK6      ;
    Real bShift6  ;
    Real bF0      ;
    Real bAlpha   ;
} PairwiseInteractionABFSFactors ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . ABFS check distances. */
# define PairwiseInteractionABFS_CheckDistances( self, r2, s, s2 ) \
    { \
        if ( r2 > self.r2Off  ) continue ; \
        else if ( r2 < self.r2Damp ) { s2 = 0.0e+00 ; s = 0.0e+00 ; } \
        else { s2 = 1.0e+00 / r2 ; s = sqrt ( s2 ) ; } \
    }

/* . ABFS electrostatic interaction. */
# define PairwiseInteractionABFS_ElectrostaticTerm( self, r2, s, qij, f, dF ) \
    { \
        if ( r2 > self.r2On ) \
        { \
            f   = qij * ( s * ( self.a - r2 * ( self.b + r2 * ( self.c + self.d * r2 ) ) ) + self.qShift2 ) ; \
            dF += - qij * 0.5e+00 * s * ( self.a + r2 * ( self.b + r2 * ( 3.0e+00 * self.c + 5.0e+00 * self.d * r2 ) ) ) / r2 ; \
        } \
        else if ( r2 > self.r2Damp ) \
        { \
            f   = qij * ( s + self.qShift1 ) ; \
            dF += - qij * 0.5e+00 * s / r2 ; \
        } \
        else \
        { \
            f   = qij * ( self.qF0 - self.qAlpha * r2 ) ; \
            dF += - qij * self.qAlpha ; \
        } \
    }

/* . ABFS Lennard-Jones interactions - A and B. */
# define PairwiseInteractionABFS_LennardJonesTerm( self, r2, s, s2, Aij, Bij, f, dF ) \
    { \
        auto Real s6 = s2 * s2 * s2 ; \
        if ( r2 > self.r2On ) \
        { \
            auto Real l1 = s6 - self.aF6, l2 = ( s / r2 ) - self.bF3 ; \
            f   = Aij * self.aK12 * pow ( l1, 2 ) - Bij * self.bK6 * pow ( l2, 2 ) ; \
            dF += - 3.0e+00 * s6 * ( 2.0e+00 * Aij * self.aK12 * l1 / r2 - Bij * self.bK6 * l2 / s ) ; \
        } \
        else if ( r2 > self.r2Damp ) \
        { \
            f   = Aij * ( s6 * s6 - self.aShift12 ) - Bij * ( s6 - self.bShift6 ) ; \
            dF += - 3.0e+00 * s6 * ( 2.0e+00 * Aij * s6 - Bij ) / r2 ; \
        } \
        else \
        { \
            f   = Aij * ( self.aF0 - self.aAlpha * r2 ) - Bij * ( self.bF0 - self.bAlpha * r2 ) ; \
            dF += - Aij * self.aAlpha + Bij * self.bAlpha ; \
        } \
    }

/* . ABFS Lennard-Jones A interaction. */
# define PairwiseInteractionABFS_LennardJonesATerm( self, r2, s, s2, Aij, f, dF ) \
    { \
        auto Real s6 = s2 * s2 * s2 ; \
        if ( r2 > self.r2On ) \
        { \
            auto Real l1 = s6 - self.aF6 ; \
            f   =   Aij * self.aK12 * pow ( l1, 2 ) ; \
            dF += - Aij * self.aK12 * 6.0e+00 * s6 * l1 / r2 ; \
        } \
        else if ( r2 > self.r2Damp ) \
        { \
            f   = Aij * ( s6 * s6 - self.aShift12 ) ; \
            dF += - 6.0e+00 * s6 * Aij * s6 / r2 ; \
        } \
        else \
        { \
            f   = Aij * ( self.aF0 - self.aAlpha * r2 ) ; \
            dF += - Aij * self.aAlpha ; \
        } \
    }

/* . ABFS Lennard-Jones B interaction. */
# define PairwiseInteractionABFS_LennardJonesBTerm( self, r2, s, s2, Bij, f, dF ) \
    { \
        auto Real s6 = s2 * s2 * s2 ; \
        if ( r2 > self.r2On ) \
        { \
            auto Real l2 = ( s / r2 ) - self.bF3 ; \
            f   = - Bij * self.bK6 * pow ( l2, 2 ) ; \
            dF +=   Bij * self.bK6 * 3.0e+00 * s6 * l2 / s ; \
        } \
        else if ( r2 > self.r2Damp ) \
        { \
            f   = - Bij * ( s6 - self.bShift6 ) ; \
            dF +=   Bij * 3.0e+00 * s6 / r2 ; \
        } \
        else \
        { \
            f   = - Bij * ( self.bF0 - self.bAlpha * r2 ) ; \
            dF +=   Bij * self.bAlpha ; \
        } \
    }

/* . Calculate the number of spline points. */
# define PairwiseInteractionABFS_NumberOfSplinePoints( self, numberOfPoints ) \
    { numberOfPoints = Maximum ( Round ( ( Real ) self->splinePointDensity * self->outerCutoff ) + 1, 2 ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Pairwise interaction - basic procedures. */
extern PairwiseInteractionABFS *PairwiseInteractionABFS_Allocate    ( void ) ;
extern PairwiseInteractionABFS *PairwiseInteractionABFS_Clone       ( const PairwiseInteractionABFS  *self, Status *status ) ;
extern void                     PairwiseInteractionABFS_Deallocate  (       PairwiseInteractionABFS **self ) ;
extern void                     PairwiseInteractionABFS_MakeFactors ( const PairwiseInteractionABFS  *self, PairwiseInteractionABFSFactors *factors ) ;

/* . Spline procedures. */
extern CubicSpline *PairwiseInteractionABFS_MakeElectrostaticSpline ( const PairwiseInteractionABFS  *self, const Boolean useAtomicUnits, Status *status ) ;
extern CubicSpline *PairwiseInteractionABFS_MakeLennardJonesASpline ( const PairwiseInteractionABFS  *self, Status *status ) ;
extern CubicSpline *PairwiseInteractionABFS_MakeLennardJonesBSpline ( const PairwiseInteractionABFS  *self, Status *status ) ;

/* . Pairwise interaction - interaction procedures. */
extern void PairwiseInteractionABFS_MMMMEnergy     ( const PairwiseInteractionABFS *self               ,
                                                     const MMAtomContainer         *mmAtoms            ,
                                                     const LJParameterContainer    *ljParameters       ,
                                                           PairList                *pairList           ,
                                                     const Real                     electrostaticScale ,
                                                     const Real                     lennardJonesScale  ,
                                                     const Coordinates3            *crd1               ,
                                                     const Coordinates3            *crd2               ,
                                                           Real                    *eElectrostatic     ,
                                                           Real                    *eLennardJones      ,
# ifdef USEOPENMP
                                                     const Integer                  numberOfThreads    ,
                                                           Coordinates3           **grds1              ,
                                                           Coordinates3           **grds2              ) ;
# else
                                                           Coordinates3            *grd1               ,
                                                           Coordinates3            *grd2               ) ;
# endif
extern void PairwiseInteractionABFS_QCMMGradients  ( const PairwiseInteractionABFS *self               ,
                                                     const Real1DArray             *mmCharges          ,
                                                     const QCAtomContainer         *qcAtoms            ,
                                                           PairList                *pairList           ,
                                                     const Real                     electrostaticScale ,
                                                     const Coordinates3            *crd1               ,
                                                     const Coordinates3            *crd2               ,
                                                     const Real1DArray             *qcCharges          ,
                                                           Coordinates3            *grd1               ,
                                                           Coordinates3            *grd2               ) ;
extern void PairwiseInteractionABFS_QCMMPotentials ( const PairwiseInteractionABFS *self               ,
                                                     const Real1DArray             *mmCharges          ,
                                                     const QCAtomContainer         *qcAtoms            ,
                                                           PairList                *pairList           ,
                                                     const Real                     electrostaticScale ,
                                                     const Coordinates3            *crd1               ,
                                                     const Coordinates3            *crd2               ,
                                                     const Real1DArray             *potentials         ) ;
extern void PairwiseInteractionABFS_QCQCGradients  ( const PairwiseInteractionABFS *self               ,
                                                           PairList                *pairList           ,
                                                     const Real                     electrostaticScale ,
                                                     const Coordinates3            *crd1               ,
                                                     const Coordinates3            *crd2               ,
                                                     const Real1DArray             *qcCharges          ,
                                                           Coordinates3            *grd1               ,
                                                           Coordinates3            *grd2               ) ;
extern void PairwiseInteractionABFS_QCQCPotentials ( const PairwiseInteractionABFS *self               ,
                                                           PairList                *pairList           ,
                                                     const Real                     electrostaticScale ,
                                                     const Coordinates3            *crd1               ,
                                                     const Coordinates3            *crd2               ,
                                                           SymmetricMatrix         *potentials         ) ;

# endif
