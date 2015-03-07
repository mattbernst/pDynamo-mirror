/*------------------------------------------------------------------------------
! . File      : QCAtomContainer.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _QCATOMCONTAINER
# define _QCATOMCONTAINER

# include "Boolean.h"
# include "Coordinates3.h"
# include "GridFunctionDataBlock.h"
# include "Integer.h"
# include "Integer1DArray.h"
# include "MMAtomContainer.h"
# include "PairList.h"
# include "QCParameters.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Selection.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The qcatom type. */
typedef struct {
    Boolean    QBOUNDARY    ;
    Integer    atomicNumber ;
    Integer    center       ;
    Integer    index        ;
    Integer    qcpartner    ;
    Real       linkfactor   ;
    Real       widthe       ;
    Real       widthn       ;
    /* . Basis function indices. */
    /* . Proper representation. */
    Integer    dstart       ;
    Integer    fstart       ;
    Integer    ostart       ;
    Integer    pstart       ;
    Integer    ndbasis      ;
    Integer    nfbasis      ;
    Integer    nobasis      ;
    Integer    npbasis      ;
    /* . Working representation. */
    Integer    dstartw      ;
    Integer    fstartw      ;
    Integer    ostartw      ;
    Integer    pstartw      ;
    Integer    ndbasisw     ;
    Integer    nfbasisw     ;
    Integer    nobasisw     ;
    Integer    npbasisw     ;
    /* . MM partners. */
    Integer    nmmpartners  ;
    Integer   *mmpartners   ;
} QCAtom ;

/* . The qcatomcontainer type. */
typedef struct {
    Boolean    QLINKRATIO     ;
    Boolean    QTOSPHERICAL   ; /* . Flag to indicate whether working in Cartesians or sphericals - if one of the basis sets is spherical harmonic. */
    Integer    natoms         ;
    Integer    nboundary      ;
    Integer    nuclearCharge  ; /* . The integral nuclear charge (equivalent to the number of electrons in the neutral system). */
    Real       energybaseline ;
    /* . Basis function counters.*/
    /* . Proper representation.*/
    Integer    ndbasis        ;
    Integer    nfbasis        ;
    Integer    nobasis        ;
    Integer    npbasis        ;
    /* . Working representation.*/
    Integer    ndbasisw       ;
    Integer    nfbasisw       ;
    Integer    nobasisw       ;
    Integer    npbasisw       ;
    /* . QC atom data. */
    QCAtom    *data           ;
    /* . Representations. */
    Selection *baselection    ;
    Selection *mmpselection   ;
} QCAtomContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern QCAtomContainer *QCAtomContainer_Allocate                           ( const Integer natoms ) ;
extern QCAtomContainer *QCAtomContainer_Clone                              ( const QCAtomContainer  *self ) ;
extern void             QCAtomContainer_Deallocate                         (       QCAtomContainer **self ) ;
extern void             QCAtomContainer_GetAtomCoordinates3                ( const QCAtomContainer *self, const Integer q, const Coordinates3 *coordinates3, Real *x, Real *y, Real *z ) ;
extern Status           QCAtomContainer_GetCoordinates3                    ( const QCAtomContainer *self, const Coordinates3 *coordinates3, const Boolean QBOHR, Coordinates3 **qccoordinates3 ) ;
extern Real1DArray     *QCAtomContainer_GetMMCharges                       (       QCAtomContainer *self, const MMAtomContainer *mmAtoms, const Boolean QRCCOUPLING ) ;
extern void             QCAtomContainer_GetMMCoordinates3                  ( const QCAtomContainer *self, const Coordinates3 *coordinates3, Coordinates3 *mmcoordinates3 ) ;
extern void             QCAtomContainer_GetNuclearCharges                  ( const QCAtomContainer *self, const QCParameter *qcParameters, Real1DArray *qcCharges ) ;
extern void             QCAtomContainer_MakeBoundaryAtomSelection          (       QCAtomContainer  *self ) ;
extern Selection       *QCAtomContainer_MakeFullSelection                  ( const QCAtomContainer  *self ) ;
extern PairList        *QCAtomContainer_MakeMMPairList                     (       QCAtomContainer  *self, const MMAtomContainer *mmAtoms, PairList *pairlist ) ;
extern void             QCAtomContainer_MakePureMMSelection                (       QCAtomContainer  *self, const Integer upperBound ) ;
extern Selection       *QCAtomContainer_MakePureSelection                  ( const QCAtomContainer  *self ) ;
extern QCAtomContainer *QCAtomContainer_Merge                              ( const QCAtomContainer  *self, const QCAtomContainer *other, const Integer atomtypeincrement ) ;
extern Integer1DArray  *QCAtomContainer_OrbitalBasisAtomIndices            ( const QCAtomContainer  *self, Status *status ) ;
extern void QCAtomContainer_OrbitalBasisGridPointValues                    ( const QCAtomContainer  *self,
                                                                             const QCParameter           *qcParameters     ,
                                                                             const Coordinates3          *qcCoordinates3   ,
                                                                             const Coordinates3          *gridCoordinates3 ,
                                                                             const Boolean                resize           ,
                                                                             const Real                  *tolerance        ,
                                                                                   GridFunctionDataBlock *data             ,
                                                                                   Status                *status           ) ;
extern QCAtomContainer *QCAtomContainer_Prune                              ( const QCAtomContainer  *self, Selection *selection ) ;
extern void             QCAtomContainer_SetAtomGradients3                  ( const QCAtomContainer *self, const Integer q, const Coordinates3 *coordinates3,
                                                                             const Real gx, const Real gy, const Real gz, Coordinates3 *gradients3 ) ;
extern Status           QCAtomContainer_SetGradients3                      ( const QCAtomContainer *self, const Coordinates3 *coordinates3, const Coordinates3 *qcgradients3, const Boolean QCONVERT, Coordinates3 **gradients3 ) ;
extern void             QCAtomContainer_SetMMGradients3                    ( const QCAtomContainer *self, const Coordinates3 *mmgradients3, Coordinates3 *gradients3 ) ;
extern Integer          QCAtomContainer_Size                               ( const QCAtomContainer  *self ) ;

# endif
