#-------------------------------------------------------------------------------
# . File      : pMolecule.QCAtomContainer.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions      cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3      cimport CCoordinates3, Coordinates3
from pCore.Memory            cimport Memory_Allocate_Array_Integer
from pCore.Selection         cimport Selection, CSelection
from pCore.Status            cimport Status
from pMolecule.GaussianBasis cimport CGaussianBasis
from pMolecule.QCParameters  cimport QCParameters, CQCParameter

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "QCAtomContainer.h":

    ctypedef struct CQCAtom "QCAtom":
        Boolean   QBOUNDARY
        Real linkfactor
        Real widthe
        Real widthn
        Integer    atomicNumber
        Integer    center
        Integer    index
        Integer    qcpartner
        Integer    dstart
        Integer    fstart
        Integer    ostart
        Integer    pstart
        Integer    ndbasis
        Integer    nfbasis
        Integer    nobasis
        Integer    npbasis
        Integer    dstartw
        Integer    fstartw
        Integer    ostartw
        Integer    pstartw
        Integer    ndbasisw
        Integer    nfbasisw
        Integer    nobasisw
        Integer    npbasisw
        Integer    nmmpartners
        Integer   *mmpartners

    ctypedef struct CQCAtomContainer "QCAtomContainer":
        Boolean   QLINKRATIO
        Boolean   QTOSPHERICAL
        Real energybaseline
        Integer    natoms
        Integer    nboundary
        Integer    nuclearCharge
        Integer    ndbasis
        Integer    nfbasis
        Integer    nobasis
        Integer    npbasis
        Integer    ndbasisw
        Integer    nfbasisw
        Integer    nobasisw
        Integer    npbasisw
        CQCAtom    *data
        CSelection *baselection
        CSelection *mmpselection

    cdef CQCAtomContainer *QCAtomContainer_Allocate                         ( Integer natoms )
    cdef CQCAtomContainer *QCAtomContainer_Clone                            ( CQCAtomContainer  *self )
    cdef void              QCAtomContainer_Deallocate                       ( CQCAtomContainer **self )
    cdef Status            QCAtomContainer_GetCoordinates3                  ( CQCAtomContainer  *self, CCoordinates3 *coordinates3, Boolean toBohrs, CCoordinates3 **qccoordinates3 )
    cdef CSelection       *QCAtomContainer_MakeFullSelection                ( CQCAtomContainer  *self )
    cdef Integer           QCAtomContainer_Size                             ( CQCAtomContainer  *self )
    cdef Status            QCAtomContainer_SetGradients3                    ( CQCAtomContainer  *self, CCoordinates3 *coordinates3, CCoordinates3 *qcgradients3, Boolean toInternalUnits, CCoordinates3 **gradients3 )

#===============================================================================
# . Class.
#===============================================================================
cdef class QCAtomContainer:

    cdef CQCAtomContainer *cObject
    cdef public object     isOwner

