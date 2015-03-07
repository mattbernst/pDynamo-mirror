#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelORCAState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions                    cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3                    cimport CCoordinates3, Coordinates3, Coordinates3_GetItem, Coordinates3_SetItem
from pCore.PairList                        cimport CPairList
from pCore.Selection                       cimport CSelection
from pCore.Real1DArray                     cimport CReal1DArray, Real1DArray_GetItem
from pMolecule.LJParameterContainer        cimport CLJParameterContainer
from pMolecule.MMAtomContainer             cimport CMMAtomContainer
from pMolecule.QCAtomContainer             cimport CQCAtomContainer
from pMolecule.QCMMLinkAtomCouplingOptions cimport QCMMLinkAtomCoupling, QCMMLinkAtomCoupling_MM, QCMMLinkAtomCoupling_RC, QCMMLinkAtomCoupling_RD

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "NBModelORCAState.h":

    ctypedef struct CNBModelOrcaState "NBModelOrcaState":
        QCMMLinkAtomCoupling   qcmmcoupling
        Real                   emmel
        Real                   emmel14
        Real                   emmlj
        Real                   emmlj14
        Real                   eqcmmlj
        Real                   eqcmmlj14
        Boolean               *QE14
        Boolean               *QFREE
        Boolean               *QINCL
        Boolean               *QMM
        CCoordinates3         *coordinates3
        CPairList             *exclusions
        CCoordinates3         *gradients3
        CLJParameterContainer *ljParameters
        CLJParameterContainer *ljParameters14
        CPairList             *interactions14
        CMMAtomContainer      *mmAtoms
        CQCAtomContainer      *qcAtoms
        CCoordinates3         *mmcoordinates3
        CCoordinates3         *mmgradients3
        CPairList             *mmexclusions
        CPairList             *mminteractions14
        CReal1DArray          *mmCharges

    cdef CNBModelOrcaState *NBModelORCAState_Allocate   ( Integer n )
    cdef void               NBModelORCAState_Deallocate ( CNBModelOrcaState **self )
    cdef void               NBModelORCAState_Finalize   ( CNBModelOrcaState  *self )
    cdef void               NBModelORCAState_Initialize ( CNBModelOrcaState  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef CNBModelOrcaState *NBModelORCAState_Setup      ( CMMAtomContainer *mmAtoms, CQCAtomContainer *qcAtoms, CSelection *fixedatoms, CPairList *exclusions, CPairList *interactions14, \
                                                                                                              CLJParameterContainer *ljParameters, CLJParameterContainer *ljParameters14, \
                                                                                                                                                       QCMMLinkAtomCoupling qcmmcoupling  )

#===============================================================================
# . Class.
#===============================================================================
cdef class NBModelORCAState:

    cdef CNBModelOrcaState *cObject
    cdef public object      isOwner
