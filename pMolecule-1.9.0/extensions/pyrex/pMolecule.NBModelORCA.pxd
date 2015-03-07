#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelORCA.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions                    cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3                    cimport CCoordinates3, Coordinates3
from pCore.PairList                        cimport CPairList, SelfPairList
from pCore.Selection                       cimport Selection, CSelection
from pMolecule.LJParameterContainer        cimport LJParameterContainer, CLJParameterContainer
from pMolecule.MMAtomContainer             cimport MMAtomContainer, CMMAtomContainer
from pMolecule.NBModel                     cimport NBModel
from pMolecule.NBModelORCAState            cimport NBModelORCAState, NBModelORCAState_Finalize, NBModelORCAState_Initialize, NBModelORCAState_Setup, CNBModelOrcaState
from pMolecule.QCAtomContainer             cimport QCAtomContainer, CQCAtomContainer
from pMolecule.QCMMLinkAtomCouplingOptions cimport QCMMLinkAtomCoupling, QCMMLinkAtomCoupling_MM, QCMMLinkAtomCoupling_RC, QCMMLinkAtomCoupling_RD
from pMolecule.QCMMInteractionState        cimport QCMMInteractionState, CQCMMInteractionState

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "NBModelORCA.h":

    ctypedef struct CNBModelOrca "NBModelOrca":
        Real               electrostaticscale14
        QCMMLinkAtomCoupling qcmmcoupling

    cdef CNBModelOrca *NBModelORCA_Allocate     ( )
    cdef CNBModelOrca *NBModelORCA_Clone        ( CNBModelOrca  *self )
    cdef void          NBModelORCA_Deallocate   ( CNBModelOrca **self )
    cdef void          NBModelORCA_MMMMEnergy   ( CNBModelOrca  *self, CNBModelOrcaState *nbState )
    cdef void          NBModelORCA_QCMMEnergyLJ ( CNBModelOrca  *self, CNBModelOrcaState *nbState )

#===============================================================================
# . Class.
#===============================================================================
cdef class NBModelORCA ( NBModel ):

    cdef CNBModelOrca *cObject
    cdef public object label
    cdef public object isOwner
