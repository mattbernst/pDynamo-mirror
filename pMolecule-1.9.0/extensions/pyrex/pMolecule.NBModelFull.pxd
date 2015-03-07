#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelFull.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions                    cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3                    cimport CCoordinates3, Coordinates3
from pCore.PairList                        cimport CPairList, SelfPairList
from pCore.Real1DArray                     cimport CReal1DArray
from pCore.Selection                       cimport Selection, CSelection
from pMolecule.LJParameterContainer        cimport LJParameterContainer, CLJParameterContainer
from pMolecule.MMAtomContainer             cimport MMAtomContainer, CMMAtomContainer
from pMolecule.NBModel                     cimport NBModel
from pMolecule.NBModelFullState            cimport NBModelFullState, NBModelFullState_Initialize, NBModelFullState_Setup, CNBModelFullState
from pMolecule.QCAtomContainer             cimport QCAtomContainer, CQCAtomContainer
from pMolecule.QCMMLinkAtomCouplingOptions cimport QCMMLinkAtomCoupling, QCMMLinkAtomCoupling_MM, QCMMLinkAtomCoupling_RC, QCMMLinkAtomCoupling_RD
from pMolecule.QCMMInteractionState        cimport QCMMInteractionState, CQCMMInteractionState

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "NBModelFull.h":

    ctypedef struct CNBModelFull "NBModelFull":
        Real                 dielectric
        Real                 electrostaticscale14
        QCMMLinkAtomCoupling qcmmcoupling

    cdef CNBModelFull *NBModelFull_Allocate       ( )
    cdef CNBModelFull *NBModelFull_Clone          ( CNBModelFull  *self )
    cdef void          NBModelFull_Deallocate     ( CNBModelFull **self )
    cdef void          NBModelFull_MMMMEnergy     ( CNBModelFull  *self, CNBModelFullState *nbState )
    cdef void          NBModelFull_QCMMEnergyLJ   ( CNBModelFull  *self, CNBModelFullState *nbState )
    cdef void          NBModelFull_QCMMGradients  ( CNBModelFull  *self, CNBModelFullState *nbState )
    cdef void          NBModelFull_QCMMPotentials ( CNBModelFull  *self, CNBModelFullState *nbState )

#===============================================================================
# . Class.
#===============================================================================
cdef class NBModelFull ( NBModel ):

    cdef CNBModelFull *cObject
    cdef public object label
    cdef public object isOwner
