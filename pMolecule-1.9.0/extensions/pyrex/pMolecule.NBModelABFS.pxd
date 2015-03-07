#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelABFS.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions                    cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3                    cimport CCoordinates3, Coordinates3
from pCore.PairList                        cimport CPairList, SelfPairList
from pCore.PairListGenerator               cimport CPairListGenerator, PairListGenerator
from pCore.Real1DArray                     cimport CReal1DArray
from pCore.Selection                       cimport Selection, CSelection
from pCore.Status                          cimport Status, Status_Continue
from pCore.SymmetricMatrix                 cimport CSymmetricMatrix
from pCore.Transformation3Container        cimport Transformation3Container, CTransformation3Container
from pMolecule.LJParameterContainer        cimport LJParameterContainer, CLJParameterContainer
from pMolecule.MMAtomContainer             cimport MMAtomContainer, CMMAtomContainer
from pMolecule.NBModel                     cimport NBModel
from pMolecule.NBModelABFSState            cimport CNBModelABFSState, NBModelABFSState, NBModelABFSState_Initialize, NBModelABFSState_SetUp, NBModelABFSState_SetUpCentering
from pMolecule.PairwiseInteraction         cimport CPairwiseInteractionABFS, PairwiseInteractionABFS
from pMolecule.QCAtomContainer             cimport QCAtomContainer, CQCAtomContainer
from pMolecule.QCMMLinkAtomCouplingOptions cimport QCMMLinkAtomCoupling, QCMMLinkAtomCoupling_MM, QCMMLinkAtomCoupling_RC, QCMMLinkAtomCoupling_RD
from pMolecule.QCMMInteractionState        cimport QCMMInteractionState, CQCMMInteractionState
from pMolecule.SymmetryParameterGradients  cimport SymmetryParameterGradients, CSymmetryParameterGradients
from pMolecule.SymmetryParameters          cimport SymmetryParameters, CSymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "NBModelABFS.h":

    ctypedef struct CNBModelABFS "NBModelABFS":
        Boolean              checkForInverses     
        Boolean              useCentering         
        Integer              imageExpandFactor    
        Real                 dampingCutoff        
        Real                 dielectric           
        Real                 electrostaticScale14 
        Real                 innerCutoff          
        Real                 listCutoff           
        Real                 outerCutoff          
        QCMMLinkAtomCoupling qcmmCoupling

    cdef CNBModelABFS *NBModelABFS_Allocate       ( )
    cdef CNBModelABFS *NBModelABFS_Clone          ( CNBModelABFS  *self )
    cdef void          NBModelABFS_Deallocate     ( CNBModelABFS **self )
    cdef void          NBModelABFS_MMMMEnergy     ( CNBModelABFS  *self, CPairwiseInteractionABFS *mmmmPairwiseInteraction, CNBModelABFSState *nbState )
    cdef void          NBModelABFS_QCMMEnergyLJ   ( CNBModelABFS  *self, CPairwiseInteractionABFS *mmmmPairwiseInteraction, CNBModelABFSState *nbState )
    cdef void          NBModelABFS_QCMMGradients  ( CNBModelABFS  *self, CPairwiseInteractionABFS *qcmmPairwiseInteraction, CPairwiseInteractionABFS *qcqcPairwiseInteraction, CNBModelABFSState *nbState )
    cdef void          NBModelABFS_QCMMPotentials ( CNBModelABFS  *self, CPairwiseInteractionABFS *qcmmPairwiseInteraction, CPairwiseInteractionABFS *qcqcPairwiseInteraction, CNBModelABFSState *nbState )
    cdef Boolean       NBModelABFS_Update         ( CNBModelABFS  *self, CPairListGenerator *generator, CNBModelABFSState *nbState, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelABFS ( NBModel ):

    cdef CNBModelABFS                  *cObject
    cdef public PairListGenerator       generator
    cdef public PairwiseInteractionABFS mmmmPairwiseInteraction
    cdef public PairwiseInteractionABFS qcmmPairwiseInteraction
    cdef public PairwiseInteractionABFS qcqcPairwiseInteraction
    cdef public object                  label
    cdef public object                  isOwner
