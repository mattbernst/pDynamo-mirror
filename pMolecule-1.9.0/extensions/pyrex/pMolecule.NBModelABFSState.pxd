#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelABFSState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions                    cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3                    cimport CCoordinates3, Coordinates3
from pCore.PairList                        cimport CPairList
from pCore.Real1DArray                     cimport CReal1DArray
from pCore.Selection                       cimport CSelection
from pCore.SelectionContainer              cimport CSelectionContainer
from pCore.Status                          cimport Status, Status_Continue
from pCore.SymmetricMatrix                 cimport CSymmetricMatrix
from pCore.Transformation3Container        cimport CTransformation3Container
from pMolecule.ImageList                   cimport ImageList_NumberOfImages, ImageList_NumberOfPairs, CImageList
from pMolecule.LJParameterContainer        cimport CLJParameterContainer
from pMolecule.MMAtomContainer             cimport CMMAtomContainer
from pMolecule.QCAtomContainer             cimport CQCAtomContainer
from pMolecule.QCMMLinkAtomCouplingOptions cimport QCMMLinkAtomCoupling, QCMMLinkAtomCoupling_MM, QCMMLinkAtomCoupling_RC, QCMMLinkAtomCoupling_RD
from pMolecule.SymmetryParameters          cimport CSymmetryParameters
from pMolecule.SymmetryParameterGradients  cimport CSymmetryParameterGradients

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "NBModelABFSState.h":

    ctypedef struct CNBModelABFSState "NBModelABFSState":
        Boolean                      isNew
        Boolean                      useCentering
        Boolean                      useGridSearch
        Integer                      numberOfCalls
        Integer                      numberOfUpdates
        QCMMLinkAtomCoupling         qcmmCoupling
        Real                         listCutoff
        Real                         outerCutoff
        Real                         eimmmel
        Real                         eimmmlj
        Real                         eimqcmmlj
        Real                         eimqcqclj
        Real                         emmel
        Real                         emmel14
        Real                         emmlj
        Real                         emmlj14
        Real                         eqcmmlj
        Real                         eqcmmlj14
        Real                         numberOfMMMMPairs
        Real                         numberOfQCMMElPairs
        Real                         numberOfQCMMLJPairs
        Real                         numberOfMMMMImageImages
        Real                         numberOfMMMMImagePairs
        Real                         numberOfQCMMElImageImages
        Real                         numberOfQCMMElImagePairs
        Real                         numberOfQCMMLJImageImages
        Real                         numberOfQCMMLJImagePairs
        Real                         numberOfQCQCElImageImages
        Real                         numberOfQCQCElImagePairs
        Real                         numberOfQCQCLJImageImages
        Real                         numberOfQCQCLJImagePairs 
        CCoordinates3               *coordinates3
        CCoordinates3               *mmCoordinates3
        CCoordinates3               *qcCoordinates3
        CCoordinates3               *referenceCoordinates3
        CSelection                  *freeSelection
        CSelection                  *mmSelection
        CSelection                  *qcbSelection
        CSelection                  *qcpSelection
#        CSelection                  *tocenter
        CSymmetryParameters         *referenceSymmetryParameters
        CReal1DArray                *mmCharges
        CImageList                  *inbmmmm
        CImageList                  *inbqcmmel
        CImageList                  *inbqcmmlj
        CImageList                  *inbqcqcel
        CImageList                  *inbqcqclj
        CPairList                   *nbmmmm
        CPairList                   *nbmmmm14
        CPairList                   *nbqcmmel
        CPairList                   *nbqcmmel14
        CPairList                   *nbqcmmlj
        CPairList                   *nbqcmmlj14
        CPairList                   *exclusions
        CCoordinates3               *gradients3
        CCoordinates3               *mmGradients3
        CPairList                   *interactions14
        CLJParameterContainer       *ljParameters
        CLJParameterContainer       *ljParameters14
        CMMAtomContainer            *mmAtoms
        CQCAtomContainer            *qcAtoms
        CReal1DArray                *qcCharges
        CReal1DArray                *qcmmPotentials
        CSelection                  *fixedAtoms
#        CSelectionContainer         *mmisolates
        CSymmetricMatrix            *qcqcPotentials
        CSymmetryParameters         *symmetryParameters
        CSymmetryParameterGradients *symmetryParameterGradients
        CTransformation3Container   *transformations

    cdef CNBModelABFSState *NBModelABFSState_Allocate             ( Integer                      n                          )
    cdef void               NBModelABFSState_Deallocate           ( CNBModelABFSState          **self                       )
    cdef void               NBModelABFSState_Initialize           ( CNBModelABFSState           *self                       ,
                                                                    CCoordinates3               *coordinates3               ,
                                                                    CSymmetryParameters         *symmetryParameters         ,
                                                                    CCoordinates3               *gradients3                 ,
                                                                    CSymmetryParameterGradients *symmetryParameterGradients )
    cdef CNBModelABFSState *NBModelABFSState_SetUp                ( CMMAtomContainer            *mmAtoms                    ,
                                                                    CQCAtomContainer            *qcAtoms                    ,
                                                                    CSelection                  *fixedAtoms                 ,
                                                                    CPairList                   *exclusions                 ,
                                                                    CPairList                   *interactions14             ,
                                                                    CLJParameterContainer       *ljParameters               ,
                                                                    CLJParameterContainer       *ljParameters14             ,
                                                                    CReal1DArray                *qcCharges                  ,
                                                                    CReal1DArray                *qcmmPotentials             ,
                                                                    CSymmetricMatrix            *qcqcPotentials             ,
                                                                    CTransformation3Container   *transformations            ,
                                                                    QCMMLinkAtomCoupling         qcmmCoupling               )
    cdef void               NBModelABFSState_SetUpCentering       ( CNBModelABFSState           *self                       ,
                                                                    Boolean                      useCentering               ,
                                                                    Status                      *status                     )
    cdef void               NBModelABFSState_StatisticsInitialize ( CNBModelABFSState           *self                       )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelABFSState:

    cdef CNBModelABFSState *cObject
    cdef public object      isOwner
