#-------------------------------------------------------------------------------
# . File      : pMolecule.SSBPModelState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions        cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Integer2DArray      cimport CInteger2DArray, Integer2DArray_Length
from pCore.Coordinates3        cimport CCoordinates3, Coordinates3
from pCore.Real1DArray         cimport CReal1DArray
from pCore.Selection           cimport CSelection, Selection_Size
from pCore.Status              cimport Status
from pCore.SymmetricMatrix     cimport CSymmetricMatrix

from pMolecule.MMAtomContainer cimport CMMAtomContainer, MMAtomContainer_NumberOfActiveAtoms
from pMolecule.QCAtomContainer cimport CQCAtomContainer, QCAtomContainer_Size

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "SSBPModelState.h":

    ctypedef struct CSSBPModelState "SSBPModelState":
        Boolean           atomOutsideCavity
        Integer           numberOfAtoms
        Integer           particleIndex
        Real              eAngular
        Real              eCavity
        Real              eEmpiricalCorrection
        Real              eHardSphere
        Real              eKirkwood
        Real              eKirkwoodCheck
        Real              eTotal
        Real              qTotal
        Real              radius
        CCoordinates3    *gradients3
        CMMAtomContainer *mmAtoms
        CQCAtomContainer *qcAtoms
        CInteger2DArray  *waterAtomIndices
        CSelection       *cavitySelection
        CSelection       *radiusSelection

    cdef CSSBPModelState *SSBPModelState_Allocate   ( )
    cdef void             SSBPModelState_Deallocate ( CSSBPModelState **self )
    cdef void             SSBPModelState_Initialize ( CSSBPModelState  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef CSSBPModelState *SSBPModelState_Setup      ( Boolean doKirkwood, Integer maximumL, CMMAtomContainer *mmAtoms,
                                                               CQCAtomContainer *qcAtoms, CSelection *cavitySelection,
                                                       CSelection *radiusSelection, CInteger2DArray *waterAtomIndices,
                                                                CReal1DArray *qcCharges, CReal1DArray *qcmmPotentials,
                                                                    CSymmetricMatrix *qcqcPotentials, Status *status )

#===============================================================================
# . Class.
#===============================================================================
cdef class SSBPModelState:

    cdef CSSBPModelState *cObject
    cdef public object    isOwner
