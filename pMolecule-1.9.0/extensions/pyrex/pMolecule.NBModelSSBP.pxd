#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelSSBP.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions             cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3             cimport CCoordinates3, Coordinates3
from pCore.Integer2DArray           cimport CInteger2DArray, Integer2DArray, Integer2DArray_Clone
from pCore.PairList                 cimport SelfPairList
from pCore.Real1DArray              cimport CReal1DArray
from pCore.Selection                cimport CSelection, Selection, Selection_Clone
from pCore.SymmetricMatrix          cimport CSymmetricMatrix
from pMolecule.LJParameterContainer cimport LJParameterContainer
from pMolecule.MMAtomContainer      cimport MMAtomContainer, CMMAtomContainer
from pMolecule.NBModelFull          cimport NBModelFull, NBModelFull_Deallocate
from pMolecule.QCAtomContainer      cimport CQCAtomContainer, QCAtomContainer
from pMolecule.QCMMInteractionState cimport CQCMMInteractionState, QCMMInteractionState, QCMMInteractionState_AllocateQCQCPotentials
from pMolecule.SSBPModelState       cimport CSSBPModelState, SSBPModelState, SSBPModelState_Initialize, SSBPModelState_Setup

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "SSBPModel.h":

    ctypedef struct CSSBPModel "SSBPModel":
        Boolean doAngularPotential
        Boolean doCavityPotential
        Boolean doEmpiricalCorrection
        Boolean doHardSphereRestriction
        Boolean doKirkwood
        Boolean fixCavityRadius
        Integer maximumL
        Real    cavityRadius
        Real    cavityRadiusIncrement
        Real    dielectricInside
        Real    dielectricOutside
        Real    empirical1
        Real    empirical2
        Real    kirkwoodRadiusIncrement
        Real    pressure
        Real    surfaceTension

    cdef CSSBPModel *SSBPModel_Allocate   ( )
    cdef CSSBPModel *SSBPModel_Clone      ( CSSBPModel  *self )
    cdef void        SSBPModel_Deallocate ( CSSBPModel **self )
    cdef void        SSBPModel_Energy     ( CSSBPModel  *self, Boolean doGradients, CSSBPModelState *state )
    cdef Real        SSBPModel_Kirkwood   ( CSSBPModel  *self, Boolean doGradients, CSSBPModelState *state )

#===============================================================================
# . Class.
#===============================================================================
cdef class NBModelSSBP ( NBModelFull ):

    cdef CSSBPModel    *cObject2
    cdef public object  doNonBondingTerms
    cdef public object  isOwner2

