#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelMonteCarlo.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions               cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3               cimport CCoordinates3, Coordinates3
from pCore.PairList                   cimport SelfPairList
from pCore.Selection                  cimport Selection, CSelection
from pCore.SelectionContainer         cimport SelectionContainer, CSelectionContainer
from pCore.Vector3                    cimport CVector3
from pMolecule.LJParameterContainer   cimport LJParameterContainer, CLJParameterContainer
from pMolecule.MMAtomContainer        cimport MMAtomContainer, CMMAtomContainer
from pMolecule.NBModel                cimport NBModel
from pMolecule.NBModelMonteCarloState cimport NBModelMonteCarloState, NBModelMonteCarloState_Initialize, NBModelMonteCarloState_Setup, CNBModelMonteCarloState
from pMolecule.QCAtomContainer        cimport QCAtomContainer
from pMolecule.SymmetryParameters     cimport SymmetryParameters, SymmetryParameters_IsMinimumImageConventionSatisfied, CSymmetryParameters

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "NBModelMonteCarlo.h":

    ctypedef struct CNBModelMonteCarlo "NBModelMonteCarlo":
        Real buffer
        Real cutoff
        Real dielectric
        Real underflowel
        Real underflowlj

    cdef CNBModelMonteCarlo *NBModelMonteCarlo_Allocate         ( )
    cdef CNBModelMonteCarlo *NBModelMonteCarlo_Clone            ( CNBModelMonteCarlo  *self )
    cdef void                NBModelMonteCarlo_Deallocate       ( CNBModelMonteCarlo **self )
    cdef Real                NBModelMonteCarlo_MMMMEnergyFull   ( CNBModelMonteCarlo  *self, CNBModelMonteCarloState *nbState )
    cdef Real                NBModelMonteCarlo_MMMMEnergySingle ( CNBModelMonteCarlo  *self, Integer isolate, CNBModelMonteCarloState *nbState )

#===============================================================================
# . Class.
#===============================================================================
cdef class NBModelMonteCarlo ( NBModel ):

    cdef CNBModelMonteCarlo *cObject
    cdef public object       label
    cdef public object       isOwner
