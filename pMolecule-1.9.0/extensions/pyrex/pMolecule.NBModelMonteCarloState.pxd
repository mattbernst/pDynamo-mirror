#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelMonteCarloState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions             cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3             cimport CCoordinates3, Coordinates3
from pCore.Selection                cimport CSelection
from pCore.SelectionContainer       cimport CSelectionContainer
from pCore.Vector3                  cimport CVector3
from pMolecule.LJParameterContainer cimport CLJParameterContainer
from pMolecule.MMAtomContainer      cimport CMMAtomContainer
from pMolecule.SymmetryParameters   cimport SymmetryParameters, CSymmetryParameters

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "NBModelMonteCarloState.h":

    ctypedef struct CNBModelMonteCarloState "NBModelMonteCarloState":
        Integer                nisolates
        Integer                isolatescale
        Real                   chargeScale
        Real                   epsilonScale
        Real                   sigmaScale
        Real                   efmmel
        Real                   e1mmel
        Real                   efmmlj
        Real                   e1mmlj
        Boolean               *QFREE
        CVector3             **centers
        CCoordinates3         *coordinates3
        CLJParameterContainer *ljParameters
        CMMAtomContainer      *mmAtoms
        CSelectionContainer   *isolates
        CSymmetryParameters   *symmetryParameters

    cdef CNBModelMonteCarloState *NBModelMonteCarloState_Allocate                          ( Integer n )
    cdef void                     NBModelMonteCarloState_Deallocate                        ( CNBModelMonteCarloState **self )
    cdef void                     NBModelMonteCarloState_Initialize                        ( CNBModelMonteCarloState  *self, CCoordinates3 *coordinates3, CSymmetryParameters *symmetryParameters )
    cdef void                     NBModelMonteCarloState_ScaleIsolateInteractionParameters ( CNBModelMonteCarloState  *self, Integer isolate, Real chargeScale, Real epsilonScale, Real sigmaScale )
    cdef CNBModelMonteCarloState *NBModelMonteCarloState_Setup                             ( CSelectionContainer *isolates, CMMAtomContainer *mmAtoms, CLJParameterContainer *ljParameters, CSelection *fixedatoms )

#===============================================================================
# . Class.
#===============================================================================
cdef class NBModelMonteCarloState:

    cdef CNBModelMonteCarloState *cObject
    cdef public object            isOwner
