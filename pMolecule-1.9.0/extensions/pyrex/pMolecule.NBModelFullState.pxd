#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelFullState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions             cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3             cimport CCoordinates3, Coordinates3
from pCore.PairList                 cimport CPairList
from pCore.Real1DArray              cimport CReal1DArray
from pCore.Selection                cimport CSelection
from pMolecule.LJParameterContainer cimport CLJParameterContainer
from pMolecule.MMAtomContainer      cimport CMMAtomContainer
from pMolecule.QCAtomContainer      cimport CQCAtomContainer

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "NBModelFullState.h":

    ctypedef struct CNBModelFullState "NBModelFullState":
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
        Integer               *baindex
        CCoordinates3         *coordinates3
        CPairList             *exclusions
        CCoordinates3         *gradients3
        CPairList             *interactions14
        CLJParameterContainer *ljParameters
        CLJParameterContainer *ljParameters14
        CMMAtomContainer      *mmAtoms
        CQCAtomContainer      *qcAtoms
        CReal1DArray          *qcCharges
        CReal1DArray          *qcmmPotentials

    cdef CNBModelFullState *NBModelFullState_Allocate   ( Integer n )
    cdef void               NBModelFullState_Deallocate ( CNBModelFullState **self )
    cdef void               NBModelFullState_Initialize ( CNBModelFullState  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef CNBModelFullState *NBModelFullState_Setup      ( CMMAtomContainer *mmAtoms, CQCAtomContainer *qcAtoms, CSelection *fixedatoms, CPairList *exclusions, CPairList *interactions14,
                                                                                                              CLJParameterContainer *ljParameters, CLJParameterContainer *ljParameters14,
                                                                                                                                  CReal1DArray *qcCharges, CReal1DArray *qcmmPotentials )

#===============================================================================
# . Class.
#===============================================================================
cdef class NBModelFullState:

    cdef CNBModelFullState *cObject
    cdef public object          isOwner
