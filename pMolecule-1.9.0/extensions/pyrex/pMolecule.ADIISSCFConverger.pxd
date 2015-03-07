#-------------------------------------------------------------------------------
# . File      : pMolecule.ADIISSCFConverger.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions               cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Real2DArray                cimport CReal2DArray, Real2DArray
from pCore.SymmetricMatrix            cimport SymmetricMatrix, CSymmetricMatrix
from pMolecule.ADIISSCFConvergerState cimport ADIISSCFConvergerState, ADIISSCFConvergerState_Deallocate, ADIISSCFConvergerState_SetUp, CADIISSCFConvergerState
from pMolecule.QCOnePDM               cimport CQCOnePDM, QCOnePDM

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "ADIISSCFConverger.h":

    ctypedef struct CADIISSCFConverger "ADIISSCFConverger":
        Boolean useEDIIS
        Boolean useODA
        Integer maximumHistory
        Integer maximumSCFCycles
        Real    densityTolerance
        Real    diisOff
        Real    diisOn
        Real    energyTolerance
        Real    minimumCoefficient

    cdef CADIISSCFConverger *ADIISSCFConverger_Allocate     ( )
    cdef CADIISSCFConverger *ADIISSCFConverger_Clone        ( CADIISSCFConverger  *self )
    cdef Boolean             ADIISSCFConverger_Continue     ( CADIISSCFConverger  *self, CADIISSCFConvergerState *state )
    cdef void                ADIISSCFConverger_Deallocate   ( CADIISSCFConverger **self )
    cdef Boolean             ADIISSCFConverger_IsConverged  ( CADIISSCFConverger  *self, CADIISSCFConvergerState *state )
    cdef void                ADIISSCFConverger_IterateStart ( CADIISSCFConverger  *self, CADIISSCFConvergerState *state, Real energy )
    cdef void                ADIISSCFConverger_IterateStop  ( CADIISSCFConverger  *self, CADIISSCFConvergerState *state )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ADIISSCFConverger:

    cdef CADIISSCFConverger *cObject
    cdef public object       debugLog
    cdef public object       isOwner
    cdef public object       minimizationsPerCycle
    cdef public object       optimizer
