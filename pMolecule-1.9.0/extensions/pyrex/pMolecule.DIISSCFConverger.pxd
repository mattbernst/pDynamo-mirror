#-------------------------------------------------------------------------------
# . File      : pMolecule.DIISSCFConverger.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions              cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Real2DArray               cimport CReal2DArray, Real2DArray
from pCore.SymmetricMatrix           cimport SymmetricMatrix, CSymmetricMatrix
from pMolecule.DIISSCFConvergerState cimport DIISSCFConvergerState, DIISSCFConvergerState_Deallocate, DIISSCFConvergerState_SetUp, CDIISSCFConvergerState
from pMolecule.QCOnePDM              cimport CQCOnePDM, QCOnePDM

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DIISSCFConverger.h":

    ctypedef struct CDIISSCFConverger "DIISSCFConverger":
        Boolean QUSERCA
        Integer maximumSCFCycles
        Integer ndiis
        Real    densityTolerance
        Real    diisonset
        Real    energytolerance

    cdef CDIISSCFConverger *DIISSCFConverger_Allocate                ( )
    cdef CDIISSCFConverger *DIISSCFConverger_Clone                   ( CDIISSCFConverger  *self )
    cdef Boolean            DIISSCFConverger_Continue                ( CDIISSCFConverger  *self, CDIISSCFConvergerState *convergerstate )
    cdef Boolean            DIISSCFConverger_Converged               ( CDIISSCFConverger  *self, CDIISSCFConvergerState *convergerstate )
    cdef void               DIISSCFConverger_Deallocate              ( CDIISSCFConverger **self )
    cdef void               DIISSCFConverger_IterateWithoutDensities ( CDIISSCFConverger  *self, CDIISSCFConvergerState *convergerstate, Real energy )
    cdef void               DIISSCFConverger_MakeDensities           ( CDIISSCFConverger  *self, CDIISSCFConvergerState *convergerstate )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DIISSCFConverger:

    cdef CDIISSCFConverger *cObject
    cdef public object      isOwner
