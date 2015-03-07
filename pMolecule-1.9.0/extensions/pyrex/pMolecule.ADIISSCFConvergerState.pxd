#-------------------------------------------------------------------------------
# . File      : pMolecule.ADIISSCFConvergerState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions        cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.AntisymmetricMatrix cimport CAntisymmetricMatrix
from pCore.Integer1DArray      cimport CInteger1DArray, Integer1DArray
from pCore.Real1DArray         cimport CReal1DArray, Real1DArray
from pCore.Real2DArray         cimport CReal2DArray, Real2DArray
from pCore.Status              cimport Status
from pCore.SymmetricMatrix     cimport SymmetricMatrix, CSymmetricMatrix
from pMolecule.QCOnePDM        cimport CQCOnePDM, QCOnePDM

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "ADIISSCFConvergerState.h":

    ctypedef enum ADIISSCFConvergerIterationType:
        ADIISSCFConvergerIterationType_ADIIS = 0
        ADIISSCFConvergerIterationType_DIIS  = 1
        ADIISSCFConvergerIterationType_Null  = 2
        ADIISSCFConvergerIterationType_ODA   = 3

    ctypedef struct CADIISSCFConvergerOnePDM "ADIISSCFConvergerOnePDM":
        pass

    ctypedef struct CADIISSCFConvergerState "ADIISSCFConvergerState":
        ADIISSCFConvergerIterationType  iterationType
        Boolean                         isConverged
        Integer                         diisActive
        Integer                         history
        Integer                         iteration
        Integer                         maximumHistory
        Real                            diisError
        Real                            energyChange
        Real                            odaOldEnergy
        Real                            oldEnergy
        Real                            odaMu
        Real                            rmsDifference
        CQCOnePDM                      *densityP
        CQCOnePDM                      *densityQ
        CReal2DArray                   *orthogonalizer
        CSymmetricMatrix               *overlap
        CADIISSCFConvergerOnePDM      **densitiesP
        CADIISSCFConvergerOnePDM      **densitiesQ
        CADIISSCFConvergerOnePDM       *odaDensityP
        CADIISSCFConvergerOnePDM       *odaDensityQ
        CAntisymmetricMatrix           *asmWork
        CInteger1DArray                *historyIndices
        CReal1DArray                   *adiisAlphas
        CReal1DArray                   *adiisGradients
        CReal1DArray                   *energies
        CReal2DArray                   *adiisTracesP
        CReal2DArray                   *adiisTracesQ
        CReal2DArray                   *diisTraces
        CSymmetricMatrix               *adiisHessian

    cdef CADIISSCFConvergerState *ADIISSCFConvergerState_Allocate   ( )
    cdef void                     ADIISSCFConvergerState_Deallocate ( CADIISSCFConvergerState **self )
    cdef CADIISSCFConvergerState *ADIISSCFConvergerState_SetUp      ( Boolean QUSERCA, Integer maximumHistory, CQCOnePDM *densityp, CQCOnePDM *densityq, \
                                                                                 CSymmetricMatrix *overlap, CReal2DArray *orthogonalizer, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ADIISSCFConvergerState:

    cdef CADIISSCFConvergerState *cObject
    cdef public object            isOwner
    cdef public object            table
