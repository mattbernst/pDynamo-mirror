#-------------------------------------------------------------------------------
# . File      : pMolecule.DIISSCFConvergerState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions             cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.AntisymmetricMatrix cimport CAntisymmetricMatrix
from pCore.Real2DArray         cimport CReal2DArray, Real2DArray
from pCore.SymmetricMatrix     cimport SymmetricMatrix, CSymmetricMatrix
from pMolecule.QCOnePDM        cimport CQCOnePDM, QCOnePDM

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DIISSCFConvergerState.h":

    ctypedef struct CDIISSCFConvergerState "DIISSCFConvergerState":
        Boolean               isConverged
        Boolean               QDAMP
        Boolean               QDIIS
        Boolean               QRCA
        Boolean               QRCAPREDICTED
        Integer               iteration
        Integer               matnum
        Integer               ndiis
        Real                  damp
        Real                  deavg
        Real                  deltaeold
        Real                  diiserror
        Real                  dmpsav
        Real                  eold
        Real                  rcamu
        Real                  rmsdifference
        Integer               *matind
        CReal2DArray          *bcoeff
        CReal2DArray          *orthogonalizer
        CAntisymmetricMatrix  *asmwork
        CAntisymmetricMatrix **errorp
        CAntisymmetricMatrix **errorq
        CQCOnePDM             *densityp
        CQCOnePDM             *densityq
        CSymmetricMatrix      *dfockp1
        CSymmetricMatrix      *dfockp2
        CSymmetricMatrix      *dfockq1
        CSymmetricMatrix      *dfockq2
        CSymmetricMatrix      *overlap
        CSymmetricMatrix      *rdensityp
        CSymmetricMatrix      *rfockp
        CSymmetricMatrix      *rdensityq
        CSymmetricMatrix      *rfockq
        CSymmetricMatrix     **xdensityp
        CSymmetricMatrix     **xdensityq
        CSymmetricMatrix     **xfockp
        CSymmetricMatrix     **xfockq

    cdef CDIISSCFConvergerState *DIISSCFConvergerState_Allocate   ( )
    cdef void                    DIISSCFConvergerState_Deallocate ( CDIISSCFConvergerState **self )
    cdef CDIISSCFConvergerState *DIISSCFConvergerState_SetUp      ( Boolean QUSERCA, Integer ndiis, CQCOnePDM *densityp, CQCOnePDM *densityq, CSymmetricMatrix *overlap, CReal2DArray *orthogonalizer )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DIISSCFConvergerState:

    cdef CDIISSCFConvergerState *cObject
    cdef public object           isOwner
    cdef public object           table
