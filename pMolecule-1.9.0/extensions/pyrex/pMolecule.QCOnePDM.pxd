#-------------------------------------------------------------------------------
# . File      : pMolecule.QCOnePDM.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions      cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Real1DArray       cimport CReal1DArray, Real1DArray, Real1DArray_Allocate, Real1DArray_Clone, Real1DArray_Deallocate
from pCore.Real2DArray       cimport CReal2DArray, Real2DArray, Real2DArray_Allocate, Real2DArray_Clone, Real2DArray_Deallocate
from pCore.Status            cimport Status, Status_Continue
from pCore.SymmetricMatrix   cimport CSymmetricMatrix, SymmetricMatrix, SymmetricMatrix_Allocate, SymmetricMatrix_Clone, SymmetricMatrix_Deallocate

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "QCOnePDM.h":

    ctypedef enum QCOnePDMDensityType:
        QCOnePDMDensityType_Alpha = 0,
        QCOnePDMDensityType_Beta  = 1,
        QCOnePDMDensityType_Spin  = 2,
        QCOnePDMDensityType_Total = 3

    ctypedef enum QCOnePDMOccupancyType:
        QCOnePDMOccupancyType_Cardinal           = 0,
        QCOnePDMOccupancyType_FractionalFixed    = 1,
        QCOnePDMOccupancyType_FractionalVariable = 2

    ctypedef struct CQCOnePDM "QCOnePDM":
       Boolean                isValid
       Integer                numberOccupied
       Integer                numberOrbitals
       Real                   fermiBroadening
       Real                   fermiEnergy
       Real                   occupancyEnergy
       Real                   occupancyFactor
       Real                   totalCharge
       QCOnePDMDensityType    densityType
       QCOnePDMOccupancyType  occupancyType
       CReal1DArray          *energies
       CReal1DArray          *occupancies
       CReal2DArray          *orbitals
       CSymmetricMatrix      *density
       CSymmetricMatrix      *fock

    cdef CQCOnePDM        *QCOnePDM_Allocate                 ( Integer length, Status *status )
    cdef void              QCOnePDM_AllocateOrbitalData      ( CQCOnePDM  *self, CReal2DArray *orthogonalizer, Status *status )
    cdef CQCOnePDM        *QCOnePDM_Clone                    ( CQCOnePDM  *self, Status *status )
    cdef void              QCOnePDM_AlphaBetaFromTotalSpin   ( CQCOnePDM  *self, CQCOnePDM *other )
    cdef void              QCOnePDM_AlphaBetaToTotalSpin     ( CQCOnePDM  *self, CQCOnePDM *other )
    cdef void              QCOnePDM_Deallocate               ( CQCOnePDM **self )
    cdef void              QCOnePDM_DeallocateOrbitalData    ( CQCOnePDM  *self )
    cdef CQCOnePDM        *QCOnePDM_FromDiagonalGuess        ( QCOnePDMDensityType densityType, QCOnePDMOccupancyType occupancyType, Integer length, Real totalCharge, Integer numberFractionalHOOs, Integer numberFractionalLUOs, Status *status )
    cdef Integer           QCOnePDM_HOMOIndex                ( CQCOnePDM  *self )
    cdef Integer           QCOnePDM_LUMOIndex                ( CQCOnePDM  *self )
    cdef Real              QCOnePDM_Make                     ( CQCOnePDM  *self, Status *status )
    cdef void              QCOnePDM_MakeElementary           ( CSymmetricMatrix *self, Integer noccupied, CReal1DArray *occupancies, CReal2DArray *orbitals, Status *status )
    cdef Real              QCOnePDM_MakeFromFock             ( CQCOnePDM  *self, CReal2DArray *orthogonalizer, Status *status )
    cdef CSymmetricMatrix *QCOnePDM_MakeWeighted             ( CQCOnePDM  *self, CQCOnePDM *other, Status *status )
    cdef void              QCOnePDM_Reorthogonalize          ( CQCOnePDM  *self, CReal2DArray *orthogonalizer, Status *status )
    cdef Boolean           QCOnePDM_ResetDensityFromDensity  ( CQCOnePDM  *self, CSymmetricMatrix *density, CSymmetricMatrix *overlap, Status *status )
    cdef Boolean           QCOnePDM_ResetDensityFromOrbitals ( CQCOnePDM  *self, CReal2DArray *orbitals, Status *status )
    cdef void              QCOnePDM_SpinExpectationValues    ( CQCOnePDM  *self, CQCOnePDM *other, CSymmetricMatrix *overlap, Real *Sz, Real *S2 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCOnePDM:

    cdef CQCOnePDM     *cObject
    cdef public object  isOwner
    cdef public object    owner
