from pCore.cDefinitions           cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3           cimport CCoordinates3
from pCore.Selection              cimport CSelection
from pMolecule.SymmetryParameters cimport CSymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "fastio.h":
    pass

cdef extern from "DCDHandle.h":

    ctypedef enum CDCDStatus "DCDStatus":
        CDCDStatus_Normal                  "DCDStatus_Normal"
        CDCDStatus_AtomNumberMismatch      "DCDStatus_AtomNumberMismatch"
        CDCDStatus_BadFormat               "DCDStatus_BadFormat"
        CDCDStatus_BadRead                 "DCDStatus_BadRead"
        CDCDStatus_BadSeek                 "DCDStatus_BadSeek"
        CDCDStatus_BadWrite                "DCDStatus_BadWrite"
        CDCDStatus_FileAccessFailure       "DCDStatus_FileAccessFailure"
        CDCDStatus_InvalidDataObject       "DCDStatus_InvalidDataObject"
        CDCDStatus_InvalidFrameIndex       "DCDStatus_InvalidFrameIndex"
        CDCDStatus_MemoryAllocationFailure "DCDStatus_MemoryAllocationFailure"
        CDCDStatus_OpenFailed              "DCDStatus_OpenFailed"

    ctypedef struct CDCDHandle "DCDHandle":
        Boolean hasUnitCell
        Boolean isXPLOR
        Boolean useVelocityHeader
        Integer currentFrame
        Integer numberOfAtomIndices
        Integer numberOfAtoms
        Integer numberOfFrames
        Integer saveFrequency
        Integer startingFrame
        Real    timeStep

    cdef CDCDHandle *DCDHandle_Allocate              ( )
    cdef CDCDStatus  DCDHandle_AllocateQW            ( CDCDHandle  *self )
    cdef CDCDStatus  DCDHandle_CheckNumberOfAtoms    ( CDCDHandle  *self, Integer numberOfAtoms )
    cdef Integer     DCDHandle_CurrentFrame          ( CDCDHandle  *self )
    cdef void        DCDHandle_Deallocate            ( CDCDHandle **self )
    cdef Integer     DCDHandle_NumberOfFrames        ( CDCDHandle  *self )
    cdef CDCDStatus  DCDHandle_SetAtomIndices        ( CDCDHandle  *self, CSelection *selection )
    cdef CDCDStatus  DCDHandle_SetData3              ( CDCDHandle  *self, CCoordinates3 *data3 )
    cdef CDCDStatus  DCDHandle_SetSymmetryParameters ( CDCDHandle  *self, CSymmetryParameters *symmetryParameters )

# . Functions.
cdef DCDStatus_Check ( CDCDStatus status )
