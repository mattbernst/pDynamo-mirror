from pBabel.DCDHandle             cimport CDCDHandle                         , \
                                          CDCDStatus                         , \
                                          CDCDStatus_MemoryAllocationFailure , \
                                          DCDHandle_Allocate                 , \
                                          DCDHandle_AllocateQW               , \
                                          DCDHandle_CheckNumberOfAtoms       , \
                                          DCDHandle_CurrentFrame             , \
                                          DCDHandle_NumberOfFrames           , \
                                          DCDHandle_SetAtomIndices           , \
                                          DCDHandle_SetData3                 , \
                                          DCDHandle_SetSymmetryParameters    , \
                                          DCDStatus_Check
from pCore.cDefinitions           cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3           cimport Coordinates3
from pCore.Selection              cimport Selection
from pMolecule.SymmetryParameters cimport SymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DCDRead.h":

    cdef void       DCDRead_Close     ( CDCDHandle **self )
    cdef CDCDStatus DCDRead_Frame     ( CDCDHandle  *self )
    cdef CDCDStatus DCDRead_GotoFrame ( CDCDHandle  *self, Integer f )
    cdef CDCDStatus DCDRead_Header    ( CDCDHandle  *self )
    cdef CDCDStatus DCDRead_Open      ( CDCDHandle  *self, char *path  )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DCDTrajectoryFileReader:

    cdef CDCDHandle    *cObject
#    cdef public object  currentFrame
    cdef public object  isOpen
    cdef public object  owner
    cdef public object  path
#    cdef public object  title
