from pBabel.DCDHandle             cimport CDCDHandle                         , \
                                          CDCDStatus                         , \
                                          CDCDStatus_MemoryAllocationFailure , \
                                          DCDHandle_Allocate                 , \
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
cdef extern from "DCDWrite.h":

    cdef void       DCDWrite_Close  ( CDCDHandle **self )
    cdef CDCDStatus DCDWrite_Frame  ( CDCDHandle  *self )
    cdef void       DCDWrite_Header ( CDCDHandle  *self, char *title )
    cdef CDCDStatus DCDWrite_Open   ( CDCDHandle  *self, char *path  )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DCDTrajectoryFileWriter:

    cdef CDCDHandle    *cObject
    cdef public object  isOpen
    cdef public object  owner
    cdef public object  path
    cdef public object  title
