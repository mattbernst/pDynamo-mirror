"""DCD handling functions."""

from pCore import CLibraryError

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Default message.
_DefaultMessage = "DCD trajectory file error."

# . Defined messages.
_Messages = { CDCDStatus_AtomNumberMismatch      : "Atom number mismatch."           ,
              CDCDStatus_BadFormat               : "Bad format."                     ,
              CDCDStatus_BadRead                 : "Read error."                     ,
              CDCDStatus_BadSeek                 : "File positioning error."         ,
              CDCDStatus_BadWrite                : "Write error."                    ,
              CDCDStatus_FileAccessFailure       : "Unable to access file."          ,
              CDCDStatus_InvalidDataObject       : "Invalid or missing data object." ,
              CDCDStatus_InvalidFrameIndex       : "Invalid frame  index."           ,
              CDCDStatus_MemoryAllocationFailure : "Memory allocation failure."      ,
              CDCDStatus_OpenFailed              : "File opening error."             }

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
cdef DCDStatus_Check ( CDCDStatus status ):
    """Check the status flag and raise an error if it is not normal."""
    if status != CDCDStatus_Normal:
        raise CLibraryError ( _Messages.get ( status, _DefaultMessage ) )
