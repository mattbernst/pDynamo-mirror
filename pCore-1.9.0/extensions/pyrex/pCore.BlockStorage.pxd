#-------------------------------------------------------------------------------
# . File      : pCore.BlockStorage.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.List         cimport CList

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "BlockStorage.h":

#    ctypedef struct CBlock "Block":
#        Integer    ndata
#        Integer16 *indices16
#        Integer32 *indices32
#        Real      *data

    ctypedef struct CBlockStorage "BlockStorage":
        Boolean  QUNDERFLOW
        Integer  blocksize
        Integer  ndata
        Integer  nindices16
        Integer  nindices32
        Real     underflow
        CList   *blocks

    cdef void BlockStorage_Deallocate ( CBlockStorage **blockstorage )
    cdef Real BlockStorage_Size       ( CBlockStorage  *self         )

#===============================================================================
# . Class.
#===============================================================================
cdef class BlockStorage:

    cdef CBlockStorage *cObject
    cdef public object  isOwner
