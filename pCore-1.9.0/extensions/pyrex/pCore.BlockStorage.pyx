#-------------------------------------------------------------------------------
# . File      : pCore.BlockStorage.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle blockstorage."""

#===============================================================================
# . Class.
#===============================================================================
cdef class BlockStorage:

    # . Public methods.
    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner: BlockStorage_Deallocate ( &self.cObject )

    def __init__ ( self ):
        """Constructor."""
        self.cObject = NULL
        self.isOwner = True
