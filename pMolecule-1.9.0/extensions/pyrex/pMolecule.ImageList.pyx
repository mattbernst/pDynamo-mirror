#-------------------------------------------------------------------------------
# . File      : pMolecule.ImageList.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Interface to the image list C type."""

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ImageList:
    """Define an image list."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            ImageList_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.ImageList"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate   ( )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = ImageList_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    @classmethod
    def Raw ( selfClass ):
        """Constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( )

    # . Properties.
    property numberOfImages:
        def __get__ ( self ):
            return ImageList_NumberOfImages ( self.cObject )

    property numberOfPairs:
        def __get__ ( self ):
            return ImageList_NumberOfPairs ( self.cObject )
