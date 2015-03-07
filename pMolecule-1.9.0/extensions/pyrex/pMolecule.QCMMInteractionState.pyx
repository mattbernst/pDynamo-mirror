#-------------------------------------------------------------------------------
# . File      : pMolecule.QCMMInteractionState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A container for quantities needed during calculation of the QC/MM interactions."""

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCMMInteractionState:

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            QCMMInteractionState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.QCMMInteractionState"

    def __init__ ( self, extent, includeQCQC = False ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( extent, includeQCQC = includeQCQC )

    def _Allocate ( self, Integer extent, includeQCQC = False ):
        """Allocation."""
        cdef Boolean cIncludeQCQC
        if includeQCQC: cIncludeQCQC = CTrue
        else:           cIncludeQCQC = CFalse
        self.cObject = QCMMInteractionState_Allocate ( extent, cIncludeQCQC, NULL )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = True

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @classmethod
    def WithExtent ( selfClass, extent, includeQCQC = False ):
        """Constructor."""
        return selfClass ( extent, includeQCQC = includeQCQC )

    # . Properties.
    property qcCharges:
        def __get__ ( self ):
            cdef Real1DArray item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.qcCharges != NULL ):
                item = Real1DArray.Raw ( )
                item.cObject = Real1DArray_Clone ( self.cObject.qcCharges, NULL )
                item.isOwner = True
            return item
    property qcmmPotentials:
        def __get__ ( self ):
            cdef Real1DArray item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.qcmmPotentials != NULL ):
                item = Real1DArray.Raw ( )
                item.cObject = Real1DArray_Clone ( self.cObject.qcmmPotentials, NULL )
                item.isOwner = True
            return item
