#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelFullState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the state for a full NB model."""

from pCore import logFile, LogFileActive

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelFullState:

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            NBModelFullState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModelFullState"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def _Allocate ( self, Integer extent ):
        """Allocation."""
        self.cObject = NBModelFullState_Allocate ( extent )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def GetEnergies ( self, energies ):
        """Append the energies for the state to energies (a simple non-zero check only is used)."""
        if ( self.cObject.emmel     != 0.0 ): energies.append ( ( "MM/MM Elect."     , self.cObject.emmel     ) )
        if ( self.cObject.emmel14   != 0.0 ): energies.append ( ( "MM/MM 1-4 Elect." , self.cObject.emmel14   ) )
        if ( self.cObject.emmlj     != 0.0 ): energies.append ( ( "MM/MM LJ"         , self.cObject.emmlj     ) )
        if ( self.cObject.emmlj14   != 0.0 ): energies.append ( ( "MM/MM 1-4 LJ"     , self.cObject.emmlj14   ) )
        if ( self.cObject.eqcmmlj   != 0.0 ): energies.append ( ( "QC/MM LJ"         , self.cObject.eqcmmlj   ) )
        if ( self.cObject.eqcmmlj14 != 0.0 ): energies.append ( ( "QC/MM 1-4 LJ"     , self.cObject.eqcmmlj14 ) )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ): log.Paragraph ( "Full NB Model State set up." )
