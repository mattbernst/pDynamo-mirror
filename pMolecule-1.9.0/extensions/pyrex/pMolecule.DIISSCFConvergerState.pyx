#-------------------------------------------------------------------------------
# . File      : pMolecule.DIISSCFConvergerState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines classes for the DIIS SCF converger state."""

from pCore import logFile, LogFileActive

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DIISSCFConvergerState:
    """Class for the DIIS SCF converger state."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            DIISSCFConvergerState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.DIISSCFConvergerState"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = DIISSCFConvergerState_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def LogIteration ( self, energy ):
        """Log an iteration."""
        if self.table is not None:
            self.table.Entry ( "{:d}".format ( self.cObject.iteration ) )
            self.table.Entry ( "{:20.8f}".format ( energy                     ) )
            self.table.Entry ( "{:20.8f}".format ( self.cObject.deltaeold     ) )
            self.table.Entry ( "{:20.8f}".format ( self.cObject.rmsdifference ) )
            self.table.Entry ( "{:20.8f}".format ( self.cObject.diiserror     ) )
            if self.cObject.iteration > 0:
                if self.cObject.QDIIS == CTrue:
                    self.table.Entry ( "DIIS" )
                    self.table.Entry ( "{:d}".format ( self.cObject.matnum ) )
                elif self.cObject.QRCA == CTrue:
                    self.table.Entry ( "ODA" )
                    self.table.Entry ( "{:12.4g}".format ( self.cObject.rcamu ) )
                elif self.cObject.QDAMP == CTrue:
                    self.table.Entry ( "Damping" )
                    self.table.Entry ( "{:12.4g}".format ( self.cObject.damp  ) )
            else:
                self.table.EndRow ( )

    def LogStart ( self, log = logFile ):
        """Start logging."""
        if LogFileActive ( log ):
            table = log.GetTable ( columns = [ 10, 20, 20, 20, 20, 8, 12 ] )
            table.Start   ( )
            table.Heading ( "Cycle"              )
            table.Heading ( "Energy"             )
            table.Heading ( "Energy Change"      )
            table.Heading ( "RMS Density Change" )
            table.Heading ( "Max. DIIS Error"    )
            table.Heading ( "Operation", columnSpan = 2 )
            self.table = table

    def LogStop ( self ):
        """Stop logging."""
        if self.table is not None:
            self.table.Stop ( )
            self.table = None

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self
