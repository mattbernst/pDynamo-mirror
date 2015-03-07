#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelORCAState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the state for a full NB model compatible with the ORCA program."""

from pCore import logFile, LogFileActive, UNITS_LENGTH_ANGSTROMS_TO_BOHRS, UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelORCAState:

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            NBModelORCAState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModelORCAState"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def _Allocate ( self, Integer extent ):
        """Allocation."""
        self.cObject = NBModelORCAState_Allocate ( extent )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def GetEnergies ( self, energies ):
        """Append the energies for the state to energies (a simple non-zero check only is used)."""
        if ( self.cObject.emmel     != 0.0 ): energies.append ( ( "MM/MM Elect.",     self.cObject.emmel     ) )
        if ( self.cObject.emmel14   != 0.0 ): energies.append ( ( "MM/MM 1-4 Elect.", self.cObject.emmel14   ) )
        if ( self.cObject.emmlj     != 0.0 ): energies.append ( ( "MM/MM LJ",         self.cObject.emmlj     ) )
        if ( self.cObject.emmlj14   != 0.0 ): energies.append ( ( "MM/MM 1-4 LJ",     self.cObject.emmlj14   ) )
        if ( self.cObject.eqcmmlj   != 0.0 ): energies.append ( ( "QC/MM LJ",         self.cObject.eqcmmlj   ) )
        if ( self.cObject.eqcmmlj14 != 0.0 ): energies.append ( ( "QC/MM 1-4 LJ",     self.cObject.eqcmmlj14 ) )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def ReadPCgradFile ( self, pcgrad_filename ):
        """Read MM data from an pcgrad file in atomic units and convert to pMolecule units."""
        cdef Real d, g
        cdef Integer    i, j
        if self.cObject.mmCharges != NULL:
            pc_file = open ( pcgrad_filename, 'r' )
            # . Skip the number of point charges.
            next ( pc_file )
            factor = UNITS_LENGTH_ANGSTROMS_TO_BOHRS * UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE
            for i from 0 <= i < self.cObject.mmCharges.length:
                line = next ( pc_file ).split ( )
                for j from 0 <= j < 3:
                    d = float ( line[j] ) * factor
                    g = Coordinates3_GetItem ( self.cObject.mmgradients3, i, j, NULL )
                    g = g + d
                    Coordinates3_SetItem ( self.cObject.mmgradients3, i, j, g, NULL )
            pc_file.close ( )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ): log.Paragraph ( "ORCA NB Model State set up." )

    def WriteInputFile ( self, pc_filename ):
        """Write MM data to an external pointcharge file."""
        cdef Real x, y, z
        cdef Integer    i
        cdef Real   q
        if self.cObject.mmCharges != NULL:
            pc_file = open ( pc_filename, 'w' )
            pc_file.write ( "{:10d}\n".format ( self.cObject.mmCharges.length ) )
            for i from 0 <= i < self.cObject.mmCharges.length:
                q = Real1DArray_GetItem  ( self.cObject.mmCharges, i, NULL )
                x = Coordinates3_GetItem ( self.cObject.mmcoordinates3, i, 0, NULL )
                y = Coordinates3_GetItem ( self.cObject.mmcoordinates3, i, 1, NULL )
                z = Coordinates3_GetItem ( self.cObject.mmcoordinates3, i, 2, NULL )
                pc_file.write ( " {:10.5f}{:20.10f}{:20.10f}{:20.10f}\n".format ( q, x, y, z ) )
            pc_file.close ( )
