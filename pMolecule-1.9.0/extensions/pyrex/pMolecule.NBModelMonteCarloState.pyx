#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelMonteCarloState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the state for a Monte Carlo NB model."""

from pCore import logFile, LogFileActive

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelMonteCarloState:

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            NBModelMonteCarloState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModelMonteCarloState"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def _Allocate ( self, Integer extent ):
        """Allocation."""
        self.cObject = NBModelMonteCarloState_Allocate ( extent )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def GetEnergies ( self, energies ):
        """Append the energies for the state to energies (a simple non-zero check only is used)."""
        if ( self.cObject.efmmel != 0.0 ): energies.append ( ( "MM/MM Elect.", self.cObject.efmmel ) )
        if ( self.cObject.efmmlj != 0.0 ): energies.append ( ( "MM/MM LJ",     self.cObject.efmmlj ) )

    def GetInteractionEnergies ( self ):
        """Return the interaction energies for the state."""
        return ( ( "MM/MM Elect.", self.cObject.e1mmel ), ( "MM/MM LJ", self.cObject.e1mmlj ) )

    def Initialize ( self, configuration ):
        """Initialize the state."""
        cdef Coordinates3              coordinates3
        cdef SymmetryParameters        symmetryParameters
        cdef CCoordinates3        *ccoordinates3
        cdef CSymmetryParameters  *csymmetryParameters
        if configuration is not None:
            if hasattr ( configuration, "coordinates3"       ): coordinates3       = configuration.coordinates3
            else:                                               coordinates3       = None
            if hasattr ( configuration, "symmetryParameters" ): symmetryParameters = configuration.symmetryParameters
            else:                                               symmetryParameters = None
            if coordinates3       is None: ccoordinates3       = NULL
            else:                          ccoordinates3       = coordinates3.cObject
            if symmetryParameters is None: csymmetryParameters = NULL
            else:                          csymmetryParameters = symmetryParameters.cObject
            NBModelMonteCarloState_Initialize ( self.cObject, ccoordinates3, csymmetryParameters )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def ScaleIsolateInteractionParameters ( self, Integer isolate, Real chargeScale, Real epsilonScale, Real sigmaScale, log = logFile ):
        """Scale the interaction parameters for an isolate."""
        if ( isolate < 0 ) or ( isolate >= self.cObject.nisolates ): raise IndexError ( "Isolate index ({:d}) out of range [0,{:d}).".format ( isolate, self.cObject.nisolates ) )
        NBModelMonteCarloState_ScaleIsolateInteractionParameters ( self.cObject, isolate, chargeScale, epsilonScale, sigmaScale )
        if LogFileActive ( log ): log.Paragraph ( "Monte Carlo isolate interaction parameters scaled." )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ): log.Paragraph ( "Monte Carlo NB Model State set up." )
