#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModel.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the quantum chemical model.

The model is defined so as to be independent of any system data.
"""

import math

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_CardinalOccupancyTolerance = 1.0e-12

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModel:
    """The base class for quantum chemical models.

    This class should not be used directly.
    """

    def __copy__ ( self ):
        """Copying."""
        options = self.__getstate__ ( )
        new     = self.__class__ ( **options )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        pass

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.QCModel"

    def __getstate__ ( self ):
        """Return the state."""
        return { }

    def __init__ ( self, *arguments, **options ):
        """Constructor with options."""
        self._Initialize ( )
        self._Allocate   ( )
        self.ProcessArguments ( *arguments )
        self.SetOptions       ( **options  )
        self.MakeLabel ( )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        self.SetOptions ( **state )

    def _Allocate ( self ):
        """Allocation."""
        pass

    def _Initialize ( self ):
        """Initialization."""
        self.converger = None
        self.label     = None

    def AtomicCharges ( self, configuration, chargeModel = None, spinDensities = False ):
        """Atomic charges."""
        return None

    def BondOrders ( self, configuration, chargeModel = None ):
        """Bond orders."""
        return None

    def Clear ( self, configuration ):
        """Clear up temporary data."""
        if configuration is not None:
            if hasattr ( configuration, "qcState" ): delattr ( configuration, "qcState" )

    def CrossCheckModelAndState ( self, charge, multiplicity, numberElectrons, numberOrbitals, permitFractionalRestrictedNonSinglet = False ):
        """Cross-check the model electronic state options versus those of a system."""
        # . Get basic variables.
        isSpinRestricted = self.GetOption ( "isSpinRestricted" )
        occupancyType    = self.GetOption ( "occupancyType"    )

        # . Initial charge estimation.
        # . Ultimately charge could be non-integer.
        electronicCharge = float ( numberElectrons ) - float ( charge )
        alphaCharge      = ( electronicCharge + ( float ( multiplicity ) - 1.0 ) ) / 2.0
        betaCharge       = ( electronicCharge - ( float ( multiplicity ) - 1.0 ) ) / 2.0

        # . Cardinal occupancy.
        if occupancyType == "Cardinal":
            if ( math.fabs ( round ( electronicCharge ) - electronicCharge ) > _CardinalOccupancyTolerance ): raise ValueError ( "Cardinal occupancy requires integral charge." )
            elif ( isSpinRestricted and ( multiplicity > 1 ) ): raise ValueError ( "A spin-unrestricted calculation is required for non-singlet states with cardinal occupancy." )
        # . Fractional occupancy.
        else:
            if isSpinRestricted:
                if multiplicity > 1:
                    if permitFractionalRestrictedNonSinglet:
                        alphaCharge = electronicCharge / 2.0
                        betaCharge  = electronicCharge / 2.0
                    else: raise ValueError ( "A spin-unrestricted calculation is required for non-singlet states with fractional occupancy." )

        # . Check the alpha and beta charge.
        numberOccupiedAlpha = int ( math.ceil ( alphaCharge ) )
        numberOccupiedBeta  = int ( math.ceil ( betaCharge  ) )
        if ( ( alphaCharge < 0.0 ) or ( numberOccupiedAlpha > numberOrbitals ) or \
             ( betaCharge  < 0.0 ) or ( numberOccupiedBeta  > numberOrbitals ) ): raise ValueError ( "Invalid alpha or beta electronic charge (< 0 or > number of available orbitals)." )

        # . Fixed fractional occupancy.
        if occupancyType == "Fixed Fractional":
            if isSpinRestricted: data = [ ( ""     , numberOccupiedAlpha ) ]
            else:                data = [ ( "Alpha", numberOccupiedAlpha ), ( "Beta", numberOccupiedBeta ) ]
            isOK = True
            for ( tag, numberOccupied ) in data:
                numberFractionalHOOs = self.GetOption ( "numberFractional" + tag + "HOOs" )
                numberFractionalLUOs = self.GetOption ( "numberFractional" + tag + "LUOs" )
                isOK = ( numberFractionalHOOs >= 1 ) and ( numberFractionalHOOs <=                    numberOccupied   ) and \
                       ( numberFractionalLUOs >= 0 ) and ( numberFractionalLUOs <= ( numberOrbitals - numberOccupied ) )
                if not isOK: break
            if not isOK: raise ValueError ( "Invalid fixed fractional occupancy orbital specification." )

        # . Finish up.
        return ( alphaCharge, betaCharge )

    def DefineParameters ( self, anlist, log = logFile ):
        """Define parameters."""
        return None

    def DipoleMoment ( self, configuration, Vector3 center = None ):
        """Dipole Moment."""
        return None

    def Energy ( self, QCAtomContainer qcAtoms, QCParameters qcParameters, electronicState, configuration, log = logFile ):
        """Calculate the quantum chemical energy using a simple iterative procedure."""
        # . Initialization.
        results = []
        # . Start the calculation.
        if ( qcAtoms is not None ) and ( qcParameters is not None ) and ( configuration is not None ) and hasattr ( configuration, "coordinates3" ) and ( configuration.coordinates3 is not None ):
            # . Check for a state which is assumed to be OK if present.
            if not hasattr ( configuration, "qcState" ): self.SetUp ( qcAtoms, qcParameters, electronicState, configuration )
            qcState = configuration.qcState
            # . Initialization - enuclear and integrals.
            self.EnergyInitialize ( configuration )
            # . Get the converger.
            if self.converger is None: self.converger = DIISSCFConverger ( )
            converger = self.converger
            # . Initialization.
            convergerState = converger.Initialize ( qcState.densityp, qcState.densityq, qcState.overlap, qcState.orthogonalizer, log = log )
            # . Iterate.
            while converger.Continue ( convergerState ):
                # . Build the Fock matrices and calculate the electronic energy.
                eelectronic = self.EnergyFock ( configuration )
                # . Apply the converger.
                converger.Iterate ( convergerState, eelectronic )
            # . Finalization.
            ( hasConverged, ncycles ) = converger.Finalize ( convergerState )
            # . Set convergence data.
            qcState.SetConvergenceData ( hasConverged, ncycles )
            # . Post-processing.
            if hasConverged: self.EnergyPostProcessing ( configuration )
            # . Print a summary of the results.
            qcState.Summary ( log )
            # . Throw an error if the SCF procedure has not converged.
            if not hasConverged: raise ValueError ( "SCF procedure not converged." )
            # . Get the energies to be returned.
            results = qcState.GetEnergies ( )
            # . Calculate the gradients.
            self.EnergyGradients ( configuration )
            # . Finish up.
            self.EnergyFinalize  ( configuration )
        # . Return the results.
        return results

    def EnergyFinalize ( self, configuration ):
        """Finalization."""
        pass

    def EnergyFock ( self, configuration ):
        """Fock matrix construction."""
        return 0.0

    def EnergyGradients ( self, configuration ):
        """Gradients."""
        pass

    def EnergyInitialize ( self, configuration ):
        """Integrals and other setup."""
        pass

    def EnergyPostProcessing ( self, configuration ):
        """Post-processing operations after the energy but before the gradients."""
        pass

    # . Test method.
    def EnergyWithTimings ( self, QCAtomContainer qcAtoms, QCParameters qcParameters, electronicState, configuration, cpuTimer, timings, log = logFile ):
        """Calculate the quantum chemical energy using a simple iterative procedure."""
        # . Initialization.
        results = []
        # . Start the calculation.
        if ( qcAtoms is not None ) and ( qcParameters is not None ) and ( configuration is not None ) and hasattr ( configuration, "coordinates3" ) and ( configuration.coordinates3 is not None ):
            tStart = cpuTimer.Current ( )
            # . Check for a state which is assumed to be OK if present.
            if not hasattr ( configuration, "qcState" ):
                self.SetUp ( qcAtoms, qcParameters, electronicState, configuration )
                tStop = cpuTimer.Current ( )
                timings["QC Set Up"] += ( tStop - tStart )
                tStart = tStop
            qcState = configuration.qcState
            # . Initialization - enuclear and integrals.
            self.EnergyInitialize ( configuration )
            tStop = cpuTimer.Current ( )
            timings["QC Initialize"] += ( tStop - tStart )
            tStart = tStop
            # . Get the converger.
            if self.converger is None: self.converger = DIISSCFConverger ( )
            converger = self.converger
            # . Initialization.
            convergerState = converger.Initialize ( qcState.densityp, qcState.densityq, qcState.overlap, qcState.orthogonalizer, log = log )
            tStop = cpuTimer.Current ( )
            timings["QC Converger Initialize"] += ( tStop - tStart )
            tStart = tStop
            # . Iterate.
            while converger.Continue ( convergerState ):
                # . Build the Fock matrices and calculate the electronic energy.
                eelectronic = self.EnergyFock ( configuration )
                tStop = cpuTimer.Current ( )
                timings["QC Fock"] += ( tStop - tStart )
                # . Apply the converger.
                converger.IterateWithTimings ( convergerState, eelectronic, cpuTimer, timings )
                tStart = cpuTimer.Current ( )
            # . Finalization.
            ( hasConverged, ncycles ) = converger.Finalize ( convergerState )
            # . Set convergence data.
            qcState.SetConvergenceData ( hasConverged, ncycles )
            tStop = cpuTimer.Current ( )
            timings["QC Converger Finalize"] += ( tStop - tStart )
            tStart = tStop
            # . Post-processing.
            if hasConverged:
                self.EnergyPostProcessing ( configuration )
                tStop = cpuTimer.Current ( )
                timings["QC Postprocessing"] += ( tStop - tStart )
                tStart = tStop
            # . Print a summary of the results.
            qcState.Summary ( log )
            # . Throw an error if the SCF procedure has not converged.
            if not hasConverged: raise ValueError ( "SCF procedure not converged." )
            # . Get the energies to be returned.
            results = qcState.GetEnergies ( )
            # . Calculate the gradients.
            tStart = cpuTimer.Current ( )
            self.EnergyGradients ( configuration )
            tStop  = cpuTimer.Current ( )
            timings["QC Gradients"] += ( tStop - tStart )
            tStart = tStop
            # . Finish up.
            self.EnergyFinalize  ( configuration )
            tStop = cpuTimer.Current ( )
            timings["QC Finalization"] += ( tStop - tStart )
            tStart = tStop
        # . Return the results.
        return results

    @classmethod
    def FromOptions ( selfClass, **options ):
        """Constructor from options."""
        return selfClass ( **options )

    def GetOption ( self, optionName ):
        """Get the value of an option."""
        return None

    def GridPointDensities ( self, configuration, Coordinates3 gridpoints, spinDensities = False ):
        """Densities at grid points."""
        return None

    def GridPointOrbitals ( self, configuration, Coordinates3 gridpoints, orbitals = None, useDensityP = True ):
        """Orbitals at grid points."""
        return None

    def GridPointPotentials ( self, configuration, Coordinates3 gridpoints ):
        """Electrostatic potentials at grid points."""
        return None

    def IsSpinRestricted ( self ):
        """Is the model spin restricted?"""
        return False

    def MakeLabel ( self ):
        """Construct a model label."""
        if ( self.label is not None ) and ( not isinstance ( self.label, basestring ) ): self.label = None

    def OrbitalEnergies ( self, configuration, useDensityP = True ):
        """Get the orbital energies and HOMO and LUMO indices."""
        return ( None, None, None )

    def ProcessArguments ( self, *arguments ):
        """Process constructor arguments."""
        pass

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        pass

    def SetUp ( self, QCAtomContainer qcAtoms, QCParameters qcParameters, electronicState, configuration ):
        """Set up a state for the model."""
        pass

    def Summary ( self, log = logFile ):
        """Summary."""
        pass
