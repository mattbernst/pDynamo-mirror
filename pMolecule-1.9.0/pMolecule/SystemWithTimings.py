#-------------------------------------------------------------------------------
# . File      : SystemWithTimings.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""System class with energy timings."""

from pCore                      import Coordinates3, CPUTime, logFile, LogFileActive, Timings
from EnergyTerms                import EnergyTerms
from SymmetryParameterGradients import SymmetryParameterGradients
from System                     import System

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SystemWithTimings ( System ):
    """A system with a timing analysis of its energy calculation."""

    def Energy ( self, log = logFile, doGradients = False ):
        """Calculate the energy and, optionally, the gradients for a system."""
        # . Initialization.
        tStartEnergy = self.cpuTimer.Current ( )
        # . Do nothing if coordinates and an energy model are absent.
        if ( self.coordinates3 is None ) or ( self.energyModel is None ): return
        # . Set up configuration.
        self.configuration.ClearTemporaryAttributes ( )
        self.configuration.SetTemporaryAttribute ( "energyTerms", EnergyTerms ( ) )
        # . Make this clearer here (no need to reallocate each time).
        if doGradients:
            g = Coordinates3.WithExtent ( len ( self.atoms ) )
            g.Set ( 0.0 )
            self.configuration.SetTemporaryAttribute ( "gradients3", g )
            if self.symmetry is not None:
                self.configuration.SetTemporaryAttribute ( "symmetryParameterGradients", SymmetryParameterGradients ( ) )
        # . Alias various data structures.
        em  = self.energyModel
        et  = self.configuration.energyTerms
        fa  = None
        crd = self.configuration.coordinates3
        grd = getattr ( self.configuration, "gradients3", None )
        if self.hardConstraints is not None: fa = self.hardConstraints.fixedAtoms
        self.timings["Initialization"] += ( self.cpuTimer.Current ( ) - tStartEnergy )
        # . Calculate MM terms.
        if ( em.mmModel is not None ) and ( em.mmTerms is not None ):
            for mmterm in em.mmTerms:
                tStart  = self.cpuTimer.Current ( )
                results = mmterm.Energy ( crd, grd )
                et.Append ( results )
                self.timings[results[0]] += ( self.cpuTimer.Current ( ) - tStart )
        # . Calculate NB terms.
        if em.nbModel is not None:
            tStart  = self.cpuTimer.Current ( )
            em.nbModel.SetUp ( em.mmAtoms, em.qcAtoms, em.ljParameters, em.ljParameters14, fa, em.interactions14, em.exclusions, self.symmetry, self.connectivity.isolates, self.configuration, log = log )
            tStop   = self.cpuTimer.Current ( )
            self.timings["NB Set Up"    ] += ( tStop - tStart )
            et.Extend ( em.nbModel.Energy ( self.configuration ) )
            self.timings["NB Evaluation"] += ( self.cpuTimer.Current ( ) - tStop )
        # . Calculate QC and QC/MM electrostatic terms.
        if em.qcModel is not None:
            tStart = self.cpuTimer.Current ( )
            if ( em.nbModel is not None ):
                em.nbModel.QCMMPotentials ( self.configuration )
                tStop  = self.cpuTimer.Current ( )
                self.timings["QC/MM Potentials"] += ( tStop - tStart )
                tStart = tStop
            et.Extend ( em.qcModel.EnergyWithTimings ( em.qcAtoms, em.qcParameters, self.electronicState, self.configuration, self.cpuTimer, self.timings, log = log ) )
            tStop = self.cpuTimer.Current ( )
            self.timings["QC Evaluation"] += ( tStop - tStart )
            if ( em.nbModel is not None ) and doGradients:
                tStart = tStop
                em.nbModel.QCMMGradients  ( self.configuration )
                self.timings["QC/MM Gradients"] += ( self.cpuTimer.Current ( ) - tStart )
        # . Calculate the soft constraint terms.
        if em.softConstraints is not None:
            tStart = self.cpuTimer.Current ( )
            ( scEnergy, scState ) = em.softConstraints.Energy ( crd, grd )
            self.configuration.SetTemporaryAttribute ( "scState", scState )
            et.Append ( scEnergy )
            self.timings["Soft Constraint"] += ( self.cpuTimer.Current ( ) - tStart )
        # . Calculate the total energy.
        tStart = self.cpuTimer.Current ( )
        et.Total ( )
        # . Zero out unnecessary gradient terms.
        if ( fa is not None ) and doGradients: grd.SetRowSelection ( fa, 0.0 )
        # . Print the energy terms if necessary.
        if LogFileActive ( log ):
            if doGradients: et.RMSGradient ( grd.RMSValue ( ) )
            et.Summary ( log = log )
        # . Finish up.
        tStop = self.cpuTimer.Current ( )
        self.timings["Energy"      ] += ( tStop - tStartEnergy )
        self.timings["Finalization"] += ( tStop - tStart       )
        self.numberOfEnergyCalls     += 1
        return et.potentialEnergy

    @classmethod
    def FromSystem ( selfClass, system ):
        """Constructor from system."""
        self = selfClass ( )
        self.__dict__.update ( system.__dict__ )
        return self

    def TimingStart ( self ):
        """Start timing."""
        self.cpuTimer            = CPUTime ( )
        self.numberOfEnergyCalls = 0
        self.timings             = Timings ( )
        self.timings.Clear ( )

    def TimingStop ( self ):
        """Stop timing."""
        self.timings["Total"] += self.cpuTimer.Current ( )

    def TimingSummary ( self, log = logFile, orderByMagnitude = False ):
        """Timing summary."""
        self.TimingStop ( )
        self.timings.Summary ( log = log, orderByMagnitude = orderByMagnitude )
        if LogFileActive ( log ): log.Paragraph ( "Number of Energy Calls = {:d}".format ( self.numberOfEnergyCalls ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
