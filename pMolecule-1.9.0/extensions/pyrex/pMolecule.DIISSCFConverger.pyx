#-------------------------------------------------------------------------------
# . File      : pMolecule.DIISSCFConverger.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines classes for the DIIS SCF converger."""

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Attributes.
_defaultAttributes = { "densityTolerance" : 1.0e-8  ,
                       "diisOnset"        : 0.2e+00 ,
                       "energyTolerance"  : 2.0e-4  ,
                       "history"          : 10      ,
                       "maximumSCFCycles" : 100     ,
                       "useRCA"           : True    }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DIISSCFConverger:
    """Class for the DIIS SCF converger."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            DIISSCFConverger_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.DIISSCFConverger"

    def __getstate__ ( self ):
        """Return the state."""
        return { "densityTolerance" : self.cObject.densityTolerance ,
                 "diisOnset"        : self.cObject.diisonset        ,
                 "energyTolerance"  : self.cObject.energytolerance  ,
                 "history"          : self.cObject.ndiis            ,
                 "maximumSCFCycles" : self.cObject.maximumSCFCycles ,
                 "useRCA"           : ( self.cObject.QUSERCA == CTrue ) }

    def __init__ ( self, **options ):
        """Constructor with options."""
        self._Initialize ( )
        self._Allocate   ( )
        self.SetOptions  ( **options  )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        self.SetOptions ( **state )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = DIISSCFConverger_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def Continue ( self, DIISSCFConvergerState convergerstate ):
        """Check to see if the calculation should continue."""
        return ( DIISSCFConverger_Continue ( self.cObject, convergerstate.cObject ) == CTrue )

    def Finalize ( self, DIISSCFConvergerState convergerstate ):
        """Finalize a calculation."""
        cdef Integer iteration
        iteration  = 0
        isConverged = False
        if convergerstate is not None:
            convergerstate.LogStop ( )
            iteration  = convergerstate.cObject.iteration
            isConverged = ( convergerstate.cObject.isConverged == CTrue )
            DIISSCFConvergerState_Deallocate ( &convergerstate.cObject )
        return ( isConverged, iteration )

    @classmethod
    def FromOptions ( selfClass, **options ):
        """Constructor from options."""
        return selfClass ( **options )

    def Initialize ( self, QCOnePDM densityp, QCOnePDM densityq, SymmetricMatrix overlap, Real2DArray orthogonalizer, log = logFile ):
        """Initialization."""
        cdef DIISSCFConvergerState  convergerstate
        cdef CQCOnePDM             *cdensityp, *cdensityq
        cdef CReal2DArray          *corthogonalizer
        cdef CSymmetricMatrix  *coverlap
        # . Get pointers.
        if densityp       is None: cdensityp        = NULL
        else:                      cdensityp        = densityp.cObject
        if densityq       is None: cdensityq        = NULL
        else:                      cdensityq        = densityq.cObject
        if overlap        is None: coverlap         = NULL
        else:                      coverlap         = overlap.cObject
        if orthogonalizer is None: corthogonalizer  = NULL
        else:                      corthogonalizer  = orthogonalizer.cObject
        # . Allocate the state.
        convergerstate         = DIISSCFConvergerState.Raw ( )
        convergerstate.cObject = DIISSCFConvergerState_SetUp ( self.cObject.QUSERCA, self.cObject.ndiis, cdensityp, cdensityq, coverlap, corthogonalizer )
        convergerstate.isOwner = True
        # . Logging.
        self.Summary            ( log = log )
        convergerstate.LogStart ( log = log )
        # . Return the object.
        return convergerstate

    def Iterate ( self, DIISSCFConvergerState convergerstate, Real ecurrent ):
        """Perform an iteration."""
        if convergerstate is not None:
            convergerstate.LogIteration ( ecurrent )
            if DIISSCFConverger_Converged ( self.cObject, convergerstate.cObject ) != CTrue:
                DIISSCFConverger_IterateWithoutDensities ( self.cObject, convergerstate.cObject, ecurrent )
                DIISSCFConverger_MakeDensities           ( self.cObject, convergerstate.cObject           )

    def IterateWithTimings ( self, DIISSCFConvergerState convergerstate, Real ecurrent, cpuTimer, timings ):
        """Perform an iteration."""
        if convergerstate is not None:
            convergerstate.LogIteration ( ecurrent )
            if DIISSCFConverger_Converged ( self.cObject, convergerstate.cObject ) != CTrue:
                tStart = cpuTimer.Current ( )
                DIISSCFConverger_IterateWithoutDensities ( self.cObject, convergerstate.cObject, ecurrent )
                tStop  = cpuTimer.Current ( )
                timings["QC Converger Iterate"  ] += ( tStop - tStart )
                tStart = tStop
                DIISSCFConverger_MakeDensities ( self.cObject, convergerstate.cObject )
                tStop  = cpuTimer.Current ( )
                timings["QC Converger Densities"] += ( tStop - tStart )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        cdef DIISSCFConverger self
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        if "densityTolerance" in keywordArguments: self.cObject.densityTolerance = keywordArguments.pop ( "densityTolerance" )
        if "diisOnset"        in keywordArguments: self.cObject.diisonset        = keywordArguments.pop ( "diisOnset"        )
        if "energyTolerance"  in keywordArguments: self.cObject.energytolerance  = keywordArguments.pop ( "energyTolerance"  )
        if "history"          in keywordArguments: self.cObject.ndiis            = keywordArguments.pop ( "history"          )
        if "maximumSCFCycles" in keywordArguments: self.cObject.maximumSCFCycles = keywordArguments.pop ( "maximumSCFCycles" )
        if "useRCA"           in keywordArguments:
            if keywordArguments.pop ( "useRCA" ): self.cObject.QUSERCA = CTrue
            else:                                 self.cObject.QUSERCA = CFalse
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            if self.cObject.QUSERCA == CTrue: summary.Start ( "DIIS/ODA SCF Converger Summary" )
            else:                            summary.Start ( "DIIS/Damping SCF Converger Summary" )
            summary.Entry ( "Maximum Cycles"    , "{:d}".format ( self.cObject.maximumSCFCycles ) )
            summary.Entry ( "DIIS Matrices"     , "{:d}".format ( self.cObject.ndiis            ) )
            summary.Entry ( "Density Tolerance" , "{:14.8g}".format ( self.cObject.densityTolerance ) )
            summary.Entry ( "Energy Tolerance"  , "{:14.8g}".format ( self.cObject.energytolerance  ) )
            summary.Entry ( "DIIS Onset"        , "{:14.8g}".format ( self.cObject.diisonset        ) )
            summary.Stop ( )
