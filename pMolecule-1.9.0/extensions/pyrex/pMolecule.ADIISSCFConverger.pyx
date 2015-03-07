#-------------------------------------------------------------------------------
# . File      : pMolecule.ADIISSCFConverger.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines classes for the ADIIS SCF converger."""

import os, os.path

from pCore import LBFGSMinimizer, logFile, LogFileActive, RawObjectConstructor, TextLogFileWriter

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Attributes.
#_defaultAttributes = { "densityTolerance"   : 1.0e-8 , \
#                       "diisOff"            : 0.8e+00, \
#                       "diisOn"             : 0.2e+00, \
#                       "energyTolerance"    : 2.0e-4 , \
#                       "maximumHistory"     : 10     , \
#                       "maximumSCFCycles"   : 100    , \
#                       "minimumCoefficient" : 1.0e-2 , \
#                       "optimizer"          : None   , \
#                       "useEDIIS"           : False  , \
#                       "useODA"             : False    }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ADIISSCFConverger:
    """Class for the ADIIS SCF converger."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.debugLog is not None: self.debugLog.Close ( )
        if self.isOwner:
            ADIISSCFConverger_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.ADIISSCFConverger"

    def __getstate__ ( self ):
        """Return the state."""
        return { "densityTolerance"      : self.cObject.densityTolerance      ,
                 "diisOff"               : self.cObject.diisOff               ,
                 "diisOn"                : self.cObject.diisOn                ,
                 "energyTolerance"       : self.cObject.energyTolerance       ,
                 "maximumHistory"        : self.cObject.maximumHistory        ,
                 "maximumSCFCycles"      : self.cObject.maximumSCFCycles      ,
                 "minimizationsPerCycle" : self.minimizationsPerCycle         ,
                 "minimumCoefficient"    : self.cObject.minimumCoefficient    ,
                 "optimizer"             : self.optimizer                     ,
                 "useEDIIS"              : ( self.cObject.useEDIIS == CTrue ) ,
                 "useODA"                : ( self.cObject.useODA   == CTrue ) }

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
        self.cObject = ADIISSCFConverger_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject               = NULL
        self.debugLog              = None
        self.isOwner               = False
        self.minimizationsPerCycle = 1
        self.optimizer             = None

    def Continue ( self, ADIISSCFConvergerState state ):
        """Check to see if the calculation should continue."""
        return ( ADIISSCFConverger_Continue ( self.cObject, state.cObject ) == CTrue )

    def DefineDefaultOptimizer ( self ):
        """Define a default optimizer."""
        if ( self.cObject.useODA == CFalse ) and ( self.optimizer is None ):
            self.optimizer = LBFGSMinimizer ( history              =  10    ,
                                              logFrequency         =   1    ,
                                              maximumIterations    = 100    ,
                                              rmsGradientTolerance = 1.0e-3 )
        # . Debugging.
#        self.debugLog = TextLogFileWriter ( os.path.join ( os.getenv ( "PDYNAMO_SCRATCH" ), "ADIIS.log" ) )

    def Finalize ( self, ADIISSCFConvergerState state ):
        """Finalize a calculation."""
        cdef Integer iteration
        iteration   = 0
        isConverged = False
        if state is not None:
            state.LogStop ( )
            iteration   = state.cObject.iteration
            isConverged = ( state.cObject.isConverged == CTrue )
            ADIISSCFConvergerState_Deallocate ( &state.cObject )
        return ( isConverged, iteration )

    @classmethod
    def FromOptions ( selfClass, **options ):
        """Constructor from options."""
        return selfClass ( **options )

    def Initialize ( self, QCOnePDM densityP, QCOnePDM densityQ, SymmetricMatrix overlap, Real2DArray orthogonalizer, log = logFile ):
        """Initialization."""
        cdef ADIISSCFConvergerState  state
        cdef CQCOnePDM              *cDensityP, *cDensityQ
        cdef CReal2DArray           *cOrthogonalizer
        cdef CSymmetricMatrix       *cOverlap
        # . Define a default optimizer if one is not already present.
        self.DefineDefaultOptimizer ( )
        # . Get a state.
        # . Get pointers.
        if densityP       is None: cDensityP       = NULL
        else:                      cDensityP       = densityP.cObject
        if densityQ       is None: cDensityQ       = NULL
        else:                      cDensityQ       = densityQ.cObject
        if overlap        is None: cOverlap        = NULL
        else:                      cOverlap        = overlap.cObject
        if orthogonalizer is None: cOrthogonalizer = NULL
        else:                      cOrthogonalizer = orthogonalizer.cObject
        # . Allocate the state.
        state         = ADIISSCFConvergerState.Raw ( )
        state.cObject = ADIISSCFConvergerState_SetUp ( self.cObject.useODA, self.cObject.maximumHistory, cDensityP, cDensityQ, cOverlap, cOrthogonalizer, NULL )
        state.isOwner = True
        # . Logging.
        self.Summary   ( log = log )
        state.LogStart ( log = log )
        # . Return the object.
        return state

    def Iterate ( self, ADIISSCFConvergerState state, Real eCurrent ):
        """Perform an iteration."""
        if state is not None:
            state.LogIteration ( eCurrent )
            if ADIISSCFConverger_IsConverged ( self.cObject, state.cObject ) != CTrue:
                ADIISSCFConverger_IterateStart ( self.cObject, state.cObject, eCurrent )
                state.SolveADIISEquations ( self.optimizer, self.cObject.minimumCoefficient, self.minimizationsPerCycle, log = self.debugLog )
                ADIISSCFConverger_IterateStop  ( self.cObject, state.cObject )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        cdef ADIISSCFConverger self
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        if "densityTolerance"      in keywordArguments: self.cObject.densityTolerance      = keywordArguments.pop ( "densityTolerance"      )
        if "diisOff"               in keywordArguments: self.cObject.diisOff               = keywordArguments.pop ( "diisOff"               )
        if "diisOn"                in keywordArguments: self.cObject.diisOn                = keywordArguments.pop ( "diisOn"                )
        if "energyTolerance"       in keywordArguments: self.cObject.energyTolerance       = keywordArguments.pop ( "energyTolerance"       )
        if "maximumHistory"        in keywordArguments: self.cObject.maximumHistory        = keywordArguments.pop ( "maximumHistory"        )
        if "maximumSCFCycles"      in keywordArguments: self.cObject.maximumSCFCycles      = keywordArguments.pop ( "maximumSCFCycles"      )
        if "minimizationsPerCycle" in keywordArguments: self.minimizationsPerCycle         = keywordArguments.pop ( "minimizationsPerCycle" )
        if "minimumCoefficient"    in keywordArguments: self.cObject.minimumCoefficient    = keywordArguments.pop ( "minimumCoefficient"    )
        if "optimizer"             in keywordArguments: self.optimizer                     = keywordArguments.pop ( "optimizer"             )
        if "useEDIIS"              in keywordArguments:
            if keywordArguments.pop ( "useEDIIS" ): self.cObject.useEDIIS = CTrue
            else:                                   self.cObject.useEDIIS = CFalse
        if "useODA"                in keywordArguments:
            if keywordArguments.pop ( "useODA" ): self.cObject.useODA = CTrue
            else:                                 self.cObject.useODA = CFalse
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            if   ( self.cObject.useODA   == CTrue ): summary.Start ( "DIIS/ODA SCF Converger Summary"   )
            elif ( self.cObject.useEDIIS == CTrue ): summary.Start ( "EDIIS/DIIS SCF Converger Summary" )
            else:                                    summary.Start ( "ADIIS/DIIS SCF Converger Summary" )
            summary.Entry ( "Maximum Cycles"     , "{:d}"    .format ( self.cObject.maximumSCFCycles   ) )
            summary.Entry ( "Maximum History"    , "{:d}"    .format ( self.cObject.maximumHistory     ) )
            summary.Entry ( "Density Tolerance"  , "{:14.8g}".format ( self.cObject.densityTolerance   ) )
            summary.Entry ( "Energy Tolerance"   , "{:14.8g}".format ( self.cObject.energyTolerance    ) )
            summary.Entry ( "DIIS Off"           , "{:14.8g}".format ( self.cObject.diisOff            ) )
            summary.Entry ( "DIIS On"            , "{:14.8g}".format ( self.cObject.diisOn             ) )
            summary.Entry ( "Minimum Coefficient", "{:14.8g}".format ( self.cObject.minimumCoefficient ) )
            if ( self.cObject.useODA == CFalse ):
                summary.Entry ( "Minimizations/Cycle", "{:d}".format ( self.minimizationsPerCycle ) )
            summary.Stop ( )
