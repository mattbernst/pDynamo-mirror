#-------------------------------------------------------------------------------
# . File      : pMolecule.ADIISSCFConvergerState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines classes for the ADIIS SCF converger state."""

from pCore import logFile, LogFileActive

from ADIISObjectiveFunction import ADIISObjectiveFunction

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ADIISSCFConvergerState:
    """Class for the ADIIS SCF converger state."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            ADIISSCFConvergerState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.ADIISSCFConvergerState"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = ADIISSCFConvergerState_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = True
        self.table   = None

    def LogIteration ( self, energy ):
        """Log an iteration."""
        if self.table is not None:
            if ( self.cObject.iteration == 0 ) and ( self.cObject.densityP.isValid == CFalse ):
                self.table.Entry ( "Guess", alignment = "right" )
            else:
                self.table.Entry ( "{:d}".format ( self.cObject.iteration     ) )
            self.table.Entry ( "{:20.8f}".format ( energy                     ) )
            self.table.Entry ( "{:20.8f}".format ( self.cObject.energyChange  ) )
            self.table.Entry ( "{:20.8f}".format ( self.cObject.rmsDifference ) )
            self.table.Entry ( "{:20.8f}".format ( self.cObject.diisError     ) )
            if self.cObject.iterationType == ADIISSCFConvergerIterationType_ADIIS:
                self.table.Entry ( "ADIIS" )
                self.table.Entry ( "{:d}".format ( self.cObject.history ) )
            elif self.cObject.iterationType == ADIISSCFConvergerIterationType_DIIS:
                self.table.Entry ( "DIIS" )
                self.table.Entry ( "{:d}".format ( self.cObject.diisActive ) )
            elif self.cObject.iterationType == ADIISSCFConvergerIterationType_ODA:
                self.table.Entry ( "ODA" )
                self.table.Entry ( "{:12.4g}".format ( self.cObject.odaMu ) )
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

    def SolveADIISEquations ( self, optimizer, minimumCoefficient, numberOfMinimizations, log = None ):
        """Solve the ADIIS equations."""
        cdef Real1DArray     alphas, alphaGradients
        cdef SymmetricMatrix alphaHessian
        if ( self.cObject.iteration > 0 ) and ( self.cObject.iterationType == ADIISSCFConvergerIterationType_ADIIS ):
            # . Set up the objective function.
            alphas                 = Real1DArray.Raw     ( )
            alphaGradients         = Real1DArray.Raw     ( )
            alphaHessian           = SymmetricMatrix.Raw ( )
            alphas.cObject         = self.cObject.adiisAlphas    ; alphas.isOwner         = False
            alphaGradients.cObject = self.cObject.adiisGradients ; alphaGradients.isOwner = False
            alphaHessian.cObject   = self.cObject.adiisHessian   ; alphaHessian.isOwner   = False
            if len ( alphas ) > self.cObject.history: numberOfSets = 2
            else:                                     numberOfSets = 1
            of = ADIISObjectiveFunction.FromOptions ( alphas, alphaGradients, alphaHessian, minimumCoefficient, numberOfSets = numberOfSets )
            # . Find a solution.
            # . Just one go.
            if numberOfMinimizations <= 1:
                of.SetRandomVariables ( )
                optimizerState = optimizer.Iterate ( of, log = log )
            # . Multiple goes.
            else:
                bestF = 1.0e+20
                bestA = Real1DArray.WithExtent ( of.NumberOfVariables ( ) )
                for tries in range ( numberOfMinimizations ):
                    of.SetRandomVariables ( )
                    optimizerState = optimizer.Iterate ( of, log = log )
                    if optimizerState["Function Value"] < bestF:
                        bestF = optimizerState["Function Value"]
                        of.alphas.CopyTo ( bestA )
                bestA.CopyTo ( alphas )
