#-------------------------------------------------------------------------------
# . File      : UniDimensionalMinimizer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Base classes for unidimensional minimization algorithms."""

# . This may be split as in the multidimensional case.

from LogFileWriter     import logFile, LogFileActive, PrintPriority_Low
from ObjectiveFunction import UniDimensionalObjectiveFunction

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class UniDimensionalMinimizerState ( object ):
    """Base class for unidimensional minimization algorithms."""

    # . Default attributes.
    defaultAttributes = { "error"                 : None  ,
                          "f"                     : None  ,
                          "g"                     : None  ,
                          "isConverged"           : False ,
                          "log"                   : None  ,
                          "numberOfFunctionCalls" : 0     ,
                          "numberOfIterations"    : 0     ,
                          "objectiveFunction"     : None  ,
                          "statusMessage"         : None  ,
                          "stepType"              : ""    ,
                          "table"                 : None  ,
                          "x"                     : None  }

    # . Objective function class.
    objectiveFunctionClass = UniDimensionalObjectiveFunction

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def ExtractSurrogateData ( self, surrogate ):
        """Extract additional data from a surrogate objective function."""
        pass

    def Finalize ( self ):
        """Finalization."""
        # . Create the report.
        report = { "Converged"      : self.isConverged           ,
                   "Function Calls" : self.numberOfFunctionCalls ,
                   "Function Value" : self.f                     ,
                   "Gradient"       : self.g                     ,
                   "Iterations"     : self.numberOfIterations    ,
                   "Status Message" : self.statusMessage         ,
                   "Variable"       : self.x                     }
        if self.error is not None: report["Error"] = self.error
        return report

    @classmethod
    def FromObjectiveFunction ( selfClass, objectiveFunction ):
        """Constructor given an objective function."""
        # . Create the object.
        self = selfClass ( )
        # . Check the objective function.
        if not isinstance ( objectiveFunction, self.__class__.objectiveFunctionClass ): raise TypeError ( "Invalid objective function." )
        self.objectiveFunction = objectiveFunction
        # . Algorithm-specific set up.
        self.SetUp ( )
        return self

    def SetUp ( self ):
       """Set up the state."""
       pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class UniDimensionalMinimizer ( object ):
    """Base class for unidimensional minimization algorithms."""

    # . Default attributes.
    defaultAttributes = { "gradientTolerance" : 1.0e-3 ,
                          "logFrequency"      : 1      ,
                          "maximumIterations" : 0      }

    # . Default attribute names.
    defaultAttributeNames = { "Gradient Tolerance" : "gradientTolerance" ,
                              "Log Frequency"      : "logFrequency"      ,
                              "Maximum Iterations" : "maximumIterations" }

    # . State class.
    stateObject = UniDimensionalMinimizerState

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        self._Initialize ( )
        self.SetOptions ( **keywordArguments )
        self.CheckAttributes ( )

    def __str__ ( self ): return self.Label ( )

    def _Initialize ( self ):
        """Initialization."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def CheckAttributes ( self ):
        """Check the attributes."""
        pass

    def Continue ( self, state ):
        """Check to see if iterations should continue."""
        state.isConverged = ( state.g <= self.gradientTolerance )
        if   state.isConverged:                                  state.statusMessage = "Minimization converged."
        elif state.numberOfIterations >= self.maximumIterations: state.statusMessage = "Too many iterations."
        elif state.error is not None:                            state.statusMessage = "Minimization error: " + state.error
        else: state.statusMessage = None
        return ( state.statusMessage is None )

    def FunctionGradient ( self, state ):
        """Evaluate the function and its gradient."""
        ( state.f, state.g ) = state.objectiveFunction.FunctionGradient ( state.x )
        state.numberOfFunctionCalls += 1

    def Initialize ( self, state ):
        """Initialization before iteration."""
        try:
            self.FunctionGradient ( state )
        except Exception as error:
            state.error = error[0]
        state.stepType = "I"

    def Iterate ( self, objectiveFunction, log = logFile ):
        """Apply the algorithm to a function."""
        state = self.StateFromObjectiveFunction ( objectiveFunction )
        self.LogStart     ( state, log = log )
        self.Initialize   ( state )
        self.LogIteration ( state )
        while ( self.Continue ( state ) ):
            self.Iteration    ( state )
            self.LogIteration ( state )
        self.LogStop ( state )
        return state.Finalize ( )

    def Iteration ( self, state ):
        """Perform an iteration."""
        try:
            pass
        except Exception as error:
            state.error = error[0]
            import traceback, sys
            traceback.print_exc(file=sys.stdout)
        state.numberOfIterations += 1

    def Label ( self ): return "Uni-Dimensional Minimizer"

    def LogIteration ( self, state ):
        """Log an iteration."""
        state.objectiveFunction.LogIteration ( state.numberOfIterations )
        if ( state.table is not None ) and ( state.numberOfIterations % self.logFrequency == 0 ):
            state.table.Entry ( "{:d}".format ( state.numberOfIterations ) )
            state.table.Entry ( state.stepType )
            state.table.Entry ( "{:20.8f}".format ( state.f ) )
            state.table.Entry ( "{:20.8f}".format ( state.g ) )
            state.table.Entry ( "{:20.8f}".format ( state.x ) )

    def LogStart ( self, state, log = logFile ):
        """Start logging."""
        state.objectiveFunction.LogStart ( )
        if ( self.logFrequency > 0 ) and LogFileActive ( log, pushPriority = PrintPriority_Low ):
            state.log   = log
            state.table = log.GetTable ( columns = [ 6, 6, 20, 20, 20 ] )
            state.table.Start ( )
            state.table.Heading ( "Iteration", columnSpan = 2 )
            state.table.Heading ( "Function" )
            state.table.Heading ( "Gradient" )
            state.table.Heading ( "Variable" )

    def LogStop ( self, state ):
        """Stop logging."""
        state.objectiveFunction.LogStop ( )
        if state.log is not None:
            state.table.Stop ( )
            if state.statusMessage is not None: state.log.Paragraph ( state.statusMessage )
            summary = state.log.GetSummary ( )
            summary.Start ( self.Label ( ) + " Statistics" )
	    summary.Entry ( "Iterations",     "{:d}".format ( state.numberOfIterations    ) )
	    summary.Entry ( "Function Calls", "{:d}".format ( state.numberOfFunctionCalls ) )
	    summary.Stop ( )
            state.log.PriorityPop ( )

    def Restart ( self, state ):
        """Restart the minimization."""
        state.stepType = "R"

    def SetOptions ( self, **options ):
        """Set options from keyword arguments."""
        unknowns = set ( )
        for ( key, value ) in options.iteritems ( ):
            if key in self.__class__.defaultAttributes: setattr ( self, key, value )
            else: unknowns.add ( key )
        if len ( unknowns ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( unknowns ) ) + "." )

    def StateFromObjectiveFunction ( self, objectiveFunction ):
        """Set up the state."""
        return self.__class__.stateObject.FromObjectiveFunction ( objectiveFunction )

    def Summary ( self, log = logFile, pageWidth = 100, valueWidth = 14 ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( pageWidth = pageWidth, valueWidth = valueWidth )
            summary.Start ( self.Label ( ) + " Options" )
            keys = self.__class__.defaultAttributeNames.keys ( )
            keys.sort ( )
            for key in keys:
                value = getattr ( self, self.__class__.defaultAttributeNames[key] )
                if   isinstance ( value, bool  ):
                    if value: valuestring = "True"
                    else:     valuestring = "False"
                elif isinstance ( value, float      ): valuestring = "{:g}".format ( value )
                elif isinstance ( value, basestring ): valuestring =  value
                else:                                  valuestring = str ( value )
                summary.Entry ( key, valuestring )
            summary.Stop ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
