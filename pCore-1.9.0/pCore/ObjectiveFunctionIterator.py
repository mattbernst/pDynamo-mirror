#-------------------------------------------------------------------------------
# . File      : ObjectiveFunctionIterator.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Objective function iterator."""

from LogFileWriter     import logFile, LogFileActive, PrintPriority_Low
from ObjectiveFunction import ObjectiveFunction
from Real1DArray       import Real1DArray

#
# . Notes:
#
#   It might be preferable to have linear constraints at this level rather than
#   in the objective function. Thus, the iterator applies the LCs and the iterator
#   state stores them. One could apply them normally directly after gradient
#   evaluation. Alternatively, explicitly call ApplyLCs (maybe better if LCs
#   are not vectors) when required rather than implicitly by OF.
#
#   It might also be preferable to have Iterate use a state instead of an
#   objective function with the stete set up explicitly beforehand.
#

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ObjectiveFunctionIteratorState ( object ):
    """Base class for algorithm states that require objective functions."""

    # . Default attributes.
    defaultAttributes = { "error"              : None ,
                          "f"                  : None ,
                          "log"                : None ,
                          "numberOfIterations" : 0    ,
                          "numberOfVariables"  : 0    ,
                          "objectiveFunction"  : None ,
                          "statusMessage"      : None ,
                          "table"              : None ,
                          "x"                  : None }

    # . Objective function class.
    objectiveFunctionClass = ObjectiveFunction

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def Finalize ( self ):
        """Finalization."""
        # . Create the report.
        report = { "Function Value" : self.f                  ,
                   "Iterations"     : self.numberOfIterations }
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
        # . Set up variables.
        self.numberOfVariables = objectiveFunction.NumberOfVariables ( )
        self.x                 = objectiveFunction.VariablesAllocate ( )
        objectiveFunction.VariablesGet ( self.x )
        # . Finish up.
        self.SetUp ( )
        return self

    def SetUp ( self ):
       """Set up the state."""
       pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ObjectiveFunctionIterator ( object ):
    """Base class for algorithms that require objective functions."""

    # . Default attributes.
    defaultAttributes = { "logFrequency"      :  0 ,
                          "maximumIterations" : -1 }

    # . Default attribute names.
    defaultAttributeNames = { "Log Frequency"      : "logFrequency"      ,
                              "Maximum Iterations" : "maximumIterations" }

    # . State class.
    stateObject = ObjectiveFunctionIteratorState

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
        if   state.error is not None:                            state.statusMessage = "Iterator error: " + state.error
        elif state.numberOfIterations >= self.maximumIterations: state.statusMessage = "Iterator iterations terminated."
        else: state.statusMessage = None
        return ( state.statusMessage is None )

    def Initialize ( self, state ):
        """Initialization before iteration."""
        pass

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
        state.numberOfIterations += 1

    def Label ( self ): return "Objective Function Iterator"

    def LogIteration ( self, state ):
        """Log an iteration."""
        if ( state.log is not None ) and ( state.numberOfIterations % self.logFrequency == 0 ):
            pass

    def LogStart ( self, state, log = logFile ):
        """Start logging."""
        if ( self.logFrequency > 0 ) and LogFileActive ( log, pushPriority = PrintPriority_Low ):
            state.log = log

    def LogStop ( self, state ):
        """Stop logging."""
        if ( state.log is not None ) and ( state.statusMessage is not None ):
            state.log.Paragraph ( state.statusMessage )

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
                    if value: valueString = "True"
                    else:     valueString = "False"
                elif isinstance ( value, float      ): valueString = "{:g}".format ( value )
                elif isinstance ( value, basestring ): valueString =  value
                else:                                  valueString = str ( value )
                summary.Entry ( key, valueString )
            summary.Stop ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
