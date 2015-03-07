#-------------------------------------------------------------------------------
# . File      : FIREMinimizer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Fast inertial relaxation engine (FIRE) minimizer."""

from MultiDimensionalMinimizer import MultiDimensionalMinimizer, MultiDimensionalMinimizerState
from Real1DArray               import Real1DArray

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class FIREMinimizerState ( MultiDimensionalMinimizerState ):
    """FIRE minimizer state."""

    # . Default attributes.
    defaultAttributes = { "alpha"          : 0.0  ,
                          "gToAFactor"     : 1.0  ,
                          "numberPositive" : 0    ,
                          "timeStep"       : 0.0  ,
                          "v"              : None }
    defaultAttributes.update ( MultiDimensionalMinimizerState.defaultAttributes )

    def ExtractSurrogateData ( self, surrogate ):
        """Extract additional data from a surrogate objective function."""
        if surrogate is not None:
            self.gToAFactor = surrogate.AccelerationConversionFactor ( )

    @classmethod
    def FromObjectiveFunction ( selfClass, objectiveFunction ):
        """Constructor given an objective function."""
        self            = super ( FIREMinimizerState, selfClass ).FromObjectiveFunction ( objectiveFunction )
        self.gToAFactor = self.objectiveFunction.AccelerationConversionFactor ( )
        return self

    def SetUp ( self ):
       """Set up the state."""
       super ( FIREMinimizerState, self ).SetUp ( )
       self.v = Real1DArray.WithExtent ( self.numberOfVariables ) ; self.v.Set ( 0.0 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class FIREMinimizer ( MultiDimensionalMinimizer ):
    """FIRE multidimensional minimizer."""

    # . Default attributes.
    # . These values are those that are mentioned in the original paper.
    defaultAttributes = { "alpha"              : 0.1  ,
                          "alphaDownFactor"    : 0.99 , 
                          "maximumStep"        : 0.01 , # . RMS step.
                          "maximumTimeStep"    : 1.0 ,
                          "stepLatencyFactor"  : 5    ,
                          "timeStep"           : 0.1  ,
                          "timeStepDownFactor" : 0.5  ,
                          "timeStepUpFactor"   : 1.1  }
    defaultAttributes.update ( MultiDimensionalMinimizer.defaultAttributes )

    # . Default attribute names.
    defaultAttributeNames = { "Alpha"                 : "alpha"              ,
                              "Alpha Down Factor"     : "alphaDownFactor"    , 
                              "Maximum Step"          : "maximumStep"        ,
                              "Maximum Time Step"     : "maximumTimeStep"    ,
                              "Step Latency Factor"   : "stepLatencyFactor"  ,
                              "Time Step"             : "timeStep"           ,
                              "Time Step Down Factor" : "timeStepDownFactor" ,
                              "Time Step Up Factor"   : "timeStepUpFactor"   }
    defaultAttributeNames.update ( MultiDimensionalMinimizer.defaultAttributeNames )

    # . State class.
    stateObject = FIREMinimizerState

    def DetermineStep ( self, state ):
        """Get the step."""
        # . A velocity-Verlet algorithm is used to solve for x and v
        # . as Euler integration does not appear to be sufficient.
        # . Complete the calculation of v from the previous iteration.
        if state.numberOfIterations > 0:
            dt = state.timeStep
            state.v.AddScaledArray ( 0.5 * dt * state.gToAFactor , state.g )
            # . Perform the FIRE procedure.
            # . Negative or zero power.
            if state.g.Dot ( state.v ) >= 0.0:
                state.alpha          = self.alpha
                state.numberPositive = 0
                state.stepType       = "-"
                state.timeStep      *= self.timeStepDownFactor
                state.v.Set ( 0.0 )
            # . Positive power.
            else:
                # . Modify v.
                gNorm2 = state.g.Norm2 ( )
                vNorm2 = state.v.Norm2 ( )
                state.v.Scale ( 1.0 - state.alpha )
                if gNorm2 != 0.0: state.v.AddScaledArray ( - state.alpha * vNorm2 / gNorm2, state.g )
                # . Adjust factors.
                if state.numberPositive > self.stepLatencyFactor:
                    state.alpha   *= self.alphaDownFactor
                    state.timeStep = min ( state.timeStep * self.timeStepUpFactor, self.maximumTimeStep )
                state.numberPositive += 1
                state.stepType        = "+"
        else: state.stepType = "0"
        # . Calculate d and v at the half-step.
        dt = state.timeStep
        state.d.Set ( 0.0 )
        state.d.AddScaledArray (       dt                       , state.v )
        state.d.AddScaledArray ( 0.5 * dt**2 * state.gToAFactor , state.g )
        state.v.AddScaledArray ( 0.5 * dt    * state.gToAFactor , state.g )
        # . Apply linear constraints - unnecessary if already applied to g.
#        state.objectiveFunction.ApplyLinearConstraints ( state.d )
        # . Check the step length - this appears to be important for this optimizer.
        step = state.d.RootMeanSquare ( )
        if step > self.maximumStep: state.d.Scale ( self.maximumStep / step )
        # . Increment the variables.
        state.x.AddScaledArray ( 1.0, state.d )

    def Initialize ( self, state  ):
        """Initialization before iteration."""
        self.Restart ( state )
        super ( FIREMinimizer, self ).Initialize ( state )

    def Label ( self ): return "FIRE Minimizer"

    def Restart ( self, state ):
        """Restart the minimization."""
        state.alpha          = self.alpha
        state.numberPositive = 0
        state.stepType       = "R"
        state.timeStep       = self.timeStep
        state.v.Set ( 0.0 )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
