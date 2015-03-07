#-------------------------------------------------------------------------------
# . File      : LBFGSMinimizer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A very basic L-BFGS minimizer."""

from MultiDimensionalMinimizer import MultiDimensionalMinimizer, MultiDimensionalMinimizerState
from Real1DArray               import Real1DArray

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LBFGSMinimizerState ( MultiDimensionalMinimizerState ):
    """LBFGS minimizer state."""

    # . Default attributes.
    defaultAttributes = { "aux"      : None ,
                          "deltaG"   : None ,
                          "deltaX"   : None ,
                          "gOld"     : None ,
                          "history"  :    0 ,
                          "iHistory" :    0 ,
                          "rho"      : None ,
                          "xOld"     : None }
    defaultAttributes.update ( MultiDimensionalMinimizerState.defaultAttributes )

    @classmethod
    def FromObjectiveFunction ( selfClass, objectiveFunction, history = 0 ):
        """Constructor given an objective function."""
        self         = super ( MultiDimensionalMinimizerState, selfClass ).FromObjectiveFunction ( objectiveFunction )
        self.history = history
        self.SetUp ( )
        return self

    @classmethod
    def FromVariableArray ( selfClass, x, d = None, g = None, history = 0, surrogateObjectiveFunction = None ):
        """Constructor given a variable array."""
        self = selfClass ( )
        self.history           = history
        self.numberOfVariables = len ( x )
        self.d = d
        self.g = g
        self.x = x
        self.SetUp ( )
        self.ExtractSurrogateData ( surrogateObjectiveFunction )
        return self

    def SetUp ( self ):
       """Set up the state."""
       super ( LBFGSMinimizerState, self ).SetUp ( )
       self.aux    =   Real1DArray.WithExtent ( self.history           ) ; self.aux.Set  ( 0.0 )
       self.gOld   =   Real1DArray.WithExtent ( self.numberOfVariables ) ; self.gOld.Set ( 0.0 )
       self.xOld   =   Real1DArray.WithExtent ( self.numberOfVariables ) ; self.xOld.Set ( 0.0 )
       self.deltaG = [ Real1DArray.WithExtent ( self.numberOfVariables ) for i in range ( self.history ) ]
       self.deltaX = [ Real1DArray.WithExtent ( self.numberOfVariables ) for i in range ( self.history ) ]
       self.rho    = [ 0.0 for i in range ( self.history ) ]
       for ( g, x ) in zip ( self.deltaG, self.deltaX ):
           g.Set ( 0.0 )
           x.Set ( 0.0 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LBFGSMinimizer ( MultiDimensionalMinimizer ):
    """L-BFGS minimization."""

    # . Option attributes.
    defaultAttributes = { "history"     : 10   ,
                          "maximumStep" : 0.01 }
    defaultAttributes.update ( MultiDimensionalMinimizer.defaultAttributes )

    # . Option attribute names.
    defaultAttributeNames = { "History"      : "history"     ,
                              "Maximum Step" : "maximumStep" } # . RMS step.
    defaultAttributeNames.update ( MultiDimensionalMinimizer.defaultAttributeNames )

    # . State class.
    stateObject = LBFGSMinimizerState

    def DetermineStep ( self, state ):
        """Determine the step."""
        # . Reshuffle the data.
        if state.iHistory > self.history:
            temp = state.deltaG.pop ( 0 ) ; state.deltaG.append ( temp )
            temp = state.deltaX.pop ( 0 ) ; state.deltaX.append ( temp )
            temp = state.rho.pop    ( 0 ) ; state.rho.append    ( temp )
        # . Calculate new correction data.
        if state.iHistory > 0:
            i = min ( state.iHistory, self.history ) - 1
            state.g.CopyTo ( state.deltaG[i] ) ; state.deltaG[i].AddScaledArray ( -1.0, state.gOld )
            state.x.CopyTo ( state.deltaX[i] ) ; state.deltaX[i].AddScaledArray ( -1.0, state.xOld )
            hGG = state.deltaG[i].Dot ( state.deltaG[i] )
            hGX = state.deltaG[i].Dot ( state.deltaX[i] )
            state.rho[i] = 1.0 / hGX
            hScale       = hGX / hGG
        # . Save the old data.
        state.g.CopyTo ( state.gOld )
        state.x.CopyTo ( state.xOld )
        # . Calculate the step.
        state.g.CopyTo ( state.d )
        state.d.Scale ( -1.0 )
        if state.iHistory == 0:
            state.d.Scale ( 1.0 / state.g.Dot ( state.g ) )
        else:
            for i in reversed ( range ( min ( state.iHistory, state.history ) ) ):
                a = state.rho[i] * state.d.Dot ( state.deltaX[i] )
                state.d.AddScaledArray ( - a, state.deltaG[i] )
                state.aux[i] = a
            state.d.Scale ( hScale )
            for i in range ( min ( state.iHistory, state.history ) ):
                a = state.aux[i] - state.rho[i] * state.d.Dot ( state.deltaG[i] )
                state.d.AddScaledArray ( a, state.deltaX[i] )
                state.aux[i] = a
        # . Check the step length.
        step = state.d.RootMeanSquare ( )
        if step > self.maximumStep: state.d.Scale ( self.maximumStep / step )
        # . Increment the variables.
        state.x.AddScaledArray ( 1.0, state.d )
        # . Finish up.
        state.stepType  = "S{:d}".format ( min ( state.iHistory, self.history ) )
        state.iHistory += 1

    def Initialize ( self, state  ):
        """Initialization before iteration."""
        self.Restart ( state )
        super ( LBFGSMinimizer, self ).Initialize ( state )

    def Label ( self ): return "L-BFGS Minimizer"

    def Restart ( self, state ):
        """Restart."""
        state.iHistory = 0

    def StateFromObjectiveFunction ( self, objectiveFunction ):
        """Set up the state."""
        return self.__class__.stateObject.FromObjectiveFunction ( objectiveFunction, history = self.history )

    def StateFromVariableArray ( self, x, d = None, g = None, surrogateObjectiveFunction = None ):
        """Set up the state."""
        return self.__class__.stateObject.FromVariableArray ( x, d = d, g = g, history = self.history, surrogateObjectiveFunction = surrogateObjectiveFunction )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
