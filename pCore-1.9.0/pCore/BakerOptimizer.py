#-------------------------------------------------------------------------------
# . File      : BakerOptimizer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Find stationary points of multidimensional functions using the Baker algorithm.

Details of Baker's algorithm may be found in: J. Baker. 'An Algorithm for the Location of Transition States'. J. Comput. Chem. 4, 385-95 (1986).

The determination of the step length uses modifications described in: D. J. Wales. 'Rearrangements of 55-Atom Lennard-Jones and (C60)55 Clusters'. J. Chem. Phys. 101, 3750-3762 (1994).
"""

import math

from MultiDimensionalMinimizer import MultiDimensionalMinimizer, MultiDimensionalMinimizerState
from Real1DArray               import Real1DArray
from Real2DArray               import Real2DArray
from SymmetricMatrix           import SymmetricMatrix

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class BakerOptimizerState ( MultiDimensionalMinimizerState ):
    """Baker optimizer state."""

    # . Default attributes.
    defaultAttributes = { "currentMode"           :   -1 ,
                          "eigenValues"           : None ,
                          "eigenVectors"          : None ,
                          "gOld"                  : None ,
                          "h "                    : None ,
                          "hessian"               : None ,
                          "hW"                    : None ,
                          "numberOfNegativeModes" :    0 ,
                          "w"                     : None }
    defaultAttributes.update ( MultiDimensionalMinimizerState.defaultAttributes )

    def SetUp ( self ):
        """Set up the state."""
        super ( BakerOptimizerState, self ).SetUp ( )
        n = self.numberOfVariables
        self.eigenValues  = Real1DArray.WithExtent     ( n    ) ; self.eigenValues.Set  ( 0.0 )
        self.gOld         = Real1DArray.WithExtent     ( n    ) ; self.gOld.Set         ( 0.0 )
        self.h            = Real1DArray.WithExtent     ( n    ) ; self.h.Set            ( 0.0 )
        self.w            = Real1DArray.WithExtent     ( n    ) ; self.w.Set            ( 0.0 )
        self.eigenVectors = Real2DArray.WithExtents    ( n, n ) ; self.eigenVectors.Set ( 0.0 )
        self.hW           = SymmetricMatrix.WithExtent ( n    ) ; self.hW.Set           ( 0.0 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class BakerOptimizer ( MultiDimensionalMinimizer ):
    """Baker optimizer."""

    # . Default attributes.
    defaultAttributes = { "analyticFrequency"     :  0      ,
                          "followMode"            : -1      ,
                          "hessianUpdatingOption" : "BFGS"  ,
                          "maximumEigenvalue"     :  1.0e+5 ,
                          "maximumStep"           :  0.3    ,
                          "minimumEigenvalue"     :  1.0e-3 ,
                          "minimumStep"           :  1.0e-5 ,
                          "numberOfNegativeModes" :  0      ,
                          "resetStep"             :  1.0e-1 ,
                          "useNewtonRaphsonStep"  :  False  }
    defaultAttributes.update ( MultiDimensionalMinimizer.defaultAttributes )

    # . Default attribute names.
    defaultAttributeNames = { "Analytic Frequency"   : "analyticFrequency"     ,
                              "Follow Mode"          : "followMode"            ,
                              "Hessian Updating"     : "hessianUpdatingOption" ,
                              "Maximum Eigenvalue"   : "maximumEigenvalue"     ,
                              "Maximum Step"         : "maximumStep"           ,
                              "Minimum Eigenvalue"   : "minimumEigenvalue"     ,
                              "Minimum Step"         : "minimumStep"           ,
                              "Negative Eigenvalues" : "numberOfNegativeModes" ,
                              "Newton-Raphson Step"  : "useNewtonRaphsonStep"  ,
                              "Reset Step"           : "resetStep"             }
    defaultAttributeNames.update ( MultiDimensionalMinimizer.defaultAttributeNames )

    # . State class.
    stateObject = BakerOptimizerState

    def Continue ( self, state ):
        """Check to see if iterations should continue."""
        state.isConverged = ( state.rmsGradient <= self.rmsGradientTolerance ) and ( state.numberOfNegativeModes == self.numberOfNegativeModes )
        if   state.isConverged:                                  state.statusMessage = "Minimization converged."
        elif state.numberOfIterations >= self.maximumIterations: state.statusMessage = "Too many iterations."
        elif state.error is not None:                            state.statusMessage = "Minimization error: " + state.error
        else: state.statusMessage = None
        return ( state.statusMessage is None )

    def DetermineStep ( self, state ):
        """Get the step vector."""
        # . Transform the gradient vector to the local Hessian modes.
        state.eigenVectors.VectorMultiply ( state.g, state.h, transpose = True )
        # . Perform a Newton-Raphson step if the number of negative modes is OK and useNewtonRaphsonStep is on.
        if self.useNewtonRaphsonStep and ( state.numberOfNegativeModes == self.numberOfNegativeModes ):
            state.h.Divide ( state.eigenValues )
            state.eigenVectors.VectorMultiply ( state.h, state.d )
        # . Determine the step in the usual way.
        else:
	    # . Determine which mode to follow for a saddle point search using a maximum overlap criterion with the old step.
	    # . Only reset the current mode if the lowest mode is not the one being followed and this is not the first step.
            if ( self.numberOfNegativeModes > 0 ) and ( state.currentMode > 0 ) and ( state.numberOfIterations > 1 ):
                state.eigenVectors.VectorMultiply ( state.d, state.w, transpose = True )
                state.currentMode = state.w.AbsoluteMaximumIndex ( )
            # . Check for convergence to a stationary point of the wrong type.
            if ( state.rmsGradient <= self.rmsGradientTolerance ) and ( state.numberOfNegativeModes != self.numberOfNegativeModes ):
                # . Choose the direction as an appropriate eigenVector.
                v = max ( 0, state.currentMode )
                for i in range ( state.numberOfVariables ):
                    state.d[i] = self.resetStep * state.eigenVectors[i,v]
	    # . Determine the step sizes using Wales's formula.
            else:
	        # . Assume that a minimization is to be done along all modes.
                for ( i, ( e, g ) ) in enumerate ( zip ( state.eigenValues, state.h ) ):
                    state.w[i] = - 2.0 * g / ( math.fabs ( e ) + math.sqrt ( e * e + 4.0 * g * g ) )
                # . Make the step negative for the mode along which a maximization is being done.
                if state.currentMode >= 0: state.w[state.currentMode] *= -1.0
                # . Finally convert the step back to real space from mode space.
                state.eigenVectors.VectorMultiply ( state.w, state.d )
        # . Apply linear constraints.
        state.objectiveFunction.ApplyLinearConstraints ( state.d )
        # . Adjust the step size if it is too big or too small.
        dSize = state.d.Norm2 ( )
        if   dSize < self.minimumStep: state.d.Scale ( self.minimumStep / dSize )
        elif dSize > self.maximumStep: state.d.Scale ( self.maximumStep / dSize )
        # . Determine the new point.
        state.stepType = "N"
        state.x.AddScaledArray ( 1.0, state.d )

    def FunctionGradients ( self, state ):
        """Evaluate the function, gradients and hessian at the current point."""
        # . Analytic calculation.
        if self.UseAnalyticHessian ( state ):
            state.f = state.objectiveFunction.FunctionGradientsHessian ( state.x, state.g, state.hessian )
        # . Update the Hessian.
        else:
            state.g.CopyTo ( state.gOld )
            state.f = state.objectiveFunction.FunctionGradients ( state.x, state.g )
            state.g.CopyTo ( state.w )
            state.w.AddScaledArray ( -1.0, state.gOld )
            state.hessian.Update ( state.d, state.w, option = self.hessianUpdatingOption )
        # . Update some counters.
        state.numberOfFunctionCalls += 1
        state.rmsGradient = state.g.RootMeanSquare ( )
        # . Treat the Hessian.
        self.TreatHessian ( state )

    def Initialize ( self, state ):
        """Initialization before iteration."""
        # . Get the starting hessian.
        if not hasattr ( state.objectiveFunction, "StartingHessian" ): raise ValueError ( "Objective function has no hessian evaluator." )
        state.hessian = state.objectiveFunction.StartingHessian ( state.x )
        # . Initial function value.
        super ( BakerOptimizer, self ).Initialize ( state )
        # . Other initialization.
        state.currentMode = self.followMode
        if state.currentMode >= state.numberOfVariables: state.currentMode = -1

    def Label ( self ): return "Baker Optimizer"

    def TreatHessian ( self, state ):
        """Treat the Hessian."""
        # . Get the eigenValues and eigenVectors.
        state.hessian.CopyTo ( state.hW )
        state.hW.Diagonalize ( state.eigenValues, eigenVectors = state.eigenVectors )
        # . Count the number of negative eigenValues and reset excessively small or large values.
        state.numberOfNegativeModes = 0
        for ( i, e ) in enumerate ( state.eigenValues ):
            if e < 0.0:
                state.numberOfNegativeModes +=  1
                sign = -1.0
            else:
                sign =  1.0
            if   math.fabs ( e ) > self.maximumEigenvalue: state.eigenValues[i] = self.maximumEigenvalue * sign
            elif math.fabs ( e ) < self.minimumEigenvalue: state.eigenValues[i] = self.minimumEigenvalue * sign
        # . Update the step type.
        state.stepType += "{:d}".format ( state.numberOfNegativeModes )

    def UseAnalyticHessian ( self, state ):
        """Decide whether to use an analytic Hessian."""
        return hasattr ( state.objectiveFunction, "FunctionGradientsHessian" ) and ( self.analyticFrequency > 0 ) and ( state.numberOfIterations % self.analyticFrequency == 0 )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
