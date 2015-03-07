#-------------------------------------------------------------------------------
# . File      : QuasiNewtonMinimizer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes for performing quasi-Newton multidimensional minimization."""

import math

from MultiDimensionalMinimizer import MultiDimensionalMinimizer, MultiDimensionalMinimizerState
from Real1DArray               import Real1DArray
from Real2DArray               import Real2DArray
from SymmetricMatrix           import SymmetricMatrix

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QuasiNewtonMinimizerState ( MultiDimensionalMinimizerState ):
    """Quasi-Newton minimizer state."""

    # . Default attributes.
    defaultAttributes = { "eigenValues"  : None ,
                          "eigenVectors" : None ,
                          "fOld"         : None ,
                          "gOld"         : None ,
                          "hessian"      : None ,
                          "rfoMatrix"    : None ,
                          "rfoScale"     : None ,
                          "trustRadius"  : None ,
                          "w"            : None }
    defaultAttributes.update ( MultiDimensionalMinimizerState.defaultAttributes )

    def SetUp ( self ):
        """Set up the state."""
        super ( QuasiNewtonMinimizerState, self ).SetUp ( )
        # . Allocation.
        n = self.numberOfVariables
        p = n + 1
        self.eigenValues  = Real1DArray.WithExtent     ( p    ) ; self.eigenValues.Set ( 0.0 )
        self.gOld         = Real1DArray.WithExtent     ( n    ) ; self.gOld.Set        ( 0.0 )
        self.w            = Real1DArray.WithExtent     ( n    ) ; self.w.Set           ( 0.0 )
        self.eigenVectors = Real2DArray.WithExtents    ( p, p )
        self.rfoMatrix    = SymmetricMatrix.WithExtent ( p    ) ; self.rfoMatrix.Set   ( 0.0 )
        # . Initialization.
        self.rfoScale = 1.0e+00 / math.sqrt ( float ( n ) ) # . The RFO scale factor.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QuasiNewtonMinimizer ( MultiDimensionalMinimizer ):
    """Quasi-Newton minimization."""

    defaultAttributes = { "initialStep"          : 0.3e+00  ,
                          "maximumTrust"         : 1.0e+00  ,
                          "minimumTrust"         : 1.0e-02  ,
                          "trustScaleDownBound"  : 0.25e+00 ,
                          "trustScaleFactor"     : 2.0e+00  ,
                          "trustScaleUpBound"    : 0.75e+00 ,
                          "trustUpdateTolerance" : 1.0e-03  }
    defaultAttributes.update ( MultiDimensionalMinimizer.defaultAttributes )

    defaultAttributeNames = { "Initial Step"           : "initialStep"          ,
                              "Maximum Trust Radius"   : "maximumTrust"         ,
                              "Minimum Trust Radius"   : "minimumTrust"         ,
                              "Trust Scale Down Bound" : "trustScaleDownBound"  ,
                              "Trust Scale Factor"     : "trustScaleFactor"     ,
                              "Trust Scale Up Bound"   : "trustScaleUpBound"    ,
                              "Trust Update Tolerance" : "trustUpdateTolerance" }
    defaultAttributeNames.update ( MultiDimensionalMinimizer.defaultAttributeNames )

    stateObject = QuasiNewtonMinimizerState

    def DetermineStep ( self, state ):
        """Get the step."""
        if state.numberOfIterations > 0:
            # . Calculate the predicted change in the function value.
            state.hessian.VectorMultiply ( state.d, state.w )
            state.w.Scale ( 0.5 )
            state.w.AddScaledArray ( 1.0, state.gOld )
            dFPredicted = state.d.Dot ( state.w )
            # . Calculate the change in the function value.
            dFActual = state.f - state.fOld
            # . Update the hessian.
            state.g.CopyTo ( state.w )
            state.w.AddScaledArray ( -1.0, state.gOld )
            state.hessian.Update   ( state.d, state.w )
            # . Update the trust radius.
            self.UpdateTrustRadius ( state, dFActual, dFPredicted )
        # . Set up and solve the RFO equations.
        # . Should use slices in the following.
        scale = 1.0 / math.sqrt ( state.rfoScale )
        state.rfoMatrix.Copy_DB_From_DB         ( 0, state.hessian, 0, state.numberOfVariables )
        state.rfoMatrix.Scale_DB                ( 0, scale**2     ,    state.numberOfVariables )
        state.rfoMatrix.Copy_Column_From_Vector ( state.numberOfVariables, 0, state.g, 0, state.numberOfVariables )
        state.rfoMatrix.Scale_OB                ( state.numberOfVariables, 0, scale  , 1, state.numberOfVariables )
        state.rfoMatrix[state.numberOfVariables,state.numberOfVariables] = 0.0
        state.rfoMatrix.Diagonalize ( state.eigenValues, eigenVectors = state.eigenVectors )
        # . The step corresponds to the appropriate elements of the intermediately-normalized eigenVector of smallest eigenValue.
        scale /= state.eigenVectors[state.numberOfVariables,0]
        for i in range ( state.numberOfVariables ): state.d[i] = state.eigenVectors[i,0]
        state.d.Scale ( scale )
        # . Apply linear constraints.
        state.objectiveFunction.ApplyLinearConstraints ( state.d )
        # . Check the step length.
        stepSize = state.d.Norm2 ( )
        if stepSize > state.trustRadius:
            state.d.Scale ( state.trustRadius / stepSize )
            stepSize = state.trustRadius
        # . Save various old values.
        state.fOld = state.f
        state.g.CopyTo ( state.gOld )
        # . Determine the new point.
        state.x.AddScaledArray ( 1.0, state.d )
        # . Finish up.
        state.stepType = "N"

    def Initialize ( self, state ):
        """Initialization before iteration."""
        super ( QuasiNewtonMinimizer, self ).Initialize ( state )
        # . Get the starting hessian.
        if not hasattr ( state.objectiveFunction, "StartingHessian" ): raise ValueError ( "Objective function has no hessian evaluator." )
        state.hessian     = state.objectiveFunction.StartingHessian ( state.x )
        # . Other initialization.
        state.fOld        = state.f
        state.trustRadius = self.initialStep

    def Label ( self ): return "Quasi-Newton Minimizer"

    def UpdateTrustRadius ( self, state, dFActual, dFPredicted ):
        """Update the trust radius."""
        if ( math.fabs ( dFActual ) > self.trustUpdateTolerance ) and ( math.fabs ( dFPredicted ) > self.trustUpdateTolerance ):
            ratio = dFActual / dFPredicted
            if   ratio < self.trustScaleDownBound: state.trustRadius /=             self.trustScaleFactor
            elif ratio > self.trustScaleUpBound  : state.trustRadius *= math.sqrt ( self.trustScaleFactor )
        if state.trustRadius < self.minimumTrust: state.trustRadius = self.minimumTrust
        if state.trustRadius > self.maximumTrust: state.trustRadius = self.maximumTrust

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
