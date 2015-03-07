#-------------------------------------------------------------------------------
# . File      : ObjectiveFunction.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines objective function classes for ObjectiveFunctionIterator algorithms."""

import math

from LogFileWriter   import logFile, LogFileActive
from Real1DArray     import Real1DArray
from Serialization   import RawObjectConstructor
from SymmetricMatrix import SymmetricMatrix

#===================================================================================================================================
# . Class.
#===================================================================================================================================
# . Base class for functions with derivatives.
class ObjectiveFunction ( object ):
    """The base class for objective functions."""

    # . Attributes.
    defaultAttributes = { "doLogging"    : False ,
                          "trajectories" : None  }

    def __getstate__ ( self ):
        return dict ( self.__dict__ )

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        self.__dict__.update ( state )

    def _Initialize ( self ):
        """Initialization."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        self.trajectories = []
 
    def ApplyLinearConstraints ( self, array ):
        """Apply linear constraints to an one-dimensional array."""
        pass

    def DefineTrajectory ( self, trajectory, logFrequency ):
        """Define a trajectory and logFrequency."""
        if ( trajectory is not None ) and ( logFrequency > 0 ):
            self.trajectories.append ( ( logFrequency, trajectory ) )

    def Function ( self, variables ):
        """Evaluate the function."""
        return 0.0

    def FunctionGradients ( self, variables, gradients ):
        """Evaluate the function and gradients."""
        return 0.0

    def FunctionGradientsHessian ( self, variables, gradients, hessian ):
        """Evaluate the function, gradients and hessian."""
        return 0.0

    def LogIteration ( self, iteration, identifier = None ):
        """Log an iteration of an iterator."""
        if self.doLogging:
            for ( logFrequency, trajectory ) in self.trajectories:
                if ( iteration % logFrequency == 0 ): trajectory.WriteOwnerData ( ) # identifier?

    def LogStart ( self, **keywordArguments ):
        """Start logging."""
        if len ( self.trajectories ) > 0:
            self.doLogging = True
            for ( logFrequency, trajectory ) in self.trajectories: trajectory.WriteHeader ( **keywordArguments )
        else:
            self.doLogging = False

    def LogStop ( self ):
        """Stop logging."""
        if self.doLogging:
            for ( logFrequency, trajectory ) in self.trajectories:
                trajectory.WriteFooter ( )
                trajectory.Close       ( )
            self.doLogging = False

    def NumberOfVariables ( self ):
        """Return the number of variables."""
        return 0

    def NumericalGradients ( self, variables, delta = 1.0e-4 ):
        """Calculate the gradients numerically.

        |delta| is the absolute step length but it should probably be modified so that it is a step length relative to the size of the variable.
        """
        n = self.NumberOfVariables ( )
        if n > 0:
            gradients = Real1DArray.WithExtent ( n )
            for i in range ( n ):
                temp         = variables[i]
                variables[i] = temp + delta
                fplus        = self.Function ( variables )
                variables[i] = temp - delta
                fminus       = self.Function ( variables )
                variables[i] = temp
                gradients[i] = ( fplus - fminus )
            gradients.Scale ( 1.0 / ( 2.0 * delta ) )
            self.ApplyLinearConstraints ( gradients )
            return gradients
        else:
            return None

    def NumericalHessian ( self, variables, delta = 1.0e-4 ):
        """Calculate the Hessian numerically."""
        n = self.NumberOfVariables ( )
        if n > 0:
            gminus  = Real1DArray.WithExtent ( n )
            gplus   = Real1DArray.WithExtent ( n )
            hessian = SymmetricMatrix.WithExtent ( n )
            hessian.Set ( 0.0 )
            for i in range ( n ):
                temp         = variables[i]
                variables[i] = temp + delta
                self.FunctionGradients ( variables, gplus  )
                variables[i] = temp - delta
                self.FunctionGradients ( variables, gminus )
                variables[i] = temp
                gplus.AddScaledArray ( -1.0, gminus )
                for j in range ( n ): hessian[i,j] += gplus[j]
            hessian.Scale ( 1.0 / ( 4.0 * delta ) )
            for i in range ( n ): hessian[i,i] *= 2.0
            return hessian
        else:
            return None
       
    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def StartingHessian ( self, variables ):
        """Get a starting hessian."""
        return self.NumericalHessian ( variables )

    def TestGradients ( self, delta = 1.0e-4, log = logFile, tolerance = 1.0e-4 ):
        """Compare analytical and numerical derivatives."""
        maxDiff    = 0.0
        nvariables = self.NumberOfVariables ( )
        if LogFileActive ( log ) and ( nvariables > 0 ):
            # . Allocate arrays.
            ga = Real1DArray.WithExtent ( nvariables ) ; ga.Set ( 0.0 )
            x  = self.VariablesAllocate ( )
            # . Analytical calculation.
            self.VariablesGet ( x )
            f  = self.FunctionGradients ( x, ga )
            # . Numerical calculation.
            gn = self.NumericalGradients ( x, delta = delta )
            if ( ga is not None ) and ( gn is not None ):
                # . Compute the maximum absolute difference between the gradients.
                for i in range ( nvariables ):
                    diff    = math.fabs ( ga[i] - gn[i] )
                    maxDiff = max ( maxDiff, diff )
                # . Output only those differences which are bigger than tolerance.
                if maxDiff > tolerance:
                    table = log.GetTable ( columns = [ 10, 20, 20, 20 ] )
                    table.Start ( )
                    table.Title ( "Analytical and Numerical Gradients" )
                    table.Heading ( "Variable"     )
                    table.Heading ( "Analytical"   )
                    table.Heading ( "Numerical"    )
                    table.Heading ( "|Difference|" )
                    for i in range ( nvariables ):
                        diff = math.fabs ( ga[i] - gn[i] )
                        if diff > tolerance:
                            table.Entry ( "{:d}".format ( i ) )
                            table.Entry ( "{:20.6f}".format ( ga[i] ) )
                            table.Entry ( "{:20.6f}".format ( gn[i] ) )
                            table.Entry ( "{:20.6f}".format ( diff  ) )
                    table.Stop ( )
                # . Output the biggest difference.
                log.Paragraph ( "Maximum absolute gradient difference = {:.6f}".format ( maxDiff ) )
        return maxDiff

    def TestHessian ( self, delta = 1.0e-4, log = logFile, tolerance = 1.0e-4 ):
        """Compare analytical and numerical hessian."""
        maxDiff    = 0.0
        nvariables = self.NumberOfVariables ( )
        if LogFileActive ( log ) and ( nvariables > 0 ):
            # . Allocate arrays.
            g  = Real1DArray.WithExtent     ( nvariables )
            ha = SymmetricMatrix.WithExtent ( nvariables )
            x  = self.VariablesAllocate ( )
            # . Analytical calculation.
            self.VariablesGet ( x )
            f  = self.FunctionGradientsHessian ( x, g, ha )
            # . Numerical calculation.
            hn = self.NumericalHessian ( x, delta = delta )
            if ( ha is not None ) and ( hn is not None ):
                # . Compute the maximum absolute difference between the gradients.
                for i in range ( nvariables ):
                    for j in range ( i+1 ):
                        diff    = math.fabs ( ha[i,j] - hn[i,j] )
                        maxDiff = max ( maxDiff, diff )
                # . Output only those differences which are bigger than tolerance.
                if maxDiff > tolerance:
                    table = log.GetTable ( columns = [ 10, 10, 20, 20, 20 ] )
                    table.Start ( )
                    table.Title ( "Analytical and Numerical Hessians" )
                    table.Heading ( "Variables", columnSpan = 2 )
                    table.Heading ( "Analytical"   )
                    table.Heading ( "Numerical"    )
                    table.Heading ( "|Difference|" )
                    for i in range ( nvariables ):
                        for j in range ( i+1 ):
                            diff = math.fabs ( ha[i,j] - hn[i,j] )
                            if diff > tolerance:
                                table.Entry ( "{:d}".format ( i ) )
                                table.Entry ( "{:d}".format ( j ) )
                                table.Entry ( "{:20.6f}".format ( ha[i,j] ) )
                                table.Entry ( "{:20.6f}".format ( hn[i,j] ) )
                                table.Entry ( "{:20.6f}".format ( diff    ) )
                    table.Stop ( )
                # . Output the biggest difference.
                log.Paragraph ( "Maximum absolute hessian difference = {:.6f}".format ( maxDiff ) )
        return maxDiff

    def VariablesAllocate ( self ):
        """Return a variables object of the correct size."""
        return None

    def VariablesGet ( self, variables ):
        """Fill the variable array."""
        pass

    def VariablesPut ( self, variables ):
        """Empty the variable array."""
        pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
# . Base class for uni-dimensional functions with derivatives.
class UniDimensionalObjectiveFunction ( object ):
    """The base class for uni-dimensional objective functions."""

    # . Attributes.
    defaultAttributes = { "doLogging"    : False ,
                          "trajectories" : None  ,
                          "variable"     : None  }

    def __getstate__ ( self ):
        return dict ( self.__dict__ )

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        self.__dict__.update ( state )

    def _Initialize ( self ):
        """Initialization."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        self.trajectories = []

    def DefineTrajectory ( self, trajectory, logFrequency ):
        """Define a trajectory and logFrequency."""
        if ( trajectory is not None ) and ( logFrequency > 0 ):
            self.trajectories.append ( ( logFrequency, trajectory ) )

    def Function ( self, variable ):
        """Evaluate the function."""
        return 0.0

    def FunctionGradient ( self, variable ):
        """Evaluate the function and gradient."""
        return ( 0.0, 0.0 )

    def FunctionGradientHessian ( self, variable ):
        """Evaluate the function, gradient and hessian."""
        return ( 0.0, 0.0, 0.0 )

    def LogIteration ( self, iteration, identifier = None ):
        """Log an iteration of an iterator."""
        if self.doLogging:
            for ( logFrequency, trajectory ) in self.trajectories:
                if ( iteration % logFrequency == 0 ): trajectory.WriteOwnerData ( )

    def LogStart ( self, **keywordArguments ):
        """Start logging."""
        if len ( self.trajectories ) > 0:
            self.doLogging = True
            for ( logFrequency, trajectory ) in self.trajectories: trajectory.WriteHeader ( **keywordArguments )
        else:
            self.doLogging = False

    def LogStop ( self ):
        """Stop logging."""
        if self.doLogging:
            for ( logFrequency, trajectory ) in self.trajectories:
                trajectory.WriteFooter ( )
                trajectory.Close       ( )
            self.doLogging = False

    def NumericalGradient ( self, variable, delta = 1.0e-4 ):
        """Calculate the gradient numerically."""
        temp     = variable
        variable = temp + delta
        fPlus    = self.Function ( variable )
        variable = temp - delta
        fMinus   = self.Function ( variable )
        variable = temp
        return ( ( fPlus - fMinus ) / ( 2.0 * delta ) )

    def NumericalHessian ( self, variable, delta = 1.0e-4 ):
        """Calculate the Hessian numerically."""
        temp               = variable
        variable           = temp + delta
        ( fPlus , gPlus  ) = self.FunctionGradient ( variable )
        variable           = temp - delta
        ( fMinus, gMinus ) = self.FunctionGradient ( variable )
        variable           = temp
        return ( ( gPlus - gMinus ) / ( 2.0 * delta ) )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def RestoreBestToCurrent ( self ):
        """Make the best point the current one."""
        pass

    def StartingHessian ( self, variable ):
        """Get a starting hessian."""
        return self.NumericalHessian ( variable )

    def StoreCurrentAsBest ( self ):
        """Save the current point as the best point."""
        pass

    def TestGradient ( self, delta = 1.0e-4, log = logFile, tolerance = 1.0e-4 ):
        """Compare analytical and numerical derivatives."""
        maxDiff = 0.0
        if LogFileActive ( log ):
            # . Analytical and numerical calculation.
            ( f, gA ) = self.FunctionGradient  ( self.variable )
            gN        = self.NumericalGradient ( self.variable, delta = delta )
            maxDiff   = math.fabs ( gA - gN )
            # . Output the difference only if it is bigger than tolerance.
            if maxDiff > tolerance:
                table = log.GetTable ( columns = [ 20, 20 ] )
                table.Start ( )
                table.Title ( "Analytical and Numerical Gradients" )
                table.Entry ( "Analytical"   ) ; table.Entry ( "{:20.6f}".format ( gA      ) )
                table.Entry ( "Numerical"    ) ; table.Entry ( "{:20.6f}".format ( gN      ) )
                table.Entry ( "|Difference|" ) ; table.Entry ( "{:20.6f}".format ( maxDiff ) )
                table.Stop ( )
        return maxDiff

    def TestHessian ( self, delta = 1.0e-4, log = logFile, tolerance = 1.0e-4 ):
        """Compare analytical and numerical hessian."""
        maxDiff = 0.0
        if LogFileActive ( log ):
            # . Analytical and numerical calculation.
            ( f, gA, hA ) = self.FunctionGradientHessian ( self.variable )
            hN            = self.NumericalHessian        ( self.variable, delta = delta )
            maxDiff       = math.fabs ( hA - hN )
            # . Output the difference only if it is bigger than tolerance.
            if maxDiff > tolerance:
                table = log.GetTable ( columns = [ 20, 20 ] )
                table.Start ( )
                table.Title ( "Analytical and Numerical Hessian" )
                table.Entry ( "Analytical"   ) ; table.Entry ( "{:20.6f}".format ( hA      ) )
                table.Entry ( "Numerical"    ) ; table.Entry ( "{:20.6f}".format ( hN      ) )
                table.Entry ( "|Difference|" ) ; table.Entry ( "{:20.6f}".format ( maxDiff ) )
                table.Stop ( )
        return maxDiff

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
