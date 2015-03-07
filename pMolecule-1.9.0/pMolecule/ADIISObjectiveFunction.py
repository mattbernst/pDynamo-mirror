#-------------------------------------------------------------------------------
# . File      : ADIISObjectiveFunction.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the objective function for the ADIIS SCF converger."""

import math

from pCore import ObjectiveFunction, RandomNumberGenerator, Real1DArray

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ADIISObjectiveFunction ( ObjectiveFunction ):
    """The objective function for the ADIIS SCF converger."""

    # . Attributes.
    defaultAttributes = { "alphas"                  : None, \
                          "gradientsAlpha"          : None, \
                          "hessianAlpha"            : None, \
                          "minimumCoefficient"      : None, \
                          "numberOfSets"            :    1, \
                          "numberOfVariables"       :    0, \
                          "numberOfVariablesPerSet" :    0, \
                          "work"                    : None, \
                          "workA"                   : None, \
                          "workS"                   : None, \
                          "workT"                   : None  }
    defaultAttributes.update ( ObjectiveFunction.defaultAttributes )

    @classmethod
    def FromOptions ( selfClass, alphas, gradientsAlpha, hessianAlpha, minimumCoefficient, numberOfSets = 1 ):
        """Constructor from options."""
        self = selfClass ( )
        self.alphas                  = alphas
        self.gradientsAlpha          = gradientsAlpha
        self.hessianAlpha            = hessianAlpha
        self.minimumCoefficient      = minimumCoefficient
        self.numberOfSets            = numberOfSets
        self.numberOfVariables       = len ( self.alphas )
        self.numberOfVariablesPerSet = self.numberOfVariables // self.numberOfSets
        self.work                    = Real1DArray.WithExtent ( self.numberOfVariables )
        self.workA                   = Real1DArray.WithExtent ( self.numberOfVariables )
        self.workS                   = Real1DArray.WithExtent ( self.numberOfVariables )
        self.workT                   = Real1DArray.WithExtent ( self.numberOfVariables )
        return self

    def Function ( self, variables ):
        """Evaluate the function."""
        self.VariablesPut ( variables )
        self.hessianAlpha.VectorMultiply ( self.alphas, self.work   )
        self.work.AddScaledArray         ( 2.0, self.gradientsAlpha )
        f = 0.5 * self.work.Dot ( self.alphas )
        return f

    def FunctionGradients ( self, variables, gradients ):
        """Evaluate the function and gradients."""
        self.VariablesPut ( variables )
        self.hessianAlpha.VectorMultiply ( self.alphas, self.work   )
        self.work.AddScaledArray         ( 1.0, self.gradientsAlpha )
        for s in range ( self.numberOfSets ):
            increment     = s * self.numberOfVariablesPerSet
            factor        = 0.0
            normalization = 0.0
            for u in range ( self.numberOfVariablesPerSet ):
                i = increment + u
                factor        += self.alphas[i] * self.work[i]
                normalization += variables[i]**2
            factor /= ( 1.0 - self.minimumCoefficient )
            scaling = 2.0 * ( 1.0 - self.minimumCoefficient ) / normalization
            for u in range ( self.numberOfVariablesPerSet ):
                i = increment + u
                gradients[i] = scaling * variables[i] * ( self.work[i] - factor )
        self.work.AddScaledArray ( 1.0, self.gradientsAlpha )
        f = 0.5 * self.work.Dot  ( self.alphas )
        return f

# . Need to check (for minimumCoefficient).
#    def FunctionGradientsHessian ( self, variables, gradients, hessian ):
#        """Evaluate the function, gradients and hessian."""
#        # . Do for the general case even though there is no coupling in the Hessian.
#        self.VariablesPut ( variables )
#        # . Calculate the set normalizations.
#        normalizations = []
#        for s in range ( self.numberOfSets ):
#            increment     = s * self.numberOfVariablesPerSet
#            normalization = 0.0
#            for u in range ( self.numberOfVariablesPerSet ):
#                i = increment + u
#                normalization += variables[i]**2
#            normalizations.append ( normalization )
#       # . H * a.
#       self.hessianAlpha.CopyTo ( hessian )
#       for s in range ( self.numberOfSets ):
#           incrementS = s * self.numberOfVariablesPerSet
#           self.workA.Set ( 0.0 )
#           for u in range ( self.numberOfVariablesPerSet ):
#               i = incrementS + u
#               self.workA[i] = self.alphas[i]
#           self.hessianAlpha.VectorMultiply ( self.workA, self.workS )
#           for t in range ( s + 1 ):
#               incrementT = t * self.numberOfVariablesPerSet
#               for v in range ( self.numberOfVariablesPerSet ):
#                   j = incrementT + v
#                   self.workA[j] = self.alphas[j]
#               self.hessianAlpha.VectorMultiply ( self.workA, self.workT )
#               factor = self.workS.Dot ( self.workT )
#               for u in range ( self.numberOfVariablesPerSet ):
#                   i = incrementS + u
#                   for v in range ( self.numberOfVariablesPerSet ):
#                       j = incrementT + v
#                       if j <= i:
#                           hessian[i,j] += ( factor - self.workT[i] - self.workS[j] )
#       # . g + H * a.
#       self.work.AddScaledArray ( 1.0, self.gradientsAlpha )
#       # . Factors.
#       factors = []
#       for s in range ( self.numberOfSets ):
#           increment = s * self.numberOfVariablesPerSet
#           factor    = 0.0
#           for u in range ( self.numberOfVariablesPerSet ):
#               i = increment + u
#               factor += ( self.alphas[i] * self.work[i] )
#           factors.append ( factor )
#       # . Diagonal blocks.
#       for s in range ( self.numberOfSets ):
#           increment = s * self.numberOfVariablesPerSet
#           factor    = 0.0
#           for u in range ( self.numberOfVariablesPerSet ):
#               i = increment + u
#               for v in range ( u + 1 ):
#                   j = increment + v
#                   hessian[i,j] -= ( self.work[i] + self.work[j] - 2.0 * factors[s] )
#       # . Scaling.
#       for s in range ( self.numberOfSets ):
#           incrementS = s * self.numberOfVariablesPerSet
#           for u in range ( self.numberOfVariablesPerSet ):
#               i = incrementS + u
#               for t in range ( s + 1 ):
#                   incrementT = t * self.numberOfVariablesPerSet
#                   for v in range ( self.numberOfVariablesPerSet ):
#                       j = incrementT + v
#                       if j <= i:
#                           hessian[i,j] *= ( 4.0 * variables[i] * variables[j] ) / ( normalizations[s] * normalizations[t] )
#       # . Diagonal terms.
#       for s in range ( self.numberOfSets ):
#           increment = s * self.numberOfVariablesPerSet
#           for u in range ( self.numberOfVariablesPerSet ):
#               i = increment + u
#               hessian[i,i] += ( 2.0 / normalizations[s] ) * ( self.work[i] - factors[s] )
#       # . Gradients and function.
#       for s in range ( self.numberOfSets ):
#           increment     = s * self.numberOfVariablesPerSet
#           factor        = 0.0
#           normalization = 0.0
#           for u in range ( self.numberOfVariablesPerSet ):
#               i = increment + u
#               factor        += self.alphas[i] * self.work[i]
#               normalization += variables[i]**2
#           scaling = 2.0 / normalization
#           for u in range ( self.numberOfVariablesPerSet ):
#               i = increment + u
#               gradients[i] = scaling * variables[i] * ( self.work[i] - factor )
#       self.work.AddScaledArray ( 1.0, self.gradientsAlpha )
#       f = 0.5 * self.work.Dot  ( self.alphas )
#       return f

    def NumberOfVariables ( self ):
        """Return the number of variables."""
        return self.numberOfVariables

    def SetRandomVariables ( self ):
        """Initialization to random variables."""
        rng = RandomNumberGenerator.WithRandomSeed ( )
        rng.NextReals     ( self.work )
        self.VariablesPut ( self.work )

    def VariablesAllocate ( self ):
        """Return an object to hold the variables."""
        variables = Real1DArray.WithExtent ( self.numberOfVariables )
        variables.Set ( 0.0 )
        return variables

    def VariablesGet ( self, variables ):
        """Get a starting set of variables."""
        for s in range ( self.numberOfSets ):
            increment = s * self.numberOfVariablesPerSet
            for u in range ( self.numberOfVariablesPerSet ):
                i = increment + u
                if u == 0: variables[i] = math.sqrt ( max ( self.alphas[i] - self.minimumCoefficient, 0.0 ) )
                else:      variables[i] = math.sqrt ( max ( self.alphas[i]                          , 0.0 ) )

    def VariablesPut ( self, variables ):
        """Put variables in the proper place."""
        for s in range ( self.numberOfSets ):
            increment     = s * self.numberOfVariablesPerSet
            normalization = 0.0
            for u in range ( self.numberOfVariablesPerSet ):
                i = increment + u
                self.alphas[i] = variables[i]**2
                normalization += self.alphas[i]
            normalization /= ( 1.0 - self.minimumCoefficient )
            check = 0.0
            for u in range ( self.numberOfVariablesPerSet ):
                i = increment + u
                self.alphas[i] /= normalization
                if u == 0: self.alphas[i] += self.minimumCoefficient
                check += self.alphas[i]
            if math.fabs ( check - 1.0 ) > 1.0e-10: raise ValueError ( "Normalization error." )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
