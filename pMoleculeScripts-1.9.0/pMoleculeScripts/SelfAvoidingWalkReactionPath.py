#-------------------------------------------------------------------------------
# . File      : SelfAvoidingWalkReactionPath.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Self-avoiding walk method for finding reaction paths."""

import math

from pCore     import LBFGSMinimizer, logFile, LogFileActive, ObjectiveFunction, Real1DArray
from pMolecule import SystemGeometryObjectiveFunction

# . Allow alternative minimizers.

#===================================================================================================================================
# . Parameter definitions.
#===================================================================================================================================
# . Default SAW objective function parameters.
_DefaultGamma =  100.0
_DefaultKappa =    2.0
_DefaultRho   = 5000.0

# . Default minimizer options.
_DefaultMinimizerOptions = { "logFrequency"         :   1     ,
                             "maximumIterations"    : 250     ,
                             "maximumStep"          :   0.001 ,
                             "rmsGradientTolerance" :   0.1   }

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
class _SAWObjectiveFunction ( ObjectiveFunction ):
    """An object function class for the SAW method."""

    # . Attributes.
    attributes = { "distances1"                : None          ,   
                   "distances2"                : None          ,   
                   "energies"                  : None          ,   
                   "g"                         : None          ,   
                   "gd"                        : None          ,   
                   "gamma"                     : _DefaultGamma ,
                   "kappa"                     : _DefaultKappa , 
                   "objectiveFunction"         : None          , 
                   "removeRotationTranslation" : True          ,          
                   "rho"                       : _DefaultRho   ,
                   "trajectories"              : None          ,
                   "trajectory"                : None          ,   
                   "x"                         : None          ,   
                   "x0"                        : None          ,   
                   "xn"                        : None          ,   
                   "y"                         : None          }
    attributes.update ( ObjectiveFunction.defaultAttributes )

    def __init__ ( self, system, trajectory, **keywordArguments ):
        """Constructor."""
        # . Set defaults for all options.
        for ( key, value ) in self.__class__.attributes.iteritems ( ): setattr ( self, key, value )
        # . Set values of keyword arguments.
        for ( key, value ) in keywordArguments.iteritems ( ): setattr ( self, key, value )
        # . Define the system geometry object function.
        self.objectiveFunction = SystemGeometryObjectiveFunction.FromSystem ( system )
        # . Define the trajectory - a check should be made here to ensure that the trajectory is a direct access one with coordinates!
        # . Reading and writing the header/footer is not required but maybe should be included for completeness?
        self.trajectory = trajectory
        # . Get some scratch space.
        self.g  = Real1DArray.WithExtent ( self.objectiveFunction.NumberOfVariables ( ) )
        self.gd = Real1DArray.WithExtent ( self.objectiveFunction.NumberOfVariables ( ) + self.NumberOfVariables ( ) )
        self.x  = Real1DArray.WithExtent ( self.objectiveFunction.NumberOfVariables ( ) )
        self.x0 = Real1DArray.WithExtent ( self.objectiveFunction.NumberOfVariables ( ) )
        self.xn = Real1DArray.WithExtent ( self.objectiveFunction.NumberOfVariables ( ) )
        self.y  = Real1DArray.WithExtent ( self.objectiveFunction.NumberOfVariables ( ) )
        self.g.Set  ( 0.0 )
        self.gd.Set ( 0.0 )
        self.x.Set  ( 0.0 )
        self.x0.Set ( 0.0 )
        self.xn.Set ( 0.0 )
        self.y.Set  ( 0.0 )
        # . Get the first structure on the trajectory.
        self.trajectory.RestoreOwnerData ( index = 0 )
        self.objectiveFunction.VariablesGet ( self.x0 )
        # . Make the first structure the reference structure.
        # . It is assumed that the remaining structures on the trajectory are already
        # . orientated with respect to this one. This is the case for trajectories
        # . from linear-interpolation or expansion or previous SAW calculations but
        # . may not be otherwise.
        if self.removeRotationTranslation: self.objectiveFunction.RemoveRotationTranslation ( reference = self.objectiveFunction.system.coordinates3 )
        # . Get the last structure on the trajectory.
        self.trajectory.RestoreOwnerData ( index = self.trajectory.frames - 1 )
        self.objectiveFunction.VariablesGet ( self.xn )
        self.trajectories = []

    def ApplyLinearConstraints ( self, vector ):
        """Apply linear constraints to a vector."""
        if self.objectiveFunction.linearVectors is not None:
            nframes    = len ( self.trajectory )
            nvariables = self.objectiveFunction.NumberOfVariables ( )
            for f in range ( 1, nframes - 1 ):
                inc = ( f - 1 ) * nvariables
                for i in range ( nvariables ): self.x[i] = vector[inc+i]
                self.objectiveFunction.linearVectors.ProjectOutOfArray ( self.x )
                for i in range ( nvariables ): vector[inc+i] = self.x[i]

    def ConstraintTerms ( self, variables, gradients = None ):
        """Evaluate the constraint terms."""
        # . Set some options.
        doGradients = ( gradients is not None )
        nframes    = len ( self.trajectory )
        nvariables = self.objectiveFunction.NumberOfVariables ( )
        # . Initialization.
        of = 0.0
        if doGradients:
            dofda = 0.0
            gradients.Set ( 0.0 )
        # . Check that terms need to be done.
        if ( self.gamma != 0.0 ) or ( self.rho != 0.0 ):
            # . Get the interpoint distances.
            self.distances1 = []
            for f in range ( 1, nframes ):
                if f == 1: a = self.x0
                else:
                    inc = ( f - 2 ) * nvariables
                    a   = self.x
                    for i in range ( nvariables ): a[i] = variables[i+inc]
                if f == ( nframes - 1 ): b = self.xn
                else:
                    inc = ( f - 1 ) * nvariables
                    b   = self.y
                    for i in range ( nvariables ): b[i] = variables[i+inc]
                d1 = math.sqrt ( self.objectiveFunction.DistanceSquared ( a, b, self.g ) )
                self.distances1.append ( d1 )
                if doGradients:
                    self.g.Scale ( 0.5 / d1 )
                    inc = ( f - 1 ) * nvariables
                    for i in range ( nvariables ): self.gd[i+inc] = self.g[i]
            # . Get the average distance.
            da = sum ( self.distances1 ) / float ( nframes - 1 )
            # . Get the gamma term.
            for f in range ( 1, nframes ):
                d1  = self.distances1[f-1]
                d1a = ( d1 - da )
                df  = self.gamma * d1a
                of += df * d1a
                if doGradients:
                    df    *= 2.0
                    dofda -= df
                    incg   = ( f - 1 ) * nvariables
                    if f != 1:
                        inc = ( f - 2 ) * nvariables
                        for i in range ( nvariables ): gradients[i+inc] += df * self.gd[i+incg]
                    if f != ( nframes - 1 ):
                        inc = ( f - 1 ) * nvariables
                        for i in range ( nvariables ): gradients[i+inc] -= df * self.gd[i+incg]
            # . Get the rho/kappa term.
            self.distances2 = []
            for f in range ( 2, nframes ):
                if f == 2: a = self.x0
                else:
                    inc = ( f - 3 ) * nvariables
                    a   = self.x
                    for i in range ( nvariables ): a[i] = variables[i+inc]
                if f == ( nframes - 1 ): b = self.xn
                else:
                    inc = ( f - 1 ) * nvariables
                    b   = self.y
                    for i in range ( nvariables ): b[i] = variables[i+inc]
                d2  = self.objectiveFunction.DistanceSquared ( a, b, self.g )
                df  = self.rho * math.exp ( - self.kappa * d2 / ( da**2 ) )
                of += df / self.kappa
                self.distances2.append ( d2 )
                if doGradients:
                    dofda +=   2.0 * d2 * df / ( da**3 )
                    df    *= - 1.0 / ( da**2 )
                    self.g.Scale ( df )
                    if f != 2:
                        inc = ( f - 3 ) * nvariables
                        for i in range ( nvariables ): gradients[i+inc] += self.g[i]
                    if f != ( nframes - 1 ):
                        inc = ( f - 1 )* nvariables
                        for i in range ( nvariables ): gradients[i+inc] -= self.g[i]
            # . Finish the gradients.
            if doGradients:
                # . dE/da terms.
                dofda /= float ( nframes - 1 )
                for f in range ( 1, nframes ):
                    incg = ( f - 1 ) * nvariables
                    if f != 1:
                        inc = ( f - 2 ) * nvariables
                        for i in range ( nvariables ): gradients[i+inc] += dofda * self.gd[i+incg]
                    if f != ( nframes - 1 ):
                        inc = ( f - 1 ) * nvariables
                        for i in range ( nvariables ): gradients[i+inc] -= dofda * self.gd[i+incg]
                # . This is inefficient as linear projection is done twice.
                self.ApplyLinearConstraints ( gradients )
        return of

    def Function ( self, variables ):
        """Evaluate the function."""
        nframes    = len ( self.trajectory )
        nvariables = self.objectiveFunction.NumberOfVariables ( )
        # . Constraints.
        of = self.ConstraintTerms ( variables )
        # . Energies.
        self.energies = []
        for f in range ( 1, nframes - 1 ):
            inc = ( f - 1 ) * nvariables
            for i in range ( nvariables ): self.x[i] = variables[inc+i]
            e   = self.objectiveFunction.Function ( self.x )
            of += e
            self.energies.append ( e )
        return of

    def FunctionGradients ( self, variables, gradients ):
        """Evaluate the function and gradients."""
        nframes    = len ( self.trajectory )
        nvariables = self.objectiveFunction.NumberOfVariables ( )
        # . Constraints.
        of = self.ConstraintTerms ( variables, gradients = gradients )
        # . Energies.
        for f in range ( 1, nframes - 1 ):
            inc = ( f - 1 ) * nvariables
            for i in range ( nvariables ): self.x[i] = variables[inc+i]
            of += self.objectiveFunction.FunctionGradients ( self.x, self.g )
            for i in range ( nvariables ): gradients[inc+i] += self.g[i]
        return of

    def NumberOfVariables ( self ):
        """Return the number of variables."""
        return self.objectiveFunction.NumberOfVariables ( ) * ( len ( self.trajectory ) - 2 )

    def PathSummary ( self, log = logFile ):
        """Output a path summary."""
        if LogFileActive ( log ):
            # . Calculate all data for the current path.
            variables = Real1DArray.WithExtent ( self.NumberOfVariables ( ) )
            self.VariablesGet ( variables )
            self.Function ( variables )
            # . Get the energies of the first and last structures.
            self.energies[0:0] = [ self.objectiveFunction.Function ( self.x0 ) ]
            self.energies.append ( self.objectiveFunction.Function ( self.xn ) )
            # . Modify the distance lists.
            self.distances1.append ( None )
            self.distances2.extend ( [ None, None ] )
            # . Output.
            table  = log.GetTable ( columns = [ 10, 20, 20, 20 ] )
            table.Start   ( )
            table.Title   ( "Path Summary" )
            table.Heading ( "Structure"    )
            table.Heading ( "Energy"       )
            table.Heading ( "Dist(i,i+1)"  )
            table.Heading ( "Dist(i,i+2)"  )
            for ( i, ( e, d1, d2 ) ) in enumerate ( zip ( self.energies, self.distances1, self.distances2 ) ):
                table.Entry ( "{:d}"  .format ( i ) )
                table.Entry ( "{:.3f}".format ( e ) )
                if d1 is None: table.Entry ( "" )
                else:          table.Entry ( "{:.3f}".format ( d1 ) )
                if d2 is None: table.Entry ( "" )
                else:          table.Entry ( "{:.3f}".format ( d2 ) )
            table.Stop ( )

    def Summary ( self, log = logFile, pageWidth = 100, valueWidth = 14 ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( pageWidth = pageWidth, valueWidth = valueWidth )
            summary.Start ( "Self-Avoiding Walk Options" )
            summary.Entry ( "Gamma" , "{:.3g}".format ( self.gamma ) )
            summary.Entry ( "Kappa" , "{:.3g}".format ( self.kappa ) )
            summary.Entry ( "Rho"   , "{:.3g}".format ( self.rho   ) )
            summary.Entry ( "Remove Rotation Translation" , "{!r}".format ( self.removeRotationTranslation ) )
            summary.Stop ( )

    def VariablesAllocate ( self ):
        """Return a variables object of the correct size."""
        if self.NumberOfVariables ( ) > 0:
            variables = Real1DArray.WithExtent ( self.NumberOfVariables ( ) )
            variables.Set ( 0.0 )
            return variables
        else:
            return None

    def VariablesGet ( self, variables ):
        """Fill the variable array."""
        nframes    = len ( self.trajectory )
        nvariables = self.objectiveFunction.NumberOfVariables ( )
        for f in range ( 1, nframes - 1 ):
            self.trajectory.RestoreOwnerData ( index = f )
            self.objectiveFunction.VariablesGet ( self.x )
            inc = ( f - 1 ) * nvariables
            for i in range ( nvariables ): variables[inc+i] = self.x[i]

    def VariablesPut ( self, variables ):
        """Empty the variable array."""
        nframes    = len ( self.trajectory )
        nvariables = self.objectiveFunction.NumberOfVariables ( )
        for f in range ( 1, nframes - 1 ):
            inc = ( f - 1 ) * nvariables
            for i in range ( nvariables ): self.x[i] = variables[inc+i]
            self.objectiveFunction.VariablesPut ( self.x )
            self.trajectory.WriteOwnerData ( index = f )

#===================================================================================================================================
# . SAW path optimization.
#===================================================================================================================================
def SAWOptimize_SystemGeometry ( system, trajectory,
                                 gamma                     = _DefaultGamma ,
                                 kappa                     = _DefaultKappa ,
                                 rho                       = _DefaultRho   ,
                                 removeRotationTranslation = True          ,
                                 minimizer                 = None          ,
                                 log                       = logFile       ):
    """Self-avoiding walk path optimization."""
    # . Create an object function.
    of = _SAWObjectiveFunction ( system, trajectory, gamma = gamma, kappa = kappa, removeRotationTranslation = removeRotationTranslation, rho = rho )
    of.Summary ( log = log )
    # . Set up the minimizer.
    if minimizer is None:
        minimizer = LBFGSMinimizer ( **_DefaultMinimizerOptions )
        minimizer.Summary ( log = log )
    # . Optimization.
    report = minimizer.Iterate ( of, log = log )
    # . Finish up.
    of.PathSummary ( log = log )
    return report

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
