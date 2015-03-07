#-------------------------------------------------------------------------------
# . File      : SEAMObjectiveFunction.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
#===================================================================================================================================
# . Contributed by Anirban Bhattacharjee (2013).
#===================================================================================================================================

"""Defines an objective function for SEAM calculations."""

import math

from pCore                           import Clone, logFile, Real1DArray
from SystemGeometryObjectiveFunction import SystemGeometryObjectiveFunction

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Penalty function constants.
_PenaltyFunctionConstant = 5.0

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SEAMObjectiveFunction ( SystemGeometryObjectiveFunction ):
    """SEAM objective function."""

    # . Attributes.
    defaultAttributes = { "c1"      : _PenaltyFunctionConstant , # . Parameters for the penalty function.
                          "c2"      : _PenaltyFunctionConstant ,
                          "g1"      : None                     , # . Spare gradient arrays.
                          "g2"      : None                     ,
                          "method"  : "GP"                     , # . GP (gradient projection) or PF (penalty function).
                          "system1" : None                     ,
                          "system2" : None                     }
    defaultAttributes.update ( SystemGeometryObjectiveFunction.defaultAttributes )

    def __len__ ( self ): return self.nvariables

    def Energies ( self, doGradients = False, log = logFile ):
        """Calculate the energies of the two systems."""
        f1 = self.system1.Energy ( doGradients = doGradients, log = log )
        f2 = self.system2.Energy ( doGradients = doGradients, log = log )
        return ( f1, f2 )

    @classmethod
    def FromSystem ( selfClass, system, electronicState1, electronicState2, method = "GP" ):
        """Constructor given a system."""
        # . Basic setup.
        system.electronicState          = electronicState1
        self                            = super ( SEAMObjectiveFunction, selfClass ).FromSystem ( system )
        self.method                     = method
        # . Define the first system.
        self.system1                    = self.system
        # . Define the second system without coordinates.
        coordinates3                    = self.system1.coordinates3
        symmetryParameters              = self.system1.symmetryParameters
        self.system1.coordinates3       = None
        self.system1.symmetryParameters = None
        self.system2                    = Clone ( self.system1 )
        self.system2.electronicState    = electronicState2
        # . Reset the coordinates for both systems so that they are the same.
        self.system1.coordinates3       = coordinates3
        self.system2.coordinates3       = coordinates3
        if symmetryParameters is not None:
            self.system1.symmetryParameters = symmetryParameters
            self.system2.symmetryParameters = symmetryParameters
        # . Allocate space.
        self.g1 = Real1DArray.WithExtent ( len ( self ) )
        self.g2 = Real1DArray.WithExtent ( len ( self ) )
        # . Finish up.
        return self

    def Function ( self, variables ):
        """Evaluate the function."""
        self.VariablesPut ( variables )
        f1 = self.system1.Energy ( log = None )
        f2 = self.system2.Energy ( log = None )
        if self.method == "GP": f = math.fabs ( f1 - f2 )
        else:                   f = self.PenaltyFunction ( f1, f2, None, None, None )
        return f

    def FunctionGradients ( self, variables, gradients ):
        """Evaluate the function and gradients."""
        # . Function.
        self.VariablesPut ( variables )
        f1 = self.system1.Energy ( doGradients = True, log = None )
        f2 = self.system2.Energy ( doGradients = True, log = None )
        # . Gradients.
        self.GradientsGet ( self.g1 ) ; self.system = self.system2
        self.GradientsGet ( self.g2 ) ; self.system = self.system1
        # . Finish up.
        if self.method == "GP":
            f = math.fabs ( f1 - f2 )
            self.GradientProjection  ( f1, f2, self.g1, self.g2, gradients )
        else:
            f = self.PenaltyFunction ( f1, f2, self.g1, self.g2, gradients )
        return f

    def GradientProjection ( self, f1, f2, g1, g2, gradients ):
        """Harvey's gradient projection method."""
        # . ( g1 - g2 ) in g1.
        g1.AddScaledArray ( -1.0, g2 )
        # . First part of g in gradients.
        g2.CopyTo ( gradients )
        # . f in gradients.
        gradients.AddScaledArray ( ( f1 - f2 ), g1 )
        # . Normalize g1.
        g1.Normalize ( )
        # . Term to project out.
        factor = g1.Dot ( g2 )
        # . Remainder of g.
        gradients.AddScaledArray ( -factor, g1 )

    def PenaltyFunction ( self, f1, f2, g1, g2, gradients ):
        """Penalty function calculation."""
        deltaF    = f1 - f2
        argument  = 1.0 + ( deltaF / self.c2 )**2
        prefactor = self.c1 * self.c2**2
        f         = 0.5 * ( f1 + f2 ) + prefactor * math.log ( argument )
        if gradients is not None:
            c1 = 0.5 + ( 2.0 * self.c1 * deltaF / argument )
            c2 = 1.0 - c1
            g1.Scale ( c1 )
            g1.CopyTo ( gradients )
            gradients.AddScaledArray ( c2, g2 )
        return f

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
