#-------------------------------------------------------------------------------
# . File      : LeapFrogIntegrator.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Define classes for leapfrog Verlet dynamics."""

import math

from MultiDimensionalDynamics import MultiDimensionalDynamics, MultiDimensionalDynamicsState

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LeapFrogIntegratorState ( MultiDimensionalDynamicsState ):
    """State for the leap frog integrator."""

    defaultAttributes = { "pressure" : 0.0 ,
                          "volume"   : 0.0 }
    defaultAttributes.update ( MultiDimensionalDynamicsState.defaultAttributes )

    def DefineAccumulables ( self ):
        """Define the quantities to accumulate."""
        super ( LeapFrogIntegratorState, self ).DefineAccumulables ( )
        self.accumulableAttributes.extend ( [ "pressure" , "volume" ] )
        self.accumulableLabels.extend     ( [ "Pressure" , "Volume" ] )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LeapFrogIntegrator ( MultiDimensionalDynamics ):
    """Class for leap frog dynamics."""

    # . Default attributes.
    defaultAttributes = { "pressure"            :     1.0 ,
                          "pressureControl"     :   False ,
                          "pressureCoupling"    :  2000.0 ,
                          "temperature"         :   300.0 ,
                          "temperatureControl"  :   False ,
                          "temperatureCoupling" :     0.1 }
    defaultAttributes.update ( MultiDimensionalDynamics.defaultAttributes )

    # . Default attribute names.
    defaultAttributeNames = { "Pressure"             : "pressure"            ,
                              "Pressure Control"     : "pressureControl"     ,
                              "Pressure Coupling"    : "pressureCoupling"    ,
                              "Temperature"          : "temperature"         ,
                              "Temperature Control"  : "temperatureControl"  ,
                              "Temperature Coupling" : "temperatureCoupling" }
    defaultAttributeNames.update ( MultiDimensionalDynamics.defaultAttributeNames )

    # . State class.
    stateObject = LeapFrogIntegratorState

    def CheckAttributes ( self ):
        """Check the attributes."""
        super ( LeapFrogIntegrator, self ).CheckAttributes ( )
        self.PressureHandlingOptions ( )

    def Initialize ( self, state ):
        """Initialization before iteration."""
        if self.pressureControl: state.objectiveFunction.DefinePressure ( )
        super ( LeapFrogIntegrator, self ).Initialize ( state )
        ( state.pressure, state.volume ) = state.objectiveFunction.Pressure ( state.kineticEnergy )

    def Iteration ( self, state ):
        """Perform one dynamics step."""
        state.numberOfIterations += 1
        state.time               += self.timeStep
        # . Calculate the new velocities.
        state.v.AddScaledArray ( self.timeStep, state.a )
        # . Temperature scaling.
        self.TemperatureScale ( state )
        # . Calculate the new variables.
        state.x.AddScaledArray ( self.timeStep, state.v )
        # . Pressure scaling.
        self.PressureScale    ( state )
        # . Calculate the function and accelerations at the new point.
        state.f = state.objectiveFunction.Accelerations ( state.x, state.a )
        # . Calculate the temperature, pressure, etc.
        ( state.kineticEnergy, state.temperature  ) = state.objectiveFunction.Temperature ( state.v             )
        ( state.pressure     , state.volume       ) = state.objectiveFunction.Pressure    ( state.kineticEnergy )
        state.totalEnergy = state.f + state.kineticEnergy

    def Label ( self ): return "Leapfrog Verlet Integrator"

    def PressureHandlingOptions ( self ):
        """Set up the pressure-handling options."""
        self.pressureControl = self.pressureControl and ( self.pressure > 0.0 ) and ( self.pressureCoupling > 0.0 )

    def PressureScale ( self, state ):
        """Do pressure scaling."""
        if self.pressureControl:
            zetaP = 1.0 - ( self.timeStep / ( 3.0 * self.pressureCoupling ) ) * ( self.pressure - state.pressure ) # . Use first-order formula.
            state.objectiveFunction.VolumeScale ( zetaP )
            state.x.Scale                       ( zetaP )

    def TemperatureHandlingOptions ( self ):
        """Set up the temperature-handling options."""
        self.temperatureControl = self.temperatureControl and ( self.temperature > 0.0 ) and ( self.temperatureCoupling > 0.0 )

    def TemperatureScale ( self, state ):
        """Do temperature scaling."""
        if self.temperatureControl:
            zetaT = math.sqrt ( 1.0 + ( self.timeStep / self.temperatureCoupling ) * ( ( self.temperature / state.temperature ) - 1.0 ) )
            state.v.Scale ( zetaT )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
