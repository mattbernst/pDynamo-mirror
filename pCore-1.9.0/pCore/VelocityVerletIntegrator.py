#-------------------------------------------------------------------------------
# . File      : VelocityVerletIntegrator.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Define classes for velocity Verlet dynamics."""

import math

from MultiDimensionalDynamics import MultiDimensionalDynamics, MultiDimensionalDynamicsState

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Temperature-handling options.
_TemperatureHandlingOptions = ( "Constant", "Exponential", "Linear" )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class VelocityVerletIntegratorState ( MultiDimensionalDynamicsState ):
    """State for the velocity Verlet integrator."""
    pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class VelocityVerletIntegrator ( MultiDimensionalDynamics ):
    """Class for velocity Verlet dynamics."""

    # . Default attributes.
    defaultAttributes = { "facX1"                     :  0.0 ,
                          "facX2"                     :  0.0 ,
                          "facV"                      :  0.0 ,
                          "temperatureScale"          :  0.0 ,
                          "temperatureScaleFrequency" :    0 ,
                          "temperatureScaleOption"    : None ,
                          "temperatureStart"          : None ,
                          "temperatureStop"           : None }
    defaultAttributes.update ( MultiDimensionalDynamics.defaultAttributes )

    # . Default attribute names.
    defaultAttributeNames = { "Temperature Handling"        : "temperatureScaleOption"    ,
                              "Temperature Scale"           : "temperatureScale"          ,
                              "Temperature Scale Frequency" : "temperatureScaleFrequency" ,
                              "Temperature Start"           : "temperatureStart"          ,
                              "Temperature Stop"            : "temperatureStop"           }
    defaultAttributeNames.update ( MultiDimensionalDynamics.defaultAttributeNames )

    # . State class.
    stateObject = VelocityVerletIntegratorState

    def CalculateIntegrationConstants ( self, state ):
        """Calculate constants for the integration."""
        self.facX1 =       self.timeStep
        self.facX2 = 0.5 * self.timeStep**2
        self.facV  = 0.5 * self.timeStep

    def Iteration ( self, state ):
        """Perform one dynamics step."""
        state.numberOfIterations += 1
        state.time               += self.timeStep
        # . Calculate the new variables.
        state.x.AddScaledArray ( self.facX1, state.v )
        state.x.AddScaledArray ( self.facX2, state.a )
        # . Calculate the half-velocities.
        state.v.AddScaledArray ( self.facV, state.a )
        # . Calculate the function and accelerations at the new point.
        state.f = state.objectiveFunction.Accelerations ( state.x, state.a )
        # . Complete the velocity calculation.
        state.v.AddScaledArray ( self.facV, state.a )
        # . Calculate other quantities and perform temperature scaling if necessary.
        ( state.kineticEnergy, state.temperature ) = state.objectiveFunction.Temperature ( state.v )
        if ( self.temperatureScaleOption is not None ) and ( state.numberOfIterations % self.temperatureScaleFrequency == 0 ):
            scale = self.TargetTemperature ( state.time ) / state.temperature
            state.v.Scale ( math.sqrt ( scale ) )
            state.kineticEnergy *= scale
            state.temperature   *= scale
        # . Finish up.
        state.totalEnergy = state.f + state.kineticEnergy

    def Label ( self ): return "Velocity Verlet Integrator"

    def TargetTemperature ( self, time ):
        """Get the target temperature for temperature scaling."""
        if   self.temperatureScaleOption == "Constant"    : tNeeded = self.temperatureStart
        elif self.temperatureScaleOption == "Exponential" : tNeeded = self.temperatureStart * math.exp ( self.temperatureScale * time )
        elif self.temperatureScaleOption == "Linear"      : tNeeded = self.temperatureStart +            self.temperatureScale * time
        return max ( 0.0, tNeeded )

    def TemperatureHandlingOptions ( self ):
        """Set up the temperature-handling options."""
        # . Convert the option to capitalized form.
        if isinstance ( self.temperatureScaleOption, basestring ): self.temperatureScaleOption = self.temperatureScaleOption.capitalize ( )
        # . No handling.
        if ( self.temperatureScaleFrequency <= 0 ) or ( self.temperatureScaleFrequency >= self.maximumIterations ) or \
           ( self.temperatureScaleOption is None ) or ( self.temperatureScaleOption not in _TemperatureHandlingOptions ):
            self.temperatureScaleOption = None
        # . Constant temperature.
        elif self.temperatureScaleOption == "Constant":
            self.temperatureStop = self.temperatureStart
        # . Exponential changes.
        elif self.temperatureScaleOption == "Exponential":
            self.temperatureScale = math.log ( self.temperatureStop / self.temperatureStart ) / ( float ( self.maximumIterations ) * self.timeStep )
        # . Linear changes.
        elif self.temperatureScaleOption == "Linear":
            self.temperatureScale =          ( self.temperatureStop - self.temperatureStart ) / ( float ( self.maximumIterations ) * self.timeStep )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
