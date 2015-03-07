#-------------------------------------------------------------------------------
# . File      : ChainOfStatesObjectiveFunction.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Chain-of-states objective function."""

# . This is a misnomer as an objective function cannot be defined for many of these methods.
# . Nevertheless, the object serves as an appropriate interface to the system being optimized.

from pMolecule import MultiLayerSystemGeometryObjectiveFunction

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ChainOfStatesObjectiveFunction ( MultiLayerSystemGeometryObjectiveFunction ):
    """The chain-of-states objective function."""

    # . Attributes.
    defaultAttributes = { "imageTrajectory" : None ,
                          "numberOfImages"  : 0    }
    defaultAttributes.update ( MultiLayerSystemGeometryObjectiveFunction.defaultAttributes )

    def DumpImage ( self, image ):
        """Dump an image."""
        self.imageTrajectory.WriteOwnerData ( index = image )

    def FinalizeImages ( self ):
        """Finalize image data."""
        self.imageTrajectory.Close ( )

    def GetGradients ( self, gradients ):
        """Get the gradients (without applying linear constaints)."""
        self.system.configuration.gradients3.CopyToToVector ( gradients, selection = self.freeAtoms )
        if self.variableWeights is not None: gradients.Divide ( self.variableWeights )

    def InitializeImages ( self, imageTrajectory ):
        """Initialize image data given a trajectory."""
        self.imageTrajectory = imageTrajectory
        self.numberOfImages  = len ( imageTrajectory )
        if ( self.NumberOfImages ( ) <= 2 ): raise ValueError ( "Invalid number of images on trajectory: {:d}.".format ( self.NumberOfImages ( ) ) )

    def LoadImage ( self, image ):
        """Load an image."""
        self.imageTrajectory.RestoreOwnerData ( index = image )

    def NumberOfImages ( self ):
        """Return the number of images."""
        return self.numberOfImages

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
