#-------------------------------------------------------------------------------
# . File      : AmberTrajectoryFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for reading Amber Trajectory files."""

from pCore                    import TextFileReader
from SystemGeometryTrajectory import SystemGeometryTrajectory

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The number and the widths of symmetry elements on a line.
_SYMMETRYNUMBER = 3
_SYMMETRYWIDTH  = 8

# . The number and the widths of XYZ elements on a line.
_XYZNUMBER = 10
_XYZWIDTH  =  8

#===================================================================================================================================
# . AmberTrajectory file reader class.
#===================================================================================================================================
class AmberTrajectoryFileReader ( TextFileReader ):
    """Class for reading Amber Trajectory files."""

    def __init__ ( self, filename, owner ):
        """Constructor."""
        super ( AmberTrajectoryFileReader, self ).__init__ ( filename )
        self.Open ( )
        # . Owner checks.
        self.frames = 0
        self.log    = None
        self.owner  = owner
        self.QBOX   = False
        if hasattr ( owner, "symmetry" ) and ( owner.symmetry is not None ):
            crystalClass = getattr ( owner.symmetry, "crystalClass", None )
            self.QBOX    = ( crystalClass is not None ) and crystalClass.IsOrthogonal ( )

    def __len__ ( self ):
        return self.frames

    def ReadFooter ( self ):
        """Read the trajectory footer."""
        pass

    def ReadHeader ( self ):
        """Read the trajectory header."""
        # . Activate the trajectory if necessary.
        if not self.QACTIVE:
            self.frames = 0
            self.Open ( )
        # . Get the title.
        self.title = self.GetLine ( )

    def RestoreOwnerData ( self ):
        """Restore data from a frame to the owner."""
        try:
            natoms = len ( self.owner.atoms )
            xyz    = self.owner.coordinates3
            items  = self.GetFixedFormatArray ( 3 * natoms, _XYZNUMBER, _XYZWIDTH, converter = float, default = 0.0 )
            for n in range ( natoms ):
                for i in range ( 3 ): xyz[n,i] = items[3*n+i]
            if self.QBOX:
                items = self.GetFixedFormatArray ( 3, _SYMMETRYNUMBER, _SYMMETRYWIDTH, converter = float, default = 0.0 )
                self.owner.symmetryParameters.SetCrystalParameters ( items[0], items[1], items[2], 90.0, 90.0, 90.0 )
            self.frames += 1
            return True
        except EOFError:
            return False

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def AmberTrajectory_ToSystemGeometryTrajectory ( inPath, outPath, system ):
    """Convert an AMBER trajectory to a SystemGeometryTrajectory."""
    # . Define the input and output trajectories.
    inTrajectory  = AmberTrajectoryFileReader ( inPath , system )
    outTrajectory = SystemGeometryTrajectory  ( outPath, system, mode = "w" )
    # . Loop over frames.
    inTrajectory.ReadHeader ( )
    while inTrajectory.RestoreOwnerData ( ): outTrajectory.WriteOwnerData ( )
    inTrajectory.Close  ( )
    outTrajectory.Close ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
