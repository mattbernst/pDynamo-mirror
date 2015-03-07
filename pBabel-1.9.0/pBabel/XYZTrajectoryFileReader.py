#-------------------------------------------------------------------------------
# . File      : XYZTrajectoryFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for reading XYZ trajectory files."""

# . All frames must contain the same system.

from pCore                    import TextFileReader
from pMolecule                import PeriodicTable
from SystemGeometryTrajectory import SystemGeometryTrajectory

#===================================================================================================================================
# . XYZTrajectory file reader class.
#===================================================================================================================================
class XYZTrajectoryFileReader ( TextFileReader ):
    """Class for reading XYZ trajectory files."""

    def __init__ ( self, filename, owner ):
        """Constructor."""
        super ( XYZTrajectoryFileReader, self ).__init__ ( filename )
        self.Open ( )
        # . Owner checks.
        self.frames     = 0
        self.isComplete = True
        self.log        = None
        self.owner      = owner

    def __len__ ( self ):
        return self.frames

    def _NextFrame ( self ):
        """Get the next frame on the file."""
        # . Number of atoms.
        items           = self.GetTokens ( converters = [ int ] )
        numberOfAtoms   = items[0]
        self.isComplete = False
        # . Title line.
        self.title    = self.GetLine ( )
        # . XYZ lines.
        atoms = self.owner.atoms
        xyz   = self.owner.coordinates3
        isOK  = ( numberOfAtoms == len ( atoms ) )
        if isOK:
            for i in range ( numberOfAtoms ):
                items = self.GetTokens ( converters = [ PeriodicTable.AtomicNumber, float, float, float ] )
                if items[0] != atoms[i].atomicNumber:
                    isOK = False
                    break
                xyz[i,0] = items[1]
                xyz[i,1] = items[2]
                xyz[i,2] = items[3]
            self.isComplete = True
        if not isOK: raise IOError ( "Frame data does not match owner definition." )

    def ReadFooter ( self ):
        """Read the trajectory footer."""
        pass

    def ReadHeader ( self ):
        """Read the trajectory header."""
        # . Activate the trajectory if necessary.
        if not self.QACTIVE:
            self.frames = 0
            self.Open ( )

    def RestoreOwnerData ( self ):
        """Restore data from a frame to the owner."""
        try:
            self._NextFrame ( )
            self.frames += 1
            return True
        except EOFError:
            if self.isComplete: return False
            else: raise IOError ( "Incomplete frame." )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def XYZTrajectory_ToSystemGeometryTrajectory ( inPath, outPath, system ):
    """Convert an XYZ trajectory to a SystemGeometryTrajectory."""
    # . Define the input and output trajectories.
    inTrajectory  = XYZTrajectoryFileReader  ( inPath , system )
    outTrajectory = SystemGeometryTrajectory ( outPath, system, mode = "w" )
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
