#-------------------------------------------------------------------------------
# . File      : AmberTrajectoryFileWriter.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for writing Amber Trajectory files.

Taken from the Amber web page:

AMBER trajectory (coordinate or velocity) file specification
This file is optionally written during dynamics in SANDER or GIBBS.

FORMAT(20A4) ITITL
  ITITL  : the title of the current run, from the AMBER
           parameter/topology file

The following snapshot is written every NTWX steps in the trajectory
(specified in the control input file):

FORMAT(10F8.3)  (X(i), Y(i), Z(i), i=1,NATOM)
  X,Y,Z  : coordinates or velocities (velocity units: Angstroms per 1/20.455 ps)

If constant pressure (in 4.1, also constant volume)
For each snapshot:

FORMAT(10F8.3)  BOX(1), BOX(2), BOX(3)
  BOX    : size of periodic box
"""

from pCore                    import TextFileWriter
from SystemGeometryTrajectory import SystemGeometryTrajectory

#===================================================================================================================================
# . AmberTrajectory file writer class.
#===================================================================================================================================
class AmberTrajectoryFileWriter ( TextFileWriter ):
    """Class for writing Amber Trajectory files."""

    def __init__ ( self, filename, owner, title = None ):
        """Constructor."""
        super ( AmberTrajectoryFileWriter, self ).__init__ ( filename )
        self.Open ( )
        if ( title is None ):
            self.file.write ( "Amber Trajectory File\n" )
        else:
            self.file.write ( title[0:min ( 80, len ( title ) )] + "\n" )
        # . Owner checks.
        self.frames = 0
        self.QBOX   = False
        self.owner  = owner
        if hasattr ( owner, "symmetry" ) and ( owner.symmetry is not None ):
            crystalClass = getattr ( owner.symmetry, "crystalClass", None )
            self.QBOX = ( crystalClass is not None ) and crystalClass.IsOrthogonal ( )

    def __len__ ( self ):
        return self.frames

    def WriteFooter ( self ):
        """Write the trajectory footer."""
        pass

    def WriteHeader ( self, **keywordArguments ):
        """Write the trajectory header."""
        pass

    def WriteOwnerData ( self ):
        """Write data from the owner to a frame."""
        xyz = self.owner.coordinates3
        n   = 0
        for i in range ( xyz.rows ):
            for j in range ( 3 ):
                self.file.write ( "{:8.3f}".format ( xyz[i,j] ) )
                n += 1
                if ( n >= 10 ):
                    n = 0
                    self.file.write ( "\n" )
        if n != 0: self.file.write ( "\n" )
        if self.QBOX:
            sp = self.owner.symmetryParameters
            self.file.write ( "{:8.3f}{:8.3f}{:8.3f}\n".format ( sp.a, sp.b, sp.c ) )
        self.frames += 1

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def AmberTrajectory_FromSystemGeometryTrajectory ( outPath, inPath, system ):
    """Convert a SystemGeometryTrajectory to an AMBER trajectory."""
    # . Define the input and output trajectories.
    inTrajectory  = SystemGeometryTrajectory  ( inPath , system, mode = "r" )
    outTrajectory = AmberTrajectoryFileWriter ( outPath, system )
    # . Loop over frames.
    while inTrajectory.RestoreOwnerData ( ): outTrajectory.WriteOwnerData ( )
    inTrajectory.Close  ( )
    outTrajectory.Close ( )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
