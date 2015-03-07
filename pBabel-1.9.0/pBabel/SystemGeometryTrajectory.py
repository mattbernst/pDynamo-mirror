#-------------------------------------------------------------------------------
# . File      : SystemGeometryTrajectory.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for manipulating system geometry trajectories.

Trajectory opening modes (same as C's fopen from stdio.h):

r      Open for reading.
r+     Open for reading and writing.
w      Create a new trajectory for writing. Existing
       trajectories are cleared.
w+     Open for reading and writing. Existing trajectories
       are cleared.
a      Open for writing. Create a new trajectory if one does
       not exist.
a+     Open for reading and writing. Create a new trajectory
       if one does not exist.

This file contains base classes only.
"""

import glob, os, os.path

from pCore     import Clone, Coordinates3, Pickle, PickleFileExtension, Unpickle
from pMolecule import SymmetryParameters

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Frame naming.
_FramePrefix  = "frame"
_FramePostfix = PickleFileExtension

# . Allowed modes.
_ReadableModes   = [ "a+", "r", "r+", "w", "w+" ]
_TrajectoryModes = [ "a", "a+", "r", "r+", "w", "w+" ]
_WritableModes   = [ "a", "a+", "r+", "w", "w+" ]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SystemGeometryTrajectory ( object ):
    """Class for system geometry trajectories."""

    defaultattributes = { "frames"                :     0 ,
                          "isLocked"              : False ,
                          "isReadable"            : False ,
                          "isWritable"            : False ,
                          "hasSymmetryParameters" : False ,
                          "mode"                  :  None ,
                          "owner"                 :  None ,
                          "path"                  :  None ,
                          "position"              :    -1 }

    def __getitem__ ( self, index ):
        """Get an item."""
        if ( index < 0 ) or ( index >= len ( self ) ): raise IndexError
        else:                                          return self.ReadFrame ( frame = index )

    def __init__ ( self, path, owner, mode = None ):
        """Constructor."""
        self.__dict__.update ( self.__class__.defaultattributes )
        # . Arguments.
        self.mode  = mode
        self.path  = path
        self.owner = owner
        # . Open the trajectory.
        self.Open ( )

    def __len__ ( self ):
        return self.frames

    def Close ( self ):
        """Close the trajectory."""
        pass

    def LinearlyExpand ( self, ninsert ):
        """Expand a trajectory by inserting points determined by linear interpolation.

        |ninsert| is the number of points to insert between each pair of existing points.
        """
        if not self.isWritable: raise IOError ( "Writing to trajectory that is not writeable." )
        if ( ninsert > 0 ) and ( self.frames > 1 ):
            # . Check for fixed atoms and get the first frame if necessary.
            hasFixedAtoms = ( self.owner.hardConstraints is not None ) and ( self.owner.hardConstraints.NumberOfFixedAtoms ( ) > 0 )
            if not hasFixedAtoms: reference = self.ReadFrame ( frame = 0, instance = Coordinates3 )
            # . Calculate the number of new frames.
            nold = self.frames
            nnew = nold + ( nold - 1 ) * ninsert
            # . Get space for the step.
            dxyz = Coordinates3.WithExtent ( len ( self.owner.atoms ) )
            dxyz.Set ( 0.0 )
            # . Loop over the frames in reverse order.
            xyz = self.ReadFrame ( frame = nold - 1, instance = Coordinates3 )
            if not hasFixedAtoms: xyz.Superimpose ( reference )
            inew = nnew
            for n in range ( nold - 2, -1, -1 ):
                # . Get the next reorientated frame.
                start = self.ReadFrame ( frame = n, instance = Coordinates3 )
                if not hasFixedAtoms: start.Superimpose ( reference )
                # . Get the step.
                start.CopyTo ( dxyz )
                dxyz.AddScaledMatrix ( -1.0, xyz )
                dxyz.Scale ( 1.0 / float ( ninsert + 1 ) )
                # . Write out the frames.
                inew -= 1
                self.WriteFrame ( xyz, frame = inew )
                for i in range ( ninsert ):
                    inew -= 1
                    xyz.AddScaledMatrix ( 1.0, dxyz )
                    self.WriteFrame ( xyz, frame = inew )
#                    print i, inew
                # . Make start the new stop.
                xyz = start
            # . Don't need to write out the last (first) frame as this is unchanged.
            # . Reset the number of frames.
            self.frames = nnew

    @classmethod
    def LinearlyInterpolate ( selfClass, path, system, npoints, point0, pointn ):
        """Generate structures on a trajectory by linear interpolation between two end points."""
        # . Check the number of points.
        npoints = max ( npoints, 2 )
        # . Check for fixed atoms.
        hasFixedAtoms = ( system.hardConstraints is not None ) and ( system.hardConstraints.NumberOfFixedAtoms ( ) > 0 )
        # . Create the trajectory.
        self = selfClass ( path, system, mode = "w" )
        # . Create an intermediate array.
        xyz  = Clone ( pointn )
        if not hasFixedAtoms: xyz.Superimpose ( point0 )
        # . Find the step.
        dxyz = Clone ( xyz )
        dxyz.AddScaledMatrix ( -1.0, point0 )
        dxyz.Scale ( 1.0 / float ( npoints - 1 ) )
        # . First point.
        self.WriteFrame ( point0 )
        # . Intermediate and last points.
        point0.CopyTo ( xyz )
        for i in range ( npoints - 1 ):
            xyz.AddScaledMatrix ( 1.0, dxyz )
            self.WriteFrame ( xyz )
        return self

    def Lock ( self ):
        """Lock the trajectory so that it can only be accessed by one process at a time.

        Locks are of two types - reading or writing. Reading locks mean that other processes
        that only read the trajectory can also read it. Writing locks block both reading and
        writing.
        """
        pass

    def Open ( self ):
        """Open the trajectory."""
        # . Basic mode checks.
        if self.mode not in [ None ] + _TrajectoryModes: raise IOError ( "Unrecognized mode: " + self.mode + "." )
        # . Check the path and assign a default mode.
        pathExists = os.access ( self.path, os.F_OK )
        if pathExists:
            if not os.path.isdir ( self.path ): raise IOError ( "Trajectory exists that is not a directory." )
            if self.mode is None: self.mode = "r"
        else:
            if   self.mode is None: self.mode = "w"
            elif self.mode in ( "r", "r+" ): raise IOError ( "Read mode specified for a trajectory that does not exist." )
            os.mkdir ( self.path )
        # . Check for readability and writeability.
        isReadable = os.access ( self.path, os.R_OK )
        isWritable = os.access ( self.path, os.W_OK )
        # . Set the mode flags.
        self.isReadable = ( self.mode in _ReadableModes  )
        self.isWritable = ( self.mode in _WritableModes )
        # . Check for consistency with the modes.
        if self.isReadable and not isReadable: raise IOError ( "Read mode specified for a trajectory that is unreadable."   )
        if self.isWritable and not isWritable: raise IOError ( "Write mode specified for a trajectory that is unwriteable." )
        # . Check the contents of an existing trajectory.
        if pathExists:
            frames = glob.glob ( os.path.join ( self.path, _FramePrefix + "*" + _FramePostfix ) )
            # . Clear the files for write modes.
            if self.mode in [ "w", "w+" ]:
                for frame in frames: os.remove ( frame )
            # . Set the number of frames.
            else: self.frames = len ( frames )
        # . Check to see if the owner has symmetry parameters.
        self.hasSymmetryParameters = hasattr ( self.owner, "symmetry" ) and ( self.owner.symmetry is not None )

    def ReadFooter ( self ):
        """Read the trajectory footer."""
        return {}

    def ReadFrame ( self, frame = None, identifier = None, instance = None ):
        """Read a frame."""
        if self.isReadable:
            data = Unpickle ( os.path.join ( self.path, _FramePrefix + "{:d}".format ( frame ) + _FramePostfix ) )
            if ( instance is None ) or isinstance ( data, instance ):
                return data
            else:
                for item in data:
                    if isinstance ( item, instance ): return item
        else: raise IOError ( "Reading from trajectory that is not readable." )

    def ReadHeader ( self ):
        """Read the trajectory header."""
        return {}

    def RestoreOwnerData ( self, index = -1 ):
        """Restore data from a frame to the owner."""
        # . Initialization.
        data = None
        # . Restore data from a specific frame.
        if index >= 0:
            data = self.ReadFrame ( frame = index )
            self.position = index + 1
        # . Restore data from the next frame in the sequence.
        else:
            self.position += 1
            if self.position >= self.frames:
                self.position = -1
            else:
                data = self.ReadFrame ( frame = self.position )
        # . Restore any data.
        if data is None:
            return False
        else:
            if isinstance ( data, Coordinates3 ):
                self.owner.coordinates3       = data
            elif isinstance ( data, tuple ):
                self.owner.coordinates3       = data[0]
                self.owner.symmetryParameters = data[1]
            return True

    def Superimpose ( self, reference3 = None, weights = None ):
        """Superimpose the trajectory on a reference structure using the first trajectory frame by default."""
        # . Initialization.
        coordinates3 = self.owner.coordinates3
        # . Save the owner coordinates to restore later.
        saved3 = Clone ( coordinates3 )
        # . Get the reference coordinates.
        if reference3 is coordinates3:
            reference3 = saved3
        elif reference3 is None:
            self.RestoreOwnerData ( index = 0 )
            reference3 = Clone ( coordinates3 )
        # . Loop over trajectory frames.
        results = []
        for index in range ( self.frames ):
            self.RestoreOwnerData ( index = index )
            rms0 = coordinates3.RMSDeviation ( reference3, weights = weights )
            coordinates3.Superimpose ( reference3, weights = weights )
            rms1 = coordinates3.RMSDeviation ( reference3, weights = weights )
            self.WriteOwnerData   ( index = index )
            results.append ( ( rms0, rms1 ) )
        # . Finish up.
        saved3.CopyTo ( coordinates3 )
        return results

    def WriteFooter ( self ):
        """Write the trajectory footer."""
        pass

    def WriteFrame ( self, data, frame = None, identifier = None ):
        """Write a frame."""
        if self.isWritable:
            if frame is None: f = self.frames
            else:             f = frame
            Pickle ( os.path.join ( self.path, _FramePrefix + "{:d}".format ( f ) + _FramePostfix ), data )
            if frame is None: self.frames += 1
        else: raise IOError ( "Writing to trajectory that is not writeable." )

    def WriteHeader ( self, **keywordArguments ):
        """Write the trajectory header."""
        pass

    def WriteOwnerData ( self, index = -1 ):
        """Write data from the owner to a frame."""
        if index >= 0 : frame = index
        else:           frame = None
        if self.hasSymmetryParameters: data = ( self.owner.coordinates3, self.owner.symmetryParameters )
        else:                          data = self.owner.coordinates3
        self.WriteFrame ( data, frame = frame )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
