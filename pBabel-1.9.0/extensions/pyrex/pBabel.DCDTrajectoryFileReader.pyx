"""Classes and functions for reading DCD trajectory files."""

from pCore                    import logFile, LogFileActive
from SystemGeometryTrajectory import SystemGeometryTrajectory

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DCDTrajectoryFileReader:
    """DCD trajectory file reader."""

    def __dealloc__ ( self ):
        """Finalization."""
        self.Close ( )

    def __getmodule__ ( self ): return "pBabel.DCDTrajectoryFileReader"

    def __init__ ( self, path, owner ):
        """Constructor."""
        self._Initialize     ( )
        self._Allocate       ( )
        self.path  = path
        self.owner = owner
        self.AssignOwnerData ( )
        self.Open            ( )

    def __len__ ( self ): return self.numberOfFrames

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = DCDHandle_Allocate ( )
        if self.cObject == NULL: DCDStatus_Check ( CDCDStatus_MemoryAllocationFailure )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOpen  = False
        self.owner   = None
        self.path    = None

    def AssignOwnerData ( self ):
        """Assign owner data to the trajectory."""
        cdef Coordinates3       data3
        cdef SymmetryParameters symmetryParameters
        # . Get objects.
        data3              = self.owner.coordinates3
        symmetryParameters = self.owner.symmetryParameters
        # . Assignment.
        if data3              is not None: DCDStatus_Check ( DCDHandle_SetData3              ( self.cObject, data3.cObject              ) )
        if symmetryParameters is not None: DCDStatus_Check ( DCDHandle_SetSymmetryParameters ( self.cObject, symmetryParameters.cObject ) )

    def Close ( self ):
        """Close the file."""
        if self.isOpen:
            DCDRead_Close ( &self.cObject )
            self.isOpen = False

    def Open ( self ):
        """Open the file."""
        cdef char *path
        if not self.isOpen:
            path = self.path
            DCDStatus_Check ( DCDRead_Open ( self.cObject, path ) )
            self.isOpen = True

# . Useful if owner data sources change?
#    def ResetOwnerData ( self ):
#        """Reset owner data."""
#        pass

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ) and ( self.cObject != NULL ) and self.isOpen:
            summary = log.GetSummary ( )
            summary.Start ( "DCD Trajectory File Reader" )
            summary.Entry ( "Atoms"        , "{:d}".format ( self.cObject.numberOfAtoms        ) )
            summary.Entry ( "Frames"       , "{:d}".format ( self.cObject.numberOfFrames       ) )
            summary.Entry ( "Atom Indices" , "{:d}".format ( self.cObject.numberOfAtomIndices  ) )
            summary.Entry ( "Has Symmetry" , "{!r}".format ( self.cObject.hasUnitCell == CTrue ) )
            summary.Stop ( )

    def ReadFooter ( self ):
        """Read the trajectory footer."""
        pass

    def ReadHeader ( self, **keywordArguments ):
        """Read the trajectory header."""
        DCDStatus_Check ( DCDRead_Header               ( self.cObject ) )
        DCDStatus_Check ( DCDHandle_CheckNumberOfAtoms ( self.cObject, len ( self.owner.atoms ) ) )
        DCDStatus_Check ( DCDHandle_AllocateQW         ( self.cObject ) )

    def RestoreOwnerData ( self ):
        """Restore data from a frame to the owner."""
        try:
            if self.currentFrame >= self.numberOfFrames: raise EOFError
            DCDStatus_Check ( DCDRead_Frame ( self.cObject ) )
            return True
        except EOFError:
            DCDStatus_Check ( DCDRead_GotoFrame ( self.cObject, 0 ) )
            return False

    property currentFrame:
        def __get__ ( self ): return DCDHandle_CurrentFrame ( self.cObject )

    property numberOfFrames:
        def __get__ ( self ): return DCDHandle_NumberOfFrames ( self.cObject )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def DCDTrajectory_ToSystemGeometryTrajectory ( inPath, outPath, system ):
    """Convert a DCD trajectory to a SystemGeometryTrajectory."""
    # . Define the input and output trajectories.
    inTrajectory  = DCDTrajectoryFileReader  ( inPath , system )
    outTrajectory = SystemGeometryTrajectory ( outPath, system, mode = "w" )
    # . Loop over frames.
    inTrajectory.ReadHeader ( )
    while inTrajectory.RestoreOwnerData ( ): outTrajectory.WriteOwnerData ( )
    inTrajectory.Close  ( )
    outTrajectory.Close ( )
