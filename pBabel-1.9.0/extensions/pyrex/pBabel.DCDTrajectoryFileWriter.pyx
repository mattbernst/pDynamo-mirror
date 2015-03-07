"""Classes and functions for writing DCD trajectory files."""

from pCore                    import logFile, LogFileActive
from SystemGeometryTrajectory import SystemGeometryTrajectory

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DCDTrajectoryFileWriter:
    """DCD trajectory file writer."""

    def __dealloc__ ( self ):
        """Finalization."""
        self.Close ( )

    def __getmodule__ ( self ): return "pBabel.DCDTrajectoryFileWriter"

    def __init__ ( self, path, owner, title = "Created by pDynamo" ):
        """Constructor."""
        self._Initialize     ( )
        self._Allocate       ( )
        self.path  = path
        self.owner = owner
        self.title = title
        self.AssignOwnerData ( )
        self.Open            ( )

    def __len__ ( self ):
        return DCDHandle_NumberOfFrames ( self.cObject )

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
        self.title   = None

    def AssignOwnerData ( self ):
        """Assign owner data to the trajectory."""
        cdef Coordinates3       data3
        cdef Selection          atomIndices
        cdef SymmetryParameters symmetryParameters
        # . Get objects.
        atomIndices = None
        if self.owner.hardConstraints is not None:
            fixedAtoms = getattr ( self.owner.hardConstraints, "fixedAtoms", None )
            if fixedAtoms is not None: atomIndices = fixedAtoms.Complement ( upperBound = len ( self.owner.atoms ) )
        data3              = self.owner.coordinates3
        symmetryParameters = self.owner.symmetryParameters
        # . Assignment.
        if atomIndices        is not None: DCDStatus_Check ( DCDHandle_SetAtomIndices        ( self.cObject, atomIndices.cObject        ) )
        if data3              is not None: DCDStatus_Check ( DCDHandle_SetData3              ( self.cObject, data3.cObject              ) )
        if symmetryParameters is not None: DCDStatus_Check ( DCDHandle_SetSymmetryParameters ( self.cObject, symmetryParameters.cObject ) )

    def Close ( self ):
        """Close the file."""
        if self.isOpen:
            DCDWrite_Close ( &self.cObject )
            self.isOpen = False

    def Open ( self ):
        """Open the file."""
        cdef char *path
        if not self.isOpen:
            path = self.path
            DCDStatus_Check ( DCDWrite_Open ( self.cObject, path ) )
            self.isOpen = True

# . Useful if owner data sources change?
#    def ResetOwnerData ( self ):
#        """Reset owner data."""
#        pass

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ) and ( self.cObject != NULL ) and self.isOpen:
            summary = log.GetSummary ( )
            summary.Start ( "DCD Trajectory File Writer" )
            summary.Entry ( "Atoms"        , "{:d}".format ( self.cObject.numberOfAtoms        ) )
            summary.Entry ( "Frames"       , "{:d}".format ( self.cObject.numberOfFrames       ) )
            summary.Entry ( "Atom Indices" , "{:d}".format ( self.cObject.numberOfAtomIndices  ) )
            summary.Entry ( "Has Symmetry" , "{!r}".format ( self.cObject.hasUnitCell == CTrue ) )
            summary.Stop ( )

    def WriteFooter ( self ):
        """Write the trajectory footer."""
        pass

    def WriteHeader ( self, **keywordArguments ):
        """Write the trajectory header."""
        cdef char *title
        title = self.title
        DCDWrite_Header ( self.cObject, title )

    def WriteOwnerData ( self ):
        """Write data from the owner to a frame."""
        self.AssignOwnerData ( )
        DCDStatus_Check ( DCDWrite_Frame ( self.cObject ) )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def DCDTrajectory_FromSystemGeometryTrajectory ( outPath, inPath, system ):
    """Convert a SystemGeometryTrajectory to a DCD trajectory."""
    # . Define the input and output trajectories.
    inTrajectory  = SystemGeometryTrajectory ( inPath , system, mode = "r" )
    outTrajectory = DCDTrajectoryFileWriter  ( outPath, system )
    # . Loop over frames.
    inTrajectory.ReadHeader   ( )
    outTrajectory.WriteHeader ( )
    while inTrajectory.RestoreOwnerData ( ): outTrajectory.WriteOwnerData  ( )
    inTrajectory.Close  ( )
    outTrajectory.Close ( )
