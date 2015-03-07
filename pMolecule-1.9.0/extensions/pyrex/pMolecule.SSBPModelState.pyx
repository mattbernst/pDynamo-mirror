#-------------------------------------------------------------------------------
# . File      : pMolecule.SSBPModelState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the state for a SSBP NB model."""

from pCore import logFile, LogFileActive

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SSBPModelState:

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            SSBPModelState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.SSBPModelState"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = SSBPModelState_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = True

    def GetEnergies ( self, energies ):
        """Append the energies for the state to energies (a simple non-zero check only is used)."""
        if self.cObject != NULL:
            if ( self.cObject.eTotal != 0.0 ): energies.append ( ( "SSBP", self.cObject.eTotal ) )
            if ( self.cObject.atomOutsideCavity == CTrue ): raise ValueError ( "Atom outside cavity in SSBP calculation." )

    def GetResults ( self ):
        """Return all results."""
        results = {}
        if self.cObject != NULL:
            results = { "Atom Outside Cavity"         : ( self.cObject.atomOutsideCavity == CTrue ), \
                        "Cavity Potential Atoms"      : Selection_Size                      ( self.cObject.cavitySelection     ), \
			"MM Atoms"                    : MMAtomContainer_NumberOfActiveAtoms ( self.cObject.mmAtoms             ), \
			"QC Atoms"                    : QCAtomContainer_Size                ( self.cObject.qcAtoms             ), \
                        "Radius Determination Atoms"  : Selection_Size                      ( self.cObject.radiusSelection     ), \
                        "Water Molecules"             : Integer2DArray_Length               ( self.cObject.waterAtomIndices, 0 ), \
                        "Cavity Radius"               : self.cObject.radius              , \
                        "Angular Energy"              : self.cObject.eAngular            , \
                        "Cavity Energy"               : self.cObject.eCavity             , \
                        "Empirical Correction Energy" : self.cObject.eEmpiricalCorrection, \
                        "HardSphere Energy"           : self.cObject.eHardSphere         , \
                        "Kirkwood Energy"             : self.cObject.eKirkwood           , \
                        "Kirkwood Energy (Check)"     : self.cObject.eKirkwoodCheck      , \
                        "Total Charge"                : self.cObject.qTotal                }
        return results

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ): log.Paragraph ( "SSBP NB Model State set up." )
