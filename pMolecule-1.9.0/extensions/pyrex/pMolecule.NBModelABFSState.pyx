#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelABFSState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the state for an ABFS NB model."""

from pCore import logFile, LogFileActive

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelABFSState:

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            NBModelABFSState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModelABFSState"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def _Allocate ( self, Integer extent ):
        """Allocation."""
        self.cObject = NBModelABFSState_Allocate ( extent )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def GetEnergies ( self, energies ):
        """Append the energies for the state to energies (a non-NULL list check is used)."""
        if ( self.cObject.nbmmmm     != NULL ):
            energies.append ( ( "MM/MM Elect.",       self.cObject.emmel     ) )
            energies.append ( ( "MM/MM LJ",           self.cObject.emmlj     ) )
        if ( self.cObject.nbmmmm14   != NULL ):
            energies.append ( ( "MM/MM 1-4 Elect.",   self.cObject.emmel14   ) )
            energies.append ( ( "MM/MM 1-4 LJ",       self.cObject.emmlj14   ) )
        if ( self.cObject.nbqcmmlj   != NULL ):
            energies.append ( ( "QC/MM LJ",           self.cObject.eqcmmlj   ) )
        if ( self.cObject.nbqcmmlj14 != NULL ):
            energies.append ( ( "QC/MM 1-4 LJ",       self.cObject.eqcmmlj14 ) )
        if ( self.cObject.inbmmmm    != NULL ):
            energies.append ( ( "MM/MM Image Elect.", self.cObject.eimmmel   ) )
            energies.append ( ( "MM/MM Image LJ",     self.cObject.eimmmlj   ) )
        if ( self.cObject.inbqcmmlj  != NULL ):
            energies.append ( ( "QC/MM Image LJ",     self.cObject.eimqcmmlj ) )
        if ( self.cObject.inbqcqclj  != NULL ):
            energies.append ( ( "QC/QC Image LJ",     self.cObject.eimqcqclj ) )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def StatisticsSummary ( self, log = logFile ):
        """Statistics Summary."""
        cdef Real n
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "ABFS NB Model State Statistics Summary" )
            summary.Entry ( "Number of Calls"   , "{:d}".format ( self.cObject.numberOfCalls   ) )
            summary.Entry ( "Number of Updates" , "{:d}".format ( self.cObject.numberOfUpdates ) )
            if self.cObject.numberOfUpdates > 0:
                n = self.cObject.numberOfUpdates
                if self.cObject.numberOfMMMMPairs         != 0.0 : summary.Entry ( "<MM/MM Pairs>"            , "{:.1f}".format ( self.cObject.numberOfMMMMPairs         / n ) )
                if self.cObject.nbmmmm14                  != NULL: summary.Entry ( "MM/MM 1-4 Pairs"          , "{:d}"  .format ( self.cObject.nbmmmm14.npairs               ) )
                if self.cObject.numberOfQCMMElPairs       != 0.0 : summary.Entry ( "<QC/MM Elect. Pairs>"     , "{:.1f}".format ( self.cObject.numberOfQCMMElPairs       / n ) )
                if self.cObject.nbqcmmel14                != NULL: summary.Entry ( "QC/MM Elect. 1-4 Pairs"   , "{:d}"  .format ( self.cObject.nbqcmmel14.npairs             ) )
                if self.cObject.numberOfQCMMLJPairs       != 0.0 : summary.Entry ( "<QC/MM LJ Pairs>"         , "{:.1f}".format ( self.cObject.numberOfQCMMLJPairs       / n ) )
                if self.cObject.nbqcmmlj14                != NULL: summary.Entry ( "QC/MM LJ 1-4 Pairs"       , "{:d}"  .format ( self.cObject.nbqcmmlj14.npairs             ) )
                if self.cObject.numberOfMMMMImageImages   != 0.0 : summary.Entry ( "<MM/MM Image Images>"     , "{:.1f}".format ( self.cObject.numberOfMMMMImageImages   / n ) )
                if self.cObject.numberOfMMMMImagePairs    != 0.0 : summary.Entry ( "<MM/MM Image Pairs>"      , "{:.1f}".format ( self.cObject.numberOfMMMMImagePairs    / n ) )
                if self.cObject.numberOfQCMMElImageImages != 0.0 : summary.Entry ( "<QC/MM El. Image Images>" , "{:.1f}".format ( self.cObject.numberOfQCMMElImageImages / n ) )
                if self.cObject.numberOfQCMMElImagePairs  != 0.0 : summary.Entry ( "<QC/MM El. Image Pairs>"  , "{:.1f}".format ( self.cObject.numberOfQCMMElImagePairs  / n ) )
                if self.cObject.numberOfQCMMLJImageImages != 0.0 : summary.Entry ( "<QC/MM LJ Image Images>"  , "{:.1f}".format ( self.cObject.numberOfQCMMLJImageImages / n ) )
                if self.cObject.numberOfQCMMLJImagePairs  != 0.0 : summary.Entry ( "<QC/MM LJ Image Pairs>"   , "{:.1f}".format ( self.cObject.numberOfQCMMLJImagePairs  / n ) )
                if self.cObject.numberOfQCQCElImageImages != 0.0 : summary.Entry ( "<QC/QC El. Image Images>" , "{:.1f}".format ( self.cObject.numberOfQCQCElImageImages / n ) )
                if self.cObject.numberOfQCQCElImagePairs  != 0.0 : summary.Entry ( "<QC/QC El. Image Pairs>"  , "{:.1f}".format ( self.cObject.numberOfQCQCElImagePairs  / n ) )
                if self.cObject.numberOfQCQCLJImageImages != 0.0 : summary.Entry ( "<QC/QC LJ Image Images>"  , "{:.1f}".format ( self.cObject.numberOfQCQCLJImageImages / n ) )
                if self.cObject.numberOfQCQCLJImagePairs  != 0.0 : summary.Entry ( "<QC/QC LJ Image Pairs>"   , "{:.1f}".format ( self.cObject.numberOfQCQCLJImagePairs  / n ) )
                summary.Entry ( "Calls per Update" , "{:.1f}".format ( float ( self.cObject.numberOfCalls ) / n, ) )
                summary.Entry ( "Grid Updating"    , "{!r}".format ( self.cObject.useGridSearch == CTrue ) )
            summary.Stop ( )
            NBModelABFSState_StatisticsInitialize ( self.cObject )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "ABFS NB Model State Summary" )
            # . Regular lists.
            if ( self.cObject.nbmmmm     != NULL ): summary.Entry ( "MM/MM Pairs"            , "{:d}".format ( self.cObject.nbmmmm.npairs     ) )
            if ( self.cObject.nbqcmmlj   != NULL ): summary.Entry ( "QC/MM LJ Pairs"         , "{:d}".format ( self.cObject.nbqcmmlj.npairs   ) )
            if ( self.cObject.nbqcmmel   != NULL ): summary.Entry ( "QC/MM Elect. Pairs"     , "{:d}".format ( self.cObject.nbqcmmel.npairs   ) )
            if ( self.cObject.nbmmmm14   != NULL ): summary.Entry ( "MM/MM 1-4 Pairs"        , "{:d}".format ( self.cObject.nbmmmm14.npairs   ) )
            if ( self.cObject.nbqcmmlj14 != NULL ): summary.Entry ( "QC/MM LJ 1-4 Pairs"     , "{:d}".format ( self.cObject.nbqcmmlj14.npairs ) )
            if ( self.cObject.nbqcmmel14 != NULL ): summary.Entry ( "QC/MM Elect. 1-4 Pairs" , "{:d}".format ( self.cObject.nbqcmmel14.npairs ) )
            # . Image lists.
            if ( self.cObject.inbmmmm   != NULL ):
                summary.Entry ( "MM/MM Image Images"     , "{:d}".format ( ImageList_NumberOfImages ( self.cObject.inbmmmm   ) ) )
                summary.Entry ( "MM/MM Image Pairs"      , "{:d}".format ( ImageList_NumberOfPairs  ( self.cObject.inbmmmm   ) ) )
            if ( self.cObject.inbqcmmel != NULL ):
                summary.Entry ( "QC/MM El. Image Images" , "{:d}".format ( ImageList_NumberOfImages ( self.cObject.inbqcmmel ) ) )
                summary.Entry ( "QC/MM El. Image Pairs"  , "{:d}".format ( ImageList_NumberOfPairs  ( self.cObject.inbqcmmel ) ) )
            if ( self.cObject.inbqcmmlj != NULL ):
                summary.Entry ( "QC/MM LJ Image Images"  , "{:d}".format ( ImageList_NumberOfImages ( self.cObject.inbqcmmlj ) ) )
                summary.Entry ( "QC/MM LJ Image Pairs"   , "{:d}".format ( ImageList_NumberOfPairs  ( self.cObject.inbqcmmlj ) ) )
            if ( self.cObject.inbqcqcel != NULL ):
                summary.Entry ( "QC/QC El. Image Images" , "{:d}".format ( ImageList_NumberOfImages ( self.cObject.inbqcqcel ) ) )
                summary.Entry ( "QC/QC El. Image Pairs"  , "{:d}".format ( ImageList_NumberOfPairs  ( self.cObject.inbqcqcel ) ) )
            if ( self.cObject.inbqcqclj != NULL ):
                summary.Entry ( "QC/QC LJ Image Images"  , "{:d}".format ( ImageList_NumberOfImages ( self.cObject.inbqcqclj ) ) )
                summary.Entry ( "QC/QC LJ Image Pairs"   , "{:d}".format ( ImageList_NumberOfPairs  ( self.cObject.inbqcqclj ) ) )
            summary.Stop ( )
