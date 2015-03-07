#-------------------------------------------------------------------------------
# . File      : pCore.PairListGenerator.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Pairlist generation class and functions."""

from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairListGenerator ( object ):
    """A pairlist generator."""

    def __copy__ ( self ):
        """Copying."""
        options = self.__getstate__ ( )
        new     = self.__class__ ( **options )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner: PairListGenerator_Deallocate ( &self.cObject )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.PairListGenerator"

    def __getstate__ ( self ):
        """Return the state."""
        return { "sortIndices"          : self.cObject.sortIndices   == CTrue ,
                 "useGridByCell"        : self.cObject.useGridByCell == CTrue ,
                 "minimumCellExtent"    : self.cObject.minimumCellExtent      ,
                 "minimumPoints"        : self.cObject.minimumPoints          ,
                 "cellSize"             : self.cObject.cellSize               ,
                 "cutoff"               : self.cObject.cutoff                 ,
                 "cutoffCellSizeFactor" : self.cObject.cutoffCellSizeFactor   ,
                 "minimumCellSize"      : self.cObject.minimumCellSize        ,
                 "minimumExtentFactor"  : self.cObject.minimumExtentFactor    }

    def __init__ ( self, **options ):
        """Constructor with options."""
        self._Initialize ( )
        self._Allocate   ( )
        self.SetOptions ( **options )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        self.SetOptions ( **state )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = PairListGenerator_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    @classmethod
    def FromOptions ( selfClass, **options ):
        """Constructor from options."""
        return selfClass ( **options )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        if "sortIndices"          in keywordArguments:
            value = keywordArguments.pop ( "sortIndices"   )
            if value: self.cObject.sortIndices   = CTrue
            else:     self.cObject.sortIndices   = CFalse
        if "useGridByCell"        in keywordArguments:
            value = keywordArguments.pop ( "useGridByCell" )
            if value: self.cObject.useGridByCell = CTrue
            else:     self.cObject.useGridByCell = CFalse
        if "minimumCellExtent"    in keywordArguments: self.cObject.minimumCellExtent    = keywordArguments.pop ( "minimumCellExtent"    )
        if "minimumPoints"        in keywordArguments: self.cObject.minimumPoints        = keywordArguments.pop ( "minimumPoints"        )
        if "cellSize"             in keywordArguments: self.cObject.cellSize             = keywordArguments.pop ( "cellSize"             )
        if "cutoff"               in keywordArguments: self.cObject.cutoff               = keywordArguments.pop ( "cutoff"               )
        if "cutoffCellSizeFactor" in keywordArguments: self.cObject.cutoffCellSizeFactor = keywordArguments.pop ( "cutoffCellSizeFactor" )
        if "minimumCellSize"      in keywordArguments: self.cObject.minimumCellSize      = keywordArguments.pop ( "minimumCellSize"      )
        if "minimumExtentFactor"  in keywordArguments: self.cObject.minimumExtentFactor  = keywordArguments.pop ( "minimumExtentFactor"  )
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )
        # . Ensure the cell size is correct.
        self.cObject.cellSize = self.cObject.cutoffCellSizeFactor * self.cObject.cutoff

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "Pairlist Generator Summary" )
            summary.Entry ( "Cell Size"               , "{:.3f}".format ( self.cObject.cellSize               ) )
            summary.Entry ( "Cutoff"                  , "{:.3f}".format ( self.cObject.cutoff                 ) )
            summary.Entry ( "Cutoff Cell Size Factor" , "{:.3f}".format ( self.cObject.cutoffCellSizeFactor   ) )
            summary.Entry ( "Grid Cell/Cell Method"   , "{!r}"  .format ( self.cObject.useGridByCell == CTrue ) )
            summary.Entry ( "Minimum Cell Extent"     , "{:d}"  .format ( self.cObject.minimumCellExtent      ) )
            summary.Entry ( "Minimum Cell Size"       , "{:.3f}".format ( self.cObject.minimumCellSize        ) )
            summary.Entry ( "Minimum Extent Factor"   , "{:.3f}".format ( self.cObject.minimumExtentFactor    ) )
            summary.Entry ( "Minimum Points"          , "{:d}"  .format ( self.cObject.minimumPoints          ) )
            summary.Entry ( "Sort Indices"            , "{!r}"  .format ( self.cObject.sortIndices == CTrue   ) )
            summary.Stop ( )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def CrossPairList_FromDoubleCoordinates3 ( Coordinates3  coordinates31     ,
                                           Coordinates3  coordinates32     ,
                                           CrossPairList exclusions = None ,
                                           Real1DArray   radii1     = None , 
                                           Real1DArray   radii2     = None ,
                                           Real          safety     = 0.0 ):
    """Create a cross-pairlist from two sets of coordinates."""
    cdef CrossPairList      pairList
    cdef CPairList         *cExclusions
    cdef CReal1DArray      *cRadii1, *cRadii2
    cdef PairListGenerator  self
    pairList = CrossPairList.Raw ( )
    if exclusions is None: cExclusions = NULL
    else:                  cExclusions = exclusions.cObject
    if radii1     is None: cRadii1     = NULL
    else:                  cRadii1     = radii1.cObject
    if radii2     is None: cRadii2     = NULL
    else:                  cRadii2     = radii2.cObject
    self = PairListGenerator.FromOptions ( cutoff = safety )
    pairList.cObject = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( self.cObject           ,
                                                                               coordinates31.cObject  ,
                                                                               coordinates32.cObject  ,
                                                                               cRadii1                ,
                                                                               cRadii2                ,
                                                                               NULL, NULL, NULL, NULL ,
                                                                               cExclusions            ,
                                                                               NULL, NULL, NULL       )
    pairList.isOwner = True
    return pairList

def CrossPairList_FromSingleCoordinates3 ( Coordinates3  coordinates3      ,
                                           CrossPairList exclusions = None ,
                                           Real1DArray   radii      = None ,
                                           Selection     selection1 = None ,
                                           Selection     selection2 = None ,
                                           Real          safety     = 0.0 ):
    """Create a cross-pairlist from a single set of coordinates."""
    cdef CrossPairList  pairList
    cdef CPairList         *cExclusions
    cdef CSelection        *cSelection1, *cSelection2
    cdef CReal1DArray      *cRadii
    cdef PairListGenerator  self
    pairList = CrossPairList.Raw ( )
    if exclusions is None: cExclusions = NULL
    else:                  cExclusions = exclusions.cObject
    if radii      is None: cRadii      = NULL
    else:                  cRadii      = radii.cObject
    if selection1 is None: cSelection1 = NULL
    else:                  cSelection1 = selection1.cObject
    if selection2 is None: cSelection2 = NULL
    else:                  cSelection2 = selection2.cObject
    self = PairListGenerator.FromOptions ( cutoff = safety )
    pairList.cObject = PairListGenerator_CrossPairListFromSingleCoordinates3 ( self.cObject            ,
                                                                               coordinates3.cObject    ,
                                                                               cRadii                  ,
                                                                               cSelection1             ,
                                                                               cSelection2             ,
                                                                               NULL                    ,
                                                                               cExclusions             ,
                                                                               CTrue, NULL, NULL, NULL )
    pairList.isOwner = True
    return pairList

def SelfPairList_FromCoordinates3 ( Coordinates3 coordinates3      ,
                                    SelfPairList exclusions = None ,
                                    Real1DArray  radii      = None ,
                                    Real         safety     = 0.0 ):
    """Create a self-pairlist from a set of coordinates."""
    cdef CPairList         *cExclusions
    cdef CReal1DArray      *cRadii
    cdef PairListGenerator  self
    cdef SelfPairList       pairList
    pairList = SelfPairList.Raw ( )
    if exclusions is None: cExclusions = NULL
    else:                  cExclusions = exclusions.cObject
    if radii      is None: cRadii      = NULL
    else:                  cRadii      = radii.cObject
    self = PairListGenerator.FromOptions ( cutoff = safety )
    pairList.cObject = PairListGenerator_SelfPairListFromCoordinates3 ( self.cObject            ,
                                                                        coordinates3.cObject    , 
                                                                        cRadii                  ,
                                                                        NULL, NULL              ,
                                                                        cExclusions             ,
                                                                        NULL, NULL, NULL        )
    pairList.isOwner = True
    return pairList
