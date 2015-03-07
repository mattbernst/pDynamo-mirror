#-------------------------------------------------------------------------------
# . File      : pCore.PairListGenerator.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions         cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3         cimport CCoordinates3, Coordinates3
from pCore.PairList             cimport CPairList, CrossPairList, SelfPairList
from pCore.Real1DArray          cimport CReal1DArray, Real1DArray
from pCore.RegularGrid          cimport CRegularGrid
from pCore.RegularGridOccupancy cimport CRegularGridOccupancy
from pCore.Selection            cimport CSelection, Selection
from pCore.Status               cimport Status

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PairListGenerator.h":

    ctypedef struct CPairListGenerator "PairListGenerator":
        Boolean sortIndices
        Boolean useGridByCell
        Integer minimumCellExtent
        Integer minimumPoints
        Real    cellSize
        Real    cutoff
        Real    cutoffCellSizeFactor
        Real    minimumCellSize
        Real    minimumExtentFactor

    cdef CPairListGenerator *PairListGenerator_Allocate        ( )
    cdef CPairListGenerator *PairListGenerator_Clone           ( CPairListGenerator  *self )
    cdef void                PairListGenerator_Deallocate      ( CPairListGenerator **self )
    cdef Boolean             PairListGenerator_DetermineMethod ( CPairListGenerator  *self, CCoordinates3 *coordinates3, CSelection *andSelection )

    cdef CPairList *PairListGenerator_CrossPairListFromDoubleCoordinates3 ( CPairListGenerator    *self          ,
                                                                            CCoordinates3         *coordinates31 ,
                                                                            CCoordinates3         *coordinates32 ,
                                                                            CReal1DArray          *radii1        ,
                                                                            CReal1DArray          *radii2        ,
                                                                            CSelection            *andSelection1 ,
                                                                            CSelection            *andSelection2 ,
                                                                            CSelection            *orSelection1  ,
                                                                            CSelection            *orSelection2  ,
                                                                            CPairList             *exclusions    ,
                                                                            CRegularGrid          *grid1         ,
                                                                            CRegularGridOccupancy *occupancy1    ,
                                                                            Status                *status        )
    cdef CPairList *PairListGenerator_CrossPairListFromSingleCoordinates3 ( CPairListGenerator    *self          ,
                                                                            CCoordinates3         *coordinates3  ,
                                                                            CReal1DArray          *radii         ,
                                                                            CSelection            *andSelection1 ,
                                                                            CSelection            *andSelection2 ,
                                                                            CSelection            *orSelection   ,
                                                                            CPairList             *exclusions    ,
                                                                            Boolean                excludeSelf   ,
                                                                            CRegularGrid          *grid1         ,
                                                                            CRegularGridOccupancy *occupancy1    ,
                                                                            Status                *status        )
    cdef CPairList *PairListGenerator_SelfPairListFromCoordinates3        ( CPairListGenerator    *self          ,
                                                                            CCoordinates3         *coordinates3  ,
                                                                            CReal1DArray          *radii         ,
                                                                            CSelection            *andSelection  ,
                                                                            CSelection            *orSelection   ,
                                                                            CPairList             *exclusions    ,
                                                                            CRegularGrid          *grid          ,
                                                                            CRegularGridOccupancy *occupancy     ,
                                                                            Status                *status        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairListGenerator ( object ):

    cdef CPairListGenerator *cObject
    cdef public object       isOwner
