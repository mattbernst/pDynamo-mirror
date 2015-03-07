#-------------------------------------------------------------------------------
# . File      : pMolecule.PairwiseInteraction.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.CubicSpline  cimport CCubicSpline, CubicSpline_Deallocate, CubicSpline_MakeFromReal1DArrays
from pCore.Real1DArray  cimport CReal1DArray, Real1DArray
from pCore.Status       cimport Status

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PairwiseInteraction.h":

    ctypedef struct CPairwiseInteractionABFS "PairwiseInteractionABFS":
        Boolean       useAnalyticForm   
        Integer       splinePointDensity
        Real          dampingCutoff     
        Real          innerCutoff       
        Real          outerCutoff       
        CCubicSpline *electrostaticSpline
        CCubicSpline *lennardJonesASpline
        CCubicSpline *lennardJonesBSpline

    cdef CPairwiseInteractionABFS *PairwiseInteractionABFS_Allocate   ( )
    cdef CPairwiseInteractionABFS *PairwiseInteractionABFS_Clone      ( CPairwiseInteractionABFS  *self, Status *status )
    cdef void                      PairwiseInteractionABFS_Deallocate ( CPairwiseInteractionABFS **self )

    cdef CCubicSpline *PairwiseInteractionABFS_MakeElectrostaticSpline ( CPairwiseInteractionABFS  *self, Boolean useAtomicUnits, Status *status )
    cdef CCubicSpline *PairwiseInteractionABFS_MakeLennardJonesASpline ( CPairwiseInteractionABFS  *self, Status *status )
    cdef CCubicSpline *PairwiseInteractionABFS_MakeLennardJonesBSpline ( CPairwiseInteractionABFS  *self, Status *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionABFS ( object ):

    cdef CPairwiseInteractionABFS *cObject
    cdef public object             electrostaticModel
    cdef public object             isOwner
    cdef public object             width1
    cdef public object             width2
