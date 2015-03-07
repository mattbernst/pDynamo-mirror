#-------------------------------------------------------------------------------
# . File      : pCore.Selection.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "Selection.h":

    ctypedef struct CSelection "Selection":
        Boolean  QSORTED
        Integer  nflags
        Integer  nindices
        Integer  npositions
        Boolean *flags
        Integer *indices
        Integer *positions

    cdef CSelection *Selection_Allocate             ( Integer nindices )
    cdef void        Selection_ClearFlags           ( CSelection  *self )
    cdef void        Selection_ClearPositions       ( CSelection  *self )
    cdef void        Selection_ClearRepresentations ( CSelection  *self )
    cdef CSelection *Selection_Clone                ( CSelection  *self )
    cdef Integer     Selection_Compare              ( CSelection  *self, CSelection *other )
    cdef CSelection *Selection_Complement           ( CSelection  *self, Integer upperBound )
    cdef void        Selection_Deallocate           ( CSelection **self )
    cdef CSelection *Selection_FromFlags            ( Integer nflags, Boolean *flags )
    cdef CSelection *Selection_FromIntegerArray     ( Integer nindices, Integer *indices )
    cdef CSelection *Selection_Intersection         ( Integer nselections, ... )
    cdef void        Selection_MakeFlags            ( CSelection  *self, Integer upperBound )
    cdef void        Selection_MakePositions        ( CSelection  *self, Integer upperBound )
    cdef CSelection *Selection_Merge                ( CSelection  *self, CSelection *other, Integer indexincrement )
    cdef CSelection *Selection_Prune                ( CSelection  *self, CSelection *selection )
    cdef Integer     Selection_Size                 ( CSelection  *self )
    cdef void        Selection_Sort                 ( CSelection  *self )
    cdef CSelection *Selection_Union                ( Integer nselections, ... )
    cdef Integer     Selection_UpperBound           ( CSelection  *self )

#===============================================================================
# . Class.
#===============================================================================
cdef class Selection:

    cdef CSelection   *cObject
    cdef public object isOwner
