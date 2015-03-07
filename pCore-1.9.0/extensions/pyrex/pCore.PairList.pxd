#-------------------------------------------------------------------------------
# . File      : pCore.PairList.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions         cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3         cimport CCoordinates3, Coordinates3
from pCore.IndexedSelection     cimport CIndexedSelection
from pCore.List                 cimport CList, List_Iterate_Initialize
from pCore.Memory               cimport Memory_Allocate_Array_Integer, Memory_Deallocate_Integer
from pCore.Real1DArray          cimport CReal1DArray, Real1DArray
from pCore.RegularGrid          cimport CRegularGrid
from pCore.RegularGridOccupancy cimport CRegularGridOccupancy
from pCore.Selection            cimport CSelection, Selection
from pCore.SelectionContainer   cimport CSelectionContainer, SelectionContainer
from pCore.Status               cimport Status

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PairList.h":

    ctypedef struct CPairList "PairList":
        Boolean  QSELF
        Integer  nconnectionsi
        Integer  npairs
        Integer  upperBoundi
        Integer  upperBoundj
        Integer *connectionsi
        Integer *connectionsj
        CList   *pairs

    cdef CPairList           *PairList_Allocate                        ( Boolean QSELF )
    cdef void                 PairList_ClearRepresentations            ( CPairList  *self )
    cdef void                 PairList_Deallocate                      ( CPairList **self )
    cdef CPairList           *PairList_FromIntegerPairArray            ( Boolean QSELF, Integer npairs, Integer *indices )
    cdef CIndexedSelection   *PairList_Iterate                         ( CPairList  *self )
    cdef Integer              PairList_Length                          ( CPairList  *self )
    cdef Integer             *PairList_ToIntegerPairArray              ( CPairList  *self )
    cdef CPairList           *SelfPairList_FromSelfPairList            ( CPairList  *self, CSelection *andselection, CSelection *orselection, Boolean QRENUMBER )
    cdef CSelectionContainer *SelfPairList_ToIsolateSelectionContainer ( CPairList  *self, Integer upperBound )
    cdef Status               SelfPairList_MakeConnections             ( CPairList  *self, Integer upperBound )
    cdef Integer              SelfPairList_UpperBound                  ( CPairList  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CrossPairList:

    cdef CPairList     *cObject
    cdef public object  label
    cdef public object  isOwner


#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CrossPairListIterator:

    cdef Integer            position
    cdef CPairList         *cObject
    cdef CIndexedSelection *indexedSelection

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SelfPairList:

    cdef CPairList     *cObject
    cdef public object  label
    cdef public object  isOwner

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SelfPairListIterator:

    cdef Integer            position
    cdef CPairList         *cObject
    cdef CIndexedSelection *indexedSelection
