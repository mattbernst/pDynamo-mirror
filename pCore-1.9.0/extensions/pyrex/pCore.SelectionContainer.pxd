#-------------------------------------------------------------------------------
# . File      : pCore.SelectionContainer.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Selection    cimport CSelection, Selection, Selection_Allocate, Selection_Sort
from pCore.Status       cimport Status     

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "SelectionContainer.h":

    ctypedef struct CSelectionContainer "SelectionContainer":
        Boolean          isOwner
        Integer          nitems
        Integer          upperBound
        CSelection **items

    cdef CSelectionContainer *SelectionContainer_Allocate      ( Integer nitems )
    cdef CSelectionContainer *SelectionContainer_Clone         ( CSelectionContainer  *self )
    cdef void                 SelectionContainer_Deallocate    ( CSelectionContainer **self )
    cdef CSelection          *SelectionContainer_GetAllIndices ( CSelectionContainer  *self, CSelection *indices )
    cdef Status               SelectionContainer_MergeIsolates ( CSelectionContainer  *self, CSelection *tomerge )

#===============================================================================
# . Class.
#===============================================================================
cdef class SelectionContainer:

    cdef CSelectionContainer *cObject
    cdef public object        itemName
    cdef public object        isOwner
