#-------------------------------------------------------------------------------
# . File      : pMolecule.CMAPDihedralContainer.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions  cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.BicubicSpline cimport BicubicSpline_MakeFromReal2DArray, BicubicSplineType_Periodic, CBicubicSpline
from pCore.Coordinates3  cimport CCoordinates3, Coordinates3
from pCore.Real1DArray   cimport CReal1DArray, Real1DArray_Allocate, Real1DArray_GetItem, Real1DArray_SetItem
from pCore.Real2DArray   cimport CReal2DArray, Real2DArray_Allocate, Real2DArray_GetItem, Real2DArray_SetItem
from pCore.Selection     cimport Selection, CSelection
from pMolecule.MMTerm    cimport MMTerm

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CMAPDihedralContainer.h":

    ctypedef struct CCMAPDihedral "CMAPDihedral":
        Boolean isActive
        Integer atom1
        Integer atom2
        Integer atom3
        Integer atom4
        Integer atom5
        Integer atom6
        Integer atom7
        Integer atom8
        Integer type

    ctypedef struct CCMAPDihedralContainer "CMAPDihedralContainer":
        Boolean          isSorted
        Integer          nparameters
        Integer          nterms
        CCMAPDihedral   *terms
        CBicubicSpline **parameters

    cdef void                    CMAPDihedralContainer_ActivateTerms            ( CCMAPDihedralContainer  *self )
    cdef CCMAPDihedralContainer *CMAPDihedralContainer_Allocate                 ( Integer nterms, Integer nparameters )
    cdef CCMAPDihedralContainer *CMAPDihedralContainer_Clone                    ( CCMAPDihedralContainer  *self )
    cdef void                    CMAPDihedralContainer_DeactivateFixedAtomTerms ( CCMAPDihedralContainer  *self, CSelection *fixedatoms )
    cdef void                    CMAPDihedralContainer_DeactivateQCAtomTerms    ( CCMAPDihedralContainer  *self, CSelection *qcAtoms, CSelection *boundaryatoms )
    cdef void                    CMAPDihedralContainer_Deallocate               ( CCMAPDihedralContainer **self )
    cdef Real                    CMAPDihedralContainer_Energy                   ( CCMAPDihedralContainer  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef CCMAPDihedralContainer *CMAPDihedralContainer_Merge                    ( CCMAPDihedralContainer  *self, CCMAPDihedralContainer *other, Integer atomincrement )
    cdef Integer                 CMAPDihedralContainer_NumberOfInactiveTerms    ( CCMAPDihedralContainer  *self )
    cdef CCMAPDihedralContainer *CMAPDihedralContainer_Prune                    ( CCMAPDihedralContainer  *self, CSelection *selection )
    cdef void                    CMAPDihedralContainer_Sort                     ( CCMAPDihedralContainer  *self )
    cdef Integer                 CMAPDihedralContainer_UpperBound               ( CCMAPDihedralContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CMAPDihedralContainer ( MMTerm ):

    cdef CCMAPDihedralContainer *cObject
    cdef public object           isOwner
    cdef public object           label
    cdef public object           parameterKeys
