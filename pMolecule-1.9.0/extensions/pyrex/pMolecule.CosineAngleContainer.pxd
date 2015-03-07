#-------------------------------------------------------------------------------
# . File      : pMolecule.CosineAngleContainer.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3 cimport CCoordinates3, Coordinates3
from pCore.Selection    cimport Selection, CSelection
from pMolecule.MMTerm   cimport MMTerm

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CosineAngleContainer.h":

    ctypedef struct CCosineAngle "CosineAngle":
        Boolean QACTIVE
        Integer  atom1
        Integer  atom2
        Integer  atom3
        Integer  type

    ctypedef struct CCosineAngleParameter "CosineAngleParameter":
        Integer     npowers
        Integer     nterms
        Integer    *periods
        Real *coefficients
        Real *powercoefficients

    ctypedef struct CCosineAngleContainer "CosineAngleContainer":
        Boolean                   QSORTED
        Integer                    nperiods
        Integer                    nparameters
        Integer                    nterms
        CCosineAngle          *terms
        CCosineAngleParameter *parameters

    cdef void                   CosineAngleContainer_ActivateTerms            ( CCosineAngleContainer  *self )
    cdef CCosineAngleContainer *CosineAngleContainer_Allocate                 ( Integer nterms, Integer nparameters )
    cdef CCosineAngleContainer *CosineAngleContainer_Clone                    ( CCosineAngleContainer  *self )
    cdef void                   CosineAngleContainer_DeactivateFixedAtomTerms ( CCosineAngleContainer  *self, CSelection *fixedatoms )
    cdef void                   CosineAngleContainer_DeactivateQCAtomTerms    ( CCosineAngleContainer  *self, CSelection *qcAtoms, CSelection *boundaryatoms )
    cdef void                   CosineAngleContainer_Deallocate               ( CCosineAngleContainer **self )
    cdef Real                   CosineAngleContainer_Energy                   ( CCosineAngleContainer  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef void                   CosineAngleContainer_FindMaximumPeriod        ( CCosineAngleContainer  *self )
    cdef void                   CosineAngleContainer_MakePowers               ( CCosineAngleContainer  *self )
    cdef CCosineAngleContainer *CosineAngleContainer_Merge                    ( CCosineAngleContainer  *self, CCosineAngleContainer *other, Integer atomincrement )
    cdef Integer                CosineAngleContainer_NumberOfInactiveTerms    ( CCosineAngleContainer  *self )
    cdef CCosineAngleContainer *CosineAngleContainer_Prune                    ( CCosineAngleContainer  *self, CSelection *selection )
    cdef void                   CosineAngleContainer_Sort                     ( CCosineAngleContainer  *self )
    cdef Integer                CosineAngleContainer_UpperBound               ( CCosineAngleContainer  *self )

    cdef void CosineAngleParameter_Allocate ( CCosineAngleParameter *self, Integer nterms )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CosineAngleContainer ( MMTerm ):

    cdef CCosineAngleContainer *cObject
    cdef public object          isOwner
    cdef public object          label
    cdef public object          parameterKeys
