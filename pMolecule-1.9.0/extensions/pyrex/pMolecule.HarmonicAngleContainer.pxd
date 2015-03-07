#-------------------------------------------------------------------------------
# . File      : pMolecule.HarmonicAngleContainer.pxd
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
cdef extern from "HarmonicAngleContainer.h":

    ctypedef struct CHarmonicAngle "HarmonicAngle":
        Boolean QACTIVE
        Integer  atom1
        Integer  atom2
        Integer  atom3
        Integer  type

    ctypedef struct CHarmonicAngleParameter "HarmonicAngleParameter":
        Real eq
        Real fc

    ctypedef struct CHarmonicAngleContainer "HarmonicAngleContainer":
       Integer                      nparameters
       Integer                      nterms
       CHarmonicAngle          *terms
       CHarmonicAngleParameter *parameters

    cdef void                     HarmonicAngleContainer_ActivateTerms            ( CHarmonicAngleContainer  *self )
    cdef CHarmonicAngleContainer *HarmonicAngleContainer_Allocate                 ( Integer nterms, Integer nparameters )
    cdef CHarmonicAngleContainer *HarmonicAngleContainer_Clone                    ( CHarmonicAngleContainer  *self )
    cdef void                     HarmonicAngleContainer_DeactivateFixedAtomTerms ( CHarmonicAngleContainer  *self, CSelection *fixedatoms )
    cdef void                     HarmonicAngleContainer_DeactivateQCAtomTerms    ( CHarmonicAngleContainer  *self, CSelection *qcAtoms, CSelection *boundaryatoms )
    cdef void                     HarmonicAngleContainer_Deallocate               ( CHarmonicAngleContainer **self )
    cdef Real                     HarmonicAngleContainer_Energy                   ( CHarmonicAngleContainer  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef CHarmonicAngleContainer *HarmonicAngleContainer_Merge                    ( CHarmonicAngleContainer  *self, CHarmonicAngleContainer *other, Integer atomincrement )
    cdef Integer                  HarmonicAngleContainer_NumberOfInactiveTerms    ( CHarmonicAngleContainer  *self )
    cdef CHarmonicAngleContainer *HarmonicAngleContainer_Prune                    ( CHarmonicAngleContainer  *self, CSelection *selection )
    cdef void                     HarmonicAngleContainer_Sort                     ( CHarmonicAngleContainer  *self )
    cdef Integer                  HarmonicAngleContainer_UpperBound               ( CHarmonicAngleContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class HarmonicAngleContainer ( MMTerm ):

    cdef CHarmonicAngleContainer *cObject
    cdef public object            isOwner
    cdef public object            label
    cdef public object            parameterKeys
