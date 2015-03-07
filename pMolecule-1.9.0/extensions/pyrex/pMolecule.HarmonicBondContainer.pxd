#-------------------------------------------------------------------------------
# . File      : pMolecule.HarmonicBondContainer.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3 cimport CCoordinates3, Coordinates3
from pCore.Memory       cimport Memory_Deallocate_Integer
from pCore.Selection    cimport Selection, Selection_MakeFlags, CSelection
from pMolecule.MMTerm   cimport MMTerm

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "HarmonicBondContainer.h":

    ctypedef struct CHarmonicBond "HarmonicBond":
        Boolean QACTIVE
        Integer  atom1
        Integer  atom2
        Integer  type

    ctypedef struct CHarmonicBondParameter "HarmonicBondParameter":
        Real eq
        Real fc

    ctypedef struct CHarmonicBondContainer "HarmonicBondContainer":
       Integer                     nparameters
       Integer                     nterms
       CHarmonicBond          *terms
       CHarmonicBondParameter *parameters

    cdef void                    HarmonicBondContainer_ActivateTerms            ( CHarmonicBondContainer  *self )
    cdef CHarmonicBondContainer *HarmonicBondContainer_Allocate                 ( Integer nterms, Integer nparameters )
    cdef CHarmonicBondContainer *HarmonicBondContainer_Clone                    ( CHarmonicBondContainer  *self )
    cdef void                    HarmonicBondContainer_DeactivateFixedAtomTerms ( CHarmonicBondContainer  *self, CSelection *fixedatoms )
    cdef void                    HarmonicBondContainer_DeactivateQCAtomTerms    ( CHarmonicBondContainer  *self, CSelection *qcAtoms, CSelection *boundaryatoms )
    cdef void                    HarmonicBondContainer_Deallocate               ( CHarmonicBondContainer **self )
    cdef Real                    HarmonicBondContainer_Energy                   ( CHarmonicBondContainer  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef Integer                 HarmonicBondContainer_IdentifyBoundaryAtoms    ( CHarmonicBondContainer  *self, CSelection *qcAtoms, Integer **mmboundary, Integer **qcpartners )
    cdef CHarmonicBondContainer *HarmonicBondContainer_Merge                    ( CHarmonicBondContainer  *self, CHarmonicBondContainer *other, Integer atomincrement )
    cdef Integer                 HarmonicBondContainer_NumberOfInactiveTerms    ( CHarmonicBondContainer  *self )
    cdef CHarmonicBondContainer *HarmonicBondContainer_Prune                    ( CHarmonicBondContainer  *self, CSelection *selection )
    cdef void                    HarmonicBondContainer_Sort                     ( CHarmonicBondContainer  *self )
    cdef Integer                 HarmonicBondContainer_UpperBound               ( CHarmonicBondContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class HarmonicBondContainer ( MMTerm ):

    cdef CHarmonicBondContainer *cObject
    cdef public object           isOwner
    cdef public object           is12Interaction
    cdef public object           label
    cdef public object           parameterKeys
