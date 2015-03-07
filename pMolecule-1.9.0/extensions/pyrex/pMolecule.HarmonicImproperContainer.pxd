#-------------------------------------------------------------------------------
# . File      : pMolecule.HarmonicImproperContainer.pxd
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
cdef extern from "HarmonicImproperContainer.h":

    ctypedef struct CHarmonicImproper "HarmonicImproper":
        Boolean QACTIVE
        Integer  atom1
        Integer  atom2
        Integer  atom3
        Integer  atom4
        Integer  type

    ctypedef struct CHarmonicImproperParameter "HarmonicImproperParameter":
         Real eq
         Real fc
         Real coseq
         Real sineq

    ctypedef struct CHarmonicImproperContainer "HarmonicImproperContainer":
       Integer                         nparameters
       Integer                         nterms
       CHarmonicImproper          *terms
       CHarmonicImproperParameter *parameters

    cdef void                        HarmonicImproperContainer_ActivateTerms            ( CHarmonicImproperContainer  *self )
    cdef CHarmonicImproperContainer *HarmonicImproperContainer_Allocate                 ( Integer nterms, Integer nparameters )
    cdef CHarmonicImproperContainer *HarmonicImproperContainer_Clone                    ( CHarmonicImproperContainer  *self )
    cdef void                        HarmonicImproperContainer_DeactivateFixedAtomTerms ( CHarmonicImproperContainer  *self, CSelection *fixedatoms )
    cdef void                        HarmonicImproperContainer_DeactivateQCAtomTerms    ( CHarmonicImproperContainer  *self, CSelection *qcAtoms, CSelection *boundaryatoms )
    cdef void                        HarmonicImproperContainer_Deallocate               ( CHarmonicImproperContainer **self )
    cdef void                        HarmonicImproperContainer_FillCosSinValues         ( CHarmonicImproperContainer  *self )
    cdef Real                        HarmonicImproperContainer_Energy                   ( CHarmonicImproperContainer  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef CHarmonicImproperContainer *HarmonicImproperContainer_Merge                    ( CHarmonicImproperContainer  *self, CHarmonicImproperContainer *other, Integer atomincrement )
    cdef Integer                     HarmonicImproperContainer_NumberOfInactiveTerms    ( CHarmonicImproperContainer  *self )
    cdef CHarmonicImproperContainer *HarmonicImproperContainer_Prune                    ( CHarmonicImproperContainer  *self, CSelection *selection )
    cdef void                        HarmonicImproperContainer_Sort                     ( CHarmonicImproperContainer  *self )
    cdef Integer                     HarmonicImproperContainer_UpperBound               ( CHarmonicImproperContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class HarmonicImproperContainer ( MMTerm ):

    cdef CHarmonicImproperContainer *cObject
    cdef public object               isOwner
    cdef public object               label
    cdef public object               parameterKeys
