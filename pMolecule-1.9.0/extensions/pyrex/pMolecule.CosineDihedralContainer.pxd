#-------------------------------------------------------------------------------
# . File      : pMolecule.CosineDihedralContainer.pxd
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
cdef extern from "CosineDihedralContainer.h":

    ctypedef struct CCosineDihedral "CosineDihedral":
        Boolean QACTIVE
        Integer atom1
        Integer atom2
        Integer atom3
        Integer atom4
        Integer type

    ctypedef struct CCosineDihedralParameter "CosineDihedralParameter":
        Integer  npowers
        Integer  nterms
        Integer *periods
        Real    *coefficients
        Real    *powercoefficients

    ctypedef struct CCosineDihedralContainer "CosineDihedralContainer":
        Boolean                   QSORTED
        Integer                   nperiods
        Integer                   nparameters
        Integer                   nterms
        CCosineDihedral          *terms
        CCosineDihedralParameter *parameters

    cdef void                      CosineDihedralContainer_ActivateTerms            ( CCosineDihedralContainer  *self )
    cdef CCosineDihedralContainer *CosineDihedralContainer_Allocate                 ( Integer nterms, Integer nparameters )
    cdef CCosineDihedralContainer *CosineDihedralContainer_Clone                    ( CCosineDihedralContainer  *self )
    cdef void                      CosineDihedralContainer_DeactivateFixedAtomTerms ( CCosineDihedralContainer  *self, CSelection *fixedatoms )
    cdef void                      CosineDihedralContainer_DeactivateQCAtomTerms    ( CCosineDihedralContainer  *self, CSelection *qcAtoms, CSelection *boundaryatoms )
    cdef void                      CosineDihedralContainer_Deallocate               ( CCosineDihedralContainer **self )
    cdef Real                      CosineDihedralContainer_Energy                   ( CCosineDihedralContainer  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef void                      CosineDihedralContainer_FindMaximumPeriod        ( CCosineDihedralContainer  *self )
    cdef void                      CosineDihedralContainer_MakePowers               ( CCosineDihedralContainer  *self )
    cdef CCosineDihedralContainer *CosineDihedralContainer_Merge                    ( CCosineDihedralContainer  *self, CCosineDihedralContainer *other, Integer atomincrement )
    cdef Integer                   CosineDihedralContainer_NumberOfInactiveTerms    ( CCosineDihedralContainer  *self )
    cdef CCosineDihedralContainer *CosineDihedralContainer_Prune                    ( CCosineDihedralContainer  *self, CSelection *selection )
    cdef void                      CosineDihedralContainer_Sort                     ( CCosineDihedralContainer  *self )
    cdef Integer                   CosineDihedralContainer_UpperBound               ( CCosineDihedralContainer  *self )

    cdef void CosineDihedralParameter_Allocate ( CCosineDihedralParameter *self, Integer nterms )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CosineDihedralContainer ( MMTerm ):

    cdef CCosineDihedralContainer *cObject
    cdef public object             isOwner
    cdef public object             label
    cdef public object             parameterKeys
