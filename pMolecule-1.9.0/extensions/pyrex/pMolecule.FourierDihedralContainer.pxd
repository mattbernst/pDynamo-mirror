#-------------------------------------------------------------------------------
# . File      : pMolecule.FourierDihedralContainer.pxd
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
cdef extern from "FourierDihedralContainer.h":

    ctypedef struct CFourierDihedral "FourierDihedral":
        Boolean QACTIVE
        Integer  atom1
        Integer  atom2
        Integer  atom3
        Integer  atom4
        Integer  type

    ctypedef struct CFourierDihedralParameter "FourierDihedralParameter":
         Real fc
         Real phase
         Real cosphase
         Real sinphase
         Integer    period

    ctypedef struct CFourierDihedralContainer "FourierDihedralContainer":
       Integer                        nparameters
       Integer                        nterms
       CFourierDihedral          *terms
       CFourierDihedralParameter *parameters

    cdef void                       FourierDihedralContainer_ActivateTerms            ( CFourierDihedralContainer  *self )
    cdef CFourierDihedralContainer *FourierDihedralContainer_Allocate                 ( Integer nterms, Integer nparameters )
    cdef CFourierDihedralContainer *FourierDihedralContainer_Clone                    ( CFourierDihedralContainer  *self )
    cdef void                       FourierDihedralContainer_DeactivateFixedAtomTerms ( CFourierDihedralContainer  *self, CSelection *fixedatoms )
    cdef void                       FourierDihedralContainer_DeactivateQCAtomTerms    ( CFourierDihedralContainer  *self, CSelection *qcAtoms, CSelection *boundaryatoms )
    cdef void                       FourierDihedralContainer_Deallocate               ( CFourierDihedralContainer **self )
    cdef void                       FourierDihedralContainer_FillCosSinPhases         ( CFourierDihedralContainer  *self )
    cdef Real                       FourierDihedralContainer_Energy                   ( CFourierDihedralContainer  *self, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )
    cdef CFourierDihedralContainer *FourierDihedralContainer_Merge                    ( CFourierDihedralContainer  *self, CFourierDihedralContainer *other, Integer atomincrement )
    cdef Integer                    FourierDihedralContainer_NumberOfInactiveTerms    ( CFourierDihedralContainer  *self )
    cdef CFourierDihedralContainer *FourierDihedralContainer_Prune                    ( CFourierDihedralContainer  *self, CSelection *selection )
    cdef void                       FourierDihedralContainer_Sort                     ( CFourierDihedralContainer  *self )
    cdef Integer                    FourierDihedralContainer_UpperBound               ( CFourierDihedralContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class FourierDihedralContainer ( MMTerm ):

    cdef CFourierDihedralContainer *cObject
    cdef public object              isOwner
    cdef public object              label
    cdef public object              parameterKeys
