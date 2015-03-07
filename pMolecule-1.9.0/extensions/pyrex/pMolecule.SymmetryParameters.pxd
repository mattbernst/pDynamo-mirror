#-------------------------------------------------------------------------------
# . File      : pMolecule.SymmetryParameters.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions       cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3       cimport CCoordinates3, Coordinates3
from pCore.Matrix33           cimport Matrix33, CMatrix33
from pCore.Selection          cimport Selection, CSelection
from pCore.SelectionContainer cimport SelectionContainer, CSelectionContainer
from pCore.Status             cimport Status
from pCore.Vector3            cimport CVector3, Vector3

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "SymmetryParameters.h":

    ctypedef struct CSymmetryParameters "SymmetryParameters":
        Boolean   QM
        Real a
        Real b
        Real c
        Real alpha
        Real beta
        Real gamma
        CMatrix33 *inverseM
        CMatrix33 *M

    cdef CSymmetryParameters *SymmetryParameters_Allocate                          ( )
    cdef Status               SymmetryParameters_CenterCoordinates3ByIndex         ( CSymmetryParameters  *self, CSelection *selection, CCoordinates3 *coordinates3 )
    cdef Status               SymmetryParameters_CenterCoordinates3ByIsolate       ( CSymmetryParameters  *self, CSelectionContainer *isolates, CSelection *selection, CCoordinates3 *coordinates3 )
    cdef void                 SymmetryParameters_Deallocate                        ( CSymmetryParameters **self )
    cdef Boolean              SymmetryParameters_IsMinimumImageConventionSatisfied ( CSymmetryParameters  *self, Real length )
    cdef void                 SymmetryParameters_MakeMinimumImageVector3           ( CSymmetryParameters  *self, CVector3 *r, CVector3 *dr )
    cdef void                 SymmetryParameters_SetCrystalParameters              ( CSymmetryParameters  *self, Real a, Real b, Real c, Real alpha, Real beta, Real gamma )
    cdef Real                 SymmetryParameters_Volume                            ( CSymmetryParameters  *self )

#===============================================================================
# . Class.
#===============================================================================
cdef class SymmetryParameters:

    cdef CSymmetryParameters *cObject
    cdef public object        isOwner
    cdef public object        M
    cdef public object        inverseM
