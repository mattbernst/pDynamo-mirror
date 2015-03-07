#-------------------------------------------------------------------------------
# . File      : pMolecule.SymmetryParameterGradients.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions           cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3           cimport CCoordinates3, Coordinates3
from pCore.Matrix33               cimport Matrix33, CMatrix33
from pMolecule.SymmetryParameters cimport SymmetryParameters, CSymmetryParameters

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "SymmetryParameterGradients.h":

    ctypedef struct CSymmetryParameterGradients "SymmetryParameterGradients":
        Real dEda
        Real dEdb
        Real dEdc
        Real dEdalpha
        Real dEdbeta
        Real dEdgamma
        CMatrix33 *dEdM

    cdef CSymmetryParameterGradients *SymmetryParameterGradients_Allocate              ( )
    cdef void                         SymmetryParameterGradients_Deallocate            ( CSymmetryParameterGradients **self )
    cdef void                         SymmetryParameterGradients_CrystalDerivatives    ( CSymmetryParameterGradients  *self, CSymmetryParameters *symmetryParameters )
    cdef void                         SymmetryParameterGradients_FractionalDerivatives ( CSymmetryParameterGradients  *self, CSymmetryParameters *symmetryParameters, CCoordinates3 *coordinates3, CCoordinates3 *gradients3 )

#===============================================================================
# . Class.
#===============================================================================
cdef class SymmetryParameterGradients:

    cdef CSymmetryParameterGradients *cObject
    cdef public object                isOwner
    cdef public object                dEdM
