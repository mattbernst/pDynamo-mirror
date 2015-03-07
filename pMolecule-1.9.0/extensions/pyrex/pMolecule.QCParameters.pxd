#-------------------------------------------------------------------------------
# . File      : pMolecule.QCParameters.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions       cimport Boolean, CFalse, CTrue, Integer, Real
from pMolecule.GaussianBasis  cimport GaussianBasis, GaussianBasis_Deallocate, GaussianBasisType_Coulomb, GaussianBasisType_Orbital, GaussianBasisType_Poisson, GaussianBasis_ScaleExponents, CGaussianBasis
from pMolecule.MNDOParameters cimport MNDOParameters, CMNDOParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "QCParameters.h":

    ctypedef struct CQCCenter "QCCenter":
        Integer                  atomicNumber
        CGaussianBasis  *densitybasis
        CGaussianBasis  *orbitalbasis
        CGaussianBasis  *poissonbasis
        CMNDOParameters *mndoparameters

    ctypedef struct CQCParameter "QCParameter":
        Integer           ncenters
        CQCCenter *centers

    cdef CQCParameter *QCParameters_Allocate   ( Integer ncenters )
    cdef CQCParameter *QCParameters_Clone      ( CQCParameter  *self )
    cdef void          QCParameters_Deallocate ( CQCParameter **self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCParameters:

    cdef CQCParameter  *cObject
    cdef public object  isOwner
