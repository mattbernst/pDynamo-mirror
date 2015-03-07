#-------------------------------------------------------------------------------
# . File      : pMolecule.GaussianBasis.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Real2DArray  cimport CReal2DArray

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "GaussianBasis.h":

    ctypedef enum GaussianBasisType:
        GaussianBasisType_Coulomb = 0,
        GaussianBasisType_Orbital = 1,
        GaussianBasisType_Poisson = 2

    ctypedef struct CPrimitive "Primitive":
        Real *ccbf
        Real *coefficients
        Real *coefficients0
        Real  exponent
        Real  exponent0

    ctypedef struct CShellDefinition "ShellDefinition":
        Integer  angularmomentum_low
        Integer  angularmomentum_high
        Integer  cbfindex
        Integer  nbasis
        Integer  ncbf

    ctypedef struct CShell "Shell":
        Integer           nbasisw
        Integer           nprimitives
        Integer           nstart
        Integer           nstartw
        CReal2DArray     *c2s
        CReal2DArray     *s2c
        CPrimitive       *primitives
        CShellDefinition *type

    ctypedef struct CGaussianBasis "GaussianBasis":
        Boolean            QNORMALIZEDPRIMITIVES
        Boolean            QSPHERICAL
        Boolean            QTOSPHERICAL
        Integer            atomicNumber
        Integer            maximum_angularmomentum
        Integer            nbasis
        Integer            nbasisw
        Integer            nshells
        GaussianBasisType  type
        CReal2DArray      *c2o
        CReal2DArray      *o2c
        CShell            *shells

    cdef CGaussianBasis *GaussianBasis_Allocate              ( Integer nshells )
    cdef void            GaussianBasis_AllocateShell         ( CGaussianBasis  *self, Integer ishell, Integer nprimitives, Integer typeindex )
    cdef CGaussianBasis *GaussianBasis_Clone                 ( CGaussianBasis  *self )
    cdef void            GaussianBasis_Deallocate            ( CGaussianBasis **self )

cdef extern from "GaussianBasisCore.h":

    cdef void            GaussianBasis_Normalize             ( CGaussianBasis  *self )
    cdef void            GaussianBasis_ScaleExponents        ( CGaussianBasis  *self, ... )
    cdef void            GaussianBasis_UnnormalizePrimitives ( CGaussianBasis  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasis:

    cdef CGaussianBasis *cObject
    cdef public object   isOwner
