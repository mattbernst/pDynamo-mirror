#-------------------------------------------------------------------------------
# . File      : pCore.RegularGrid.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Status       cimport Status

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RegularGrid.h":

    # . The regular grid dimension type.
    ctypedef struct CRegularGridDimension "RegularGridDimension":
        Boolean isPeriodic
        Integer bins
        Integer stride
        Real    binSize
        Real    lower
        Real    midPointLower
        Real    period
        Real    upper

    # . The regular grid type.
    ctypedef struct CRegularGrid "RegularGrid":
        Integer               ndimensions
        CRegularGridDimension *dimensions

    # . Functions.
    # . Regular grid dimension.
    cdef void   RegularGridDimension_CopyTo     ( CRegularGridDimension *self, CRegularGridDimension *other )
    cdef void   RegularGridDimension_Initialize ( CRegularGridDimension *self )

    # . Regular grid.
    cdef CRegularGrid *RegularGrid_Allocate           ( Integer ndimensions, Status *status )
    cdef CRegularGrid *RegularGrid_Clone              ( CRegularGrid  *self, Status *status )
    cdef void          RegularGrid_Deallocate         ( CRegularGrid **self )
    cdef Integer       RegularGrid_NumberOfGridPoints ( CRegularGrid  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RegularGrid:

    cdef CRegularGrid  *cObject
    cdef public object  isOwner
