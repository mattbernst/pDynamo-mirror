#-------------------------------------------------------------------------------
# . File      : pCore.Correlation.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Real1DArray  cimport CReal1DArray, Real1DArray
from pCore.Real2DArray  cimport CReal2DArray, Real2DArray
from pCore.Status       cimport Status, Status_Continue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Correlation.h":

    cdef CReal1DArray *Correlation_MakeDotProduct ( CReal2DArray *x, CReal2DArray *y, CReal1DArray *c, Boolean useFFT, Boolean normalize, Boolean removeMean, Integer tCorrelation, Real *tolerance, Status *status )
    cdef CReal1DArray *Correlation_MakeSimple     ( CReal1DArray *x, CReal1DArray *y, CReal1DArray *c, Boolean useFFT, Boolean normalize, Boolean removeMean, Integer tCorrelation, Real *tolerance, Status *status )
