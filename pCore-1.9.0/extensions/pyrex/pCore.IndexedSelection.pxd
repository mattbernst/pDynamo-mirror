#-------------------------------------------------------------------------------
# . File      : pCore.IndexedSelection.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Integer

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "IndexedSelection.h":

    ctypedef struct CIndexedSelection "IndexedSelection":
        Integer  index
        Integer  nindices
        Integer *indices
