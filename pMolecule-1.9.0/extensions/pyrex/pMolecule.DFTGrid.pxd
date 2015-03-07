#-------------------------------------------------------------------------------
# . File      : pMolecule.DFTGrid.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions        cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3        cimport CCoordinates3
from pMolecule.QCAtomContainer cimport CQCAtomContainer

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "DFTGridWeights.h":

    ctypedef struct CDFTGridWeights "DFTGridWeights":
        pass

cdef extern from "List.h":

    ctypedef struct CList "List":
        pass

cdef extern from "DFTGrid.h":

    ctypedef enum DFTGridAccuracy:
        DFTGridAccuracy_VeryLow  = 0 ,
        DFTGridAccuracy_Low      = 1 ,
        DFTGridAccuracy_Medium   = 2 ,
        DFTGridAccuracy_High     = 3 ,
        DFTGridAccuracy_VeryHigh = 4

    ctypedef struct CDFTGrid "DFTGrid":
        DFTGridAccuracy  accuracy
        Integer          numberOfPoints
        Real             bfTolerance
        Real             rhoTolerance
        CDFTGridWeights *weights
        CList           *points

    cdef CDFTGrid *DFTGrid_Construct              ( DFTGridAccuracy accuracy, CQCAtomContainer *qcAtoms, CCoordinates3 *qcCoordinates3 )
    cdef void      DFTGrid_Deallocate             ( CDFTGrid **self )
    cdef Real      DFTGrid_FunctionDataSize       ( CDFTGrid  *self )
    cdef Boolean   DFTGrid_HasFunctionData        ( CDFTGrid  *self )
    cdef Integer   DFTGrid_NumberOfFunctionValues ( CDFTGrid  *self )

#===============================================================================
# . Class.
#===============================================================================
#cdef class DFTGrid:

#    cdef CDFTGrid  *cObject
#    cdef public object  isOwner
