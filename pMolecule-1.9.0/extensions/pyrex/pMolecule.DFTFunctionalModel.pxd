#-------------------------------------------------------------------------------
# . File      : pMolecule.DFTFunctionalModel.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions   cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Integer1DArray cimport CInteger1DArray
from pCore.Status         cimport Status

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DFTFunctionalModel.h":

    ctypedef struct CDFTFunctionalModel "DFTFunctionalModel":
        pass

    cdef CDFTFunctionalModel *DFTFunctionalModel_Allocate    ( Integer numberOfFunctionals, Status *status )
    cdef void                 DFTFunctionalModel_Deallocate  ( CDFTFunctionalModel **self )
    cdef CDFTFunctionalModel *DFTFunctionalModel_MakeFromIDs ( CInteger1DArray *ids, Boolean isSpinRestricted, Status *status )
