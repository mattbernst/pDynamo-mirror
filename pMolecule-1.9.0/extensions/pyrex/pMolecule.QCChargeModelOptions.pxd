#-------------------------------------------------------------------------------
# . File      : pMolecule.QCChargeModelOptions.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "QCChargeModelOptions.h":

    ctypedef enum QCChargeModel:
        QCChargeModel_CoulombFitting = 0,
        QCChargeModel_Lowdin         = 1,
        QCChargeModel_Mulliken       = 2
