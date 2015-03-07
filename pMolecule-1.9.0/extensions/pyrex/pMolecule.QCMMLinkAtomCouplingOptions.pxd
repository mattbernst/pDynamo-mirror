#-------------------------------------------------------------------------------
# . File      : pMolecule.QCMMLinkAtomCouplingOptions.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "QCMMLinkAtomCouplingOptions.h":

    ctypedef enum QCMMLinkAtomCoupling:
        QCMMLinkAtomCoupling_MM = 0,
        QCMMLinkAtomCoupling_RC = 1,
        QCMMLinkAtomCoupling_RD = 2
