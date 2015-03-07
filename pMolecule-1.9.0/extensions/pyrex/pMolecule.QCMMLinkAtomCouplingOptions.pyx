#-------------------------------------------------------------------------------
# . File      : pMolecule.QCMMLinkAtomCouplingOptions.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""QC/MM link atom coupling options."""

#===================================================================================================================================
# . Dictionaries for interconverting between C and Python representations.
#===================================================================================================================================
QCMMLinkAtomCoupling_ToEnum   = { "MM Coupling" : QCMMLinkAtomCoupling_MM, "RC Coupling" : QCMMLinkAtomCoupling_RC, "RD Coupling" : QCMMLinkAtomCoupling_RD }
QCMMLinkAtomCoupling_ToString = { QCMMLinkAtomCoupling_MM : "MM Coupling", QCMMLinkAtomCoupling_RC : "RC Coupling", QCMMLinkAtomCoupling_RD : "RD Coupling" }
