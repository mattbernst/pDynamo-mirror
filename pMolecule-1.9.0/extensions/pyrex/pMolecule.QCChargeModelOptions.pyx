#-------------------------------------------------------------------------------
# . File      : pMolecule.QCChargeModelOptions.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""QC charge model options."""

#===================================================================================================================================
# . Dictionaries for interconverting between C and Python representations.
#===================================================================================================================================
QCChargeModel_ToEnum   = { "Coulomb Fitting" : QCChargeModel_CoulombFitting , "Lowdin" : QCChargeModel_Lowdin, "Mulliken" : QCChargeModel_Mulliken }
QCChargeModel_ToString = {  QCChargeModel_CoulombFitting : "Coulomb Fitting", QCChargeModel_Lowdin : "Lowdin", QCChargeModel_Mulliken : "Mulliken" }
