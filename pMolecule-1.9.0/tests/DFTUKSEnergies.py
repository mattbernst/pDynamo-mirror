"""DFT UKS tests."""

import MNDOUHFEnergies

from pMolecule import QCModelDFT

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Tolerances.
#_EnergyTolerance                = 0.1
#_GradientAbsoluteErrorTolerance = 1.0e-02

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class DFTUKSEnergiesTest ( MNDOUHFEnergies.MNDOUHFEnergiesTest ):
    """A test case for calculating DFT UKS energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( DFTUKSEnergiesTest, self ).__init__ ( *args )
        # . Options.
        self.maximumEnergyAtoms   = 10
        self.maximumGradientAtoms =  5

    def ConvergerKeywords ( self ): return { "densityTolerance" : 1.0e-8, "maximumSCFCycles" : 500 }
    def QCModelArguments  ( self ): return []
    def QCModelClass      ( self ): return QCModelDFT
    def QCModelKeywords   ( self ): return { "densityBasis" : "demon", "functional" : "lda", "orbitalBasis" : "321g" }

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = DFTUKSEnergiesTest ( )
    test.setUp ( )
    test.runTest ( )
