"""DFT RKS tests."""

import MNDORHFEnergies

from pMolecule import QCModelDFT

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class DFTRKSEnergiesTest ( MNDORHFEnergies.MNDORHFEnergiesTest ):
    """A test case for calculating DFT RKS energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( DFTRKSEnergiesTest, self ).__init__ ( *args )
        # . Options.
        self.doPrinting           = True
        self.maximumEnergyAtoms   =    9
        self.maximumEnergyTests   =   20
        self.maximumGradientAtoms =    5
        self.maximumGradientTests =    5
        self.testGradients        = True

    def QCModelClass   ( self ): return QCModelDFT
    def QCModelOptions ( self ):
        self.modelLabels = ( "lda:321g", "blyp:321g", "lda:svp", "blyp:svp" )
        return ( ( [], { "densityBasis" : "demon"  , "functional" : "lda" , "inCore" : True, "orbitalBasis" : "321g" } ) ,
                 ( [], { "densityBasis" : "demon"  , "functional" : "blyp", "inCore" : True, "orbitalBasis" : "321g" } ) ,
                 ( [], { "densityBasis" : "weigend", "functional" : "lda" , "inCore" : True, "orbitalBasis" : "svp"  } ) ,
                 ( [], { "densityBasis" : "weigend", "functional" : "blyp", "inCore" : True, "orbitalBasis" : "svp"  } ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = DFTRKSEnergiesTest ( )
    test.setUp ( )
    test.runTest ( )
