"""Crystal tests - QC."""

import CrystalQCMMEnergies

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CrystalQCTest ( CrystalQCMMEnergies.CrystalQCMMTest ):
    """A test case for calculating crystal QC energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( CrystalQCTest, self ).__init__ ( *args )
        # . Options.
        self.doQCMM = False
        self.doQCQC = True

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = CrystalQCTest ( )
    test.run ( )
