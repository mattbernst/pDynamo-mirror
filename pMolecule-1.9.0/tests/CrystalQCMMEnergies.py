"""Crystal tests - QC/MM."""

import CrystalMMEnergies

from pMolecule import DIISSCFConverger, QCModelMNDO

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CrystalQCMMTest ( CrystalMMEnergies.CrystalMMTest ):
    """A test case for calculating crystal QC/MM energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( CrystalQCMMTest, self ).__init__ ( *args )
        # . Options.
        self.checkEnergies    = False
        self.doLong           = True
        self.doQCMM           = True
        self.doQCQC           = False
        self.geometryOptimize = True
        # . QC model.
        converger             = DIISSCFConverger ( densityTolerance = 1.0e-8, maximumSCFCycles = 250 )
        self.qcModel          = QCModelMNDO   ( converger = converger )

    def MakeShort ( self ):
        self.doLong           = False
        self.geometryOptimize = False
        return True

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = CrystalQCMMTest ( )
    test.run ( )
