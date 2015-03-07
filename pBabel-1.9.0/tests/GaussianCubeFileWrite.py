"""Test for writing Gaussian cube files."""

import glob, os, os.path

from pBabel  import GaussianCubeFile_FromSystemDensity, GaussianCubeFile_FromSystemOrbitals, GaussianCubeFile_FromSystemPotential, XYZFile_ToSystem
from pCore   import TestCase
from pMolecule import ElectronicState, QCModelDFT, QCModelMNDO

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The destination for results.
_Destination = "gaussianCubeFiles"

# . The file extension.
_Extension = ".cube"

# . Grid spacing.
_GridSpacing = 0.3

# . Energy models.
_QCModels = [ ( "am1", QCModelMNDO ( keepOrbitalData = True ) ), ( "dft", QCModelDFT ( keepOrbitalData = True ) ) ]

# . The system.
_SystemLabel = "fch3cl"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GaussianCubeFileWriteTest ( TestCase ):
    """A test case for writing Gaussian cube files."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK = True

        # . Output setup.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "xyz" )
        if self.resultPath is None: outPath  = os.path.join ( os.getenv ( "PDYNAMO_SCRATCH" ), _Destination )
        else:                       outPath  = os.path.join ( self.resultPath, _Destination )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )
        log = self.GetLog ( )

        # . Loop over energy models.
        for ( qcLabel, qcModel ) in _QCModels:

            # . Read the system.
            molecule                 = XYZFile_ToSystem ( os.path.join ( dataPath, _SystemLabel + ".xyz" ) )
            molecule.electronicState = ElectronicState  ( charge = -1 )
            molecule.DefineQCModel ( qcModel )
            molecule.Summary ( log = log )
            molecule.Energy  ( log = log )

            # . Orbital energies.
            ( energies, HOMO, LUMO ) = molecule.energyModel.qcModel.OrbitalEnergies ( molecule.configuration )
            if log is not None:
                if energies is not None: energies.Print ( log = log, title = "Orbital Energies" )
                log.Paragraph ( "HOMO and LUMO indices: {:d}, {:d}.".format ( HOMO, LUMO ) )
            indices = [ HOMO, LUMO ]

            # . Write out the cube files.
            GaussianCubeFile_FromSystemDensity   ( os.path.join ( outPath, _SystemLabel + "_" + qcLabel + "_density"  + _Extension ), molecule, gridspacing = _GridSpacing, log = log )
            GaussianCubeFile_FromSystemOrbitals  ( os.path.join ( outPath, _SystemLabel + "_" + qcLabel + "_orbitals" + _Extension ), molecule, gridspacing = _GridSpacing, orbitals = indices, log = log )
        #    GaussianCubeFile_FromSystemPotential ( os.path.join ( outPath, _SystemLabel + "_" + qclabel + "_potential" + _Extension ), molecule, gridspacing = _GridSpacing, log = log )

        # . Success/failure.
        self.assertTrue ( isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = GaussianCubeFileWriteTest ( )
    test.run ( )
