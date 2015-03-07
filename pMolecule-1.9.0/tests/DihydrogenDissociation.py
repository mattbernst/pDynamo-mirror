"""Test for calculating dihydrogen dissociation curves."""

import glob, math, os, os.path

from pBabel  import XYZFile_ToSystem
from pCore   import TestCase, TestDataSet, TestReal, UNITS_ENERGY_ELECTRON_VOLTS_TO_KILOJOULES_PER_MOLE
from pMolecule import DIISSCFConverger, QCModelMNDO

# . UHF should give a lower energy than RHF. Problem with initial guess?

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Distance parameters.
_NumberIncrements = 59
_XIncrement       = 0.1
_XStart           = 0.2

# . QC models.
# . MNDO is used because it gives good results.
_RestrictedQCModel   = QCModelMNDO ( "mndo" )
_UnrestrictedQCModel = QCModelMNDO ( "mndo", isSpinRestricted = False )
_S0CIQCModel  = QCModelMNDO ( "mndo", CIMethod = "Full", activeElectrons =  2, activeOrbitals = 2, requiredRoot = 0, rootMultiplicity = 1 )
_S1CIQCModel  = QCModelMNDO ( "mndo", CIMethod = "Full", activeElectrons =  2, activeOrbitals = 2, requiredRoot = 1, rootMultiplicity = 1 )
_T1CIQCModel  = QCModelMNDO ( "mndo", CIMethod = "Full", activeElectrons =  2, activeOrbitals = 2, requiredRoot = 0, rootMultiplicity = 3 )
_QCModels = ( ( "RHF"  , _RestrictedQCModel   ), \
              ( "UHF"  , _UnrestrictedQCModel ), \
              ( "CI S0", _S0CIQCModel         ), \
              ( "CI S1", _S1CIQCModel         ), \
              ( "CI T1", _T1CIQCModel         )  )

# . Reference values.
# . The experimental values for the ground state are ~ 4.52 eV with a minimum of ~ 0.74 A.
_DissociationEnergy = 4.52
_MinimumDistance    = 0.74

# . Tolerances.
_DistanceTolerance = 0.1
_EnergyTolerance   = 0.2

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class DihydrogenDissociationTest ( TestCase ):
    """A test case for calculating dihydrogen dissociation curves."""

    def runTest ( self ):
        """The test."""

        # . Output setup.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "xyz" )
        log = self.GetLog ( )

        # . Initialize the distances.
        distances = []
        for i in range ( _NumberIncrements ):
            distance = _XStart + float ( i ) * _XIncrement
            distances.append ( distance )

        # . Loop over the QC models.
        results = {}
        for ( label, qcModel ) in _QCModels:

            # . Define the system.
            molecule = XYZFile_ToSystem ( os.path.join ( dataPath, "dihydrogen.xyz" ) )
            molecule.DefineQCModel ( qcModel )
            molecule.Summary ( log = None )

            # . Initialize the coordinates.
            molecule.coordinates3.Set ( 0.0 )

            # . Calculate energies at different distances.
            energies = []
            for distance in distances:
                molecule.coordinates3[0,1] = distance
                energies.append ( molecule.Energy ( log = None ) / UNITS_ENERGY_ELECTRON_VOLTS_TO_KILOJOULES_PER_MOLE )
            minimumEnergy = min ( energies )
            for i in range ( len ( energies ) ):
                energies[i] -= minimumEnergy

            # . Save the results.
            results[label] = energies

        # . Print the distances.
        if log is not None:
            table = log.GetTable ( columns = [ 10 ] + len ( results ) * [ 20 ] )
            table.Start ( )
            table.Title ( "Dissociation Curves (Angstroms/eV)" )
            table.Heading ( "Distance" )
            labels = results.keys ( )
            labels.sort ( )
            for label in labels:
                table.Heading ( label )
            for i in range ( len ( distances ) ):
                table.Entry ( "{:4.1f}".format ( distances[i] ) )
                for label in labels:
                    table.Entry ( "{:20.2f}".format ( results[label][i] ) )
            table.Stop ( )

        # . Get the observed data.
        energies           = results["CI S0"]
        dissociationEnergy = energies[-1]
        minimumDifference  = dissociationEnergy
        minimumDistance    = distances[0]
        minimumEnergy      = min ( energies )
        for ( i, ( distance, energy ) ) in enumerate ( zip ( distances, energies ) ):
            energyDifference = math.fabs ( energy - minimumEnergy )
            if ( energyDifference < minimumDifference ):
                minimumDifference = energyDifference
                minimumDistance   = distance
        observed = { "Dissociation Energy" : dissociationEnergy, "Minimum Distance" : minimumDistance }

        # . Generate reference data.
        referenceData = TestDataSet ( "Dihydrogen Dissociation" )
        referenceData.AddDatum ( TestReal ( "Dissociation Energy", _DissociationEnergy, referenceData, absoluteErrorTolerance = _EnergyTolerance  , toleranceFormat = "{:.2f}", valueFormat = "{:.2f}" ) )
        referenceData.AddDatum ( TestReal ( "Minimum Distance"   , _MinimumDistance   , referenceData, absoluteErrorTolerance = _DistanceTolerance, toleranceFormat = "{:.2f}", valueFormat = "{:.2f}" ) )

        # . Check for success/failure.
        results = referenceData.VerifyAgainst ( observed )
        results.Summary ( log = log, fullSummary = self.fullVerificationSummary )
        isOK    = results.WasSuccessful ( )
        self.assertTrue ( isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = DihydrogenDissociationTest ( )
    test.run ( )
