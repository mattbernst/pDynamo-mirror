"""Test for reading Gromacs files."""

import os

from pBabel    import PDBFile_FromSystem, GromacsParameters_ToParameters, GromacsDefinitions_ToSystem, GromacsCrdFile_Process, GromacsCrdFile_ToCoordinates3
from pCore     import logFile, LogFileActive, Pickle, TestCase, TestDataSet, TestReal, TextLogFileWriter, UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE, Unpickle
from pMolecule import NBModelFull, SystemGeometryObjectiveFunction, NBModelABFS

# . Energy values should be the same as those obtained with Gromacs 4.5 to within 1.0e-3 kJ/mole for covalent terms and within 0.1 kJ/mole for LJ/elect. terms
# . They differ from Charmm and Amber program values because of parameter manipulation and (back)conversions.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Molecule names.
_SystemLabels = [ "ava", "1atp_peptide" ] 

# . Force field names.
_ForceFields  = [ "CHARMM", "AMBER" ]

# . Options.
_MaximumAtoms = 100

# . NB Model options. 
ABFS_options = { "innerCutoff" :  8.0 ,
                 "outerCutoff" : 12.0 ,
                 "listCutoff"  : 14.0 }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GromacsTopCrdReadTest ( TestCase ):
    """A test case for reading Gromacs topology and coordinate files."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        observed = {}
        gromacsReference = True

	# . We have two references: One for energies obtained with Gromacs; another for pDynamo energies
        if gromacsReference:
            referenceDataPath = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", "reference/GromacsTopCrdRead_gromacsValues.pkl" )
        else:
            referenceDataPath = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", "reference/GromacsTopCrdRead_pDynamoValues.pkl" )

        # . Output setup.
        dataPath        = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", "gromacs" )
        outPath         = None
	self.resultPath = dataPath
        if self.resultPath is not None:
            outPath = os.path.join ( self.resultPath, "pdb" )
            if not os.path.exists ( outPath ): os.mkdir ( outPath )
        log = self.GetLog ( )

        # . Generate systems.
        for label in _SystemLabels:

	    # . Force field.
	    for ff in _ForceFields:

                # . Header.
                if log is not None:
                    log.Text ( "\n" + ( 120 * "=" ) + "\n" )
                    log.Text ( label + "-" + ff + "\n" )
                    log.Text ( 120 * "=" + "\n" )

                # . Get the parameters.
	        filename = os.path.join ( dataPath, label + "_" + ff )
                parameters          = GromacsParameters_ToParameters ( filename + ".top", log = log )
                system              = GromacsDefinitions_ToSystem ( filename + ".top", parameters = parameters, log = log )
                system.coordinates3 = GromacsCrdFile_ToCoordinates3 ( filename + ".gro", log = log )
                system.label        = label
	        if hasattr ( system.configuration, "symmetryParameters" ): system.DefineNBModel ( NBModelABFS ( **ABFS_options ) )
                else                                                     : system.DefineNBModel ( NBModelFull ( ) )
                system.Summary ( log = log )
                energy = system.Energy  ( log = log )
                log.Text ( "\nEnergy (kcal/mole) = {:.4f}\n".format ( energy / UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE ) )

                # . Get the dictionary of energies.
                localObserved = system.configuration.energyTerms.Terms ( asDictionary = True )
                localObserved["Potential Energy"] = energy

		# . U-B term in Gromacs includes harmonic angle contribution, but not in pDynamo. Their sum should be equal, though.
		if gromacsReference and ( ff == "CHARMM" ):
		    localObserved["Harmonic Angle + U-B"] = localObserved["Harmonic Angle"] + localObserved["Urey-Bradley"]
		    del localObserved["Harmonic Angle"]
		    del localObserved["Urey-Bradley"]

                # . Test gradients.
                if len ( system.atoms ) <= _MaximumAtoms:
                    of = SystemGeometryObjectiveFunction.FromSystem ( system )
                    of.TestGradients ( log = log )
                    localObserved["Gradient Error"] = of.TestGradients ( log = log )

                # . Write PDB file 
                if outPath is not None:
                    PDBFile_FromSystem ( os.path.join ( outPath, label + ".pdb" ), system )

                # . Accumulate current data
		dataLabel = label + "-" + ff
                observed[dataLabel] = localObserved

        # . Verify the observed data against the reference data.
        referenceData = Unpickle ( referenceDataPath )
        results       = referenceData.VerifyAgainst ( observed )
        isOK          = results.WasSuccessful ( )
        results.Summary ( log = log, fullSummary = True )

        # . Success/failure.
        self.assertTrue ( isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = GromacsTopCrdReadTest ( )
    test.run ( )
