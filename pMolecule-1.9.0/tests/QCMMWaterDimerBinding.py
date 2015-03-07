"""Test for calculating MM, QC and QC/MM water dimer binding energies."""

import glob, math, os, os.path

from pBabel           import MOLFile_ToSystem
from pCore            import logFile, LogFileActive, Selection, TestCase, TestDataSet, TestReal
from pMolecule        import MMModelOPLS, NBModelFull, QCModelMNDO
from pMoleculeScripts import PruneByAtom

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Model names.
_Models = ( "MM", "QC" )

# . Bounds on binding energies (kJ mol^-1).
_LowerBound = -35.0
_UpperBound =   5.0

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class WaterDimerBindingTest ( TestCase ):
    """A test case for calculating water dimer binding energies."""

    @staticmethod
    def MonomerEnergies ( dimer, selection, nbModel, qcModel, log = logFile ):
        """Determine the monomer energies."""
        if LogFileActive ( log ):
            log.Heading ( "Monomer Calculation", QBLANKLINE = True )
        e = {}
        # . Get the monomer.
        monomer       = PruneByAtom ( dimer, selection )
        monomer.label = "Water Monomer"
        monomer.Summary ( log = log )
        # . MM model.
        monomer.DefineNBModel ( nbModel )
        e["MM"] = monomer.Energy ( log = log )
        # . QC model.
        monomer.DefineQCModel ( qcModel )
        e["QC"] = monomer.Energy ( log = log )
        return e

    def runTest ( self ):
        """The test."""

        # . Output setup.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "mol" )
        log = self.GetLog ( )

        # . Define the MM, NB and QC models.
        mmModel = MMModelOPLS ( "bookSmallExamples" )
        nbModel = NBModelFull ( )
        qcModel = QCModelMNDO ( )

        # . Define the dimer with an MM model.
        dimer = MOLFile_ToSystem ( os.path.join ( dataPath, "waterDimer_cs.mol" ) )
        dimer.DefineMMModel ( mmModel )
        dimer.Summary ( log = log )

        # . Define the monomer selections.
        selection1 = Selection.FromIterable ( range ( 0, 3 ) )
        selection2 = Selection.FromIterable ( range ( 3, 6 ) )

        # . Get the monomer energies.
        e1 = self.MonomerEnergies ( dimer, selection1, nbModel, qcModel, log = log )
        e2 = self.MonomerEnergies ( dimer, selection2, nbModel, qcModel, log = log )

        # . Get the binding energies.
        e12 = {}
        for model1 in _Models:
            for model2 in _Models:
                key = model1 + " " + model2
                if log is not None:
                    log.Heading ( model1 + "/" + model2 + " Dimer Calculation", QBLANKLINE = True )
                # . Define the energy model.
                if   key == "QC QC": dimer.DefineQCModel ( qcModel )
                elif key == "QC MM": dimer.DefineQCModel ( qcModel, qcSelection = selection1 )
                elif key == "MM QC": dimer.DefineQCModel ( qcModel, qcSelection = selection2 )
                else:                dimer.energyModel.ClearQCModel ( dimer.configuration )
                if "MM" in key: dimer.DefineNBModel ( nbModel )
                dimer.Summary ( log = log )
                # . Store the results.
                e12[key] = dimer.Energy ( log = log ) - e1[model1] - e2[model2]

        # . Output the results.
        if log is not None:
            keys = e12.keys ( )
            keys.sort ( )
            table = log.GetTable ( columns = [ 20, 20, 20 ] )
            table.Start  ( )
            table.Title  ( "Water Dimer Binding Energies" )
            table.Heading ( "Monomer 1" )
            table.Heading ( "Monomer 2" )
            table.Heading ( "Binding Energy" )
            for key in keys:
                ( model1, model2 ) = key.split ( )
                table.Entry ( model1 )
                table.Entry ( model2 )
                table.Entry ( "{:.1f}".format ( e12[key] ) )
            table.Stop ( )

        # . Success/failure.
        isOK = True
        for e in e12.values ( ):
            if ( e < _LowerBound ) or ( e > _UpperBound ):
                isOK = False
                break
        self.assertTrue ( isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = WaterDimerBindingTest ( )
    test.run ( )
