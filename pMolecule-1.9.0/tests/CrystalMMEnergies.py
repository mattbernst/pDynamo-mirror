"""Crystal tests - MM.

The systems are small molecular crystals.

The AMBER files were set up using fDynamo and so the parameters may not be that good.
However, they are fine for testing.

All pure QC optimizations fail except three. This is due to protons in the primary image migrating to the secondary images.

Care must be taken when choosing a QCregion that the total MM charge is zero (or integer).
Likewise some selections (in this case only for GLYGLY) lead to exploded molecules.

For MNDO - results OK with existing parameters.
For DFT  - gradients worse. Require 1.0e-4 step for DFT Cartesian part but this gives bad crystal derivatives.
"""

import math, os.path

from pBabel           import AmberTopologyFile_ToSystem, AmberCrdFile_ToCoordinates3
from pCore            import LBFGSMinimizer, logFile, LogFileActive, Selection, TestCase, TestDataSet, TestReal, Transformation3_FromSymmetryOperationString, Transformation3Container
from pMolecule        import CrystalClassOrthorhombic, CrystalClassTetragonal, CrystalClassHexagonal, CrystalClassMonoclinic, CrystalClassRhombohedral, CrystalClassTriclinic, \
                             ElectronicState, MMModelOPLS, NBModelABFS, SystemGeometryObjectiveFunction
from pMoleculeScripts import CrystalAnalyzeTransformations, CrystalExpandToP1

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . Paths.
_dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "molecularCrystals" )

# . Transformation sets as strings.
_operations = { "C2"      : [ "(x,y,z)", "(-x,y,-z)"],
                "I4"      : [ "(x,y,z)", "(-x,-y,z)", "(-y,x,z)", "(y,-x,z)", "(x+1/2,y+1/2,z+1/2)", "(-x+1/2,-y+1/2,z+1/2)", "(-y+1/2,x+1/2,z+1/2)", "(y+1/2,-x+1/2,z+1/2)" ],
                "P1"      : [ "(x,y,z)" ],
                "P21a"    : [ "(x,y,z)", "(-x+1/2,y+1/2,-z)", "(-x,-y,-z)", "(x+1/2,-y+1/2,z)" ],
                "P212121" : [ "(x,y,z)", "(-x+1/2,-y,z+1/2)", "(-x,y+1/2,-z+1/2)", "(x+1/2,-y+1/2,-z)" ],
                "P61"     : [ "(x,y,z)", "(-y,x-y,z+1/3)", "(-x+y,-x,z+2/3)", "(-x,-y,z+1/2)", "(y,-x+y,z+5/6)", "(x-y,x,z+1/6)" ],
                "P21c"    : [ "(x,y,z)", "(-x,y+1/2,-z+1/2)", "(-x,-y,-z)", "(x,-y+1/2,z+1/2)" ],
                "P21a"    : [ "(x,y,z)", "(-x+1/2,y+1/2,-z)", "(-x,-y,-z)", "(x+1/2,-y+1/2,z)" ],
                "R3"      : [ "(x,y,z)", "(z,x,y)", "(y,z,x)" ] }

# . Convert to transformation sets.
_transformationSets = {}
for ( key, olist ) in _operations.iteritems ( ):
    tset = []
    for o in olist:
        tset.append ( Transformation3_FromSymmetryOperationString ( o ) )
    _transformationSets[key] = Transformation3Container.WithTransformations ( tset )

# . Molecule definitions.
_keywordLabels = ( "label", "canQCOptimize", "qcSelection", "qcCharge", "vacuumEnergy", "crystalEnergy", "crystalClass", "transformations", "crystalParameters" )
_moleculeData  = ( ( "ALAALA",   False, Selection.FromIterable ( range (  6, 10 ) ),  0,  -178.7095,  -677.7298, CrystalClassTetragonal   ( ), _transformationSets["I4"],      { "a" : 17.9850, "c" : 5.1540 } ),
                   ( "ALAMET01", False, Selection.FromIterable ( range ( 16, 27 ) ),  0,  -135.0084,  -652.9681, CrystalClassMonoclinic   ( ), _transformationSets["P21c"],    { "a" : 13.089, "b" : 5.329, "c" : 15.921, "beta" :  108.57 } ),
                   ( "AQARUF",   False, Selection.FromIterable ( range (  6, 19 ) ),  0,  -144.5558,  -608.3551, CrystalClassHexagonal    ( ), _transformationSets["P61"],     { "a" : 14.3720, "c" : 9.8282 } ),
                   ( "BEVXEF01", False, Selection.FromIterable ( range ( 23, 29 ) ),  0,  -240.0604,  -878.7555, CrystalClassOrthorhombic ( ), _transformationSets["P212121"], { "a" : 9.6590, "b" : 9.6720, "c" : 10.7390 } ),
                   ( "GLYALB",   False, Selection.FromIterable ( range ( 13, 17 ) ),  0,   -44.8402,  -562.9905, CrystalClassOrthorhombic ( ), _transformationSets["P212121"], { "a" : 9.6930, "b" : 9.5240, "c" : 7.5370 } ),
                   ( "GLYGLY",   False, Selection.FromIterable ( range ( 14, 17 ) ), -1,  -184.8863,  -687.4449, CrystalClassMonoclinic   ( ), _transformationSets["P21a"],    { "a" : 7.812, "b" : 9.566, "c" : 9.410, "beta" :  124.60 } ),
                   ( "GUFQON",   False, Selection.FromIterable ( range ( 13, 18 ) ),  0,  -191.1490,  -511.8117, CrystalClassOrthorhombic ( ), _transformationSets["P212121"], { "a" : 7.2750, "b" : 9.0970, "c" : 10.5070 } ),
                   ( "HXACAN19", True,  Selection.FromIterable ( range ( 16, 20 ) ),  0,   157.7170,    24.6588, CrystalClassMonoclinic   ( ), _transformationSets["P21a"],    { "a" : 12.8720, "b" : 9.3700, "c" : 7.0850, "beta" :  115.6200 } ),
                   ( "IWANID",   False, Selection.FromIterable ( range ( 37, 51 ) ),  0,  4380.9282,  4070.8472, CrystalClassMonoclinic   ( ), _transformationSets["C2"],      { "a" : 23.091, "b" : 5.494, "c" : 17.510, "beta" :  117.88 } ),
                   ( "LCDMPP10", True,  Selection.FromIterable ( range (  4,  8 ) ),  0,   157.3277,    40.3629, CrystalClassTriclinic    ( ), _transformationSets["P1"],      { "a" : 8.067, "b" : 6.082, "c" : 5.155, "alpha" : 131.7, "beta" :  82.4, "gamma" : 106.6 } ),
                   ( "WIRYEB",   False, Selection.FromIterable ( range (  6, 16 ) ),  0,  -155.7355,  -612.6186, CrystalClassHexagonal    ( ), _transformationSets["P61"],     { "a" : 14.4240, "c" : 9.9960 } ),
                   ( "WABZOO",   True,  Selection.FromIterable ( range ( 0, 4 ) + range ( 10, 14 ) + range ( 18, 25 ) ), 0, 209.9697,   152.9967, CrystalClassRhombohedral ( ), _transformationSets["R3"], { "a" : 12.5940 , "alpha" : 118.0300 } ) )

# . Options for expanding to P1 symmetry.
_defaultRange = range ( 1 ) #range ( -1, 2 )
_aRange       = _defaultRange
_bRange       = _defaultRange
_cRange       = _defaultRange
_numberCells  = len ( _aRange ) * len ( _bRange ) * len ( _cRange )

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Geometry optimization.
_LogFrequency      =  1000
_NumberSteps       = 20000
_GradientTolerance = 0.1

# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-02

#===================================================================================================================================
# . Class for a crystal test system.
#===================================================================================================================================
class CrystalTestSystem ( object ):
    """Crystal test system."""

    def __init__ ( self, **kwargs ):
        """Constructor."""
        for ( attribute, value ) in kwargs.iteritems ( ):
            setattr ( self, attribute, value )
        self.mmModel  = MMModelOPLS ( )
        self.nbModel  = NBModelABFS ( )
        self.p1Factor = 1.0

    def GetSystem ( self, doQCMM = False, doQCQC = False, expandToP1 = False, log = logFile, qcModel = None, useSymmetry = True ):
        """Get the system with the energy model defined."""
        # . Read the molecule.
        molecule              = AmberTopologyFile_ToSystem  ( os.path.join ( _dataPath, self.label + ".top" ), mmModel = self.mmModel, log = log )
        molecule.coordinates3 = AmberCrdFile_ToCoordinates3 ( os.path.join ( _dataPath, self.label + ".crd" ), log = log )
        molecule.label        = self.label
        # . Set up symmetry.
        if useSymmetry:
            kwargs = { "crystalClass" : self.crystalClass, "transformations" : self.transformations }
            kwargs.update ( self.crystalParameters )
            molecule.DefineSymmetry ( **kwargs )
            if expandToP1:
                self.p1Factor = float ( len ( molecule.symmetry.transformations ) * _numberCells )
                molecule      = CrystalExpandToP1 ( molecule, aRange = _aRange, bRange = _bRange, cRange = _cRange )
        # . Set up the QC model.
        if qcModel is not None:
            # . QC/MM.
            if doQCMM:
                if self.qcCharge != 0: molecule.electronicState = ElectronicState ( charge = self.qcCharge )
                # . For the tests do not worry if the MM charge is not zero or integral.
                try:    molecule.DefineQCModel ( qcModel, qcSelection = self.qcSelection )
                except: pass
            else:
                molecule.DefineQCModel ( qcModel )
        # . Set up the NB model.
        if ( qcModel is None ) or doQCMM or doQCQC: molecule.DefineNBModel  ( self.nbModel )
        # . Summary.
        if LogFileActive ( log ):
            molecule.Summary ( log = log )
            log.Paragraph ( "\nFormula = " + molecule.atoms.FormulaString ( ) + "." )
        # . Finish up.
        return molecule

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CrystalMMTest ( TestCase ):
    """A test case for calculating crystal MM energies."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( CrystalMMTest, self ).__init__ ( *args )
        # . Options.
        self.checkEnergies    = True
        self.doQCMM           = False
        self.doQCQC           = False
        self.geometryOptimize = True
        self.expandToP1       = False
        self.qcModel          = None
        self.useSymmetry      = True
        self.testGradients    = True

    def GeometryOptimizationSummary ( self, optimizationResults, log = logFile ):
        """Output for geometry optimization results."""
        if LogFileActive ( log ) and ( len ( optimizationResults ) > 0 ):
            # . Parameters.
            tableEnergyWidth    = 20
            tableIntegerWidth   = 20
            tableNameWidth      = 30
            tableParameterWidth = 20
            # . System labels.
            labels = optimizationResults.keys ( )
            labels.sort ( )
            # . Energies.
            # . Header.
            table = log.GetTable ( columns = [ tableNameWidth, tableEnergyWidth, tableEnergyWidth, tableEnergyWidth, tableIntegerWidth,tableIntegerWidth ] )
            table.Start   ( )
            table.Title   ( "Optimized Energies" )
            table.Heading ( "Label"           )
            table.Heading ( "Initial Energy"  )
            table.Heading ( "Final Energy"    )
            table.Heading ( "Energy Lowering" )
            table.Heading ( "Function Calls"  )
            table.Heading ( "Convergence"     )
            # . Data.
            for label in labels:
                data = optimizationResults[label]
                table.Entry ( label, alignment = "left" )
                table.Entry ( "{:20.4f}".format ( data["Initial Energy"]                          ) )
                table.Entry ( "{:20.4f}".format ( data["Final Energy"  ]                          ) )
                table.Entry ( "{:20.4f}".format ( data["Final Energy"  ] - data["Initial Energy"] ) )
                table.Entry ( "{:d}".format ( data["Function Calls"] ) )
                if data["Converged"]: table.Entry ( "T" )
                else:                 table.Entry ( "F" )
            table.Stop ( )
            # . Symmetry parameters.
            if self.useSymmetry:
                table = log.GetTable ( columns = [ tableNameWidth, tableParameterWidth, tableEnergyWidth, tableEnergyWidth, tableEnergyWidth ] )
                table.Start  ( )
                table.Title  ( "Optimized Symmetry Parameters" )
                table.Heading ( "Label"      )
                table.Heading ( "Parameter"  )
                table.Heading ( "Initial"    )
                table.Heading ( "Final"      )
                table.Heading ( "Difference" )
                for label in labels:
                    # . Get the data.
                    data      = optimizationResults[label]
                    spInitial = data["Initial Parameters"]
                    spFinal   = data["Final Parameters"  ]
                    spLabels  = spInitial.keys ( )
                    spLabels.sort ( )
                    for ( i, spLabel ) in enumerate ( spLabels ):
                        if i == 0: table.Entry ( label, alignment = "left" )
                        else:      table.Entry ( "" )
                        table.Entry ( spLabel, alignment = "left" )
                        table.Entry ( "{:.4f}".format ( spInitial[spLabel]                    ) )
                        table.Entry ( "{:.4f}".format ( spFinal  [spLabel]                    ) )
                        table.Entry ( "{:.4f}".format ( spInitial[spLabel] - spFinal[spLabel] ) )
                table.Stop ( )

    def runTest ( self ):
        """The test."""

        # . Initialization.
        log = self.GetLog ( )
        if self.checkEnergies:
            energyDifferences   = {}
        if self.geometryOptimize:
            optimizationResults = {}
            self.optimizer.Summary ( log = log )
        if self.testGradients:
            maximumGradientDeviation = 0.0

        # . Loop over the molecules.
        numberEnergyFailures = 0
        for testSystem in self.crystalTestSystems:

            # . Get the molecule.
            molecule = testSystem.GetSystem ( doQCMM = self.doQCMM, doQCQC = self.doQCQC, log = log, qcModel = self.qcModel, useSymmetry = self.useSymmetry )

            # . Analyze the transformations.
            CrystalAnalyzeTransformations ( molecule, log = log )

            # . Calculate the energy and check it.
            try:
                energy   = molecule.Energy ( log = log, doGradients = True )
                energyOK = True
                if self.checkEnergies:
                    if self.useSymmetry: energyDifferences[molecule.label] = math.fabs ( energy - testSystem.p1Factor * testSystem.crystalEnergy )
                    else:                energyDifferences[molecule.label] = math.fabs ( energy -                       testSystem.vacuumEnergy  )
            except:
                energyOK = False
                numberEnergyFailures += 1
                if log is not None: log.Paragraph ( "Error calculating initial energy." )

            # . Skip if not OK.
            if not energyOK: continue

            # . Test the gradients.
            # . The gradients are definitely more sensitive to the finite-difference step when using fractional coordinates.
            if self.testGradients:
                of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                of.IncludeSymmetryParameters ( )
                deviation = of.TestGradients ( delta = 5.0e-6, log = log, tolerance = 3.0e-4 )
                maximumGradientDeviation = max ( maximumGradientDeviation, deviation )

            # . Geometry optimize.
            if self.geometryOptimize and ( ( self.qcModel is None ) or ( self.doQCMM ) or ( self.doQCQC ) or ( ( self.qcModel is not None ) and ( testSystem.canQCOptimize ) ) ):
                # . Save the symmetry parameters for later.
                if self.useSymmetry:
                    spInitial = molecule.symmetry.crystalClass.GetUniqueSymmetryParameters ( molecule.configuration.symmetryParameters )
                of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                if self.useSymmetry: of.IncludeSymmetryParameters ( )
                state               = self.optimizer.Iterate ( of, log = log )
                finalEnergy         = molecule.Energy ( log = log, doGradients = True )
                optimizationResults[molecule.label] = { "Initial Energy" : energy, "Final Energy" : finalEnergy, "Function Calls" : state["Function Calls"], "Converged" : state["Converged"] }
                if self.useSymmetry:
                    spFinal = molecule.symmetry.crystalClass.GetUniqueSymmetryParameters ( molecule.configuration.symmetryParameters )
                    optimizationResults[molecule.label].update ( { "Initial Parameters" : spInitial, "Final Parameters" : spFinal } )

        # . Results.
        # . Energy differences.
        if self.checkEnergies and ( log is not None ):
            # . System labels.
            labels = energyDifferences.keys ( )
            labels.sort ( )
            # . Header.
            table = log.GetTable ( columns = [ 30, 20 ] )
            table.Start  ( )
            if self.useSymmetry: table.Title  ( "Differences in Crystal Energies" )
            else:                table.Title  ( "Differences in Vacuum Energies"  )
            table.Heading ( "System"     )
            table.Heading ( "Difference" )
            # . Data
            for label in labels:
                table.Entry ( label, alignment = "left" )
                table.Entry ( "{:20.4f}".format ( energyDifferences[label] ) )
            table.Stop ( )

        # . Optimizations.
        if self.geometryOptimize:
            self.GeometryOptimizationSummary ( optimizationResults, log = log )

        # . Get the observed and reference data.
        observed      = {}
        referenceData = TestDataSet ( "Crystal Energies" )
        if self.checkEnergies:
            observed.update ( energyDifferences )
            for label in energyDifferences:
                referenceData.AddDatum ( TestReal ( label, 0.0, referenceData, absoluteErrorTolerance = _EnergyTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )
        if self.testGradients:
            observed["Gradient Error"] = maximumGradientDeviation
            referenceData.AddDatum ( TestReal ( "Gradient Error", 0.0, referenceData, absoluteErrorTolerance = _GradientAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )

        # . Check for success/failure.
        if len ( observed ) > 0:
            results = referenceData.VerifyAgainst ( observed )
            results.Summary ( log = log, fullSummary = self.fullVerificationSummary )
            isOK    = results.WasSuccessful ( )
        else:
            isOK    = True
        isOK = isOK and ( numberEnergyFailures == 0 )
        self.assertTrue ( isOK )

    def setUp ( self ):
        """Set up the calculation."""
        # . Set up the systems.
        self.crystalTestSystems = []
        for values in _moleculeData:
            kwargs = { key : value for ( key, value ) in zip ( _keywordLabels, values ) }
            self.crystalTestSystems.append ( CrystalTestSystem ( **kwargs ) )
        # . Set up the optimizer.
        if self.geometryOptimize:
            self.optimizer = LBFGSMinimizer ( logFrequency         = _LogFrequency      ,
                                              maximumIterations    = _NumberSteps       ,
                                              rmsGradientTolerance = _GradientTolerance )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = CrystalMMTest ( )
    test.run ( )
