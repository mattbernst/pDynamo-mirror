"""Grid updating tests."""

import math, os

from pBabel    import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, MOLFile_ToSystem, MOL2File_ToSystem, XYZFile_ToCoordinates3
from pCore     import CPUTime, logFile, PairListGenerator, TestCase
from pMolecule import CrystalClassCubic, MMModelOPLS, NBModelABFS

#===================================================================================================================================
# . Test system.
#===================================================================================================================================
class GridUpdatingTestSystem ( object ):
    """Grid updating test system."""

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        self.__dict__.update ( keywordArguments )

    def GetSystem ( self, log = logFile ):
        """Get the system."""
        if self.setUpType == "CHARMM":
            parameters = CHARMMParameterFiles_ToParameters ( self.parameterFiles, log = log )
            system     = CHARMMPSFFile_ToSystem ( self.setUpFile, isXPLOR = True, log = log, parameters = parameters )
        elif self.setUpType == "MOL"   : system = MOLFile_ToSystem  ( self.setUpFile )
        elif self.setUpType == "MOL2"  : system = MOL2File_ToSystem ( self.setUpFile )
        if self.xyzFile is not None: system.coordinates3 = XYZFile_ToCoordinates3 ( self.xyzFile )
        if self.mmModel is not None: system.DefineMMModel ( self.mmModel )
        if self.hasSymmetry: system.DefineSymmetry ( crystalClass = CrystalClassCubic ( ), a = self.a )
        system.label = self.label
        system.Summary ( log = log )
        return system

#===================================================================================================================================
# . Test systems.
#===================================================================================================================================
# . Test system definitions.
dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "gridUpdating" )
testSystemDefinitions = ( { "hasSymmetry"    : False     ,
                            "label"          : "Crambin" ,
                            "mmModel"        : MMModelOPLS ( "protein" ) ,
                            "setUpFile"      : os.path.join ( dataPath, "crambin.mol" ) ,
                            "setUpType"      : "MOL" ,
                            "xyzFile"        : None  },
                          { "a"              : 40.0 ,
                            "hasSymmetry"    : True ,
                            "label"          : "Water Box 40x40x40" ,
                            "mmModel"        : MMModelOPLS ( "protein" ) ,
                            "setUpFile"      : os.path.join ( dataPath, "waterBox40x40x40.mol2" ) ,
                            "setUpType"      : "MOL2" ,
                            "xyzFile"        : None   },
                          { "a"              : 62.23 ,
                            "hasSymmetry"    : True  ,
                            "label"          : "DHFR Benchmark" ,
                            "mmModel"        : None ,
                            "parameterFiles" : [ os.path.join ( dataPath, "par_all22_prot.inp" ) ] ,
                            "setUpFile"      : os.path.join ( dataPath, "dhfr.psfx" ) ,
                            "setUpType"      : "CHARMM" ,
                            "xyzFile"        : os.path.join ( dataPath, "dhfr.xyz" ) } )

# . Define the test systems.
testSystems = []
for keywordArguments in testSystemDefinitions:
    testSystems.append ( GridUpdatingTestSystem ( **keywordArguments ) ) 

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GridUpdatingTest ( TestCase ):
    """A test case for calculating crystal MM energies."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        log = self.GetLog ( )

        # . Set up the generators.
        _NVALUES = 11
        generators = [ ( "Direct" , PairListGenerator ( minimumPoints = 100000 ) ) ]
        for i in range ( _NVALUES ):
            s = float ( i + 1 ) * 0.25
            generators.append ( ( "CC {:6.3f}".format ( s ), PairListGenerator ( cutoffCellSizeFactor = s, minimumCellExtent = 0, minimumPoints = 100, sortIndices = False, useGridByCell = True ) ) )
        for i in range ( _NVALUES ):
            s = float ( i + 1 ) * 0.25
            generators.append ( ( "PC {:6.3f}".format ( s ), PairListGenerator ( cutoffCellSizeFactor = s, minimumCellExtent = 0, minimumPoints = 100, sortIndices = False, useGridByCell = False ) ) )

        # . Loop over test systems.
        numberEnergyFailures = 0
        for testSystem in testSystems:

            # . Get the system.
            system = testSystem.GetSystem ( log = log )

            # . Energies.
            energies = []
            for ( label, generator ) in generators:
                nbModel = NBModelABFS ( generator = generator )
                nbModel.Summary ( log = log )
                system.DefineNBModel ( nbModel )
                cpuTime = CPUTime ( )
                e = system.Energy  ( doGradients = True, log = log )
                c = cpuTime.CurrentAsString ( )
                energies.append ( ( e, c, label ) )
                system.configuration.nbState.StatisticsSummary ( log = log )

            # . Print summary of results.
            if log is not None:
                e0 = energies[0][0]
                n  = 0
                table = log.GetTable ( columns = [ 6, 10, 16, 16, 16, 20 ] )
                table.Start  ( )
                table.Title  ( "Energies and CPU Times for {:s}:".format ( system.label ) )
                for ( i, ( e, c, l ) ) in enumerate ( energies ):
                    deltaE = math.fabs ( e - e0 )
                    if deltaE > 1.0e-2: n += 1
                    table.Entry ( "{:d}".format ( i ) )
                    table.Entry (  l  )
                    table.Entry ( "{:15.3f}".format ( e      ) )
                    table.Entry ( "{:15.3f}".format ( e0     ) )
                    table.Entry ( "{:15.3f}".format ( deltaE ) )
                    table.Entry ( "    {:s}".format ( c      ) )
                table.Stop ( )
                if log is not None:
                    if n == 0: log.Paragraph ( "All tests were successful."      )
                    else:      log.Paragraph ( "{:d} tests failed.".format ( n ) )
                numberEnergyFailures += n

        # . Check for success/failure.
        self.assertTrue ( ( numberEnergyFailures == 0 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = GridUpdatingTest ( )
    test.run ( )
