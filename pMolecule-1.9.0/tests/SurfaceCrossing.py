"""Surface crossing test."""

import glob, math, os.path

from pBabel    import MOLFile_ToCoordinates3, MOLFile_ToSystem
from pCore     import CPUTime, logFile, LogFileActive, QuasiNewtonMinimizer, TestCase
from pMolecule import DIISSCFConverger, ElectronicState, QCModelMNDO, SEAMObjectiveFunction

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SurfaceCrossingTest ( TestCase ):
    """Surface crossing optimization test."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK    = True
        log     = self.GetLog ( )
        molPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "mol" )

        # . Options.
        converger = DIISSCFConverger ( densityTolerance = 1.0e-6, maximumSCFCycles = 500 )
        qcModel   = QCModelMNDO      ( "am1", converger = converger, isSpinRestricted = False )
        singlet   = ElectronicState  ( charge = 1, multiplicity = 1 )
        triplet   = ElectronicState  ( charge = 1, multiplicity = 3 )

        # . Optimizer.
        optimizer = QuasiNewtonMinimizer ( logFrequency         =   1  ,
                                           maximumIterations    = 500  ,
                                           rmsGradientTolerance = 0.05 )
        optimizer.Summary ( log = log )

        # . Set up the system.
        system                 = MOLFile_ToSystem ( os.path.join ( molPath, "phenylCation.mol" ) )
        system.electronicState = singlet
        system.label           = "Phenyl Cation"
        system.DefineQCModel ( qcModel )
        system.Summary ( log = log )

        # . Check both methods.
        numberNotConverged = 0
        results            = {}
        for method in ( "GP", "PF" ):

            # . Reset coordinates.
            system.coordinates3 = MOLFile_ToCoordinates3 ( os.path.join ( molPath, "phenylCation.mol" ) )
            system.configuration.Clear ( )

            # . Set up the objective function.
            seamOF = SEAMObjectiveFunction.FromSystem ( system, singlet, triplet, method = method )
            #seamOF.RemoveRotationTranslation ( )

            # . Minimize.
            #seamOF.TestGradients ( delta = 1.0e-05 ) # . Works with 1.0e-10 density tolerance.
            cpu                = CPUTime ( )
            report             = optimizer.Iterate ( seamOF, log = log )
            report["CPU Time"] = cpu.CurrentAsString ( )

            # . Final energies.
            ( f1, f2 )         = seamOF.Energies ( doGradients = True, log = log )
            report["Energy 1"] = f1
            report["Energy 2"] = f2
            results[method]    = report
            if not report.get ( "Converged", False ): numberNotConverged += 1

        # . Print out a summary of the results.
        if LogFileActive ( log ):
            table = log.GetTable ( columns = [ 10, 20, 20, 10, 10, 20 ] )
            table.Start   ( )
            table.Title   ( "Surface Crossing Optimizations" )
            table.Heading ( "Method"    )
            table.Heading ( "State Energies", columnSpan = 2 )
            table.Heading ( "Converged" )
            table.Heading ( "Calls"     )
            table.Heading ( "Time"      )
            for method in ( "GP", "PF" ):
                report = results[method]
                table.Entry ( method, alignment = "left" )
                table.Entry ( "{:20.1f}".format ( report["Energy 1"      ] ) )
                table.Entry ( "{:20.1f}".format ( report["Energy 2"      ] ) )
                table.Entry ( "{!r}"    .format ( report.get ( "Converged", False ) ) )
                table.Entry ( "{:d}"    .format ( report["Function Calls"] ) )
                table.Entry (                     report["CPU Time"      ]   )
            table.Stop ( )

        # . Finish up.
        self.assertTrue ( numberNotConverged == 0 )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = SurfaceCrossingTest ( )
    test.run ( )
