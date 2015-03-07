"""pDynamo script to run package tests."""

import glob, os.path, subprocess, sys, types

from optparse import OptionParser
from pCore    import CPUTime, TestCase, TestResult

#===================================================================================================================================
# . Parameters
#===================================================================================================================================
# . Extensions.
_logExtension = ".log"
_pklExtension = ".pkl"

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def RunPackageTests ( ):
    """Main function."""

    # . Parse the command line.
    epilog = "The arguments are optional but, if present, are either package or test case names. By default short tests are run for all packages."
    usage  = "%prog [options] [arg1 [arg2 ...]]"

    # . For version 2.5.
    try:    parser = OptionParser ( epilog = epilog, usage = usage )
    except: parser = OptionParser (                  usage = usage )

    # . Options.
    parser.add_option ( "-a", "--all"     , action = "store_true", dest = "runAll"                  , default = False, help = "run all tests"                                    )
    parser.add_option ( "-f", "--full"    , action = "store_true", dest = "fullVerificationSummary" , default = False, help = "print a full verification summary"                )
    parser.add_option ( "-g", "--generate", action = "store_true", dest = "generateReferenceData"   , default = False, help = "generate reference data instead of testing"       )
    parser.add_option ( "-l", "--long"    , action = "store_true", dest = "runLong"                 , default = False, help = "run long tests"                                   )
    parser.add_option ( "-p", "--package" ,                        dest = "packageName"             , default = None , help = "test case package name"   , metavar = "PACKAGE"   )
    parser.add_option ( "-s", "--short"   , action = "store_true", dest = "runShort"                , default = True , help = "run short tests [default]"                        )

    # . Parsing.
    ( options, args ) = parser.parse_args ( )

    # . Set the all/long/short options.
    runLong  = False
    runShort = False
    if options.runAll or ( options.runLong and options.runShort ):
        runLong  = True
        runShort = True
    elif options.runLong:
        runLong  = True
    elif options.runShort:
        runShort = True

    # . Find root directory.
    rootPath = os.getenv ( "PDYNAMO_ROOT" )
    if rootPath is None: raise ValueError ( "Please define pDynamo environment variables before running the tests." )

    # . Find the packages whose tests are to be run.
    if options.packageName is not None:
        packages = glob.glob ( os.path.join ( rootPath, options.packageName + "*" ) )
    elif len ( args ) > 0:
        packages = []
        for arg in args:
            paths = glob.glob ( os.path.join ( rootPath, arg + "-*" ) )
            if len ( paths ) == 1: packages.append ( paths[0] )
    else:
        packages = glob.glob ( os.path.join ( rootPath, "p*" ) )

    # . Get the packages which have test directories.
    testPackages = []
    for package in packages:
        testPath = os.path.join ( package, "tests" )
        if os.path.exists ( testPath ): testPackages.append ( package )

    # . Get all potential testcases.
    testCases = {}
    for testPackage in testPackages:
        ( head, packageName ) = os.path.split ( testPackage )
        # . All ".py" files.
        if options.packageName is None:
            tests = glob.glob ( os.path.join ( testPackage, "tests", "*.py" ) )
        # . All files specified in args.
        else:
            tests = []
            for arg in args:
                if not arg.endswith ( ".py" ): arg += ".py"
                paths = glob.glob ( os.path.join ( testPackage, "tests", arg ) )
                if len ( paths ) == 1: tests.append ( paths[0] )
        if len ( tests ) > 0:
            tests.sort ( )
            testCases[packageName] = tests

    # . Initialization.
    totalTests = 0

    # . There are potentially test cases.
    if ( len ( testCases ) > 0 ) and ( runLong or runShort ):

        # . Initialization.
        numberSuccesses = 0
        totalTime       = 0.0

        # . Define the directory for output.
        scratchPath = os.getenv ( "PDYNAMO_SCRATCH" )
        if not os.path.exists ( scratchPath ): os.mkdir ( scratchPath )
        scratchPath = os.path.join ( scratchPath, "packageTests" )
        if not os.path.exists ( scratchPath ): os.mkdir ( scratchPath )

        # . Header.
        if options.generateReferenceData: print ( "\nGenerating reference data ..." )
        else:                             print ( "\nRunning tests ..."             )
        if runLong: print ( "\nPlease be patient ... this will take some time ..." )
        else:       print ( "\nThis may take a few minutes ..." )
        print ( "\nLogs and other files are in " + scratchPath + "." )

        # . Loop over packages.
        packageNames = testCases.keys ( )
        packageNames.sort ( )
        for packageName in packageNames:

            # . Set up the import path.
            tests = testCases[packageName]
            ( testPath, tail ) = os.path.split ( tests[0] )
            sys.path.insert ( 0, testPath )

            # . Define some paths.
            packagePath = os.path.join ( scratchPath, packageName      )
            errorPath   = os.path.join ( packagePath, "errors"         )
            logPath     = os.path.join ( packagePath, "logs"           )
            resultPath  = os.path.join ( packagePath, "generatedFiles" )
            if options.generateReferenceData: referenceDataPath = os.path.join ( packagePath                         , "reference" )
            else:                             referenceDataPath = os.path.join ( os.path.dirname ( testPath ), "data", "reference" )

             # . Gather the cases.
             # . This is equivalent to TestLoader producing a TestSuite.
            cases              = []
            maximumLabelLength = 0
            for test in tests:
                ( head, tail ) = os.path.split ( test )
                moduleName = tail[0:-3]
                module = __import__ ( moduleName, globals(), locals(), [], -1 )
                for name in dir ( module ):
                    obj = getattr ( module, name )
                    if isinstance ( obj, ( type, types.ClassType ) ) and issubclass ( obj, TestCase ) and ( obj is not TestCase ):
                        case = obj ( )
                        if ( runLong or ( runShort and case.MakeShort ( ) ) ) and ( ( options.generateReferenceData and case.GenerateReferenceData ( ) ) or ( not options.generateReferenceData ) ):
                            case.fullVerificationSummary = options.fullVerificationSummary
                            case.outputPath              = os.path.join ( logPath          , moduleName + _logExtension )
                            case.referenceDataPath       = os.path.join ( referenceDataPath, moduleName + _pklExtension )
                            case.resultPath              = resultPath
                            cases.append ( case )
                            maximumLength = max ( maximumLabelLength, len ( case.label ) )
                        break

            # . Reset the import path.
            sys.path.pop ( 0 )

            # . Run the tests.
            if len ( cases ) > 0:

                # . Ensure all directories exist.
                if not os.path.exists ( packagePath ): os.mkdir ( packagePath )
                if not os.path.exists ( errorPath   ): os.mkdir ( errorPath   )
                if not os.path.exists ( logPath     ): os.mkdir ( logPath     )
                if not os.path.exists ( resultPath  ): os.mkdir ( resultPath  )
                if options.generateReferenceData:
                    if not os.path.exists ( referenceDataPath  ): os.mkdir ( referenceDataPath  )

                # . Header.
                print ( "\n" + 80 * "-" )
                print ( "Tests for package " + packageName + ":" )
                print ( 80 * "-" + "\n" )

                # . Run the tests.
                # . Equivalent to TestRunner.
                labelWidth = max ( 50, maximumLabelLength + 10 )
                for case in cases:
                    # . Remove existing error and log files.
                    errorLog = os.path.join ( errorPath, case.label + ".err" )
                    if os.path.exists ( errorLog        ): os.remove ( errorLog        )
                    if os.path.exists ( case.outputPath ): os.remove ( case.outputPath )
                    # . Run the test.
                    result = TestResult ( )
                    cpu    = CPUTime ( )
                    case.run ( result )
                    time   = cpu.Current ( )
                    # . Get the result of the test.
                    if len ( result.errors ) == 1:
                        resultLabel = "Error"
                        efile = open ( errorLog, "w" )
                        efile.write ( "Error in " + case.label + ":\n\n" )
                        for item in result.errors:
                            efile.write ( item[1] )
                        efile.close ( )
                    elif len ( result.failures ) == 1:
                        resultLabel = "Fail"
                    else:
                        numberSuccesses += 1
                        resultLabel = "Pass"
                    # . Finish up.
                    print ( "{:s}{:<10s}{:>20s}".format ( case.label.ljust ( labelWidth ), resultLabel, CPUTime.TimeToString ( time ) ) )
                    totalTests += 1
                    totalTime  += time

        # . Print a terminating message.
        if totalTests > 0:
            print ( "\n" + 80 * "-" )
            if numberSuccesses == 1:
                if numberSuccesses == totalTests: print ( "\nThe test was successful." )
                else:                             print ( "\n{:d} test out of {:d} was successful.".format ( numberSuccesses, totalTests ) )
            else:
                if numberSuccesses == totalTests: print ( "\nAll {:d} tests were successful.".format ( totalTests ) )
                else:                             print ( "\n{:d} tests out of {:d} were successful.".format ( numberSuccesses, totalTests ) )
            print ( "\nTotal time used: " + CPUTime.TimeToString ( totalTime, compact = False ) + "." )

    # . No tests run.
    if totalTests == 0:
        print ( "\nNo test cases were found." )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    # . Run the script.
    RunPackageTests ( )
