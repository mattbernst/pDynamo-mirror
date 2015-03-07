"""Run the book examples."""

# . Some way of verifying the results of the examples is necessary (to do away with DifferenceLogs).
# . Move this and DifferenceLogs to installation.

import glob, os.path, subprocess, time

from Definitions import errorExtension, errorPath, logExtension, logPath, scratchPath, xhtmlExtension
from optparse    import OptionParser
from pCore       import CPUTime

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Additional long examples.
# . 26Q before 26LJ.
_AdditionalLongExamples = [ "Example26Q.py", "Example26LJ.py", "NewExample27.py" ]

# . The examples which take a long time.
_LongExamples = [ 18, 20, 21, 22, 25, 26, 27, 28 ]

# . Largest example index.
_MaximumIndex = 28

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def RunExamples ( ):
    """Main function."""

    # . Parse the command line.
    epilog = "The arguments are optional but, if present, should be the names of the examples to run. By default all short examples are run."
    usage  = "%prog [arg1 [arg2 ...]]"

    # . For version 2.5.
    try:    parser = OptionParser ( epilog = epilog, usage = usage )
    except: parser = OptionParser (                  usage = usage )

    # . Options.
    parser.add_option ( "-a", "--all",   action = "store_true", dest = "runAll",   default = False, help = "run all examples"             )
    parser.add_option ( "-l", "--long",  action = "store_true", dest = "runLong",  default = False, help = "run long examples"            )
    parser.add_option ( "-s", "--short", action = "store_true", dest = "runShort", default = False, help = "run short examples [default]" )

    # . Parsing.
    ( options, args ) = parser.parse_args ( )

    # . Initialization.
    toRun = []

    # . The presence of arguments overrides anything else.
    # . Just run those examples that were specified.
    if len ( args ) > 0:
        for arg in args:
            if not arg.endswith ( ".py" ): arg += ".py"
            if os.path.exists ( arg ): toRun.append ( arg )
            else: raise ValueError ( "Unknown example - " + arg + "." )
    # . Process flag options only.
    else:

        # . Get the indices of the examples to run.
        # . All examples.
        if options.runAll or ( options.runLong and options.runShort ):
            indices = range ( _MaximumIndex + 1 )
        # . Long examples.
        elif options.runLong:
            indices = _LongExamples
        # . Short examples - everything else.
        else:
            indices = list ( set ( range ( _MaximumIndex + 1 ) ).difference ( set ( _LongExamples ) ) )
            indices.sort ( )

        # . Generate the file names.
        for index in indices:
            toRun.append ( "Example{:d}.py".format ( index ) )
            xhtmlPath = "Example{:d}XHTML.py".format ( index )
            if os.path.exists ( xhtmlPath ): toRun.append ( xhtmlPath )

        # . Add additional long examples.
        if ( options.runAll or options.runLong ): toRun.extend ( _AdditionalLongExamples )

    # . Run the examples.
    if len ( toRun ) > 0:

        # . Make sure the error and log directories exist.
        if not os.path.exists ( errorPath ): os.mkdir ( errorPath )
        if not os.path.exists ( logPath   ): os.mkdir ( logPath   )

        # . Header.
        print ( "\nRunning examples ..." )
        if options.runLong: print ( "\nPlease be patient ... this will take some time ..."  )
        else:               print ( "\nThis may take a few minutes ..."                     )
        print ( "\nLogs and other files are in " + os.path.dirname ( scratchPath ) + "." )
        print ( "\n" + 80 * "-"  )
        print ( "Book Examples:" )
        print ( 80 * "-" + "\n"  )

        # . Loop over files.
        labelWidth      = 50
        numberSuccesses = 0
        totalTime       = 0
        for inFile in toRun:

            # . Get the error file name.
            errorFile = os.path.join ( errorPath, inFile[0:-3] + errorExtension )

            # . Get the log extension and log file name.
            if inFile.find ( "XHTML" ) >= 0: extension = xhtmlExtension
            else:                            extension = logExtension
            logFile = os.path.join ( logPath, inFile[0:-3] + extension )

            # . Remove existing log of this name.
            if os.path.exists ( errorFile ): os.remove ( errorFile )
            if os.path.exists ( logFile   ): os.remove ( logFile   )

            # . Run the job.
            isOK  = True
            time0 = time.time ( )
            eFD   = open ( errorFile, "w" )
            oFD   = open ( logFile  , "w" )
            try:
                process = subprocess.Popen ( [ "python", inFile ], stderr = eFD, stdout = oFD )
                process.wait ( )
            except Exception as e:
                eFD.write ( e )
                isOK = False
            elapsedTime = time.time ( ) - time0
            totalTime  += elapsedTime

            # . Close files.
            eFD.close ( )
            oFD.close ( )

            # . Check for a non-empty error file.
            if isOK:
                eFD   = open ( errorFile, "r" )
                lines = eFD.readlines ( )
                eFD.close ( )
                isOK  = ( len ( lines ) == 0 )
                if isOK: os.remove ( errorFile )

            # . Get the result.
            if isOK:
                numberSuccesses += 1
                resultLabel = "Pass"
            else:
                resultLabel = "Error"

            # . Printing.
            print ( "{:s}{:<10s}{:>20s}".format ( inFile[0:-3].ljust ( labelWidth ), resultLabel, CPUTime.TimeToString ( elapsedTime ) ) )

        # . Print a terminating message.
        print ( "\n" + 80 * "-" )
        if numberSuccesses == 1:
            if numberSuccesses == len ( toRun ): print ( "\nThe example was successful." )
            else:                                print ( "\n{:d} example out of {:d} was successful.".format ( numberSuccesses, len ( toRun ) ) )
        else:
            if numberSuccesses == len ( toRun ): print ( "\nAll {:d} examples were successful.".format ( len ( toRun ) ) )
            else:                                print ( "\n{:d} examples out of {:d} were successful.".format ( numberSuccesses, len ( toRun ) ) )
        print ( "\nTotal time used: " + CPUTime.TimeToString ( totalTime, compact = False ) + "." )

    # . No examples run.
    else:
        print ( "\nNo examples were found." )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    # . Run the script.
    RunExamples ( )
