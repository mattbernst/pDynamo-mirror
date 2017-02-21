"""Run the tutorials."""

import glob, os.path, subprocess, time, sys

from argparse import ArgumentParser
from pCore    import CPUTime, YAMLUnpickle

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Definitions.
_ErrorExtension  = ".err"
_LabelWidth      = 50
_OutputExtension = ".log"
_PythonCommand   = sys.executable
_ScriptExtension = ".py"

# . Index file name.
_IndexFileName = "index.yaml"

# . Path names.
_ErrorPath     = "errors"
_OutputPath    = "logs"
_TutorialPath  = "tutorials"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class InstallationOptions ( object ):
    """Class for determining installation options."""

    def __init__ ( self ):
        """Constructor."""
        # . Basic initialization.
        self.maximumThreads    = self.MaximumNumberOfThreads ( )
        self.possibleTutorials = self.FindPossibleTutorials  ( )

    def FindPossibleTutorials ( self ):
        """Find the possible tutorials."""
        paths     = glob.glob ( "*" )
        tutorials = []
        for path in paths:
            if os.path.isdir ( path ) and os.path.isfile ( os.path.join ( path, _IndexFileName ) ): tutorials.append ( path )
        return set ( tutorials )

    @classmethod
    def FromCommandLine ( selfClass ):
        """Constructor by getting arguments from a command line."""
        # . Initialization.
        self = selfClass ( )
        # . Parse the command line.
        parser = ArgumentParser ( epilog = "The arguments are optional but, if present, should be the names of the tutorials to run. By default all tutorials are run." )
        parser.add_argument ( "-a" , "--all"        , action = "store_true" , dest = "useSerialAtlas" , default = True , help = "run all tutorials              [default]" )
        parser.add_argument (        "--threads"    , type   = int          , dest = "threads"        , default = 1    , help = "the number of threads"                    )
        parser.add_argument (        "tutorialName" , nargs  = "*"          ,                                            help = "the name of the tutorial to run"          )
        arguments = parser.parse_args ( )
        # . Find the number of threads.
        self.numberOfThreads = min ( max ( arguments.threads, 1 ), self.maximumThreads )
        # . Tutorials have been specified.
        if len ( arguments.tutorialName ) > 0:
            specified  = set ( arguments.tutorialName )
            found      = specified & self.possibleTutorials
            notFound   = specified - self.possibleTutorials
        # . Do everything.
        else:
            found      = self.possibleTutorials
            notFound   = set ( )
        self.found     = list ( found    )
        self.notFound  = list ( notFound )
        self.found.sort    ( )
        self.notFound.sort ( )
        # . Finish up.
        return self

    def MaximumNumberOfThreads ( self ):
        """Find the maximum number of threads."""
        try:
            import multiprocessing
            return multiprocessing.cpu_count ( )
        except:
            return 1

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def RunTutorials ( ):
    """Main function."""

    # . Get the run options.
    options = InstallationOptions.FromCommandLine ( )

    # . There are errors.
    if len ( options.notFound ) > 0:
        print ( "\nUnknown tutorials were specified:" )
        for name in options.notFound:
            print ( "  {:s}".format ( name ) )

    # . There are no tutorials.
    elif len ( options.found ) == 0:
        print ( "\nNo tutorials were found." )

    # . There are tutorials.
    else:

        # . Global header.
        print ( "\nRunning tutorials ..." )
        print ( "\nPlease be patient ... this will take some time ..." )

        # . Set the root directory.
        rootPath = os.path.join ( os.getenv ( "PDYNAMO_SCRATCH" ), _TutorialPath )
        if not os.path.exists ( rootPath ): os.mkdir ( rootPath )
        print ( "\nOutputs and other files are in " + rootPath + "." )

        # . Loop over tutorials.
        totalTime = 0.0
        for name in options.found:

            # . Set the output directories.
            namePath   = os.path.join ( rootPath, name        )
            errorPath  = os.path.join ( namePath, _ErrorPath  )
            outputPath = os.path.join ( namePath, _OutputPath )
            if not os.path.exists ( namePath   ): os.mkdir ( namePath   )
            if not os.path.exists ( errorPath  ): os.mkdir ( errorPath  )
            if not os.path.exists ( outputPath ): os.mkdir ( outputPath )

            # . Get the scripts.
            scripts = YAMLUnpickle ( os.path.join ( name, _IndexFileName ) )
            if len ( scripts ) <= 0: continue

            # . Local header.
            print ( "\n" + 80 * "-"   )
            print ( "Tutorial: {:s}".format ( name ) )
            print ( 80 * "-" + "\n"   )

            # . Loop over scripts.
            numberOfSuccesses = 0
            tutorialTime      = 0.0
            for items in scripts:

                # . Script name and path.
                if isinstance ( items, list ):
                    scriptName = items[-1]
                    scriptPath = os.path.join ( *items )
                else:
                    scriptName = items
                    scriptPath = items

                # . File names.
                errorFile  = os.path.join ( errorPath , scriptName + _ErrorExtension  )
                outFile    = os.path.join ( outputPath, scriptName + _OutputExtension )
                scriptFile = os.path.join ( name      , scriptPath + _ScriptExtension )
                if     os.path.exists ( errorFile  ): os.remove ( errorFile )
                if     os.path.exists ( outFile    ): os.remove ( outFile   )

                # . Check that the script file exists.
                if os.path.exists ( scriptFile ):

                    # . Run the script.
                    isOK  = True
                    time0 = time.time ( )
                    eFD   = open ( errorFile, "w" )
                    oFD   = open ( outFile  , "w" )
                    try:
                        process = subprocess.Popen ( [ _PythonCommand, scriptFile ], stderr = eFD, stdout = oFD )
                        process.wait ( )
                    except Exception as e:
                        eFD.write ( e )
                        isOK = False
                    scriptTime = time.time ( ) - time0
                    tutorialTime  += scriptTime

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

                # . Script file not found.
                else:
                    isOK = False

                # . Get the result.
                if isOK:
                    numberOfSuccesses += 1
                    resultLabel = "Pass"
                else:
                    resultLabel = "Error"

                # . Printing.
                print ( "{:s}{:10s}{:>20s}".format ( scriptPath.ljust ( _LabelWidth ), resultLabel, CPUTime.TimeToString ( scriptTime ) ) )

            # . Print a terminating message.
            if len ( scripts ) > 1:
                resultLabel = "{:d}/{:d}".format ( numberOfSuccesses, len ( scripts ) )
                print ( "{:s}{:10s}{:>20s}".format ( "** Totals **".ljust ( _LabelWidth ), resultLabel, CPUTime.TimeToString ( tutorialTime ) ) )

            # . Finish up.
            totalTime += tutorialTime

        # . Final totals.
        print ( "\n" + 80 * "-" )
        if len ( options.found ) > 1:
            print ( "\nTotal time used: " + CPUTime.TimeToString ( totalTime, compact = False ) + "." )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    # . Run the script.
    RunTutorials ( )
