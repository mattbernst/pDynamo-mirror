"""Testing for point groups."""

import glob, math, os

from pBabel           import XYZFile_ToSystem
from pCore            import logFile, LogFileActive, TestCase, TextLogFileWriter
from pMoleculeScripts import FindSystemPointGroup

# . This test is quite rapid apart from the icosahedral point groups (I, Ih)
# . which take > 97% of the time.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_doPrinting = True

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PointGroupTest ( TestCase ):
    """A point group test case."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( PointGroupTest, self ).__init__ ( *args )
        self.doLong = True

    def MakeShort ( self ):
        self.doLong = False
        return True

    def runTest ( self ):
        """The test."""

        # . Initialization.
        log = self.GetLog ( )
        if not _doPrinting: log = None

        # . Paths.
        paths = glob.glob ( os.path.join ( os.getenv ( "PDYNAMO_ROOT" ), "molecularStructures", "pointGroupExamples", "xyz", "*.xyz" ) )
        paths.sort ( )

        # . Calculation.
        results   = []
        successes = 0
        for path in paths:

            ( head, tail ) = os.path.split ( path )
            inGroup = tail[0:-4].split ( "_" )[0]

            # . Skip icosahedral cases for short jobs.
            if ( not self.doLong ) and ( ( inGroup == "I" ) or ( tail[0:-4] == "Ih_c" ) ): continue

            system = XYZFile_ToSystem ( path )
            system.Summary ( log = log )

            outGroup = FindSystemPointGroup ( system, log = log )

            if outGroup is None:
                if LogFileActive ( log ):
                    log.Text ( "\n\nNo point group found for " + system.label + "\n" )
                outGroup = "-"
            else:
                if LogFileActive ( log ):
                    log.Text ( "\n\nPoint group " + outGroup + " for " + system.label + "\n" )
                if inGroup == outGroup: successes += 1

            results.append ( ( inGroup, outGroup, system.label ) )

        # . Printing.
        if LogFileActive ( log ):
            log.Text ( "\nSummary of Results (OK, expected, obtained, system):\n\n" )
            for ( inGroup, outGroup, system.label ) in results:
                if inGroup == outGroup: tag = ""
                else:                   tag = "*"
                log.Text ( "{:1s}   {:<10s}  {:<10s}  {:s}\n".format ( tag, inGroup, outGroup, system.label ) )

            log.Text ( "\nNumber of successes = {:d} (out of {:d}).\n".format ( successes, len ( results ) ) )

            if successes == len ( results ): log.Text ( "\nTest Succeeded.\n" )
            else:                            log.Text ( "\nTest Failed.\n"    )

            # . Close up.
            if self.outputPath is not None: log.Close ( )

        # . Success/failure.
        self.assertEqual ( successes, len ( results ), msg = "Number of point group matches = {:d} (/{:d})\n".format ( successes, len ( results ) ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = PointGroupTest ( )
    test.run ( )
