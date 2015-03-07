"""Testing for relevant cycles."""

import glob, math, os

from pCore  import LogFileActive, TestCase, TextLogFileWriter
from pGraph import *

from BiconnectedGraphExamples import connectivities, relevantCycles

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
testFigueras = True
testRelevant = True

if testFigueras:
    from pCore import FindRings

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class RelevantCycleTest ( TestCase ):
    """A relevant cycle test case."""

    def runTest ( self ):
        """The test."""

        # . Output setup.
        log = self.GetLog ( )

        # . Initialization.
        totalTests       = len ( connectivities )
        vismaraSuccesses = 0
        vismaraTests     = 0

        # . Loop over examples.
        keys = connectivities.keys ( )
        keys.sort ( )
        for key in keys:
            connectivity = connectivities[key]

            # . Create graph.
            graph = Graph ( )

            # . Nodes.
            nodes = {}
            for i in range ( len ( connectivity ) ):
                node = Node ( )
	        nodes[i] = node
	        graph.AddNode ( node )

            # . Edges.
            for ( i, connections ) in enumerate ( connectivity ):
                for j in connections:
                    if j > i: graph.AddEdge ( Edge ( nodes[i], nodes[j] ) )

            # . Various graph statistics.
            connectedComponents   =   ConnectedComponents ( graph )
            biconnectedComponents = BiconnectedComponents ( graph )

            # . Output.
            if LogFileActive ( log ):
                log.Text ( "\nGraph Example - " + key + ":\n" )
                log.Text ( "Number of nodes        = {:d}\n".format ( len ( graph.nodes           ) ) )
                log.Text ( "Number of edges        = {:d}\n".format ( len ( graph.edges           ) ) )
                log.Text ( "Connected components   = {:d}\n".format ( len ( connectedComponents   ) ) )
                log.Text ( "Biconnected components = {:d}\n".format ( len ( biconnectedComponents ) ) )
                log.Text ( "Cyclomatic number      = {:d}\n".format ( len ( graph.edges ) - len ( graph.nodes ) + len ( connectedComponents ) ) )

            # . Figueras testing.
            if testFigueras:
                table = {}
                for ( i, connections ) in enumerate ( connectivity ):
                    if len ( connections ) > 0: table[i] = connections
                SSSR = FindRings ( table )
                if len ( SSSR ) > 0:
                    SSSR.sort ( key = len )
                    lengths = []
                    for ring in SSSR: lengths.append ( "{:d}".format ( len ( ring ) ) )
                    if LogFileActive ( log ):
                        log.Text ( "Figueras cycle lengths =  " + " ".join ( lengths ) + "\n" )

            # . Cycles.
            if testRelevant:
                cycles = VismaraRelevantCycles ( graph, biconnectedComponents = biconnectedComponents )
                if len ( cycles ) > 0:
                    cycles.sort ( key = len )
                    lengths = []
                    for cycle in cycles: lengths.append ( "{:d}".format ( len ( cycle ) ) )
	            vismaraLengths = " ".join ( lengths )
                    if LogFileActive ( log ):
                        log.Text ( "Vismara cycle lengths  =  {:s} ({:d})\n".format ( " ".join ( lengths ), len ( cycles ) ) )

                    # . Check for success.
	            reference = relevantCycles.get ( key, None )
	            if reference is not None:
	                lengths = []
		        sizes   = reference.keys ( )
		        sizes.sort ( )
		        for size in sizes:
		            lengths.extend ( reference[size] * [ "{:d}".format ( size ) ] )
		        referenceLengths = " ".join ( lengths )
		        if referenceLengths == vismaraLengths: vismaraSuccesses += 1
		        vismaraTests += 1

        # . Final printing.
        log.Text ( "\nTotal tests = {:d}\n".format ( totalTests ) )
        if vismaraTests > 0:
            log.Text ( "\nNumber of successful relevant cycle tests = {:d} (/{:d})\n".format ( vismaraSuccesses, vismaraTests ) )
            if vismaraSuccesses == vismaraTests: log.Text ( "\nTest passed.\n" )
            else:                                log.Text ( "\nTest failed.\n" )

            # . Close up.
            if self.outputPath is not None: log.Close ( )

        # . Success/failure.
        self.assertEqual ( vismaraSuccesses, vismaraTests, msg = "Number of relevant cycle matches = {:d} (/{:d})\n".format ( vismaraSuccesses, vismaraTests ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = RelevantCycleTest ( )
    test.run ( )


