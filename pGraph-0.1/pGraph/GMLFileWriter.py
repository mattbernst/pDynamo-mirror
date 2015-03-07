"""Simple graph tests."""

import glob, math, os

from BiconnectedComponents import *
from ConnectedComponents   import *
from Vismara               import *

#from Dijkstra import *

from Edge  import *
from Graph import *
from Node  import *

def WriteGML ( graph, path ):
    gfile = open ( path, "w" )
    gfile.write (  "graph [\n" )
    for node in graph.nodes:
        gfile.write (  "node [\n" )
        gfile.write (  "id {:d}\n".format ( node.index ) )
        gfile.write (  "]\n" )
    for edge in graph.edges:
        gfile.write (  "edge [\n" )
        gfile.write (  "source {:d}\n".format ( edge.node1.index ) )
        gfile.write (  "target {:d}\n".format ( edge.node2.index ) )
        gfile.write (  "]\n" )
    gfile.write (  "]\n" )
    gfile.close ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    paths = ( ( 4, ( ( 0, 1 ), ( 0, 2 ), ( 0, 3 ), ( 1, 2 ), ( 1, 3 ), ( 2, 3 ) ) ),
              ( 6, ( ( 0, 1 ), ( 0, 2 ), ( 0, 3 ), ( 0, 4 ), ( 1, 2 ), ( 1, 4 ), ( 1, 5 ), ( 2, 3 ), ( 2, 5 ), ( 3, 4 ), ( 3, 5 ), ( 4, 5 ) ) ) )

    X = 0
    for ( numberNodes, edges ) in paths:

        # . Create graph.
        graph = Graph ( )
        nodes = {}
        for i in range ( numberNodes ):
            node = Node ( )
	    node.index = i
	    nodes[i] = node
	    graph.AddNode ( node )
        for ( i, j ) in edges:
            graph.AddEdge ( Edge ( nodes[i], nodes[j] ) )

        WriteGML ( graph, "g{:d}.gml".format ( X ) )
        X += 1
