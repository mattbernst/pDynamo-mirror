"""Recursive enumeration.

An implementation of REA from:

V. M. Jimenez and A. Marzal
"Computing the K Shortest Paths: A New Algorithm and an Experimental Comparison"
WAE 1999, LNCS 1668: 15-29, 1999

"""

# . Brute force computation of all paths too expensive.

# . Note that index only needs to be defined for paths in bestPaths.

import operator, sys

from collections import defaultdict

from BellmanFord import BellmanFordShortestPaths
from GraphStatus import GraphError
from PathHead    import PathHead

doPrint = True
def PrintPath ( path, k, tag, n ):
    tags = []
    for node in path.Nodes ( ): tags.append ( node.tag )
    print ( tag, k, tags, round ( path.weight, 3 ), n )

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _NextPath ( graph, source, node, k, bestPaths, candidatePaths ):
    """Calculate the next longest path for a node."""

    # . Initialization.
    bestPath = None

    # . Step B1.
    # . Initialize a set of candidate paths.
    if k == 1:
        lastBestPath = bestPaths[node][0]
        candidates   = candidatePaths[node]
        for edge in graph.adjacentEdges.get ( node, [] ): # . Should be inEdges.
            predecessor = edge.Opposite ( node )
            newPath = PathHead ( headEdge = edge, headNode = node, tail = bestPaths[predecessor][0] )
            if ( not newPath.IsIdenticalTo ( lastBestPath ) ) and newPath.IsSimplePath ( ):
                include = True
                for candidate in candidates:
                    if candidate == newPath:
                        include = False
                        break
                if include:
                    candidates.append ( newPath )
                    if doPrint: PrintPath ( newPath, k, "C1>", len ( candidates ) )

    # . Step B2.
    # . Recursion check.
    if ( k != 1 ) or ( node is not source ):

        # . Step B3.
        lastBestPath = bestPaths[node][k-1]
        predecessor  = lastBestPath.tail.headNode
        newTailIndex = lastBestPath.tail.index + 1

        # . Step B4.
        # . Compute next tail path if it does not already exist.
        tailPaths = bestPaths[predecessor]
        if newTailIndex not in tailPaths:
            newTail = _NextPath ( graph, source, predecessor, newTailIndex, bestPaths, candidatePaths )
            if newTail is not None: tailPaths[newTailIndex] = newTail

        # . Try again.
        newTail = tailPaths.get ( newTailIndex, None )

        # . Step B5.
        if newTail is not None:
            newPath = PathHead ( headEdge = lastBestPath.headEdge, headNode = node, tail = newTail )
            if ( not newPath.IsIdenticalTo ( lastBestPath ) ) and newPath.IsSimplePath ( ):
                include = True
                for candidate in candidatePaths[node]:
                    if candidate.IsIdenticalTo ( newPath ):
                        include = False
                        break
                if include:
                    candidatePaths[node].append ( newPath )
                    if doPrint: PrintPath ( newPath, k, "Cn>", len ( candidatePaths[node] ) )

    # . Step B6.
    # . Get the best path.
    candidates = candidatePaths[node]
    if len ( candidates ) > 0:
        candidates.sort ( key = operator.attrgetter ( "weight" ) )
        bestPath = candidates.pop ( 0 )
        bestPath.index = k
        if doPrint: PrintPath ( bestPath, k, "BP>", len ( candidates ) )

    # . Finish up.
    if bestPath is None and doPrint:
        print ( "NO BEST PATH>", node.tag, k )
    return bestPath

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def RecursiveEnumeration ( graph, source, target, maximumPaths = 100, maximumWeight = None ):
    """Compute the K shortest paths."""

    # . Preliminaries.
    if maximumWeight is None: maximumWeight = sys.float_info.max
    else:                     maximumWeight = max ( maximumWeight, 0.0 )

    # . Best paths and candidate paths for the nodes.
    bestPaths      = defaultdict ( dict )
    candidatePaths = defaultdict ( list )
    finalPaths     = []

    # . Step A1.
    # . All shortest paths - level 0.
    paths = BellmanFordShortestPaths ( graph, source )
    finalPaths.append ( paths[target] )
    for ( node, path ) in paths.iteritems ( ):
        path.index = 0
        bestPaths[node][0] = path

    # . Step A2.
    # . Loop over the number of paths required - levels 1 and up.
    for k in range ( 1, maximumPaths ):
        newPath = _NextPath ( graph, source, target, k, bestPaths, candidatePaths )
        if ( newPath is None ) or ( newPath.weight > maximumWeight ):
            if newPath is None: print ( "\nNO NEW PATH FOUND." )
            break
        bestPaths[target][k] = newPath
        finalPaths.append ( newPath )

    # . Finish up.
    return finalPaths

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
