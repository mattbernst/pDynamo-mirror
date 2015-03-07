"""Yen's algorithm for finding the K shortest path between two nodes in a graph."""

# . Yen's original algorithm more straightforward. Use instead?
# . Many paths could be generated. Too many?
# . Use weight restriction?

"""Dijkstra shortest paths."""

import heapq, sys # . heapq is a priority queue.

from GraphStatus import GraphError
from PathHead    import PathHead

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def DijkstraShortestPaths ( graph, source, reverse = False, target = None ):
    """Returns shortest paths in a weighted graph."""

    # . Check edge weights.
    if graph.MinimumEdgeWeight ( ) < 0: raise GraphError ( "Graph has negative edge weights." )

    # . Edges to use.
    if reverse: outEdges = graph.adjacentEdges # inEdges
    else:       outEdges = graph.adjacentEdges # outEdges

    # . Initialization.
    predecessorEdges = { source : None }
    predecessorNodes = { source : None }
    weights          = {}
    priorityQueue    = []
    seen             = { source : 0.0 }
    heapq.heappush ( priorityQueue, ( 0.0, source ) )

    # . Loop until the queue is empty.
    while len ( priorityQueue ) > 0:
        ( weight, node ) = heapq.heappop ( priorityQueue )
        if node in weights: continue
        weights[node]    = weight
        if node is target: break
        for edge in outEdges.get ( node, [] ):
	    other     = edge.Opposite ( node )
            newWeight = weight + edge.weight
            if other in weights:
                if newWeight < weights[other]:
                    raise GraphError ( "Logic error: shorter path found." )
            elif ( other not in seen ) or ( newWeight < seen[other] ):
                heapq.heappush ( priorityQueue, ( newWeight, other ) )
                predecessorEdges[other] = edge
                predecessorNodes[other] = node
                seen            [other] = newWeight

    # . Build paths.
# . Make prettier.
    paths = {}
    for node in weights:
        paths[node] = PathHead ( headNode = node )
    for node in weights:
        paths[node].UpdateTail ( headEdge = predecessorEdges[node], tail = paths.get ( predecessorNodes[node], None ) )
    for node in weights:
        paths[node].CalculateWeight ( )

    # . Finish up.
    if target is None: result = paths
    else:              result = paths[target]
    return result

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def AllShortestPaths ( graph, source, target, maximumPaths = 100, maximumWeight = None ):
    """Return all shortest paths."""
    if maximumWeight is None: maximumWeight = sys.float_info.max
    iterator = NextShortestPath ( graph, source, target )
    paths    = []
    for path in iterator:
#        print "NEW PATH> ", path, len ( path ), round ( path.weight, 3 ), PathList ( path )
        if ( path is None ) or ( path.weight > maximumWeight ): break
        paths.append ( path )
        if len ( paths ) >= maximumPaths: break
    return paths

#def PathList ( path ):
#    if path is None: return None
#    else:
#        tags = []
#        for node in path.Nodes ( ):
#            tags.append ( node.tag )
#        return tags

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def NextShortestPath ( graph, source, target ):
    """Find the next shortest path."""

    # . Initialization.
    candidateNumber    = 0
    pathCandidates     = []
    deviationNodeIndex = {}
    previousPaths      = []

    # . Checks.
    numberOfEdges = len ( graph.edges )
    numberOfNodes = len ( graph.nodes )
#    print "Initial number of edges and nodes = ", numberOfEdges, numberOfNodes

    # . Get the shortest path.
    shortestPath = DijkstraShortestPaths ( graph, source, target = target )
    if shortestPath is not None:
        deviationNodeIndex[shortestPath] = source
        heapq.heappush ( pathCandidates, ( shortestPath.weight, candidateNumber, shortestPath ) )
        candidateNumber += 1

    # . Loop until there are no more path candidates.
    while len ( pathCandidates ) > 0:

        # . Get the current path.
        ( weight, number, currentPath ) = heapq.heappop ( pathCandidates )
        yield currentPath

#        print "CurrentPath> ", PathList ( currentPath ), deviationNodeIndex[currentPath].tag

        # . Use the current path to generate more path candidates.

        # . Get edges and nodes in path in forward order and without terminating node.
        pathEdges = currentPath.Edges ( )
        pathNodes = currentPath.Nodes ( )
        pathNodes.pop ( )

        # . Remove nodes and edges along the path before the current deviation node.
	currentDeviationNode = deviationNodeIndex[currentPath]
        edgesToRemove        = set ( )
        nodesToRemove        = set ( )
        while True:
            if pathNodes[0] is currentDeviationNode: break
            node = pathNodes.pop ( 0 )
            edge = pathEdges.pop ( 0 )
            edgesToRemove.update ( graph.RemoveNode ( node ) )
            nodesToRemove.add ( node )

        # . Remove edges from the current deviation node which coincide with previous paths with identical tails.
        ( node, edge, currentPathTail ) = currentPath.Split ( currentDeviationNode )
        for previousPath in previousPaths:
            ( node, edge, previousPathTail ) = previousPath.Split ( currentDeviationNode )
            if previousPathTail is not None:
                if ( currentPathTail == previousPathTail ) and ( edge not in edgesToRemove ):
                    edgesToRemove.add ( edge )
                    graph.RemoveEdge  ( edge )
        previousPaths.append ( currentPath )

        # . Define the remaining nodes and edges along the current path that are to be removed and restored.
        elementsToRestore = []
        for ( node, edge ) in zip ( pathNodes, pathEdges ):
            graph.RemoveEdge ( edge )
            otherEdges = graph.RemoveNode ( node )
            elementsToRestore.append ( ( node, edge, otherEdges ) )
        elementsToRestore.reverse ( )

        # . Get the reverse shortest distances from the target to all other nodes in the graph.
        reversePaths = DijkstraShortestPaths ( graph, target, reverse = True )

        # . Restore the relevant nodes and edges.
        successorNode = target
        for ( currentNode, currentEdge, otherEdges ) in elementsToRestore:

            if currentEdge.Opposite ( currentNode ) is not successorNode: raise GraphError ( "Edge mismatch." )

            # . Restore the elements.
	    graph.AddNode  ( currentNode )
            graph.AddEdges ( otherEdges  )

            # . Update the path to the node using forward star form.
            subPath = _UpdateForwardPath ( graph, reversePaths, currentNode )
#            print "SubPath> ", PathList ( subPath )

	    # . Build the new path if possible.
	    if subPath is not None:

                # . Correct paths.
                _CorrectBackwardPaths ( graph, reversePaths, currentNode )


# . Make prettier? Use split first?
                newPath = PathHead ( headNode = source )
                for edge in currentPath.Edges ( ):
                    headNode = newPath.headNode
                    if headNode is currentNode: break
                    newPath = newPath.PushHead ( edge )
# . Use MakeReversed.
# newPath = subPath.MakeReversed ( )
# newPath???
                subPathEdges = subPath.Edges ( )
                subPathEdges.reverse ( )
                for edge in subPathEdges:
                    newPath = newPath.PushHead ( edge )

		# . Save the path if it is new.
# . Is this necessary and, if so, how check?
# . Some hash key for path?
# . Explicit check against all paths with same weight?
                if True:
                    heapq.heappush ( pathCandidates, ( newPath.weight, candidateNumber, newPath ) )
                    candidateNumber += 1
		    deviationNodeIndex[newPath] = currentNode

	    # . Restore the current edge.
	    graph.AddEdge ( currentEdge )

	    # . Update the reverse paths if necessary.
            currentReversePath = reversePaths.get ( currentNode, None )
            successorPath      = reversePaths[successorNode]
            if currentReversePath is None:
                currentReversePath = PathHead ( headEdge = currentEdge, headNode = currentNode, tail = successorPath )
                reversePaths[currentNode] = currentReversePath
            else:
   	        weight = currentEdge.weight + successorPath.weight
	        if currentReversePath.weight > weight:
                    currentReversePath.UpdateTail ( headEdge = currentEdge, tail = successorPath )
		    _CorrectBackwardPaths ( graph, reversePaths, currentNode )

            # . Reset successor node.
            successorNode = currentNode

	# . Restore all nodes and edges to the graph.
	graph.AddNodes ( nodesToRemove )
	graph.AddEdges ( edgesToRemove )
        if ( len ( graph.edges ) != numberOfEdges ) or ( len ( graph.nodes ) != numberOfNodes ):
#            print "Edges = ", len ( graph.edges ), numberOfEdges
#            print "Nodes = ", len ( graph.nodes ), numberOfNodes
            raise GraphError ( "Restored graph has invalid number of edges or nodes." )

    # . Finish up.
#    return None

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _CorrectBackwardPaths ( graph, paths, node ):
    """Given a node update all shortest paths that may include it."""

    # . Loop over successor nodes until 
    # . Successors are actually predecessors because the paths are reversed.
    queue = [ node ]
    while len ( queue ) > 0:
        currentNode   = queue.pop ( 0 )
        currentPath   = paths[currentNode]
        currentWeight = currentPath.weight
        for edge in graph.adjacentEdges.get ( currentNode, [] ): # . inEdges.
            successorNode   = edge.Opposite ( currentNode )
            successorPath   = paths.get ( successorNode, None )
            if successorPath is None:
                successorPath = PathHead ( headEdge = edge, headNode = successorNode, tail = currentPath )
                paths[successorNode] = successorPath
            else:
                successorWeight = successorPath.weight
                newWeight       = currentWeight + edge.weight
                if successorWeight > newWeight:
                    successorPath.UpdateTail ( headEdge = edge, tail = currentPath ) # . Why should successorNode not be on currentPath?
                    queue.append ( successorNode )

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _UpdateForwardPath ( graph, paths, node ):
    """Update the path to node from target and return it."""

    # . Get the path for the input node.
    path = paths.get ( node, None )

    # . Loop over predecessor edges to find a lower weight path.
    # . Predecessors are actually successors because the paths are reversed.
    for edge in graph.adjacentEdges.get ( node, [] ): # . outEdges.
        predecessorNode = edge.Opposite ( node )
        predecessorPath = paths.get ( predecessorNode, None )
        if predecessorPath is not None:
            predecessorWeight = predecessorPath.weight
            newWeight         = predecessorWeight + edge.weight
            if path is None:
                path = PathHead ( headEdge = edge, headNode = node, tail = predecessorPath )
            elif path.weight > newWeight:
                path.UpdateTail ( headEdge = edge, tail = predecessorPath )

    # . Update paths and return.
    if path is not None: paths[node] = path
    return path

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass




