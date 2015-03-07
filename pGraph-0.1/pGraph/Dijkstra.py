"""Solve the single source shortest path problem for a graph with non-negative edge weights."""

import heapq # . A priority queue.

from GraphStatus import GraphError

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def DijkstraSingleSource ( graph, source, cutoff = None, target = None ):
    """Returns shortest paths and lengths in a weighted graph."""

    # . Initialization.
    distances = { source : 0 }
    paths     = { source : [ source ] }

    # . Do nothing if there are no edges or the source and target are identical.
    if ( len ( graph.edges ) > 0 ) and ( source is not target ):

        # . Check edge weights.
	if graph.MinimumEdgeWeight ( ) < 0: raise GraphError ( "Graph has negative edge weights." )

        # . Initialization.
	distances = {} # . Needed to avoid immediate exit.
	queue     = []
	seen      = { source : 0 }
        heapq.heappush ( queue, ( 0, source ) )

        # . Loop until the queue is empty.
	while len ( queue ) > 0:
            ( distance, node ) = heapq.heappop ( queue )
            if node in distances: continue
            distances[node] = distance
            if node is target: break
            for edge in graph.adjacentEdges.get ( node, [] ):
	        other    = edge.Opposite ( node )
        	distance = distances[node] + edge.weight
        	if cutoff is not None:
                    if distance > cutoff: continue
        	if other in distances:
                    if distance < distances[other]:
                	raise GraphError ( "Logic error: shorter path found." )
        	elif ( other not in seen ) or ( distance < seen[other] ):
                    heapq.heappush ( queue, ( distance, other ) )
                    paths[other] = paths[node] + [other]
                    seen [other] = distance

    # . Finish up.
    return ( distances, paths )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
