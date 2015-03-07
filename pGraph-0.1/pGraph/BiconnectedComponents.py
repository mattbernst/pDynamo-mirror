"""Calculate the biconnected components of an undirected graph."""

# . Test for recursion.
# . This could cause a problem on big graphs.
#maxRecursion = 0
#recursion    = 0

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def BiconnectedComponents ( graph ):
    """Calculate the biconnected components of an undirected graph."""
    # . Find the biconnected components.
#    global maxRecursion, recursion
#    maxRecursion = 0
#    recursion    = 0
    components  = []
    count1      =  0
    count2      =  0
    current     = []
    depth       = {}
    low         = {}
    predecessor = {}
    for node in graph.nodes:
        if node not in depth:
	    count1 += 1
	    depth[node] = count1
	    current.append ( node )
	    ( count1, count2 ) = _DepthFirstSearch ( graph, node, components, depth, low, predecessor, current, count1, count2 )
#    print "Maximum recursion depth = ", maxRecursion
    # . Sort into components - largest to smallest.
    components.sort ( key = len, reverse = True )
    return components

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _DepthFirstSearch ( graph, node, components, depth, low, predecessor, current, count1, count2 ):
    """Depth first search."""
#    global maxRecursion, recursion
#    recursion   += 1
#    maxRecursion = max ( maxRecursion, recursion )
    low[node] = depth[node]
    for opposite in graph.adjacentNodes.get ( node, [] ):
        if opposite not in depth:
            count1 += 1
	    depth[opposite] = count1
	    current.append ( opposite )
	    predecessor[opposite] = node
	    ( count1, count2 ) = _DepthFirstSearch ( graph, opposite, components, depth, low, predecessor, current, count1, count2 )
            low[node] = min ( low[node], low[opposite] )
        else:
	    low[node] = min ( low[node], depth[opposite] )
    if ( node in predecessor ) and ( low[node] == depth[predecessor[node]] ):
	nodes = set ( )
        while True:
            next = current.pop ( -1 )
	    for opposite in graph.adjacentNodes.get ( next, [] ):
	        if depth[next] > depth[opposite]:
		    nodes.add ( next     )
		    nodes.add ( opposite )
	    if next is node: break
	if len ( nodes ) > 2:
	    components.append ( nodes )
            count2 += 1
#    recursion -= 1
    return ( count1, count2 )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

