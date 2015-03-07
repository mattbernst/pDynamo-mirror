"""Calculate the connected components of an undirected graph."""

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def ConnectedComponents ( graph ):
    """Calculate the connected components of an undirected graph."""
    # . Assign component indices to each node.
    count   = 0
    indices = {}
    for node in graph.nodes:
        if node not in indices:
            toScan = [ node ]
            indices[node] = count
	    while len ( toScan ) > 0:
                next = toScan.pop ( -1 )
		for opposite in graph.adjacentNodes.get ( next, [] ):
                    if opposite not in indices:
		        indices[opposite] = count
			toScan.append ( opposite )
            count += 1
    # . Sort into components - largest to smallest.
    components = [ set ( ) for i in range ( count ) ]
    for node in graph.nodes:
        components[indices[node]].add ( node )
    components.sort ( key = len, reverse = True )
    return components

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
