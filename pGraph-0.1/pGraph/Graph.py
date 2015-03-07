"""Graph class."""

from GraphStatus import GraphError
from Path        import Path

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Graph ( object ):
    """A basic graph class."""
 
    defaultAttributes = { "adjacentEdges" : None,
                          "adjacentNodes" : None,
                          "edges"         : None,
                          "nodes"         : None,
                          "nodeIndex"     : None } # . Needed? Use id ( node ) instead?

# . Currently nodeIndex won't work if nodes are removed and then others are readded to the graph.

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        self.Clear ( )

    def AddEdge ( self, edge ):
        """Add an edge to the graph."""
	iNodes = self.adjacentNodes.get ( edge.node1, set ( ) )
	jNodes = self.adjacentNodes.get ( edge.node2, set ( ) )
	if ( edge.node1 in jNodes ) and ( edge.node2 in iNodes ):
	    raise ValueError ( "Edge already present in graph." )
	else:
   	    self.edges.append ( edge )
            # . Add to adjacent nodes.
	    iNodes.add ( edge.node2 ) ; self.adjacentNodes[edge.node1] = iNodes
	    jNodes.add ( edge.node1 ) ; self.adjacentNodes[edge.node2] = jNodes
	    # . Add to adjacent edges.
    	    iEdges = self.adjacentEdges.get ( edge.node1, set ( ) )
	    jEdges = self.adjacentEdges.get ( edge.node2, set ( ) )
	    iEdges.add ( edge ) ; self.adjacentEdges[edge.node1] = iEdges
	    jEdges.add ( edge ) ; self.adjacentEdges[edge.node2] = jEdges

    def AddEdges ( self, edges ):
        """Add edges to the graph."""
        for edge in edges: self.AddEdge ( edge )

    def AddNode ( self, node ):
        """Add a node to the graph."""
	if node in self.adjacentNodes:
	    raise GraphError ( "Node already present in graph." )
	else:
	    self.nodeIndex[node] = len ( self.nodes )
    	    self.nodes.append ( node )
#            print "Added number of nodes = ", len ( self.nodes ), node.tag

    def AddNodes ( self, nodes ):
        """Add nodes to the graph."""
        for node in nodes: self.AddNode ( node )

    def Clear ( self ):
        """Clear all edges and nodes."""
        self.ClearEdges ( )
	self.nodes     = []
	self.nodeIndex = {}

    def ClearEdges ( self ):
        """Clear all edges."""
        self.adjacentEdges = {}
	self.adjacentNodes = {}
	self.edges         = []

    def MakePathFromNodes ( self, nodes ):
        """Generate a path from a sequence of nodes."""
        path = Path ( )
        for node in nodes:
            path.AddNode ( node )
        if len ( nodes ) > 1:
            tail = nodes[0]
            for head in nodes[1:]:
                found = False
                for edge in self.adjacentEdges[tail]:
                    if ( ( edge.node1 is tail ) and ( edge.node2 is head ) ) or \
                       ( ( edge.node2 is tail ) and ( edge.node1 is head ) ) :
                       found = True
                       path.AddEdge ( edge )
                       break
                if not found: raise GraphError ( "Unable to find edge between two putative path nodes." )
                tail = head
        return path

    def MakeSubgraph ( self, nodes, induced = True ):
        """Generate a subgraph with the selected nodes."""
        pruned = self.__class__ ( )
	if len ( nodes ) > 0:
	    for node in nodes:
		if node in self.nodes: pruned.AddNode ( node )
		else: raise GraphError ( "Node not present in graph." )
	    if induced:
    	        for edge in self.edges:
		    if ( edge.node1 in nodes ) and ( edge.node2 in nodes ):
                        pruned.AddEdge ( edge )
        return pruned

    def MinimumEdgeWeight ( self ):
        """Return the minimum edge weight."""
        if len ( self.edges ) > 0:
            minimumEdgeWeight = self.edges[0].weight
            for edge in self.edges[1:]:
                minimumEdgeWeight = min ( minimumEdgeWeight, edge.weight )
        else:
            minimumEdgeWeight = 0.0
        return minimumEdgeWeight

    def RemoveEdge ( self, edge ):
        """Remove an edge from the graph."""
#        print "Removing Edge> ", edge.node1.tag, edge.node2.tag
        try:    self.edges.remove ( edge )
        except: raise GraphError ( "Edge not in graph." )
        iNode = edge.node1
        jNode = edge.node2
        self.adjacentEdges[iNode].remove ( edge )
        self.adjacentEdges[jNode].remove ( edge )
# . This won't work if the graph is a multigraph.
# . if self.IsMultiGraph ( ):
#        adjacentNodes = set ( )
#        for edge in self.adjacentEdges[iNode]:
#            adjacentNodes.add ( edge.Opposite ( iNode ) )
#        self.adjacentNodes[iNode] = adjacentNodes
# . and repeat for jNode.
# . else:
        self.adjacentNodes[iNode].remove ( jNode )
        self.adjacentNodes[jNode].remove ( iNode )

    def RemoveEdges ( self, edges ):
        """Remove edges from the graph."""
        for edge in edges: self.RemoveEdge ( edge )

    def RemoveNode ( self, node ):
        """Remove a node from the graph."""
        try:    self.nodes.remove ( node )
        except: raise GraphError ( "Node not in graph." )
        del self.nodeIndex[node]
        edges = self.adjacentEdges.pop ( node, set ( ) )
        for edge in edges:
            self.edges.remove ( edge )
            other = edge.Opposite ( node )
            self.adjacentEdges[other].remove ( edge )
        for other in self.adjacentNodes.pop ( node, set ( ) ):
            self.adjacentNodes[other].remove ( node )
#        print "Removed number of nodes = ", len ( self.nodes ), node.tag
        return edges

    def RemoveNodes ( self, nodes ):
        """Remove nodes from the graph."""
        edges = set ( )
        for node in nodes:
            edges.update ( self.RemoveNode ( node ) )
        return edges

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

