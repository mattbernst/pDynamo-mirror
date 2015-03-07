"""Graph edge."""

from GraphStatus import GraphError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Edge ( object ):
    """A graph edge."""
 
    defaultAttributes = { "node1"  : None ,
                          "node2"  : None ,
                          "weight" :  1.0 }

    def __init__ ( self, node1, node2, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in keywordArguments.iteritems ( ):
            if key in self.__class__.defaultAttributes:
                setattr ( self, key, value )
	self.node1 = node1
	self.node2 = node2

    def Opposite ( self, node ):
        """Return the other node for the edge."""
	if   node is self.node1: return self.node2
	elif node is self.node2: return self.node1
	else: raise GraphError ( "Node not associated with edge." )

    def SourceTargetPairs ( self ):
        """Return a tuple of source/target pairs."""
        return ( ( self.node1, self.node2 ), ( self.node2, self.node1 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
