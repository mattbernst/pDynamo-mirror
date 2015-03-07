#-------------------------------------------------------------------------------
# . File      : RingFinder.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Ring-finding module using Figueras's SSSR algorithm.

This module finds the smallest set of smallest rings (SSSR) for a given
network based upon the algorithm of J. Figueras "Ring Perception Using
Breadth-First Search" J. Chem. Inf. Comput. Sci. 36, 986-991 (1996).

The module works for 2-D networks (either in a plane or a surface in 3-D)
and for networks that have at least some nodes with connectivities of 3
or less.

The connectivity table that is input to FindSSSR must be valid which means
that it consists of a list of N nodes (numbered from 0 to N-1) with no
nodes of other index appearing in the table. All connections are specified
twice (i.e. N1-N2 is specified as N2 in the table for N1 and vice-versa).

The number of rings should equal ( Nring_bonds - Nring_atoms + 1 ) for each
disconnected set of rings.
"""

# . This module definitely needs replacing although it works normally for simple graphs with few rings.

import bisect

from Clone import Clone

#===============================================================================
# . Find SSSR Rings.
#===============================================================================
def FindRings ( connectivitytable ):
    """Find the SSSR for a network."""

    # . Initialize SSSR.
    SSSR = [ ]

    # . Save the connectivity and sort each list.
    connectivity = Clone ( connectivitytable )
    for connections in connectivity.values ( ): connections.sort ( )

    # . Loop until all nodes have been eliminated.
    while len ( connectivity ) > 0:

        # . Create information about the current state of the connectivity.
        currentstate = []
        for ( node, connections ) in connectivity.iteritems ( ):
            currentstate.append ( ( len ( connections ), node ) )
        currentstate.sort ( ) # . Smallest to largest degree.

        # . Automatically remove all connectivity 0 nodes but process remaining nodes one-by-one.
        while len ( currentstate ) > 0:

            # . Get the lowest connectivity node.
            ( degree0, node0 ) = currentstate.pop ( 0 )

            # . Connectivity 0 - remove the node.
            if degree0 == 0: del connectivity[node0]

            # . Higher connectivities.
            else:

                # . Connectivity 1 - remove connections to node0.
                if degree0 == 1: _Trim ( node0, connectivity )

                # . Connectivity 2.
                elif degree0 == 2:

                    # . Identify a single node in each separate chain of nodes of order 2.
                    currentorder2 = _IdentifyChains ( degree0, node0, currentstate, connectivity )

                    # . Find rings for each chain of nodes of order 2.
                    for node2 in currentorder2:

                        # . Get the ring for node2.
                        ring = _GetRing ( node2, connectivity )

                        # . Save the ring if it is not already present.
                        if len ( ring ) > 0:
                            if ring not in SSSR: SSSR.append ( ring )

                    # . Break a bond for each chain of nodes of order 2.
                    for node2 in currentorder2:

                        # . Get the first connection to node2 (it does not matter which).
                        node1 = connectivity[node2][0]

                        # . Remove this connection from the table.
                        connectivity[node1].remove ( node2 )
                        connectivity[node2].remove ( node1 )

                # . Connectivity 3 or higher (although there is no guarantee that it will work for connectivities higher than 3).
                else:

                    # . Get the ring for node0.
                    ring = _GetRing ( node0, connectivity )

                    # . Save the ring if it is not already present.
                    if len ( ring ) > 0:
                        if ring not in SSSR: SSSR.append ( ring )

                    # . Eliminate an edge.
                    connectivity = _CheckEdges ( ring, connectivity )

                # . Go back.
                break

    # . Return the set of SSSR.
    return SSSR

#===============================================================================
# . Find ring sets.
# . A set of rings consists of rings that are joined together.
#
# . In principle this function needs to be modified to exclude cases where
# . rings are connected by single connections but do not form part of
# . super-rings. E.g. Ph-Ph either by itself or with closing ring sets at
# . either end.
#===============================================================================
def FindRingSets ( connectivitytable ):
    """Find the SSSR for a network but divided into ring sets."""
    # . Find the rings.
    SSSR = FindRings ( connectivitytable )
    # . Divide the ring into sets.
    nrings   = len ( SSSR )
    ringsets = []
    if nrings > 0:
#        for ring in SSSR:
#            ringsets.append ( [ ring ] )
        # . Define which nodes belong to which rings.
        noderings = {}
        for ( iring, ring ) in enumerate ( SSSR ):
            for node in ring:
                if node not in noderings: noderings[node] = set ( )
                noderings[node].add ( iring )
        # . Get the ring connections.
        ringconnections = [ set ( ) for i in range ( nrings ) ]
        for ( iring, ring ) in enumerate ( SSSR ):
            for node in ring:
                for other in connectivitytable[node]:
                    if other in noderings: ringconnections[iring].update ( noderings[other] )
        for iring in range ( nrings ): ringconnections[iring].discard ( iring )
        # . Initialization.
        QASSIGNED = [ False for i in range ( nrings ) ]
        # . Loop over unassigned rings.
        for sring in range ( nrings ):
            if not QASSIGNED[sring]:
                # . Start the new isolate.
                QASSIGNED[sring] = True
                new = [ sring ]
                # . Flag all rings in the new isolate.
                for iring in new:
                    for jring in ringconnections[iring]:
                        if not QASSIGNED[jring]:
                            QASSIGNED[jring] = True
                            new.append ( jring )
                # . Save the ringset.
                ringset = [ SSSR[i] for i in new ]
                ringsets.append ( ringset )
    return ringsets

#===============================================================================
# . Private Functions.
#===============================================================================
def _CheckEdges ( ring, connectivity ):
    """Select an optimum edge to eliminate from a network with nodes of order 3."""

    # . Initialize the set of edge nodes.
    edges = [ ]

    # . Get the nodes at each edge.
    # . The number of nodes equals the number of edges.
    for edge in range ( len ( ring ) ):
        if edge == len ( ring ) - 1: edges.append ( ( ring[edge], ring[0]      ) )
        else:                        edges.append ( ( ring[edge], ring[edge+1] ) )

    # . Initialize the set of largest rings.
    largestrings = [ ]

    # . Loop over the edges in the ring.
    for ( node1, node2 ) in edges:

        # . Remove the edge from the connectivity table.
        connectivity[node1].remove ( node2 )
        connectivity[node2].remove ( node1 )

        # . Find the rings associated with node1 and node2.
        ring1 = _GetRing ( node1, connectivity )
        ring2 = _GetRing ( node2, connectivity )

        # . Get the lengths of the rings.
        lenring1 = len ( ring1 )
        lenring2 = len ( ring2 )

        # . Save the ring of largest size.
        if lenring1 >= lenring2: largestrings.append ( lenring1 )
        else:                    largestrings.append ( lenring2 )

        # . Put back the edge into the connectivity table (maintain sorted order).
        bisect.insort ( connectivity[node1], node2 )
        bisect.insort ( connectivity[node2], node1 )

    # . Find the index of the ring of smallest size in largestrings.
    smallest = min ( largestrings )
    edge     = largestrings.index ( smallest )

    # . Remove this edge from the connectivity table.
    ( node1, node2 ) = edges[edge]
    connectivity[node1].remove ( node2 )
    connectivity[node2].remove ( node1 )
    if len ( connectivity[node1] ) <= 0: del connectivity[node1]
    if len ( connectivity[node2] ) <= 0: del connectivity[node2]

    # . Return the connectivity table.
    return connectivity

def _GetRing ( node0, connectivity ):
    """Find the smallest ring in the network that contains the node node0."""

    # . Initialize the paths for each node as a list of lists.
    path = {}
    for node in connectivity: path[node] = []

    # . Initialize path for node0.
    path[node0] = [ node0 ]

    # . Initialize the queue with data about node0.
    queue = [ ( node0, -1 ) ]

    # . Loop over the elements in the queue.
    while len ( queue ) > 0:

        # . Remove the first node and its source from the queue.
        ( frontnode, source ) = queue.pop ( 0 )

        # . Loop over the connections of frontnode.
        for m in connectivity[frontnode]:

            # . Avoid frontnode's source.
            if m == source: continue

            # . The path for m is empty.
            if len ( path[m] ) == 0:

                # . Compute a path for m.
                path[m] = path[frontnode] + [ m ]

                # . Put m on the back of the queue.
                queue.append ( ( m, frontnode ) )

            # . The path for m is not empty.
            else:

                # . Determine the intersection of the paths.
                intersection = [ ]
                for n in path[m]:
                    if n in path[frontnode]: intersection.append ( n )

                # . Check for a singleton intersection.
                if len ( intersection ) == 1:

                    # . Create and return the ring.
                    ring = _OrderPath ( path[frontnode], path[m] )
                    return ring

    # . Return a zero ring.
    return [ ]

def _IdentifyChains ( degree0, node0, currentstate, connectivity ):
    """Identify a single node in each separate chain of order 2 nodes."""

    # . Create a connectivity table for nodes of degree 2.
    connectivity2 = {}
    for ( degree, node ) in [ ( degree0, node0 ) ] + currentstate:
        if degree == 2:
            connections = []
            for m in connectivity[node]:
                if len ( connectivity[m] ) == 2: connections.append ( m )
            connectivity2[node] = connections
#        else:
#            break
#    print "AA>", connectivity2
    # . Get the list of degree 2 nodes.
    currentorder2 = [ ]
    nodes         = connectivity2.keys ( )
    nodes.sort ( )
    for node in nodes:
        connections = connectivity2[node]

        # . Save this node if it has no connections.
        if len ( connections ) == 0: currentorder2.append ( node )
        # . Delete the connections for this node.
        else:
            for m in connections: connectivity2[m].remove ( node )
            del connectivity2[node]

    # . Return the nodes.
    return currentorder2

def _OrderPath ( path1, path2 ):
    """Order the path."""

    # . Initialize newpath.
    newpath = path1

    # . Add the elements of path2 in reverse order skipping the first element.
    i = len ( path2 ) - 1
    while i > 0:
       newpath.append ( path2[i] )
       i = i - 1

    # . Find the minimum element in the list and its position.
    node0 = min ( newpath )
    pos0  = newpath.index ( node0 )

    # . Find the values of the elements either side of node0.
    if pos0 != 0:
        nodem1 = newpath[pos0-1]
    else:
        nodem1 = newpath[-1]
    if pos0 != len ( newpath ) - 1:
        nodep1 = newpath[pos0+1]
    else:
        nodep1 = newpath[0]

    # . If nodem1 is greater than nodep1 reverse the list and find the new position of node0.
    # . Note the newpath.reverse() does not seem to function very well so the reversing is
    # . done explicitly.
    if nodep1 > nodem1:
        temppath = newpath
        newpath  = [ ]
        i = len ( temppath ) - 1
        while i >= 0:
            newpath.append ( temppath[i] )
            i = i - 1
        pos0 = newpath.index ( node0 )

    # . Cyclically shift the elements in the list until node0 is the first one.
    if pos0 != 0:
        temppath = newpath
        newpath  = temppath[pos0:] + temppath[:pos0]

    # . Return the reordered path.
    return newpath

def _Trim ( node0, connectivity ):
    """Remove from the connection table all connections to node node0."""
    toremove = connectivity[node0]
    for node in toremove:
        nodes = connectivity[node]
        nodes.remove ( node0 )
        if len ( nodes ) <= 0: del connectivity[node]
    del connectivity[node0]

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
