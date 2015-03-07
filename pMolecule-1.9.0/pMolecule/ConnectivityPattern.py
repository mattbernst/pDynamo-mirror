#-------------------------------------------------------------------------------
# . File      : ConnectivityPattern.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A connectivity pattern is a set of atom and bond patterns that can be used
   to match various parts of a connectivity."""

import copy

from Bond import DoubleBond, SingleBond, TripleBond, UndefinedBond

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class _MatchObject ( object ):
    """A local class for holding match data."""

    def __init__ ( self ):
        """Constructor."""
        self.edges = {}
        self.nodes = {}

#===================================================================================================================================
# . Class.
#===================================================================================================================================
# . Can generalize to allow ranges or sets of values to match, rather than single values.

class AtomPattern ( object ):
    """A class to represent an atom pattern."""

    defaultAttributes = { "attributes" : None }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        self.attributes = {}
        for ( key, value ) in keywordArguments.iteritems ( ): self.attributes[key] = value

    @classmethod
    def FromList ( selfClass, attributeList, attributes, undefined ):
        """Constructor from list."""
        keywordArguments = {}
        for ( key, value ) in zip ( attributes, attributeList ):
            if value != undefined: keywordArguments[key] = value
        self = selfClass ( **keywordArguments )
        return self

    def Match ( self, atom ):
        """Match the pattern against an atom."""
        # . Initialization.
        isMatched = True
        # . Loop over all attributes.
        keys = self.attributes.keys ( )
        for key in keys:
            pAttribute = self.attributes[key]
            if pAttribute is not None:
                try:    attribute = getattr ( atom, key )
                except: attribute = None
#                print "ATOM>", key, pAttribute, attribute
                isMatched = ( attribute == pAttribute )
                if not isMatched: break
#        if not isMatched and self.attributes["atomicNumber"] == atom.atomicNumber:
#            print "ATOM FINAL>", isMatched, self.attributes, atom.__dict__
        return isMatched

    def ToList ( self, attributes, undefined ):
        """Get a list representation of the object."""
        selfList = []
        for attribute in attributes:
            selfList.append ( getattr ( self, attribute, undefined ) )
        return selfList

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class BondPattern ( object ):
    """A class to represent a bond pattern."""

    defaultAttributes = { "atomKey1"   :   -1 ,
                          "atomKey2"   :   -1 ,
                          "type"       : None ,
                          "typeObject" : None }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )
        self.DefineBondType ( )

    def DefineBondType ( self ):
        """Define the bond type."""
        bT = self.type
        if bT is not None:
            if   bT == "Double": self.typeObject = DoubleBond    ( )
            elif bT == "Single": self.typeObject = SingleBond    ( )
            elif bT == "Triple": self.typeObject = TripleBond    ( )
            else               : self.typeObject = UndefinedBond ( )

    @classmethod
    def FromList ( selfClass, attributeList, attributes, undefined ):
        """Constructor from list."""
        keywordArguments = {}
        for ( key, value ) in zip ( attributes, attributeList ):
            if value != undefined: keywordArguments[key] = value
        self = selfClass ( **keywordArguments )
        return self

    def Match ( self, bond ):
        """Match the pattern against a bond."""
        if self.typeObject is UndefinedBond ( ): isMatched = True
        else:  isMatched = ( self.typeObject is bond.type )
#        print "BOND MATCH> ", self.typeObject, bond.type, isMatched
        return isMatched

    def Other ( self, atomKey ):
        """Return the index of the other atom."""
        if   self.atomKey1 == atomKey: return self.atomKey2
        elif self.atomKey2 == atomKey: return self.atomKey1
        else: return -1

    def ToList ( self, attributes, undefined ):
        """Get a list representation of the object."""
        selfList = []
        for attribute in attributes:
            selfList.append ( getattr ( self, attribute, undefined ) )
        return selfList

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ConnectivityPattern ( object ):
    """A class to represent a connectivity pattern."""

    defaultAttributes = { "atomConnections" : None,
                          "atomPatterns"    : None,
                          "bondConnections" : None,
                          "bondPatterns"    : None,
                          "label"           : None,
                          "upperBound"      :  -1 }

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def ClearRepresentations ( self ):
        """Clear all representations."""
        self.atomConnections = None
        self.bondConnections = None
        self.upperBound      = -1

    def FindAllMatches ( self, connectivity, selection = None ):
        """Find all occurrences of |self| that occur in |connectivity|.

        This method uses a depth-first method to find all matches. Uniqueness
        is assured only at the end when incorporating matches into the final
        list and not when generating them. The algorithm in this method could
        undoubtedly be made more efficient.
        """
        # . Initialization.
        matches = []
        # . Treat the selection.
        if selection is None: selection = range ( len ( connectivity.atoms ) )
        # . Simple checks to avoid too much work.
        if ( len ( connectivity.atoms ) >= len ( self.atomPatterns ) ) and ( len ( selection ) > 0 ):
            # . Find matches for the first atom pattern.
            firstAtom    = self.atomPatterns[0]
            isSingleAtom = ( len ( self.atomPatterns ) == 1 )
            for iAtom in selection:
                if firstAtom.Match ( connectivity.atoms[iAtom] ):
                    # . Single atom patterns.
                    if isSingleAtom: matches.append ( [ iAtom ] )
                    # . Multiple atom patterns.
                    else:
                        mObject          = _MatchObject ( )
                        mObject.nodes[0] = iAtom
                        mObjects         = self.MatchBranches ( iAtom, 0, mObject, connectivity, selection )
                        for mObject in mObjects:
                            if ( len ( mObject.nodes ) == len ( self.atomPatterns ) ) and ( len ( mObject.edges ) == len ( self.bondPatterns ) ):
                                match = []
                                for i in range ( len ( mObject.nodes ) ): match.append ( mObject.nodes[i] )
                                if match not in matches: matches.append ( match )
        return matches

    def GetConnections ( self, i ):
        """Get the connections for an atom pattern."""
        if ( self.bondConnections is not None ):
            try:    iAtom = int ( i )
            except: raise TypeError ( "Atom index must be an integer." )
            if ( iAtom >= 0 ) and ( iAtom < self.upperBound ):
                try:    return self.bondConnections[iAtom]
                except: return []
            else: raise IndexError ( "Connection index out of range." )
        else: raise AttributeError ( "Connection representation does not exist." )

    def MakeConnections ( self, upperBound = -1 ):
        """Make the connections representation of the bond container."""
        # . Make the representation.
        if ( self.atomConnections is None ) or ( self.bondConnections is None ):
            self.ClearRepresentations ( )
            # . Initialization.
            maximumIndex = len ( self.atomPatterns )
            self.atomConnections = {}
            self.bondConnections = {}
            for ( b, bond ) in enumerate ( self.bondPatterns ):
                i = bond.atomKey1
                j = bond.atomKey2
                if i not in self.atomConnections: self.atomConnections[i] = []
                if j not in self.atomConnections: self.atomConnections[j] = []
                if i not in self.bondConnections: self.bondConnections[i] = []
                if j not in self.bondConnections: self.bondConnections[j] = []
                self.atomConnections[i].append ( j )
                self.atomConnections[j].append ( i )
                self.bondConnections[i].append ( b )
                self.bondConnections[j].append ( b )
                maximumIndex = max ( maximumIndex, i, j )
            # . Sorting.
            for value in self.atomConnections.values ( ): value.sort ( )
            for value in self.bondConnections.values ( ): value.sort ( )
            # . Set the upper bound of the representation.
            self.upperBound = max ( maximumIndex, upperBound )
        # . Check the upperBound.
        else:
            if upperBound > self.upperBound: self.upperBound = upperBound
#        print "CONNECTIONS>", len ( self.atomConnections ), len ( self.bondConnections )

    def MatchBranches ( self, cNode, pNode, mObject, connectivity, selection ):
        """Match nodes connected to a root node."""
        # . Get the edges for the current nodes.
        cEdges = connectivity.bonds.GetConnections ( cNode )
        pEdges =   copy.copy ( self.GetConnections ( pNode ) )
        # . Remove matched edges from the data.
        for pEdge in pEdges:
            if pEdge in mObject.edges: pEdges.remove ( pEdge )
        # . There are no edges to match so return the current data as is.
        if len ( pEdges ) == 0:
            mObjects = [ mObject ]
        # . Match all edge patterns and descending branches in order.
        else:
            # . Initialization.
            mObjects = []
            # . Always do full double loop without exiting so that all pattern/real edge permutations are explored.
            for pEdge in pEdges:
                qNode      = self.bondPatterns[pEdge].Other ( pNode )
                isQMatched = ( qNode in mObject.nodes )
                # . Loop over edges for cNode.
                for cEdge in cEdges:
                    dNode = connectivity.bonds[cEdge].Other ( cNode )
                    # . dNode must be in selection.
                    if dNode in selection:
                        isDMatched = ( dNode in mObject.nodes.values ( ) )
                        # . Check for an edge match.
                        if self.bondPatterns[pEdge].Match ( connectivity.bonds[cEdge] ):
                            # . Check for a node match.
                            if self.atomPatterns[qNode].Match ( connectivity.atoms[dNode] ):
                                # . Copy and update the match.
                                nObject = copy.copy ( mObject )
                                nObject.edges[pEdge] = cEdge
                                # . A previously matched pattern node.
                                if isQMatched:
                                    # . The match is OK if the real node corresponds.
                                    if mObject.nodes[qNode] == dNode: mObjects.append ( nObject )
                                # . A previously unmatched real node.
                                elif not isDMatched:
                                    nObject.nodes[qNode] = dNode
                                    # . Branches arising from the new node.
                                    tObjects = self.MatchBranches ( dNode, qNode, nObject, connectivity, selection )
                                    # . Branches arising from unmatched edges of the current node (if there are any).
                                    if len ( tObjects ) > 0:
                                        if tObjects[0] is nObject:
                                            mObjects.append ( nObject )
                                        else:
                                            for tObject in tObjects:
                                                mObjects.extend ( self.MatchBranches ( cNode, pNode, tObject, connectivity, selection ) )
        # . Finish up.
        return mObjects

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
