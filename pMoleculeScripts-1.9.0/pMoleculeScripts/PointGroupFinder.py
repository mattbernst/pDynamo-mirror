#-------------------------------------------------------------------------------
# . File      : PointGroupFinder.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Point group tests."""

import glob, math, os

from pBabel            import XYZFile_ToSystem
from pCore             import Clone, logFile, LogFileActive, Real1DArray, Real2DArray, Selection, Vector3
from PointGroup        import PointGroups_FromText
from SymmetryOperation import Identity, ImproperRotation, Inversion, ProperRotation, Reflection

# . Which ring-finding to use.
try:
    from pGraph import BiconnectedComponents, Edge, EdgeVectorToPath, Graph, Node, PathToEdgeVector, VismaraRelevantCycles
    useFigueras = False
except:
    from pCore  import FindRingSets
    useFigueras = True

#
# . This algorithm is quite slow for very high symmetry systems (e.g. icosahedral groups) due to the ring finding. It may be
# . possible to optimize it by not including every group of connections in PostulateHigherOrderRotationAxes. However, simple
# . restrictions (like checking the length of the number of connections of a given distance) do not work.
#
# . Possible additions include:
#
# - Better way of handling indeterminate possible C2 axes.
# - Flexible point group assignment (if keys don't match).
# - Tolerances higher for larger Cn.
#
# - Include a pointGroupFinderState to avoid putting results into Finder itself.
#

#
# . Irreducible representation determination is quite sensitive to the degeneracy tolerances. Too high a tolerance bunches
# . vectors or states together and then illegal combinations of IRs can occur. This might be improved by a better algorithm.
# . As high an accuracy as possible is best in the CI, SCF or normal mode calculation.
#
# . Accidental degeneracies are not yet handled.
#
# . Also problems will arise if the wavefunction does not have the symmetry of the nuclear framework (possible in some cases,
# . particularly when there are distorted geometries).
#

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_includeCompositeImproperRotations = False

# . Tolerances.
# . The results can be quite sensitive to these tolerances.
_distanceTolerance = 0.1 # . In Angstroms.
_angleTolerance    = 0.1 # . In radians.

# . Other tolerances.
_characterMatchTolerance            = 0.1    # . Dimensionless.
_electronicStateDegeneracyTolerance = 0.01   # . Hartrees.
_orbitalDegeneracyTolerance         = 1.0e-3 # . Hartrees.
_normalModeDegeneracyTolerance      = 5.0    # . cm^-1.

#===================================================================================================================================
# . Utility functions.
#===================================================================================================================================
def AllCombinations ( items, k ):
    """Return a list of all combinations of k items from items."""
    combinations = []
    n = len ( items )
    if ( k > 0 ) and ( k <= n ):
        # . First combination.
        data = [ i for i in range ( k ) ]
        combinations.append ( [ items[i] for i in data ] )
        # . Successive combinations.
        while True:
            i = k - 1
            while ( i > 0 ) and ( data[i] == n - k + i ): i -= 1
            if ( i == 0 ) and ( data[i] == ( n - k ) ): break
            data[i] += 1
            for j in range ( i, k-1 ): data[j+1] = data[j] + 1
            combinations.append ( [ items[i] for i in data ] )
    return combinations

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PointGroupFinder ( object ):
    """Class for finding the point group of an object."""

    # . Attributes.
    defaultAttributes = { "angleTolerance"                     : _angleTolerance                     ,
                          "characterMatchTolerance"            : _characterMatchTolerance            ,
                          "distanceTolerance"                  : _distanceTolerance                  ,
                          "electronicStateDegeneracyTolerance" : _electronicStateDegeneracyTolerance ,
                          "orbitalDegeneracyTolerance"         : _orbitalDegeneracyTolerance         ,
                          "normalModeDegeneracyTolerance"      : _normalModeDegeneracyTolerance      ,
                          "characterSymmetryOperations"        : None  ,
                          "indeterminateC2Axes"                : False ,
                          "pointGroup"                         : None  ,
                          "symmetryOperations"                 : None  ,
                          "system"                             : None  }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        # . Defaults.
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        # . Input keywords - tolerances only.
        for ( key, value ) in keywordArguments.iteritems ( ):
            if key.find ( "Tolerance" ) >= 0: setattr ( self, key, value )
        # . Set additional tolerances.
        self.cosineParallelTolerance      = math.fabs ( 1.0 - math.cos ( self.angleTolerance ) )
        self.cosinePerpendicularTolerance = math.fabs (       math.sin ( self.angleTolerance ) )

    #===============================================================================================================================
    def AddAdditionalCycles ( self, graph, cycles ):
        """Add cycles by fusing existing cycles."""

        # . Do something only if there are sufficient cycles.
        if len ( cycles ) > 1:

            # . Initialization.
            newCycles = []

            # . Find the cycles each node is involved in.
            maximumCycleLength = 0
            nodeCycles         = {}
            for cycle in cycles:
                maximumCycleLength = max ( maximumCycleLength, len ( cycle ) )
                for node in cycle:
                    existing = nodeCycles.get ( node, [] )
                    existing.append ( cycle )
                    nodeCycles[node] = existing

            # . Loop over nodes.
            # . Only do something if there are sufficient cycles to fuse.
            for ( node, local ) in nodeCycles.iteritems ( ):
                if len ( local ) > maximumCycleLength:

                    # . Find the indices of the edges that this node is involved in.
                    edgeIndices = []
                    for ( e, edge ) in enumerate ( graph.edges ):
                        if ( edge.node1 is node ) or ( edge.node2 is node ): edgeIndices.append ( e )

                    # . Double loop over cycles.
                    for ( c, cycle ) in enumerate ( local[:-1] ):
                        fusedCycle = PathToEdgeVector ( graph, cycle, closePath = True )
                        for other in local[c+1:]:
                            edgeVector = PathToEdgeVector ( graph, other, closePath = True )

                            # . Check for common edges at node.
                            overlap = False
                            for e in edgeIndices:
                                if edgeVector[e]:
                                    overlap = True
                                    break

                            # . Fuse cycles and then check for edges at node.
                            if overlap:
                                for i in range ( len ( edgeVector ) ):
		                    fusedCycle[i] ^= edgeVector[i]
                                found = False
                                for e in edgeIndices:
                                    if fusedCycle[e]:
                                        found = True
                                        break
                                if not found: break

                        # . Check the fused cycle.
                        # . Only include a path if it is a cycle with length > maximumCycleLength and that does not contain the node.
                        found = False
                        for e in edgeIndices:
                            if fusedCycle[e]:
                                found = True
                                break
                        if not found:
                            try:
                                newCycle = EdgeVectorToPath ( graph, fusedCycle )[:-1]
                                if len ( newCycle ) > maximumCycleLength: newCycles.append ( newCycle )
                            except:
                                pass

            # . Finish up.
            cycles.extend ( newCycles )

    #===============================================================================================================================
    def CheckForAlignment ( self, order, iRing, iCenter, jRing, jCenter, coordinates3 ):
        """Check two rings for an alignment which is appropriate for an improper rotation."""
        # . Initialization.
        areAligned    = False
        improperOrder = 0
        # . Get minimum absolute angle between a center-node vector of one ring and all those of the other ring.
        iVector = Vector3.Null ( )
        jVector = Vector3.Null ( )
        iNode   = iRing[0]
        for i in range ( 3 ): iVector[i] = coordinates3[iNode,i] - iCenter[i]
        iVector.Normalize ( )
        angles = []
        for jNode in jRing:
            for i in range ( 3 ): jVector[i] = coordinates3[jNode,i] - jCenter[i]
            jVector.Normalize ( )
	    dot = iVector.Dot ( jVector )
	    if   dot >  1.0: dot =  1.0
	    elif dot < -1.0: dot = -1.0
	    angles.append ( math.fabs ( math.acos ( dot ) ) )
        angles.sort ( )
        minimumAngle = angles[0]
        # . Check for alignment.
        areEclipsed  = ( minimumAngle <= self.angleTolerance )
        areStaggered = ( math.fabs ( minimumAngle - math.pi / float ( order ) ) <= self.angleTolerance )
        assert not ( areEclipsed and areStaggered ), "Logic error."
        # . Staggered rings imply an S2n that cannot be generated from the ring Cn and sigmah.
        if areStaggered:
            areAligned    = True
	    improperOrder = 2 * order
        # . Eclipsed rings imply existing Cn and sigmah operations.
        elif areEclipsed and _includeCompositeImproperRotations:
            areAligned    = True
	    improperOrder = order
        return ( areAligned, improperOrder )

    #===============================================================================================================================
    def CheckForPlanarity ( self, nodes, coordinates3 ):
        """Check for planarity.

        Calculate the axis parameters if the ring is planar.
        """
        # . Initialization.
        axis     = None
        center   = None
        distance = None
        isPlanar = True
        # . Determine planarity.
        n        = len ( nodes )
        target   = math.pi * ( 1.0 - 2.0 / float ( n ) )
        for j in range ( n ):
            jIndex = nodes[j]
            if j == 0:
	        iIndex = nodes[-1]
	        kIndex = nodes[ 1]
	    elif j == ( n - 1 ):
	        iIndex = nodes[j-1]
	        kIndex = nodes[0]
            else:

	        iIndex = nodes[j-1]
	        kIndex = nodes[j+1]
	    angle = math.radians ( coordinates3.Angle ( iIndex, jIndex, kIndex ) )
            if math.fabs ( angle - target ) > self.angleTolerance:
	        isPlanar = False
	        break
        # . Determine axis parameters.
        if isPlanar:
            # . Calculate the coordinates of the ring center.
            axis   = Vector3.Null ( )
	    center = Vector3.Null ( )
	    for iNode in nodes:
	        for i in range ( 3 ): center[i] += coordinates3[iNode,i]
            for i in range ( 3 ): center[i] /= float ( n )
	    distance = center.Norm2 ( )
            # . The center is at the origin so determine the ring normal.
	    # . Arbitrarily use the coordinates of the first two nodes in the ring.
	    # . Direction is also arbitrary.
	    if distance <= self.distanceTolerance:
	        iNode    = nodes[0]
	        jNode    = nodes[1]
	        axis[0]  = coordinates3[iNode,1]*coordinates3[jNode,2] - coordinates3[jNode,1]*coordinates3[iNode,2]
	        axis[1]  = coordinates3[iNode,2]*coordinates3[jNode,0] - coordinates3[jNode,2]*coordinates3[iNode,0]
	        axis[2]  = coordinates3[iNode,0]*coordinates3[jNode,1] - coordinates3[jNode,0]*coordinates3[iNode,1]
	        axis.Normalize ( )
	        distance = 0.0
	        for i in range ( 3 ): center[i] = 0.0
	    # . Normalize.
	    else:
	        for i in range ( 3 ): axis[i] = center[i]
	        axis.Scale ( 1.0 / distance )
        return ( isPlanar, distance, axis, center )

    #===============================================================================================================================
    def FindPlanarRings ( self, connections, coordinates3 ):
        """Find all planar rings given a set of connections between nodes."""

        # . Initialization.
        planarRings = {}

        # . Figueras version.
        if useFigueras:

            # . Get the connection table.
            table = {}
            for ( i, j ) in connections:
	        iconns = table.get ( i, [] )
	        jconns = table.get ( j, [] )
	        iconns.append ( j )
	        jconns.append ( i )
	        table[i] = iconns
	        table[j] = jconns

            # . Get the ring sets.
            ringSets = FindRingSets ( table )

            # . Analyze the rings to find those which are planar.
            for ringSet in ringSets:
	        for ring in ringSet:
	            ( isPlanar, distance, axis, center ) = self.CheckForPlanarity ( ring, coordinates3 )
	            if isPlanar:
  		        n = len ( ring )
		        rings = planarRings.get ( n, [] )
		        rings.append ( ( ring, distance, axis, center ) )
		        planarRings[n] = rings

        # . Relevant cycles version.
        else:

            # . Create the graph.
            graph = Graph ( )
            nodes = {}
            for connection in connections:
                for i in connection:
                    if i not in nodes:
                        node = Node ( )
                        node.index = i
                        nodes[i] = node
                        graph.AddNode ( node )
            for ( i, j ) in connections:
                graph.AddEdge ( Edge ( nodes[i], nodes[j] ) )

            # . Get the cycles.
            cycles = VismaraRelevantCycles ( graph )

            # . Add additional cycles by fusing the relevant cycles for each node.
            self.AddAdditionalCycles ( graph, cycles )

            # . Loop over cycles.
            for cycle in cycles:
                ring = []
                for node in cycle: ring.append ( node.index )
	        ( isPlanar, distance, axis, center ) = self.CheckForPlanarity ( ring, coordinates3 )
	        if isPlanar:
  	            n = len ( ring )
	            rings = planarRings.get ( n, [] )
	            rings.append ( ( ring, distance, axis, center ) )
	            planarRings[n] = rings

        # . Finish up.
        return planarRings

    #===============================================================================================================================
    def FindSystemPointGroup ( self, system, log = logFile ):
        """Find the point group of a system."""

        # . Initialization.
        symmetryOperations = {}

        # . Clone the coordinates and translate to the center of mass.
#        coordinates3 = Clone ( system.coordinates3 )
        coordinates3 = system.coordinates3
        masses       = system.atoms.GetItemAttributes ( "mass" )
        center       = coordinates3.Center ( weights = masses )
        center.Scale ( -1.0 )
        coordinates3.Translate ( center )

        # . Classify atoms into possible symmetry classes.
        elementGroups = self.PostulateSymmetryEquivalentNodes ( system.atoms, coordinates3, log = log )

        # . Check the moments of inertia for basic information about the system.
        self.ProcessMomentsOfInertia ( elementGroups, coordinates3, masses, symmetryOperations )

        # . Find higher-order operations.
        ( properRotations, improperRotations ) = self.PostulateHigherOrderRotationAxes ( elementGroups, coordinates3, log = log )
        self.ProcessHigherOrderRotationAxes ( elementGroups, coordinates3, properRotations, improperRotations, symmetryOperations )

        # . Find order-2 operations (C2 rotation, inversion, reflection).
        self.ProcessOrderTwoOperations ( elementGroups, coordinates3, symmetryOperations )

        # . Additional processing for indeterminate C2 axes.
        self.ProcessIndeterminateC2Axes ( elementGroups, coordinates3, symmetryOperations )

        # . Remove redundant S4s.
        self.RemoveRedundantS4s ( symmetryOperations )

        # . Append the identity.
        symmetryOperations["E"] = [ Identity ( ) ]

        # . Printing.
        keys = symmetryOperations.keys ( )
        keys.sort ( )
        if LogFileActive ( log ):
            if len ( symmetryOperations ) > 0:
                log.Text ( "\nSymmetry Operations:\n" )
                for key in keys:
                    log.Text ( "{:<20s} {:5d}\n".format ( key, len ( symmetryOperations[key] ) ) )
            else: log.Text ( "\nNo symmetry operations identified.\n" )

        # . Generate the operation key.
        items = []
        for key in keys:
            n = len ( symmetryOperations[key] )
            if n == 1: items.append ( key )
            else:      items.append ( repr ( n ) + "*" + key )
        operationKey = " ".join ( items )

        # . Identify point group from a point group dictionary.
        # . Maybe need to be cleverer here if there is no match by searching for best compromise group.
        pointGroups = PointGroups_FromText ( )
        group       = pointGroups.get ( operationKey, None )
        if group is None: groupName = None
        else:             groupName = group.label

        # . Save all data.
        self.elementGroups      = elementGroups
        self.pointGroup         = group
        self.symmetryOperations = symmetryOperations
        self.system             = system

        # . Finish up.
        return groupName

    #===============================================================================================================================
    def IdentifyCIStateIrreducibleRepresentations ( self, degeneracyTolerance = None, log = logFile, stateIndices = None, includeCoreOrbitals = False ):
        """Identify CI state irreducible representations."""
        # . Set up.
        self.SetUpIrreducibleRepresentationCalculation ( log = log )
        # . Get system data.
        qcModel = self.system.energyModel.qcModel
        try:
            qcState                = self.system.configuration.qcState
            ciStateEnergies        = qcState.ciStateEnergies
            numberCIConfigurations = qcState.numberCIConfigurations
        except:
            raise ValueError ( "Unable to retrieve CI state data." )
        # . Get state indices.
        if stateIndices is None:
            stateIndices = Selection.FromIterable ( range ( numberCIConfigurations ) )
        # . Initialization.
        stateIrreducibleRepresentations = []
        if len ( stateIndices ) > 0:
            # . Get characters for each state under each symmetry operation.
            characters  = Real2DArray.WithExtents ( len ( stateIndices ), len ( self.characterSymmetryOperations ) )
            eigenValues = []
            characters.Set ( 1.0 )
            # . Get state energies.
            for stateIndex in stateIndices:
                eigenValues.append ( ciStateEnergies[stateIndex] )
            # . Loop over symmetry operations.
            for ( j, operation ) in enumerate ( self.characterSymmetryOperations ):
                if operation.Label ( ) != "E":
                    localCharacters = qcModel.CIStateCharacters ( self.system.configuration, operation.transformationMatrix, operation.mapping, stateIndices, includeCoreOrbitals = includeCoreOrbitals )
                    for ( i, character ) in enumerate ( localCharacters ): characters[i,j] = character
            # . Print the character table.
            if LogFileActive ( log ): characters.Print ( title = "CI State Character Table" )
            # . Identify symmetries.
            if degeneracyTolerance is None: degeneracyTolerance = self.electronicStateDegeneracyTolerance
            stateIrreducibleRepresentations = self.IdentifyVectorIrreducibleRepresentations ( characters, eigenValues, degeneracyTolerance )
        # . Finish up.
        return stateIrreducibleRepresentations

    #===============================================================================================================================
    def IdentifyNormalModeIrreducibleRepresentations ( self, degeneracyTolerance = None, log = logFile, modeIndices = None ):
        """Identify normal mode irreducible representations."""
        # . Set up.
        self.SetUpIrreducibleRepresentationCalculation ( log = log )
        # . Get mode indices.
        if modeIndices is None:
            modeIndices = range ( 3 * len ( self.system.atoms ) )
        # . Get the normal modes.
        nmState = self.system.configuration.nmState
        if nmState is None: raise ValueError ( "Normal modes not found for system." )
        else:
            frequencies = nmState.frequencies
            modes       = nmState.modes
            weights     = nmState.weights
        # . Initialization.
        modeIrreducibleRepresentations = []
        if len ( modeIndices ) > 0:
            # . Get characters for each mode under each symmetry operation.
            characters  = Real2DArray.WithExtents ( len ( modeIndices ), len ( self.characterSymmetryOperations ) )
            eigenValues = []
            vector1     = Vector3.Uninitialized ( )
            vector2     = Vector3.Uninitialized ( )
            characters.Set ( 1.0 )
            for ( i, modeIndex ) in enumerate ( modeIndices ):
                eigenValues.append ( frequencies[modeIndex] )
                for ( j, operation ) in enumerate ( self.characterSymmetryOperations ):
                    if operation.order == 1:
                        character = 1.0
                    else:
                        character = 0.0
                        for iatom in range ( len ( self.system.atoms ) ):
                            jatom = operation.mapping[iatom]
                            # . Unweight modes to obtain normalized vectors.
                            for c in range ( 3 ): vector1[c] = modes[3*iatom+c,modeIndex] * weights[3*iatom+c]
                            for c in range ( 3 ): vector2[c] = modes[3*jatom+c,modeIndex] * weights[3*jatom+c]
                            operation.ApplyTo ( vector1 )
                            character += vector1.Dot ( vector2 )
                    characters[i,j] = character
            # . Print the character table.
            if LogFileActive ( log ): characters.Print ( title = "Normal Mode Character Table" )
            # . Identify symmetries.
            if degeneracyTolerance is None: degeneracyTolerance = self.normalModeDegeneracyTolerance
            iRs = self.IdentifyVectorIrreducibleRepresentations ( characters, eigenValues, degeneracyTolerance )
            # . Convert to lower case.
            modeIrreducibleRepresentations = [ iR.lower ( ) for iR in iRs ]
        # . Finish up.
        return modeIrreducibleRepresentations

    #===============================================================================================================================
    def IdentifyOrbitalIrreducibleRepresentations ( self, degeneracyTolerance = None, log = logFile, orbitalIndices = None, useDensityP = True ):
        """Identify orbital irreducible representations."""
        # . Set up.
        self.SetUpIrreducibleRepresentationCalculation ( log = log )
        # . Get system data.
        qcModel = self.system.energyModel.qcModel
        ( orbitalEnergies, homo, lumo ) = qcModel.OrbitalEnergies ( self.system.configuration, useDensityP = useDensityP )
        # . Get orbital indices.
        if orbitalIndices is None:
            orbitalIndices = Selection.FromIterable ( range ( self.system.energyModel.qcAtoms.nobasis ) )
        # . Initialization.
        orbitalIrreducibleRepresentations = []
        if len ( orbitalIndices ) > 0:
            # . Get characters for each orbital under each symmetry operation.
            characters  = Real2DArray.WithExtents ( len ( orbitalIndices ), len ( self.characterSymmetryOperations ) )
            eigenValues = []
            characters.Set ( 1.0 )
            # . Get orbital energies.
            for orbitalIndex in orbitalIndices:
                eigenValues.append ( orbitalEnergies[orbitalIndex] )
            # . Loop over symmetry operations.
            for ( j, operation ) in enumerate ( self.characterSymmetryOperations ):
                if operation.Label ( ) != "E":
                    localCharacters = qcModel.OrbitalCharacters ( self.system.configuration, operation.transformationMatrix, operation.mapping, orbitalIndices, useDensityP = useDensityP )
                    for ( i, character ) in enumerate ( localCharacters ): characters[i,j] = character
            # . Print the character table.
            if LogFileActive ( log ): characters.Print ( title = "Orbital Character Table" )
            # . Identify symmetries.
            if degeneracyTolerance is None: degeneracyTolerance = self.orbitalDegeneracyTolerance
            iRs = self.IdentifyVectorIrreducibleRepresentations ( characters, eigenValues, degeneracyTolerance )
            # . Convert to lower case.
            orbitalIrreducibleRepresentations = [ iR.lower ( ) for iR in iRs ]
        # . Finish up.
        return orbitalIrreducibleRepresentations

    #===============================================================================================================================
    def IdentifyVectorIrreducibleRepresentations ( self, characters, eigenValues, degeneracyTolerance ):
        """Identify irreducible representations from characters."""
        # . Basic initialization.
        assigned        = set ( )
        numberClasses   = len ( self.characterSymmetryOperations )
        numberVectors   = characters.rows
        labels          = [ "?" for i in range ( numberVectors ) ]
        # . Assign vectors first to degenerate representations.
        if self.pointGroup.maximumDegeneracy > 1:
            # . Determine degenerate groups.
            oldE   = eigenValues[0]
            groups = [ [ 0 ] ]
            for i in range ( 1, numberVectors ):
                e = eigenValues[i]
                if math.fabs ( e - oldE ) <= degeneracyTolerance:
                    groups[-1].append ( i )
                else:
                    groups.append ( [ i ] )
                oldE = e
            # . Allocate space.
            totalCharacters = Real1DArray.WithExtent ( numberClasses )
            # . Process degenerate groups.
            for group in groups:
                indices = list ( group )
                while len ( indices ) > 1:
                    n = len ( indices )
                    if n == 1:
                        break
                    else:
                        # . Loop over all possible combinations of indices in each group.
                        done = False
                        for degeneracy in range ( min ( n, self.pointGroup.maximumDegeneracy ), 1, -1 ):
                            combinations = AllCombinations ( indices, degeneracy )
                            for combination in combinations:
                                totalCharacters.Set ( 0.0 )
                                for vectorIndex in combination:
                                    for c in range ( numberClasses ): totalCharacters[c] += characters[vectorIndex,c]
                                # . Check for a match.
                                for ( r, label ) in enumerate ( self.pointGroup.irreducibleRepresentations ):
                                    matched = True
                                    for c in range ( numberClasses ):
                                        if math.fabs ( totalCharacters[c] - self.pointGroup.characterTable[r,c] ) > self.characterMatchTolerance:
                                            matched = False
                                            break
                                    if matched:
                                        for vectorIndex in combination:
                                            assigned.add   ( vectorIndex )
                                            indices.remove ( vectorIndex )
                                            labels[vectorIndex] = label
                                        done = True
                                        break
                                if done: break
                            if done: break
                        # . No match found.
                        if not done:
                            break
        # . Assign vectors to non-degenerate representations.
        if len ( assigned ) < numberVectors:
            newUnassigned = []
            unassigned    = list ( set ( range ( numberVectors ) ).difference ( assigned ) )
            unassigned.sort ( )
            for vectorIndex in unassigned:
                done = False
                for ( r, label ) in enumerate ( self.pointGroup.irreducibleRepresentations ):
                    matched = True
                    for c in range ( numberClasses ):
                        if math.fabs ( characters[vectorIndex,c] - self.pointGroup.characterTable[r,c] ) > self.characterMatchTolerance:
                            matched = False
                            break
                    if matched:
                        done                = True
                        labels[vectorIndex] = label
                        break
                if not done: newUnassigned.append ( vectorIndex )
        # . Need to handle accidental degeneracy here.
        # . Determine sum of characters then possible solutions in terms of integer combinations of irreducible representations (if any exist).
        # . Finish up.
        return labels

    #===============================================================================================================================
    def PostulateHigherOrderRotationAxes ( self, elementGroups, coordinates3, log = None ):
        """Analyze group connectivity by distance to determine possible proper and improper rotation axes.

        Only axes of order >= 3 are searched for. In addition, composite symmetry operations are omitted.
        Thus, for example all S(2n+1) can be skipped as they automatically imply a C(2n+1) axis with a
        perpendicular mirror plane.
        """

        # . Initialization.
        # . Possible symmetry operations.
        improperRotations = {}
        properRotations   = {}

        # . Process each set of groups in turn.
        for ( element, groups ) in elementGroups.iteritems ( ):
            for group in groups:

                # . Skip small groups.
	        if len ( group ) > 2:

        	    # . Get the all internode distances within the group.
		    distances = []
        	    for i in range ( 1, len ( group ) ):
	                iNode = group[i]
                        for j in range ( i ):
                	    jNode = group[j]
                	    distances.append ( ( coordinates3.Distance ( iNode, jNode ), ( iNode, jNode ) ) )

        	    # . Separate distances by value.
		    distances.sort ( )
                    ( oldr, nodes ) = distances[0]
                    connectionSets = [ [ oldr, [ nodes ] ] ]
		    if len ( distances ) > 1:
        	        for ( r, nodes ) in distances[1:]:
                	    if math.fabs ( r - oldr ) <= self.distanceTolerance:
                	        connectionSets[-1][1].append ( nodes )
                	    else:
                                if len ( connectionSets[-1][1] ) < 2: connectionSets.pop ( -1 )
                	        connectionSets.append ( [ r, [ nodes ] ] )
                	        oldr = r
                        if len ( connectionSets[-1][1] ) < 2: connectionSets.pop ( -1 )

                    # . Process connections.
		    for ( r, connections ) in connectionSets:

                        # . Look for possible S4s.
		        # . This search involves looking for perpendicular connections equidistant from the origin but on opposite sides.
                        if len ( connections ) > 1:

                            # . Get data for the connections.
			    connectionData = []
			    for ( iNode, jNode ) in connections:
			        iAxis = Vector3.Null ( )
			        for i in range ( 3 ): iAxis[i] = coordinates3[iNode,i] + coordinates3[jNode,i]
                                distance = iAxis.Norm2 ( )
			        # . Centers at the origin do not contribute to S4 operations.
			        if distance > self.distanceTolerance:
			            iAxis.Scale ( 1.0 / distance )
				    iVector = Vector3.Null ( )
				    for i in range ( 3 ): iVector[i] = coordinates3[iNode,i] - coordinates3[jNode,i]
				    iVector.Normalize ( )
                                    connectionData.append ( ( distance, iAxis, iVector ) )

                            # . Check for appropriate pairs.
                            for i in range ( 1, len ( connectionData ) ):
			        ( iDistance, iAxis, iVector ) = connectionData[i]
			        for j  in range ( i ):
			            ( jDistance, jAxis, jVector ) = connectionData[j]
				    if ( math.fabs ( iDistance - jDistance     ) <= self.distanceTolerance            ) and \
		                       ( math.fabs ( iAxis.Dot ( jAxis ) + 1.0 ) <= self.cosineParallelTolerance      ) and \
				       ( math.fabs ( iVector.Dot ( jVector )   ) <= self.cosinePerpendicularTolerance ):
				        axes     = improperRotations.get ( 4, [] )
				        isUnique = True
				        for axis in axes:
					    if ( math.fabs ( math.fabs ( iAxis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
				 	        isUnique = False
				 	        break
				        if isUnique:
					    axes.append ( iAxis )
					    improperRotations[4] = axes

                        # . Higher rotations (n > 2).
                        if len ( connections ) > 2:

                            # . Find all planar rings.
                            planarRings = self.FindPlanarRings ( connections, coordinates3 )

                            # . Proper rotations.
			    for ( order, rings ) in planarRings.iteritems ( ):
			        for ( iRing, iDistance, iAxis, iCenter ) in rings:
                                    axes     = properRotations.get ( order, [] )
				    isUnique = True
				    for axis in axes:
				        if ( math.fabs ( math.fabs ( iAxis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
					    isUnique = False
					    break
				    if isUnique:
				        axes.append ( iAxis )
				        properRotations[order] = axes

                            # . Improper rotations.
                            # . Find aligned pairs of parallel rings of the same size equidistant from the origin but on opposite sides.
                            for ( order, rings ) in planarRings.iteritems ( ):
			        for i in range ( 1, len ( rings ) ):
			            ( iRing, iDistance, iAxis, iCenter ) = rings[i]
				    if iDistance > self.distanceTolerance:
			                for j in range ( i ):
                                            ( jRing, jDistance, jAxis, jCenter ) = rings[j]
					    if ( math.fabs ( iDistance - jDistance     ) <= self.distanceTolerance ) and \
					       ( math.fabs ( iAxis.Dot ( jAxis ) + 1.0 ) <= self.cosineParallelTolerance   ):
					       ( areAligned, improperOrder ) = self.CheckForAlignment ( order, iRing, iCenter, jRing, jCenter, coordinates3 )
					       if areAligned:
					           axes     = improperRotations.get ( improperOrder, [] )
					           isUnique = True
					           for axis in axes:
					               if ( math.fabs ( math.fabs ( iAxis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
						           isUnique = False
						           break
					           if isUnique:
					               axes.append ( iAxis )
						       improperRotations[improperOrder] = axes

        # . Print results.
        if LogFileActive ( log ):
            for ( tag, rotations ) in ( ( "Proper", properRotations ), ( "Improper", improperRotations ) ):
	        if len ( rotations ) <= 0: log.Text ( "\nNo candidate " + tag.lower ( ) + " rotation axes found.\n" )
                else:
                    orders = rotations.keys ( )
		    orders.sort ( )
		    log.Text ( "\nCandidate " + tag + " Rotation Axes:\n" )
		    for order in orders:
		        log.Text ( "{:5d} {:5d}\n".format ( order, len ( rotations[order] ) ) )

        # . Finish up.
        return ( properRotations, improperRotations )

    #===============================================================================================================================
    def PostulateSymmetryEquivalentNodes ( self, atoms, coordinates3, log = None ):
        """Group atoms by element type and distance to origin."""
        elementdata   = {}
        elementGroups = {}
        for ( i, atom ) in enumerate ( atoms ):
            r = 0.0
            for j in range ( 3 ): r += coordinates3[i,j]**2
            r = math.sqrt ( r )
            data = elementdata.get ( atom.atomicNumber, [] )
            data.append ( ( r, i ) )
            elementdata[atom.atomicNumber] = data
        elements = elementdata.keys ( )
        elements.sort ( )
        for element in elements:
            data = elementdata[element]
            data.sort ( )
            ( oldr, i ) = data[0]
            groups = [ [ i ] ]
            if len ( data ) > 1:
                for ( r, i ) in data[1:]:
                    if math.fabs ( r - oldr ) <= self.distanceTolerance:
                        groups[-1].append ( i )
                    else:
                        groups.append ( [ i ] )
                        oldr = r
            elementGroups[element] = groups
        if LogFileActive ( log ):
            log.Text ( "\n\nElement Groups:\n" )
            for element in elements:
                groups = elementGroups[element]
                log.Text ( "{!r} {:d}: ".format ( element, len ( groups ) ) )
                for group in groups: log.Text ( "{:d}".format ( len ( group ) ) )
	        log.Text ( "\n" )
        return elementGroups

    #===============================================================================================================================
    def ProcessHigherOrderRotationAxes ( self, elementGroups, coordinates3, properRotations, improperRotations, symmetryOperations ):
        """Process postulated higher-order rotation axes.

        All non-unique operations of lower order are excluded.

        Thus, a S2m is redundant if there is a higher S2n axis such that (2m) * p = (2n) where p is odd.

        A Cm is redundant if there is a higher S2n axis such that m * p = (2n) where p is even.

        A Cm is also redundant if there is a higher Cn such that m * p = n where p is any integer.
        """

        # . Find the orders of the operations.
        improperOrders = improperRotations.keys ( )
        improperOrders.sort ( reverse = True )
        properOrders   = properRotations.keys   ( )
        properOrders.sort   ( reverse = True )

        # . Maximum orders.
        maximumImproperOrder = 0
        maximumProperOrder   = 0
        if len ( improperOrders ) > 0: maximumImproperOrder = improperOrders[0]
        if len (   properOrders ) > 0: maximumProperOrder   =   properOrders[0]

        # . Improper rotations.
        for order in improperOrders:
            axes       = improperRotations[order]
            operations = []

            # . Loop over axes.
            for axis in axes:
                # . Check whether to test.
                toTest = True
                toTry  = []
                i      = 3
                while ( i * order <= maximumImproperOrder ):
                    toTry.extend ( symmetryOperations.get ( "S{:d}".format ( i * order ), [] ) )
                    i += 2
                for other in toTry:
                    if ( math.fabs ( math.fabs ( other.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                        toTest = False
                        break
                if toTest:
                    operation = ImproperRotation ( axis = Clone ( axis ), order = order )
                    if operation.EstablishSymmetryRelatedPairs ( elementGroups, coordinates3, self.distanceTolerance ):
                        operations.append ( operation )
            if len ( operations ) > 0: symmetryOperations["S{:d}".format ( order )] = operations

        # . Proper rotations.
        for order in properOrders:
            axes       = properRotations[order]
            operations = []

            # . Loop over axes.
            for axis in axes:
                # . Check whether to test.
                toTest = True
                toTry  = []
                i = 2
                while ( i * order <= maximumImproperOrder ):
                    toTry.extend ( symmetryOperations.get ( "S{:d}".format ( i * order ), [] ) )
                    i += 2
                i = 2
                while ( i * order <= maximumProperOrder ):
                    toTry.extend ( symmetryOperations.get ( "C{:d}".format ( i * order ), [] ) )
                    i += 1
                for other in toTry:
                    if ( math.fabs ( math.fabs ( other.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                        toTest = False
                        break
                # . Test.
                if toTest:
                    operation = ProperRotation ( axis = Clone ( axis ), order = order )
                    if operation.EstablishSymmetryRelatedPairs ( elementGroups, coordinates3, self.distanceTolerance ):
                        operations.append ( operation )
            if len ( operations ) > 0: symmetryOperations["C{:d}".format ( order )] = operations

    #===============================================================================================================================
    def ProcessIndeterminateC2Axes ( self, elementGroups, coordinates3, symmetryOperations ):
        """Process indeterminate C2 axes."""
        # . A fudge until something better comes along.
        if self.indeterminateC2Axes:

            # . Only treat if there are more than one mirror plane.
            sigmas = symmetryOperations.get ( "sigma", [] )
            if len ( sigmas ) > 1:

                # . Gather all rotations (C2n and S4n) which can produce a C2.
                existing = []
                for ( key, operations ) in symmetryOperations.iteritems ( ):
                    if ( key.startswith ( "C" ) and ( operations[0].order % 2 == 0 ) ) or \
                       ( key.startswith ( "S" ) and ( operations[0].order % 4 == 0 ) ):
                        for operation in operations:
                            existing.append ( operation.axis )

                # . Find possible C2 axes from pairs of perpendicular mirror planes.
                axes = []
                for i in range ( 1, len ( sigmas ) ):
                    iNormal = sigmas[i].normal
                    for j in range ( i ):
                        jNormal = sigmas[j].normal
                        # . Perpendicular normals.
                        if ( math.fabs ( iNormal.Dot ( jNormal ) ) <= self.cosinePerpendicularTolerance ):
                            axis = Clone ( iNormal )
                            axis.Cross ( jNormal )
                            axis.Normalize ( )
                            # . Check for uniqueness.
                            isUnique = True
                            for other in existing:
			        if ( math.fabs ( math.fabs ( other.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
			            isUnique = False
			            break
                            if isUnique:
                                axes.append     ( axis )
                                existing.append ( axis )

                # . Process the additional axes.
                if len ( axes ) > 0:
                    rotations = symmetryOperations.get ( "C2", [] )
                    for axis in axes:
                        operation = ProperRotation ( axis = Clone ( axis ), order = 2 )
                        if operation.EstablishSymmetryRelatedPairs ( elementGroups, coordinates3, self.distanceTolerance ):
                            rotations.append ( operation )
                    if len ( rotations ) > 0: symmetryOperations["C2"] = rotations

    #===============================================================================================================================
    def ProcessMomentsOfInertia ( self, elementGroups, coordinates3, masses, symmetryOperations ):
        """Check the moments of inertia for basic information about the system."""

        # . Inertia tolerance.
        inertiaTolerance = 2.0 * sum ( masses ) * self.distanceTolerance**2

        # . Find the number of "zero" adjusted moments of inertia.
        ( moments, axes ) = coordinates3.MomentsOfInertia ( weights = masses )
        mr2    = 0.5 * sum ( moments )
        zeroes = []
        for i in range ( 3 ):
            if ( math.fabs ( mr2 - moments[i] ) <= inertiaTolerance ): zeroes.append ( i )
        nZeroes = len ( zeroes )

        # . Initialization.
        axis = Vector3.Null ( )

        # . Check for planarity.
        if nZeroes > 0:
	    for i in zeroes:
	        for j in range ( 3 ): axis[j] = axes[j,i]
	        operation = Reflection ( normal = Clone ( axis ), order = 1 )
                if operation.EstablishSymmetryRelatedPairs ( elementGroups, coordinates3, self.distanceTolerance ):
                    symmetryOperations["sigma"] = [ operation ]
	            break

        # . Check for linearity.
        if nZeroes > 1:
	    for i in range ( 3 ):
                if ( ( i + 1 ) % 3 in zeroes ) and ( ( i + 2 ) % 3 in zeroes ):
   	            for j in range ( 3 ): axis[j] = axes[j,i]
                    operation = ProperRotation ( axis = Clone ( axis ), order = 1 )
                    if operation.EstablishSymmetryRelatedPairs ( elementGroups, coordinates3, self.distanceTolerance ):
                        symmetryOperations["CInfinity"] = [ operation ]
	                break

    #===============================================================================================================================
    def ProcessOrderTwoOperations ( self, elementGroups, coordinates3, symmetryOperations ):
        """Find order-2 operations (C2 rotation, inversion, reflection)."""

        # . Initialization.
        axis          = Vector3.Null ( )
        vij           = Vector3.Null ( )
        moleculePlane = None
        rotations     = []
        reflections   = []
        if "sigma" in symmetryOperations:
            moleculePlane = symmetryOperations["sigma"][0]
            reflections.append ( moleculePlane )

        # . Gather all rotations (C2n and S4n) which can produce a C2.
        existing = []
        for ( key, operations ) in symmetryOperations.iteritems ( ):
            if ( key.startswith ( "C" ) and ( operations[0].order % 2 == 0 ) ) or \
               ( key.startswith ( "S" ) and ( operations[0].order % 4 == 0 ) ):
                existing.extend ( operations )

        # . Inversion.
        operation = Inversion ( )
        if operation.EstablishSymmetryRelatedPairs ( elementGroups, coordinates3, self.distanceTolerance ):
            symmetryOperations["i"] = [ operation ]

        # . C2 rotations and reflections.
        # . Loop over all pairs of possibly symmetry-related atoms.
        for ( element, groups ) in elementGroups.iteritems ( ):
            for group in groups:
	        n = len ( group )
                for i in range ( 1, n ):
	            inode = group[i]
		    for j in range ( i ):
		        jnode = group[j]
                        # . C2.
		        maximumValue       = 0.0
		        tryC2              = True
                        usedMolecularPlane = False
		        for s in range ( 3 ):
		            axis[s]      = coordinates3[inode,s] + coordinates3[jnode,s]
			    maximumValue = max ( maximumValue, math.fabs ( axis[s] ) )
		        if ( maximumValue < self.distanceTolerance ):
		            if ( moleculePlane is not None ):
# . Check.
                                moleculePlane.normal.CopyTo ( axis )
                                usedMolecularPlane = True
			    else:
                                self.indeterminateC2Axes = True
                                tryC2 = False
		        if tryC2:
                            axis.Normalize ( )
                            isUnique = True
                            for other in ( existing + rotations ):
			        if ( math.fabs ( math.fabs ( other.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
			            isUnique = False
			            break
   		            if isUnique:
                                operation = ProperRotation ( axis = Clone ( axis ), order = 2 )
                                if operation.EstablishSymmetryRelatedPairs ( elementGroups, coordinates3, self.distanceTolerance ):
                                    rotations.append ( operation )
                        # . Additional C2 in cases where there is a molecular plane and the (ij) mid-point passes through the origin.
                        if usedMolecularPlane:
		            for s in range ( 3 ): vij[s] = coordinates3[inode,s] - coordinates3[jnode,s]
                            normal = moleculePlane.normal
                            axis[0]  = normal[1]*vij[2] - vij[1]*normal[2]
	                    axis[1]  = normal[2]*vij[0] - vij[2]*normal[0]
	                    axis[2]  = normal[0]*vij[1] - vij[0]*normal[1]
                            axis.Normalize ( )
                            isUnique = True
                            for other in ( existing + rotations ):
			        if ( math.fabs ( math.fabs ( other.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
			            isUnique = False
			            break
   		            if isUnique:
                                operation = ProperRotation ( axis = Clone ( axis ), order = 2 )
                                if operation.EstablishSymmetryRelatedPairs ( elementGroups, coordinates3, self.distanceTolerance ):
                                    rotations.append ( operation )
                        # . Reflection.
		        for s in range ( 3 ): axis[s] = coordinates3[inode,s] - coordinates3[jnode,s]
                        axis.Normalize ( )
                        isUnique = True
                        for other in reflections:
			    if ( math.fabs ( math.fabs ( other.normal.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
			        isUnique = False
			        break
		        if isUnique:
                            operation = Reflection ( normal = Clone ( axis ), order = 2 )
                            if operation.EstablishSymmetryRelatedPairs ( elementGroups, coordinates3, self.distanceTolerance ):
                                reflections.append ( operation )

        # . Finish up.
        if len ( rotations   ) > 0: symmetryOperations["C2"   ] = rotations
        if len ( reflections ) > 0: symmetryOperations["sigma"] = reflections

    #===============================================================================================================================
    def RemoveRedundantS4s ( self, symmetryOperations ):
        """Remove S4s that can be formed by combination of a C4 and a reflection.

        This is not done previously as S4s are found before C4s and sigmas.
        """
        C4s    = symmetryOperations.get ( "C4"   , None )
        S4s    = symmetryOperations.get ( "S4"   , None )
        sigmas = symmetryOperations.get ( "sigma", None )
        if ( C4s is not None ) and ( S4s is not None ) and ( sigmas is not None ):
            newS4s = []
            for S4 in S4s:
                axis      = S4.axis
                toInclude = True
                for C4 in C4s:
                    if ( math.fabs ( math.fabs ( C4.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                        for sigma in sigmas:
                            if ( math.fabs ( math.fabs ( sigma.normal.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                                toInclude = False
                                break
                        if not toInclude: break
                if toInclude: newS4s.append ( S4 )
            if len ( newS4s ) > 0:     symmetryOperations["S4"] = newS4s
            else:                  del symmetryOperations["S4"]

    #===============================================================================================================================
    def SetUpIrreducibleRepresentationCalculation ( self, log = logFile ):
        """Set up the information necessary for determination of irreducible representations."""
        if ( self.pointGroup is None ) or ( self.system is None ) or ( self.symmetryOperations is None ):
            raise ValueError ( "Insufficient data for determining irreducible representations." )
        elif self.characterSymmetryOperations is None:
            # . Initialization.
            characterSymmetryOperations = {}
            pointGroup                  = self.pointGroup
            # . Build CInfinity operations if required.
            if ( pointGroup.cInfinityRotations is not None ) and ( len ( pointGroup.cInfinityRotations ) > 0 ):
                cInfinity = self.symmetryOperations.get ( "CInfinity", None )
                if cInfinity is None: raise ValueError ( "Unable to find CInfinity rotation for point group: " + pointGroup.label + "." )
                for rotation in pointGroup.cInfinityRotations:
                    operation = ProperRotation ( axis = Clone ( cInfinity[0].axis ), order = int ( rotation[1:] ) )
                    if operation.EstablishSymmetryRelatedPairs ( self.elementGroups, self.system.coordinates3, self.distanceTolerance ):
                        characterSymmetryOperations[rotation] = operation
                    else: raise ValueError ( "CInfinity rotation " + rotation + " for point group " + pointGroup.label + " is not a symmetry operation." )
            # . Process existing operations.
            labels = pointGroup.characterSymmetryOperations
            # . Get sigma-h.
            sigmaH = None
            sigmas = self.symmetryOperations.get ( "sigma", [] )
            if ( pointGroup.principalAxisOperation is not None ) and ( ( "sigma-h" in labels ) or ( "sigma-m" in labels ) ):
                principalAxisOperation = self.symmetryOperations[pointGroup.principalAxisOperation][0]
                for sigma in sigmas:
                    if ( math.fabs ( math.fabs ( sigma.normal.Dot ( principalAxisOperation.axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                        sigmaH = sigma
                        break
                if ( "sigma-h" in labels ):
                    if sigmaH is None: raise ValueError ( "Unable to find \"sigma-h\" operation for point group: " + pointGroup.label + "." )
                    else: characterSymmetryOperations["sigma-h"] = sigmaH
            # . Get sigma-m.
            if ( "sigma-m" in labels ):
                sigmas = list ( sigmas )
                if sigmaH is not None: sigmas.remove ( sigmaH )
                if len ( sigmas ) > 0:
                    sigmaM = sigmas[0]
                    for sigma in sigmas[1:]:
                        if sigma.selfMappings > sigmaM.selfMappings: sigmaM = sigma
                    characterSymmetryOperations["sigma-m"] = sigmaM
                else: raise ValueError ( "Unable to find \"sigma-m\" operation for point group: " + pointGroup.label + "." )
            # . Get C2-m.
            if ( "C2-m" in labels ):
                C2s = self.symmetryOperations.get ( "C2", [] )
                if len ( C2s ) > 0:
                    C2M = C2s[0]
                    for C2 in C2s[1:]:
                        if C2.selfMappings > C2M.selfMappings: C2M = C2
                    characterSymmetryOperations["C2-m"] = C2M
                else: raise ValueError ( "Unable to find \"C2-m\" operation for point group: " + pointGroup.label + "." )
            # . Existing operations.
            for operation in pointGroup.characterSymmetryOperations:
                if operation not in characterSymmetryOperations.keys ( ):
                    inputOperations = self.symmetryOperations.get ( operation, None )
                    if inputOperations is None: raise ValueError ( "Unable to find required operation \"" + operation + "\" for point group: " + pointGroup.label + "." )
                    characterSymmetryOperations[operation] = inputOperations[0]
            # . Save.
            self.characterSymmetryOperations = [ characterSymmetryOperations[key] for key in pointGroup.characterSymmetryOperations ]

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def FindSystemPointGroup ( system, **keywordArguments ):
    """Find the point group of a system."""
    log    = keywordArguments.get ( "log", logFile )
    finder = PointGroupFinder ( **keywordArguments )
    return finder.FindSystemPointGroup ( system, log = log )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

