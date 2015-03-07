#-------------------------------------------------------------------------------
# . File      : SymmetryOperation.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Symmetry operations for point group identification."""

import math

from pCore import Clone, Integer1DArray, Matrix33

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
class SymmetryElement ( object ):
    """Base class for symmetry elements."""

    defaultAttributes = { "mapping"              : None, \
                          "order"                :    1, \
                          "selfMappings"         :    0, \
                          "transformationMatrix" : None  }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ):
	    self.__dict__[key] = value
        for ( key, value ) in keywordArguments.iteritems ( ):
	    if key in self.__class__.defaultAttributes:
	        setattr ( self, key, value )
        self.MakeTransformationMatrix ( )

    def ApplyTo ( self, point ):
        """Apply the transformation to a point (in-place)."""
	self.transformationMatrix.ApplyTo ( point )

    def EstablishSymmetryRelatedPairs ( self, elementGroups, coordinates3, tolerance ):
        """Establish the symmetry related pairs for a transformation."""
        # . Initialization.
	numberNodes  = len ( coordinates3 ) // 3
        self.mapping = Integer1DArray.WithExtent ( numberNodes )
        self.mapping.Set ( -1 )
        success      = True
        # . All atoms map to themselves.
	if self.order == 1:
            for i in range ( numberNodes ): self.mapping[i] = i
            self.selfMappings = numberNodes
        # . Build the pairs.
	else:
            for ( element, groups ) in elementGroups.iteritems ( ):
        	for group in groups:
                    isUsed = set ( )
                    for iNode in group:
                	if self.mapping[iNode] < 0:
                            bestDistance = 2.0 * tolerance
                            bestJ        = -1
                            xyz          = coordinates3.GetRow ( iNode )
                            self.ApplyTo ( xyz )
                            for jNode in group:
                        	if jNode not in isUsed:
                                    distance = math.sqrt ( ( xyz[0] - coordinates3[jNode,0] )**2 + \
                                                	   ( xyz[1] - coordinates3[jNode,1] )**2 + \
                                                	   ( xyz[2] - coordinates3[jNode,2] )**2 )
                                    if distance < bestDistance:
                                	bestJ        = jNode
                                	bestDistance = distance
                            if bestDistance > tolerance:
                        	success = False
                        	break
                            else:
                        	self.mapping[iNode] = bestJ
                        	isUsed.add ( bestJ )
            self.selfMappings = 0
            if success:
                for ( i, j ) in enumerate ( self.mapping ):
                    if ( i == j ): self.selfMappings += 1
        return success

    def IsEquivalent ( self, other ):
        """Check for equivalence of two transformations."""
	return True

    def Label ( self ): return ""

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        self.transformationMatrix = Matrix33.Identity ( )

#===================================================================================================================================
class Identity ( SymmetryElement ):
    """Identity."""

    def ApplyTo ( self ): pass

    def Label ( self ): return "E"

#===================================================================================================================================
class Inversion ( SymmetryElement ):
    """Inversion."""

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        super ( Inversion, self ).__init__ ( **keywordArguments )
	self.order = 2

    def ApplyTo ( self, point ):
        """Apply the transformation to a point (in-place)."""
	point.Scale ( -1.0 )

    def Label ( self ): return "i"

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        super ( Inversion, self ).MakeTransformationMatrix ( )
        for i in range ( 3 ): self.transformationMatrix[i,i] = -1.0

#===================================================================================================================================
class ProperRotation ( SymmetryElement ):
    """Proper rotation."""

    defaultAttributes = { "angle" :  1.0, \
                          "axis"  : None  }
    defaultAttributes.update ( SymmetryElement.defaultAttributes )

    def IsEquivalent ( self, other, tolerance ):
        """Check for equivalence of two transformations."""
	self.axis.Normalize  ( )
	other.axis.Normalize ( )
	areSame = ( self.order == other.order ) and ( math.fabs ( math.fabs ( self.axis.Dot ( other.axis ) ) - 1.0 ) < tolerance )
	return areSame

    def Label ( self ):
        if self.order == 1: tag = "inf"
        else:               tag = repr ( self.order )
        return "C" + tag

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        super ( ProperRotation, self ).MakeTransformationMatrix ( )
	if self.order > 1: self.angle = ( 2.0 * math.pi ) / float ( self.order )
	self.axis.Normalize ( )
        self.transformationMatrix.RotationAboutAxis ( self.angle, self.axis )

#===================================================================================================================================
class ImproperRotation ( ProperRotation ):
    """Improper rotation."""

    def Label ( self ): return "S" + repr ( self.order )

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        super ( ImproperRotation, self ).MakeTransformationMatrix ( )
        reflection = Matrix33.Null ( )
        reflection.Reflection ( self.axis )
        self.transformationMatrix.PreMultiplyBy ( reflection )

#===================================================================================================================================
class Reflection ( SymmetryElement ):
    """Reflection."""

    defaultAttributes = { "normal" : None }
    defaultAttributes.update ( SymmetryElement.defaultAttributes )

    def IsEquivalent ( self, other, tolerance ):
        """Check for equivalence of two transformations."""
	self.normal.Normalize  ( )
	other.normal.Normalize ( )
	areSame = ( self.order == other.order ) and ( math.fabs ( math.fabs ( self.normal.Dot ( other.normal ) ) - 1.0 ) < tolerance )
	return areSame

    def Label ( self ): return "sigma"

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        super ( Reflection, self ).MakeTransformationMatrix ( )
	self.normal.Normalize ( )
        self.transformationMatrix.Reflection ( self.normal )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
