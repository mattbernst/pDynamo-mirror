#-------------------------------------------------------------------------------
# . File      : SoftConstraint.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes for manipulating soft constraints.

This should be reimplemented in C for more efficiency.
"""

import math

from pCore import Clone, Coordinates3, RawObjectConstructor, Selection, UNITS_ANGLE_RADIANS_TO_DEGREES, Vector3

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SoftConstraintEnergyModel ( object ):
    """Base class for soft constraint energy models.

    This is the most general class and can be used directly.
    Subclasses provide specialization.
    """

    defaultAttributes = [ "highEquilibriumValue" ,
                          "highForceConstant"    ,
                          "highPower"            ,
                          "isPeriodic"           ,
                          "lowEquilibriumValue"  ,
                          "lowForceConstant"     ,
                          "lowPower"             ,
                          "period"               ]

    def __getstate__ ( self ):
        state = {}
        for key in self.__class__.defaultAttributes:
            value = getattr ( self, key, None )
            if value is not None: state[key] = value
        return state

    def __init__ ( self, lowPower, lowEquilibriumValue, lowForceConstant, highPower, highEquilibriumValue, highForceConstant, period = None ):
        """Constructor."""
        QOK = isinstance ( lowPower,  int ) and isinstance ( lowEquilibriumValue,  float ) and isinstance ( lowForceConstant,  float ) and \
              isinstance ( highPower, int ) and isinstance ( highEquilibriumValue, float ) and isinstance ( highForceConstant, float ) and \
              ( ( lowPower  > 0 ) or ( lowForceConstant  == 0.0 ) ) and ( ( highPower > 0 ) or ( highForceConstant == 0.0 ) ) and \
              ( lowEquilibriumValue <= highEquilibriumValue ) and ( lowForceConstant >= 0.0 ) and ( highForceConstant >= 0.0 )
        if QOK:
            self.highEquilibriumValue = highEquilibriumValue
            self.highForceConstant    = highForceConstant
            self.highPower            = highPower
            self.lowEquilibriumValue  = lowEquilibriumValue
            self.lowForceConstant     = lowForceConstant
            self.lowPower             = lowPower
            self.period               = period
            self.isPeriodic           = ( period is not None )
        else:
            raise TypeError ( "Invalid argument to constructor." )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        for ( key, value ) in state.iteritems ( ):
            if ( key in self.__class__.defaultAttributes ) and ( value is not None ):
                setattr ( self, key, value )

    def _Initialize ( self ):
        """Initialization."""
        for key in self.__class__.defaultAttributes: setattr ( self, key, None )

    def Energy ( self, value ):
        """Calculate an energy and a derivative given a value."""
        energy = 0.0
        dedv   = 0.0
        v      = value
        if self.period is None:
            higher = self.highEquilibriumValue
            lower  = self.lowEquilibriumValue
        else:
            # . Make sure the boundaries are OK.
            higher = self.highEquilibriumValue - math.floor ( self.highEquilibriumValue / self.period ) * self.period
            lower  = self.lowEquilibriumValue  - math.floor ( self.lowEquilibriumValue  / self.period ) * self.period
            if lower > higher: lower -= self.period
            # . Move the value just to the left of lower.
            v = v - math.floor ( v / self.period ) * self.period
            while v > lower: v -= self.period
            # . Three possible cases.
            v0 = v + self.period
            if ( ( v0 >= lower ) and ( v0 <= higher ) ) or ( ( v0 > higher ) and ( ( v0 - higher ) < ( lower - v ) ) ): v = v0
        if v < lower:
            disp   = lower - v
            dedv   = self.lowForceConstant * disp**( self.lowPower-1 )
            energy = dedv * disp
            dedv  *= - float ( self.lowPower )
        elif v > higher:
            disp   = v - higher
            dedv   = self.highForceConstant * disp**( self.highPower-1 )
            energy = dedv * disp
            dedv  *= float ( self.highPower )
        return ( energy, dedv )
        
    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

#===================================================================================================================================
# . Subclasses.
#===================================================================================================================================
class SoftConstraintEnergyModelHarmonic ( SoftConstraintEnergyModel ):
    """Define a harmonic energy model."""

    def __init__ ( self, eq, fc, period = None ):
        """Constructor."""
        super ( SoftConstraintEnergyModelHarmonic, self ).__init__ ( 2, eq, fc, 2, eq, fc, period = period )

#-----------------------------------------------------------------------------------------------------------------------------------

class SoftConstraintEnergyModelHarmonicRange ( SoftConstraintEnergyModel ):
    """Define a harmonic energy model with a range."""

    def __init__ ( self, lowEquilibriumValue, highEquilibriumValue, fc, period = None ):
        """Constructor."""
        super ( SoftConstraintEnergyModelHarmonicRange, self ).__init__ ( 2, lowEquilibriumValue, fc, 2, highEquilibriumValue, fc, period = period )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SoftConstraint ( object ):
    """Base class for soft constraint types.

    Should not be used directly."""

    defaultAttributes = [ "energyModel" ]

    def __getstate__ ( self ):
        state = {}
        for key in self.__class__.defaultAttributes:
            value = getattr ( self, key, None )
            if value is not None: state[key] = value
        return state

    def __init__ ( self, energyModel ):
        """Constructor."""
        if isinstance ( energyModel, SoftConstraintEnergyModel ):
            period = self.Period ( )
            QCLONE = ( ( period is None ) and energyModel.isPeriodic ) or ( ( period is not None ) and not energyModel.isPeriodic )
            if QCLONE: self.energyModel = Clone ( energyModel )
            else:      self.energyModel =         energyModel
            self.energyModel.period     = period
            self.energyModel.isPeriodic = ( period is not None )
        else:
            raise TypeError ( "Invalid energy model." )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        for ( key, value ) in state.iteritems ( ):
            if ( key in self.__class__.defaultAttributes ) and ( value is not None ):
                setattr ( self, key, value )

    def _Initialize ( self ):
        """Initialization."""
        for key in self.__class__.defaultAttributes: setattr ( self, key, None )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        return ( 0.0, 0.0 )

    def Period ( self ):
        """Return the period of the constraint variable. None implies no periodicity."""
        return None

    def Label ( self ):
        """Return a label for the constraint."""
        return "Undefined"

    def Merge ( self, indexIncrement ):
        """Merging."""
        return None

    def Prune ( self, selection ):
        """Pruning."""
        return None
        
    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

#===================================================================================================================================
# . Subclasses.
#===================================================================================================================================
class SoftConstraintAngleDotProduct ( SoftConstraint ):
    """An angle dot product soft constraint.

    The dot product of the angle is constrained so the reference value
    must be between -1 and 1.
    """

    defaultAttributes = [ "point1", "point2", "point3" ]
    defaultAttributes.extend ( SoftConstraint.defaultAttributes )

    def __init__ ( self, point1, point2, point3, energyModel ):
        """Constructor."""
        super ( SoftConstraintAngleDotProduct, self ).__init__ ( energyModel )
        if isinstance ( point1, int ) and isinstance ( point2, int ) and isinstance ( point3, int ):
            self.point1 = point1
            self.point2 = point2
            self.point3 = point3
        else:
            raise TypeError ( "Invalid points." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        v12 = coordinates3.Displacement ( self.point1, self.point2 ) ; r12 = v12.Norm2 ( )
        v32 = coordinates3.Displacement ( self.point3, self.point2 ) ; r32 = v32.Norm2 ( )
        v12.Scale ( 1.0 / r12 )
        v32.Scale ( 1.0 / r32 )
        value = v12.Dot ( v32 )
        if   value >  1.0: value =  1.0
        elif value < -1.0: value = -1.0
        ( f, df ) = self.energyModel.Energy ( value )
        if gradients3 is not None:
            v32.AddScaledVector3 ( - value, v12 ) # . This is s12 * d(costheta)/d(p1) in d32.
            v12.Scale ( 1.0 - value**2 )
            v12.AddScaledVector3 ( - value, v32 ) # . This is s32 * d(costheta)/d(p3) in d12.
            v32.Scale ( df / r12 )
            v12.Scale ( df / r32 )
            gradients3.Increment ( self.point1, v32 )
            gradients3.Decrement ( self.point2, v12 )
            gradients3.Decrement ( self.point2, v32 )
            gradients3.Increment ( self.point3, v12 )
        return ( f, value )

    def Label ( self ): return "Angle Dot Product"

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.point1 += indexIncrement
        new.point2 += indexIncrement
        new.point3 += indexIncrement
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        if ( self.point1 in selection ) and \
           ( self.point2 in selection ) and \
           ( self.point3 in selection ):
            pruned = self.__class__ ( selection.Position ( self.point1 ), \
                                      selection.Position ( self.point2 ), \
                                      selection.Position ( self.point3 ), Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class SoftConstraintDihedral ( SoftConstraint ):
    """A dihedral angle soft constraint."""

    defaultAttributes = [ "point1", "point2", "point3", "point4" ]
    defaultAttributes.extend ( SoftConstraint.defaultAttributes )

    def __init__ ( self, point1, point2, point3, point4, energyModel ):
        """Constructor."""
        super ( SoftConstraintDihedral, self ).__init__ ( energyModel )
        if isinstance ( point1, int ) and isinstance ( point2, int ) and isinstance ( point3, int ) and isinstance ( point4, int ):
            self.point1 = point1
            self.point2 = point2
            self.point3 = point3
            self.point4 = point4
        else:
            raise TypeError ( "Invalid points." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        # . Displacements.
        v12 = coordinates3.Displacement ( self.point1, self.point2 )
        v32 = coordinates3.Displacement ( self.point3, self.point2 )
        v34 = coordinates3.Displacement ( self.point3, self.point4 )
        # . Get m and n.
        m = Clone ( v12 ) ; m.Cross ( v32 )
        n = Clone ( v32 ) ; n.Cross ( v34 )
        # . Get the sizes of m and n.
        msize = m.Norm2 ( )
        nsize = n.Norm2 ( )
        # . Normalize m and n.
        m.Scale ( 1.0 / msize )
        n.Scale ( 1.0 / nsize )
        # . Get the dot-product.
        dotfac = m.Dot ( n )
        if    dotfac > 1.0: dotfac =  1.0
        elif dotfac < -1.0: dotfac = -1.0
        # . Get the sign of the angle.
        sgnfac = 1.0
        if v12.Dot ( n ) < 0.0: sgnfac = -1.0
        # . Determine the dihedral.
        value = sgnfac * math.acos ( dotfac ) * UNITS_ANGLE_RADIANS_TO_DEGREES
        # . Get the energy.
        ( f, df ) = self.energyModel.Energy ( value )
        # . Derivatives.
        if gradients3 is not None:
            df *= UNITS_ANGLE_RADIANS_TO_DEGREES
            # . Calculate r32.
            r32 = v32.Norm2 ( )
            # . Calculate dedi and dedl in m and n respectively.
            m.Scale (   df * r32 / msize )
            n.Scale ( - df * r32 / nsize )
            # . Calculate some additional factors.
            fact12 = v12.Dot ( v32 ) / ( r32 * r32 )
            fact34 = v34.Dot ( v32 ) / ( r32 * r32 )
            # . Gradients for i and l.
            gradients3.Increment ( self.point1, m )
            gradients3.Increment ( self.point4, n )
            # . Calculate dedj and dedk in v12 and v32 respectively.
            m.CopyTo ( v12 ) ; v12.Scale ( fact12 - 1.0 ) ; v12.AddScaledVector3 ( -fact34, n )
            n.CopyTo ( v32 ) ; v32.Scale ( fact34 - 1.0 ) ; v32.AddScaledVector3 ( -fact12, m )
            # . calculate the gradients.
            gradients3.Increment ( self.point2, v12 )
            gradients3.Increment ( self.point3, v32 )
        return ( f, value )

    def Label ( self ): return "Dihedral"

    def Period ( self ): return 360.0

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.point1 += indexIncrement
        new.point2 += indexIncrement
        new.point3 += indexIncrement
        new.point4 += indexIncrement
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        if ( self.point1 in selection ) and \
           ( self.point2 in selection ) and \
           ( self.point3 in selection ) and \
           ( self.point4 in selection ):
            pruned = self.__class__ ( selection.Position ( self.point1 ), \
                                      selection.Position ( self.point2 ), \
                                      selection.Position ( self.point3 ), \
                                      selection.Position ( self.point4 ), Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class SoftConstraintDistance ( SoftConstraint ):
    """A distance soft constraint."""

    defaultAttributes = [ "point1", "point2" ]
    defaultAttributes.extend ( SoftConstraint.defaultAttributes )

    def __init__ ( self, point1, point2, energyModel ):
        """Constructor."""
        super ( SoftConstraintDistance, self ).__init__ ( energyModel )
        if isinstance ( point1, int ) and isinstance ( point2, int ):
            self.point1 = point1
            self.point2 = point2
        else:
            raise TypeError ( "Invalid points." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        v12 = coordinates3.Displacement ( self.point1, self.point2 )
        r12 = v12.Norm2 ( )
        ( f, df ) = self.energyModel.Energy ( r12 )
        if gradients3 is not None:
            v12.Scale ( df / r12 )
            gradients3.Increment ( self.point1, v12 )
            gradients3.Decrement ( self.point2, v12 )
        return ( f, r12 )

    def Label ( self ): return "Distance"

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.point1 += indexIncrement
        new.point2 += indexIncrement
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        if ( self.point1 in selection ) and \
           ( self.point2 in selection ):
            pruned = self.__class__ ( selection.Position ( self.point1 ), \
                                      selection.Position ( self.point2 ), Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class SoftConstraintMultipleDistance ( SoftConstraint ):
    """A multiple distance soft constraint."""

    defaultAttributes = [ "distances" ]
    defaultAttributes.extend ( SoftConstraint.defaultAttributes )

    def __init__ ( self, distances, energyModel ):
        """Constructor."""
        super ( SoftConstraintMultipleDistance, self ).__init__ ( energyModel )
        try:
            self.distances = []
            for ( point1, point2, weight ) in distances:
                if isinstance ( point1, int ) and isinstance ( point2, int ) and isinstance ( weight, float ):
                    self.distances.append ( ( point1, point2, weight ) )
                else:
                    raise TypeError ( "Invalid distance parameters." )
        except:
            raise TypeError ( "Invalid distance specification." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        value = 0.0
        temp  = []
        for ( point1, point2, weight ) in self.distances:
            v = coordinates3.Displacement ( point1, point2 )
            r = v.Norm2 ( )
            value += weight * r
            temp.append ( ( v, r ) )
        ( f, df ) = self.energyModel.Energy ( value )
        if gradients3 is not None:
            for ( ( point1, point2, weight ), ( v, r ) ) in zip ( self.distances, temp ):
                v.Scale ( df * weight / r )
                gradients3.Increment ( point1, v )
                gradients3.Decrement ( point2, v )
        return ( f, value )

    def Label ( self ): return "Multiple Distance"

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = None
        # . Loop over distances.
        distances = []
        for ( point1, point2, weight ) in self.distances:
            distances.append ( ( point1 + indexIncrement, point2 + indexIncrement, weight ) )
        # . Create object.
        if len ( distances ) > 0:
            new = self.__class__ ( distances, Clone ( self.energyModel ) )
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        # . Loop over distances.
        reduced = []
        for ( point1, point2, weight ) in self.distances:
            if ( point1 in selection ) and ( point2 in selection ):
                reduced.append ( ( selection.Position ( point1 ), selection.Position ( point2 ), weight ) )
        # . Create object.
        if len ( reduced ) > 0:
            pruned = self.__class__ ( reduced, Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class SoftConstraintMultipleTether ( SoftConstraint ):
    """A multiple tether soft constraint."""

    defaultAttributes = [ "reference", "selection" ]
    defaultAttributes.extend ( SoftConstraint.defaultAttributes )

    def __init__ ( self, selection, reference, energyModel ):
        """Constructor."""
        super ( SoftConstraintMultipleTether, self ).__init__ ( energyModel )
        if isinstance ( selection, Selection ) and isinstance ( reference, Coordinates3 ):
            self.selection = selection
            self.reference = reference
        else:
            raise TypeError ( "Invalid arguments." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        fTotal = 0.0
        for i in self.selection:
            v12 = coordinates3.GetRow ( i )
            v12.AddScaledVector3 ( -1.0, self.reference.GetRow ( i ) )
            r12 = v12.Norm2 ( )
            ( f, dF ) = self.energyModel.Energy ( r12 )
            fTotal += f
            if ( gradients3 is not None ) and ( r12 != 0.0 ):
                v12.Scale ( dF / r12 )
                gradients3.Increment ( i, v12 )
        return ( fTotal, None )

    def Label ( self ): return "Multiple Tether"

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.selection.Increment ( indexIncrement )
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        # . Get reduced selection.
        reduced = self.selection.Prune ( selection )
        if len ( reduced ) > 0:
            pruned = self.__class__ ( reduced, self.reference.Prune ( selection ), Clone ( self.energyModel ) )
        return pruned

#-----------------------------------------------------------------------------------------------------------------------------------

class SoftConstraintTether ( SoftConstraint ):
    """A tether soft constraint."""

    defaultAttributes = [ "origin", "point" ]
    defaultAttributes.extend ( SoftConstraint.defaultAttributes )

    def __init__ ( self, point, origin, energyModel ):
        """Constructor."""
        super ( SoftConstraintTether, self ).__init__ ( energyModel )
        if isinstance ( point, int ) and isinstance ( origin, Vector3 ):
            self.origin = origin
            self.point  = point
        else:
            raise TypeError ( "Invalid points." )

    def Energy ( self, coordinates3, gradients3 = None ):
        """Calculate the energy."""
        v12 = coordinates3.GetRow ( self.point )
        v12.AddScaledVector3 ( -1.0, self.origin )
        r12 = v12.Norm2 ( )
        ( f, df ) = self.energyModel.Energy ( r12 )
        if ( gradients3 is not None ) and ( r12 != 0.0 ):
            v12.Scale ( df / r12 )
            gradients3.Increment ( self.point, v12 )
        return ( f, r12 )

    def Label ( self ): return "Tether"

    def Merge ( self, indexIncrement ):
        """Merging."""
        new = Clone ( self )
        new.point += indexIncrement
        return new

    def Prune ( self, selection ):
        """Pruning."""
        pruned = None
        if self.point in selection:
            pruned = self.__class__ ( selection.Position ( self.point ), Clone ( self.origin ), Clone ( self.energyModel ) )
        return pruned

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
