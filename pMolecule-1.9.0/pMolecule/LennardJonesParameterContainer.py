#-------------------------------------------------------------------------------
# . File      : LennardJonesParameterContainer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Lennard-Jones parameter classes."""

from LJParameterContainer import LJParameterContainer
from pCore                import UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LennardJonesParameterContainer ( object ):
    """A container for Lennard-Jones parameters."""

    defaultAttributes = { "items"               :            None,
                          "properties"          :            None,
                          "label"               :            None,
                          "rawItems"            :            None,
                          "style"               :            None,
                          "termLabel"           : "Lennard-Jones",
                          "useStrictAssignment" :           True }

    #yaml_tag = "!LennardJonesParameterContainer"

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Columns.
        mapping["Parameter Fields"] = [ "Atom Type", "Epsilon", "Sigma" ]
        # . Rows.
        keys = self.rawItems.keys ( )
        keys.sort ( )
        rows = []
        for key in keys:
            ( epsilon, sigma ) = self.rawItems[key]
            rows.append ( [ key, epsilon, sigma ] )
        mapping["Parameter Values"] = rows
        # . Other data.
        if self.label is not None: mapping["Label"] = self.label
        if self.style is not None: mapping["Style"] = self.style
        for key in ( "Analytic Form", "Units" ):
            if key in self.properties:
                mapping[key] = self.properties[key]
        return mapping

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Basic construction.
            self.__init__ ( )
            columns         = mapping.pop ( "Parameter Fields" )
            rows            = mapping.pop ( "Parameter Values" )
            self.label      = mapping.pop ( "Label", None )
            self.style      = mapping.pop ( "Style", None )
            self.properties = dict ( mapping )
            # . Get token indices.
            l = columns.index ( "Atom Type" )
            e = columns.index ( "Epsilon"   )
            s = columns.index ( "Sigma"     )
            # . Create the raw items.
            self.rawItems = {}
            for row in rows:
                key = row[l]
                self.rawItems[key] = ( float ( row[e] ), float ( row[s] ) )
            if len ( self.rawItems ) != len ( rows ): raise
            # . Create the items that are to be used.
            self.CheckStyle      ( )
            self.ProcessRawItems ( )
#        except Exception as e:
#            print e[0]
        except:
            raise ValueError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def CheckStyle ( self ):
        """Check the style setting."""
        if ( self.style != "AMBER" ) and ( self.style != "OPLS" ):
            raise ValueError ( "Invalid Lennard-Jones parameter style: {:s}.".format ( self.style ) )

    @staticmethod
    def MakeKey ( label1, label2 ):
        """Make a key."""
        if ( label1 >= label2 ): return ( label1, label2 )
        else:                    return ( label2, label1 )

    def MakeParameterContainer ( self, atomTypeLabels ):
        """Make a LJ parameter container."""
        missingParameters = set ( )
        mm                = None
        if len ( atomTypeLabels ) > 0:
            epsilons      = []
            parameterKeys = []
	    sigmas        = []
            for key in atomTypeLabels:
                parameters = self.items.get ( key, None )
                if parameters is None:
                    if self.useStrictAssignment:
                        missingParameters.add ( ( self.termLabel, key ) )
                else:
                    ( epsilon, sigma ) = parameters
	            epsilons.append      ( epsilon )
                    parameterKeys.append ( key     )
                    sigmas.append        ( sigma   )
            numberOfParameters = len ( epsilons )
            if ( numberOfParameters > 0 ) and ( len ( missingParameters ) <= 0 ):
                state = { "label"         : self.termLabel    ,
                          "analyticForm"  : self.style        ,
                          "epsilons"      : epsilons          ,
                          "numberOfTypes" : numberOfParameters,
                          "parameterKeys" : parameterKeys     ,
                          "sigmas"        : sigmas            }
                mm = LJParameterContainer.Raw ( )
                mm.__setstate__ ( state )
        return ( mm, missingParameters )

    def ProcessRawItems ( self ):
        """Process the raw items."""
        # . Create the items that are to be used and convert to the correct units if necessary.
        units = self.properties.get ( "Units", None )
        if units is not None:
            if units.get ( "Epsilon", "" ).find ( "kcal" ) > -1:
                self.items = {}
                for ( key, ( epsilon, sigma ) ) in self.rawItems.iteritems ( ):
                    self.items[key] = ( UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE * epsilon, sigma )
            else:
                self.items = self.rawItems

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

    def ScaleEnergies ( self, scale ):
        """Scale the energy parameters."""
        if scale != 1.0:
            for ( key, ( epsilon, sigma ) ) in self.items.iteritems ( ):
                self.items[key] = ( scale * epsilon, sigma )

    def UpdateParameters ( self, other ):
        """Update the parameters with additional values from another container."""
        if self.items is None:
            self.items = {}
        if other.items is not None:
            self.items.update ( other.items )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
