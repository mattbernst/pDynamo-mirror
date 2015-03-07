#-------------------------------------------------------------------------------
# . File      : MMModel.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines MM energy model classes."""

import os, os.path

from pCore               import logFile, LogFileActive, RawObjectConstructor
from MMAtomTyper         import MMAtomTyper
from MMModelError        import MMModelError
from MMParameterAssigner import MMParameterAssigner

# . Permit here configuration file for default parameters?
# . Must be order: defaults, configuration, arguments as latter override everything else.
# . Only process if parameterSet specified.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_ForceFieldPath = "forceFields"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModel ( object ):
    """Define an MM energy model."""

    defaultAttributes = { "electrostaticScale14" :       1.0 ,
                          "fullPath"             :      None ,
                          "label"                : "Generic" ,
                          "lennardJonesScale14"  :       1.0 ,
                          "lennardJonesStyle"    :    "OPLS" ,
                          "parameterSet"         :      None ,
                          "path"                 :      None }

    def __getstate__ ( self ):
        state = {}
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ):
            state[key] = value
        return state

    def __init__ ( self, *arguments, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )
        # . Parameter set.
        if ( len ( arguments ) > 0 ) and isinstance ( arguments[0], basestring ): self.parameterSet = arguments[0]
        # . Paths, etc.
        self.SetPaths ( )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ): self.__dict__.update ( state )

    def BuildModel ( self, connectivity, energyModel, sequence, log = logFile ):
        """Build the model."""
        # . Check path.
        if ( self.fullPath is None ) or not os.path.exists ( self.fullPath ):
            raise MMModelError ( "Unable to build " + self.label + " MM Model as parameter set missing or unspecified." )
        # . Type atoms.
        typer = MMAtomTyper ( self.fullPath )
        ( atomTypes, atomCharges ) = typer.TypeAtoms ( connectivity, sequence, log )
        # . Assign parameters.
        assigner = MMParameterAssigner ( self.fullPath, lennardJonesScale14 = self.lennardJonesScale14 )
        assigner.AssignParameters ( connectivity, atomTypes, atomCharges, energyModel, log )

    def ClearModelBuildingData ( self ):
        """Clear all data associated with model building."""
        pass

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass ( )
        return self

    def SetPaths ( self ):
        """Determine the paths for finding the parameter set."""
        if self.parameterSet is not None:
            if self.path is None:
                self.path = os.path.join ( os.getenv ( "PDYNAMO_PARAMETERS" ), _ForceFieldPath, self.label.lower () )
            if self.path is not None:
                self.fullPath = os.path.join ( self.path, self.parameterSet )

    def SummaryEntry ( self, summary ):
        """Summary entry."""
        if summary is not None:
            for ( key, item ) in self.__dict__.iteritems ( ):
                if   key == "electrostaticScale14": summary.Entry ( "El. 1-4 Scaling", "{:.3f}".format ( item ) )
                elif key == "lennardJonesScale14":  summary.Entry ( "LJ 1-4 Scaling",  "{:.3f}".format ( item ) )
                elif hasattr ( item, "SummaryEntry" ): item.SummaryEntry ( summary )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModelAMBER ( MMModel ):
    """Define an AMBER MM energy model."""

    defaultAttributes = dict ( MMModel.defaultAttributes )
    defaultAttributes.update ( { "electrostaticScale14" : ( 1.0 / 1.2 ) ,
                                 "label"                :      "AMBER"  ,
                                 "lennardJonesScale14"  : ( 1.0 / 2.0 ) ,
                                 "lennardJonesStyle"    :      "AMBER"  } )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModelCHARMM ( MMModel ):
    """Define a CHARMM MM energy model."""

    defaultAttributes = dict ( MMModel.defaultAttributes )
    defaultAttributes.update ( { "label"             : "CHARMM",
                                 "lennardJonesStyle" : "AMBER" } )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModelOPLS ( MMModel ):
    """Define an OPLS MM energy model."""

    defaultAttributes = dict ( MMModel.defaultAttributes )
    defaultAttributes.update ( { "electrostaticScale14" :    0.5 ,
                                 "label"                : "OPLS" ,
                                 "lennardJonesScale14"  :    0.5 ,
                                 "lennardJonesStyle"    : "OPLS" } )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
