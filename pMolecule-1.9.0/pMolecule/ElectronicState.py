#-------------------------------------------------------------------------------
# . File      : ElectronicState.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Define the electronic state of the system.

This module currently defines only the overall charge and multiplicity of the quantum chemical part of a system.
"""

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Data.
#===================================================================================================================================
# . Allowed multiplicity names.
_MultiplicityNames = { "SINGLET" : 1, "DOUBLET" : 2, "TRIPLET" : 3, "QUARTET" : 4, "QUINTET" : 5, "SEXTET" : 6, "SEPTET" : 7 }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ElectronicState ( object ):

    defaultAttributes = { "_charge"       : 0 ,
                          "_multiplicity" : 1 }

    def __getstate__ ( self ):
        return { "charge" : self._charge, "multiplicity" : self._multiplicity }

    def __init__ ( self, charge = 0, multiplicity = 1 ):
        """Constructor."""
        # . Defaults.
        self._Initialize ( )
        # . Charge.
        try:    self.__dict__["_charge"] = int ( charge )
        except: raise TypeError ( "Invalid charge argument: {!r}.".format ( charge ) )
        # . Multiplicity.
        if isinstance ( multiplicity, basestring ): value = _MultiplicityNames.get ( multiplicity.upper ( ), multiplicity )
        else:                                       value = multiplicity
        try:    self.__dict__["_multiplicity"] = int ( value )
        except: raise TypeError ( "Invalid multiplicity argument: {!r}.".format ( multiplicity ) )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        self._charge       = state.get ( "charge"      , 0 )
        self._multiplicity = state.get ( "multiplicity", 1 )

    def _Initialize ( self ):
        """Initialization."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        
    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "Electronic State" )
            summary.Entry ( "Charge"       , "{:d}".format ( self._charge       ) )
            summary.Entry ( "Multiplicity" , "{:d}".format ( self._multiplicity ) )
            summary.Stop ( )

    # . Properties.
    @property
    def charge       ( self ): return self.__dict__["_charge"      ]
    @property
    def multiplicity ( self ): return self.__dict__["_multiplicity"]

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
