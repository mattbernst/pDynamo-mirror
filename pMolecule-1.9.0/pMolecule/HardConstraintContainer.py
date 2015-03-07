#-------------------------------------------------------------------------------
# . File      : HardConstraintContainer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines hard constraints for a system.

Hard constraints are of three types currently:
(i)   constraints that fix the positions of atoms.
(ii)  constraints that depend linearly on the atom's Cartesian coordinates -
      most notably those restricting rotations and translations.
(iii) constraints that linearly constrain the QC charges and/or spins of a system.

For simplicity, fixed atom constraints are incompatible with rotation/translation
constraints. However, other types of linear constraint may be applicable.
"""

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class HardConstraintContainer ( object ):
    """Define hard constraints."""

    defaultAttributes = { "_chargeConstraints" : None ,
                          "_fixedAtoms"        : None }

    def __getstate__ ( self ):
        state = {}
        for key in self.__class__.defaultAttributes:
            value = self.__dict__.get ( key, None )
            if value is not None:
                label        = key[1:]
                state[label] = value
        return state

    def __init__ ( self ):
        """Constructor."""
        self.__dict__.update ( self.__class__.defaultAttributes )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        for key in self.__class__.defaultAttributes:
            label = key[1:]
            value = state.get ( label, None )
            if value is not None: self.__dict__[key] = value

    def ClearChargeConstraints ( self ):
        """Clear charge constraints from the container."""
        self.__dict__.pop ( "_chargeConstraints", None )

    def ClearFixedAtoms ( self ):
        """Clear fixed atoms from the container."""
        self.__dict__.pop ( "_fixedAtoms", None )

    def DefineChargeConstraints ( self, chargeConstraints ):
        """Define charge constraints."""
        self.__dict__["_chargeConstraints"] = chargeConstraints

    def DefineFixedAtoms ( self, selection ):
        """Define fixed atoms."""
        self.__dict__["_fixedAtoms"] = selection

    def IsEmpty ( self ):
        """Is the container empty?"""
        isEmpty = True
        for key in self.__class__.defaultAttributes:
            if self.__dict__.get ( key, None ) is not None:
                isEmpty = False
                break
        return isEmpty

    def NumberOfChargeConstraints ( self ):
        """Return the number of charge constraints."""
        if self._chargeConstraints is None: return 0
        else:                               return len ( self._chargeConstraints )

    def NumberOfFixedAtoms ( self ):
        """Return the number of fixed atoms."""
        if self._fixedAtoms is None: return 0
        else:                        return len ( self._fixedAtoms )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( pageWidth = 90, valueWidth = 12 )
            summary.Start ( "Hard Constraint Container Summary" )
            summary.Entry ( "Number of Charge Constraints", "{:d}".format ( self.NumberOfChargeConstraints ( ) ) )
            summary.Entry ( "Number of Fixed Atoms"       , "{:d}".format ( self.NumberOfFixedAtoms        ( ) ) )
            summary.Stop ( )

    @property
    def chargeConstraints ( self ): return self.__dict__["_chargeConstraints"]

    @property
    def fixedAtoms ( self ): return self.__dict__["_fixedAtoms"]

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
