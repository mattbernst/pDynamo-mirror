#-------------------------------------------------------------------------------
# . File      : MMPatternContainer.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MM patterns."""

from ConnectivityPattern import AtomPattern, BondPattern, ConnectivityPattern
from MMModelError        import MMModelError

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The character indicating that a pattern attribute is not defined.
_Undefined = "."

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMPattern ( ConnectivityPattern ):
    """An MM pattern used for assigning MM atom types."""

    defaultAttributes = dict ( ConnectivityPattern.defaultAttributes )
    defaultAttributes.update ( { "atomTypeLabels" : None,
                                 "atomTypes"      : None,
                                 "charges"        : None } )

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    @classmethod
    def FromMapping ( selfClass, mapping, atomAttributes = None, bondAttributes = None ):
        """Constructor from a mapping."""
        # . Basic construction.
        self         = selfClass ( )
        atomPatterns = mapping.get ( "Atom Patterns", [] )
        bondPatterns = mapping.get ( "Bond Patterns", [] )
        self.label   = mapping.get ( "Label", None )
        # . Charges and atom type labels.
        try:
            c       = atomAttributes.index ( "charge" )
            charges = []
            for pattern in atomPatterns:
                q = pattern[c]
                if q != _Undefined:
                    charges.append ( q )
                    pattern[c] = _Undefined
            self.charges = charges
        except:
            pass
        try:
            l      = atomAttributes.index ( "type" )
            labels = []
            for pattern in atomPatterns:
                label = pattern[l]
                if label != _Undefined:
                    labels.append ( label )
                    pattern[l] = _Undefined
            self.atomTypeLabels = labels
        except:
            pass
        # . Atom keys.
        try:
            k    = atomAttributes.index ( "key" )
            keys = {}
            for ( i, pattern ) in enumerate ( atomPatterns ):
                key = pattern[k]
                if key != _Undefined:
                    keys[key]  = i
                    pattern[k] = _Undefined
            k1 = bondAttributes.index ( "atomKey1" )
            k2 = bondAttributes.index ( "atomKey2" )
            for pattern in bondPatterns:
                key1 = pattern[k1]
                key2 = pattern[k2]
                pattern[k1] = keys[key1]
                pattern[k2] = keys[key2]
        except:
            pass
        # . Atom patterns.
        patterns = []
        for patternList in atomPatterns:
            patterns.append ( AtomPattern.FromList ( patternList, atomAttributes, _Undefined ) )
        self.atomPatterns = patterns
        # . Bond patterns.
        patterns = []
        for patternList in bondPatterns:
            patterns.append ( BondPattern.FromList ( patternList, bondAttributes, _Undefined ) )
        self.bondPatterns = patterns
        return self

    def ToMapping ( self, atomAttributes = None, bondAttributes = None ):
        """Get a mapping representation of the pattern."""
        mapping = {}
        # . Atom attribute indices.
        try   : c = atomAttributes.index ( "charge" )
        except: c = -1
        try   : k = atomAttributes.index ( "key"    )
        except: k = -1
        try   : l = atomAttributes.index ( "type"   )
        except: l = -1
        # . Atom patterns.
        atomPatterns = {}
        for ( i, pattern ) in enumerate ( self.atomPatterns ):
            patternList = pattern.ToList ( atomAttributes, _Undefined )
            if c >= 0: patternList[c] = self.charges        [i]
            if k >= 0: patternList[k] = i
            if l >= 0: patternList[l] = self.atomTypesLabels[i]
            atomPatterns.append ( patternList )
        # . Bond patterns.
        bondPatterns = {}
        for pattern in self.bondPatterns:
            patternList = pattern.ToList ( bondAttributes, _Undefined )
            bondPatterns.append ( patternList )
        # . Finish up.
        mapping["Atom Patterns"] = atomPatterns
        mapping["Bond Patterns"] = bondPatterns
        if self.label is not None: mapping["Label"] = self.label
        return mapping

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMPatternContainer ( object ):
    """A container for MM patterns."""

    defaultAttributes = { "atomAttributes" :          None,
                          "atomFields"     :          None,
                          "bondAttributes" :          None,
                          "bondFields"     :          None,
                          "items"          :          None,
                          "mmAtomTypes"    :          None,
                          "properties"     :          None,
                          "label"          :          None,
                          "rawItems"       :          None,
                          "termLabel"      : "MM Pattern" }

    #yaml_tag = "!MMPatternContainer"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Keys.
        mapping["Atom Fields"] = self.atomFields
        mapping["Bond Fields"] = self.bondFields
        # . Patterns.
        patterns = []
        for pattern in self.rawItems:
            patterns.append ( pattern.ToMapping ( atomAttributes = self.atomAttributes, bondAttributes = self.bondAttributes ) )
        mapping["Patterns"] = patterns
        # . Other data.
        if self.label is not None: mapping["Label"] = self.label
        return mapping

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Basic construction.
            self.__init__ ( )
            self.atomFields = mapping.pop ( "Atom Fields"       )
            self.bondFields = mapping.pop ( "Bond Fields", None )
            patterns        = mapping.pop ( "Patterns"          )
            self.label      = mapping.pop ( "Label"      , None )
            self.properties = dict ( mapping )
            # . Convert keys to attributes.
            self.atomAttributes = self.FieldsToAttributes ( self.atomFields )
            self.bondAttributes = self.FieldsToAttributes ( self.bondFields )
            # . Create the raw items.
            self.rawItems = []
            for mapping in patterns:
                self.rawItems.append ( MMPattern.FromMapping ( mapping, atomAttributes = self.atomAttributes, bondAttributes = self.bondAttributes ) )
            # . Create the items that are to be used.
            self.ProcessRawItems ( )
#        except Exception as e:
#            print e[0]
        except:
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def FieldsToAttributes ( self, fields ):
        """Convert fields to attributes."""
        attributes = []
        if fields is not None:
            for field in fields:
                attribute = "".join ( field.split ( ) )
                attributes.append ( attribute[0:1].lower ( ) + attribute[1:] )
        return attributes

    def IndexAtomTypes ( self, mmAtomTypes ):
        """Index the atom types in the container's patterns."""
        if ( mmAtomTypes is not None ) and ( self.mmAtomTypes is None ):
            for item in self.items:
                # . Types.
                item.atomTypes = []
                for label in item.atomTypeLabels:
                    atomType = mmAtomTypes.GetItem ( label )
                    if ( atomType is None ): raise MMModelError ( "Unknown atom type in MM pattern: " + label + "." )
                    else:                    item.atomTypes.append ( atomType )
                # . Charges.
                if ( item.charges is None ) or ( len ( item.charges ) == 0 ):
                    item.charges = []
                    for atomType in item.atomTypes:
                        item.charges.append ( atomType.charge )
            # . Finish up.
            self.mmAtomTypes = mmAtomTypes

    def ProcessRawItems ( self ):
        """Process the raw items."""
        self.items = self.rawItems

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

    def TypeConnectivity ( self, connectivity, atomTypes, atomCharges, untypedAtoms ):
        """Type the atoms in a connectivity."""
        if self.mmAtomTypes is not None:
            for pattern in self.items:
                pattern.MakeConnections ( )
                matches = pattern.FindAllMatches ( connectivity, selection = untypedAtoms )
                # . Maybe here should check that there are not intersecting matches within the same call (selection takes care of the others).
                # . I.e. an atom has been matched more than once and with different atom types and/or charges - in this case which definition to take?
                # . Easiest if raise error. No problem if type/charge same.
                # . THIS NEEDS TO BE DONE.
                if ( matches is not None ) and ( len ( matches ) > 0 ):
                    matchedAtoms = []
                    for match in matches:
                        for ( atomType, charge, i ) in zip ( pattern.atomTypes, pattern.charges, match ):
                            if ( atomType is not None ) and ( atomTypes[i] is None ):
                                atomCharges[i] = charge
                                atomTypes  [i] = atomType.label
                                matchedAtoms.append ( i )
                                hydrogenType   = atomType.hydrogenType
                                if hydrogenType is not None:
                                    hCharge = hydrogenType.charge
                                    for h in connectivity.GetConnectedHydrogens ( i ):
                                        atomCharges[h]  = hCharge
                                        atomCharges[i] -= hCharge
                                        atomTypes  [h]  = hydrogenType.label
                                        matchedAtoms.append ( h )
                    untypedAtoms.Exclude ( matchedAtoms )
                    if len ( untypedAtoms ) <= 0: break

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
