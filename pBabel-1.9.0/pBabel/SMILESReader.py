#-------------------------------------------------------------------------------
# . File      : SMILESReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for reading SMILES strings.

Chirality is parsed but not handled.
"""

import string

from pCore   import logFile
from pMolecule import DoubleBond, NullBond, SingleBond, TripleBond, PeriodicTable

from SMILESUtilities import ELEMENTSAROMATIC, ELEMENTSORGANIC, SMILESAtom, SMILESBond, SMILESConnectivity, SMILESConnectivityError

#===============================================================================
# . Parameters.
#===============================================================================
# . Tokens used for parsing Smiles strings.
# . Bonds - no aromaticity distinction is made here.
_BONDTOKENS = { "."  : NullBond   ( ),
                "-"  : SingleBond ( ),
                "/"  : SingleBond ( ),
                "\\" : SingleBond ( ),
                ":"  : SingleBond ( ),
                "="  : DoubleBond ( ),
                "#"  : TripleBond ( ) }

# . Charge.
_CHARGETOKENS = { "+" :  1, "++" :  2, "+++" :  3, "++++" :  4, "-" : -1, "--" : -2, "---" : -3, "----" : -4 }

# . Chirality.
_CHIRALITYFULL    = { "AL" : "AL", "OH" : "OH", "SP" : "SP", "TB" : "TB", "TH" : "TH" }
_CHIRALITYREDUCED = { "@" : ( "??", 1 ), "@@" : ( "??", 2 ), "@@@" : ( "??", 3 ), "@@@@" : ( "??", 4 ) }

# . Element data.
# . All element tokens.
_ELEMENTTOKENS = { }
for atomicNumber in PeriodicTable.AtomicNumbers ( ):
    _ELEMENTTOKENS[PeriodicTable.Symbol ( atomicNumber )] = ( atomicNumber, False )
for atomicNumber in ELEMENTSAROMATIC:
    _ELEMENTTOKENS[string.lower ( PeriodicTable.Symbol ( atomicNumber ) )] = ( atomicNumber, True )

# . Tokens for elements that can take a reduced representation (outside of square brackets).
_REDUCEDELEMENTTOKENS = { }
for atomicNumber in ELEMENTSORGANIC:
    _REDUCEDELEMENTTOKENS[PeriodicTable.Symbol ( atomicNumber )] = ( atomicNumber, False )
    if atomicNumber in ELEMENTSAROMATIC:
        _REDUCEDELEMENTTOKENS[string.lower ( PeriodicTable.Symbol ( atomicNumber ) )] = ( atomicNumber, True )

# . Hydrogen count.
_HCOUNTTOKENS = { "H" : 1 }

# . Remove unwanted items.
del atomicNumber

#===============================================================================
# . Class for a parsable string.
#===============================================================================
class _ParsableString ( object ):

    def __init__ ( self, inputstring ):
        """Constructor."""
        # . Save the input string and set the current position.
        self.inputstring = inputstring
        self.position    = 0
        # . Convert the input string to a list of characters tagged with their positions in the input string.
        self.characters = [ ]
        for ( i, c ) in enumerate ( inputstring ): self.characters.append ( ( c, i ) )
        # . Initialize the bracket data.
        # . Bracketpairs stores the length of each bracket pair.
        self.bracketcache = [ ]
        self.bracketpairs = { }

    def CacheBrackets ( self ):
        """Add a bracket pair to the bracket cache."""
        # . Decrement self.position.
        self.position = self.position - 1
        # . Get the position of the closing bracket corresponding to the opening bracket at the current position.
        closing = self.position + self.bracketpairs[self.characters[self.position][1]]
        # . Save the unparsed bracket string.
        self.bracketcache.extend ( self.characters[self.position:closing+1] )
        # . Delete details about the bracket pair.
        del ( self.bracketpairs[self.characters[self.position][1]] )
        # . Delete the unparsed characters.
	del ( self.characters[self.position:closing+1] )

    def Finalize ( self ):
        """Finish parsing."""
        # . Check that all the string has been parsed.
        if self.position < len ( self.characters ): raise SMILESReaderError ( "There are unparsed characters.", parsestring = self )
        # . Check that the bracket cache is empty.
        if len ( self.bracketcache ) > 0: raise SMILESReaderError ( "There are unparsed brackets.", parsestring = self )

    def GetCharacter ( self, characters, default = None ):
        """Get a single character from the current position."""
        if self.position < len ( self.characters ):
            c = self.characters[self.position][0]
            if c in characters:
                self.position = self.position + 1
                return c
        return default

    def GetCharacters ( self, characters, default = None ):
        """Get a string of characters from the current position."""
        characterlist = [ ]
        while self.position < len ( self.characters ):
            c = self.characters[self.position][0]
            if c in characters:
                characterlist.append ( c )
                self.position = self.position + 1
            else:
                break
        if len ( characterlist ) > 0: return string.join ( characterlist, "" )
        else:                         return default

    def GetInteger ( self, default = None ):
        """Get an integer from the string."""
        integer = self.GetCharacters ( string.digits )
        if integer is None: return default
        else:               return int ( integer )

    def GetToken ( self, tokens, maximumtokenlength, default = None, firstcharacters = None ):
        """Get a token from the string."""
        # . Check that there are more characters to parse.
        if self.position >= len ( self.characters ): return default
        # . Check the first characters (done for speed).
        if firstcharacters is not None:
            if self.characters[self.position][0] not in firstcharacters: return default
        # . Get maximumtokenlength characters from the string (or as near as possible to it).
        characterlist = [ ]
        position      = self.position
        for length in range ( maximumtokenlength ):
            if position < len ( self.characters ): characterlist.append ( self.characters[position][0] )
            else:                                  break
            position = position + 1
        tokenstring = string.join ( characterlist, "" )
        tokenlength = len ( tokenstring )
        # . Loop over token lengths in reverse order.
        for length in range ( tokenlength, 0, -1 ):
            if tokenstring[0:length] in tokens:
                self.position = self.position + length
                return tokens[tokenstring[0:length]]
        # . Nothing found.
        return default

    def GetTokenInteger ( self, tokens, maximumtokenlength, defaultinteger = None, defaulttoken = None, firstcharacters = None ):
        """Get a token from the string with a following integer."""
        # . Save the position.
        position = self.position
        # . Get the token.
        token = self.GetToken ( tokens, maximumtokenlength, default = defaulttoken, firstcharacters = firstcharacters )
        # . Find the following integer.
        if position != self.position: integer = self.GetInteger ( default = defaultinteger )
        else:                         integer = defaultinteger
        # . Return the results.
        return ( token, integer )

    def HasCharacters ( self ):
        """Check whether there are characters remaining in the string."""
        return self.position < len ( self.characters )

    def ParseBrackets ( self ):
        """Parse brackets to check syntax and set up appropriate data structures."""
        # . Initialization.
        copen = 0
        popen = 0
        sopen = 0
        plength = 0
        psize   = [ ]
        pstart  = [ ]
        # . Loop over the characters.
        self.position = 0
        while self.position < len ( self.characters ):
            c = self.characters[self.position][0]
            # . Check for bracket characters.
            # . Parentheses.
            if c == "(":
                if ( copen == 1 ) or ( sopen == 1 ): raise SMILESReaderError ( "Cannot nest ( within [] or {} pairs.", parsestring = self )
                psize.append  ( plength       )
                pstart.append ( self.position )
                plength = 0
                popen   = popen + 1
            elif c == ")":
                if ( popen   <  1 ): raise SMILESReaderError ( "Unpaired closing ).", parsestring = self )
                if ( plength <= 0 ): raise SMILESReaderError ( "Empty () pair.",      parsestring = self )
                opening = pstart.pop ( )
                self.bracketpairs[opening] = self.position - opening
                plength = psize.pop ( )
                popen   = popen - 1
            # . Square brackets.
            elif c == "[":
                if ( copen == 1 ) or ( sopen == 1 ): raise SMILESReaderError ( "Cannot nest [ within [] or {} pairs.", parsestring = self )
                sopen  = 1
                sstart = self.position
            elif c == "]":
                if ( sopen                      != 1 ): raise SMILESReaderError ( "Unpaired closing ].", parsestring = self )
                if ( self.position - sstart - 1 <= 0 ): raise SMILESReaderError ( "Empty [] pair.",      parsestring = self )
                sopen = 0
            # . Curly brackets.
            elif c == "}":
                if ( copen == 1 ) or ( sopen == 1 ): raise SMILESReaderError ( "Cannot nest { within [] or {} pairs.", parsestring = self )
                copen  = 1
                cstart = self.position
            elif c == "{":
                if ( copen                      != 1 ): raise SMILESReaderError ( "Unpaired closing }.", parsestring = self )
                if ( self.position - cstart - 1 <= 0 ): raise SMILESReaderError ( "Empty {} pair.",      parsestring = self )
                copen = 0
            # . Increment plength.
            if ( c != "(" ) and ( c != ")" ): plength = plength + 1
            # . Increment the current position.
            self.position = self.position + 1
        # . Do some final checks.
        if ( copen != 0 ):
            self.position = cstart
            raise SMILESReaderError ( "Unpaired opening {.", parsestring = self )
        if ( popen  > 0 ):
            self.position = pstart[-1]
            raise SMILESReaderError ( "Unpaired opening (.", parsestring = self )
        if ( sopen != 0 ):
            self.position = sstart
            raise SMILESReaderError ( "Unpaired opening [.", parsestring = self )
        # . Reinitialize the current position.
        self.position = 0

    def ParsedString ( self ):
        """Return the current parsed string as a string."""
        characterlist = [ ]
        for character in self.characters: characterlist.append ( character[0] )
        return string.join ( characterlist, "" )

    def RemoveWhiteSpace ( self ):
        """Remove white space."""
        for character in self.characters:
            if character[0] in string.whitespace: self.characters.remove ( character )

    def StatusString ( self, positions = None ):
        """Return the input string and its current position."""
        length = len ( self.inputstring )
        if length > 0:
            # . Get the indices at which to place markers.
            if positions is None: indices = [ self.position ]
            else:                 indices = positions
            cindices = []
            for index in indices:
                if index >= 0: cindices.append ( self.characters[index][1] )
            # . Generate the marker string.
            cindices.sort ( )
            p = cindices[-1]
            marker = ( p + 1 ) * [ " " ]
            for index in cindices: marker[index] = "^"
            # . Create the string.
            s = " " + self.inputstring + "\n" + " " + "".join ( marker )
        else:
	    s = ""
        return s

    def UncacheBrackets ( self ):
        """Uncache any brackets."""
        # . Check for unparsed brackets.
        if len ( self.bracketcache ) > 0:
            # . Add the bracket cache at the current position.
            self.characters[self.position:self.position] = self.bracketcache
            # . Reinitialize the bracket cache.
	    self.bracketcache = [ ]

#===============================================================================
# . Error classes.
#===============================================================================
class SMILESReaderError ( Exception ):

    def __init__ ( self, *arguments, **keywordArguments ):
        """Constructor."""
        super ( SMILESReaderError, self ).__init__ ( *arguments )
        if len ( arguments ) == 0: self.args = ( "SMILES reader error.", )
        parseString = keywordArguments.get ( "parsestring", None )
        positions   = keywordArguments.get ( "positions"  , None )
        if parseString is not None:
            self.args = ( self.args[0] + "\n" + parseString.StatusString ( positions = positions ), )

#===============================================================================
# . SMILES reader class.
#===============================================================================
class SMILESReader ( object ):
    """SMILESReader is the class for reading SMILES strings."""

    def __init__ ( self, inputsmiles ):
        """Constructor."""
        self.QPARSED       = False
        self.atompositions = {}
        self.connectivity  = None
        self.parsedsmiles  = None
        self.smilesstring  = inputsmiles

    def Parse ( self ):
        """Parse a SMILES string."""
        if not self.QPARSED:

            # . Start a connectivity.
            self.connectivity = SMILESConnectivity ( )

            try:
                if len ( self.smilesstring ) > 0:

                    # . Initialization.
                    lastatom        = None
                    bracketlastatom = [ ]
                    crosslink       = { }

                    # . Initialize the input string for parsing.
                    parsestring = _ParsableString ( self.smilesstring )
                    parsestring.RemoveWhiteSpace ( )
                    parsestring.ParseBrackets    ( )

                    # . Loop over all the characters in the string.
                    while parsestring.HasCharacters ( ):

                        # . Check for parentheses.
                        c = parsestring.GetCharacter ( [ "(", ")" ], default = " " )

                        # . An opening parenthesis.
	                if c == "(":
                            # . There is a root atom so parse this bracket.
		            if lastatom is not None:
                                # . Save and reset lastatom.
                                bracketlastatom.append ( lastatom )
		                lastatom = None
                            # . There is no root atom.
		            else:
                                # . Cache this bracket pair for later.
                                parsestring.CacheBrackets ( )

                        # . A closing parenthesis.
	                elif c == ")":
                            # . Reset lastatom.
		            lastatom = bracketlastatom.pop ( )

                        # . An atom or crosslink specification.
	                else:

                            # . For the first atom in the system set a NullBond character.
                            if len ( self.connectivity ) == 0:
                                bondtype = NullBond ( )
                            # . Check for an explicit bond character.
                            else:
                                bondtype = parsestring.GetToken ( _BONDTOKENS, 1, default = SingleBond ( ) )

                            # . Check for a crosslink or a full atom specification.
                            c = parsestring.GetCharacter ( "%[" + string.digits, default = " " )

                            # . Parse the atom or crosslink string.
                            if ( c == "%" ) or ( c in string.digits ):

                                # . Check that lastatom exists.
                                if lastatom is None: raise SMILESReaderError ( "A crosslink specification must have a preceding atom.", parsestring = parsestring )

                                # . Parse the crosslink string.
                                if c == "%":
                                    link = parsestring.GetInteger ( )
                                    if link is None: raise SMILESReaderError ( "The crosslink integer is missing.", parsestring = parsestring )
                                else:
                                    link = int ( c )

                                # . Check to see if this crosslink already exists.
                                if link in crosslink:

                                    # . Get the link data.
                                    ( connectedatom, firsttype ) = crosslink[link]

                                    # . If this key is already associated with lastatom there is an error.
                                    if connectedatom == lastatom:
                                        raise SMILESReaderError ( "An atom has two identical crosslink keys.", parsestring = parsestring )

                                    # . The key is not associated with lastatom.
                                    else:

                                        # . Get the type of bond to add to the bond list.
                                        if isinstance ( bondtype, SingleBond ):
                                            bondtype = firsttype
                                        # . The bond specifications are incompatible.
                                        elif ( bondtype is not firsttype ) and ( not isinstance ( firsttype, SingleBond ) ):
                                            raise SMILESReaderError ( "A pair of crosslink keys has incompatible bond-type specifications.", parsestring = parsestring )

                                        # . Add this bond to the bond list (if it is not NullBond) and delete the associated dictionary entry.
                                        if not isinstance ( bondtype, NullBond ): self.connectivity.AddBond ( SMILESBond ( atom1 = lastatom, atom2 = connectedatom, type = bondtype ) )
                                        del crosslink[link]

                                # . The key does not exist so add it to the dictionary.
                                else:
                                    crosslink[link] = ( lastatom, bondtype )

                            # . An atom specification.
                            else:

                                # . Full atom specification.
                                if c == "[":

                                    # . Set QREDUCED.
                                    QREDUCED = False

                                    # . Parse the individual atom fields.
                                    # . Isotope.
                                    isotope = parsestring.GetInteger ( default = 0 )

                                    # . Atom symbol.
                                    smilesposition          = parsestring.position
                                    atomicNumber, isAromatic = parsestring.GetToken ( _ELEMENTTOKENS, 2, default = ( -1, False ) )
                                    if atomicNumber == -1: raise SMILESReaderError ( "Missing or unrecognized element symbol.", parsestring = parsestring )

                                    # . Chirality (full and reduced notations).
                                    ( chiralityclass, chiralitynumber ) = parsestring.GetToken ( _CHIRALITYREDUCED, 4, default = ( None, 0 ), firstcharacters = [ "@" ] )
                                    if chiralityclass is not None:
                                        if chiralitynumber == 1:
                                            ( chiralityclass, chiralitynumber ) = parsestring.GetTokenInteger ( _CHIRALITYFULL, 2, defaultinteger = 0, defaulttoken = None )
                                            if chiralityclass is None:
                                                chiralityclass  = "??"
                                                chiralitynumber = 1

                                    # . Hydrogen count.
                                    hcount, multiplier = parsestring.GetTokenInteger ( _HCOUNTTOKENS, 1, defaultinteger = 1, defaulttoken = 0 )
                                    hcount *= multiplier

                                    # . Charge.
                                    charge = parsestring.GetToken ( _CHARGETOKENS, 4, default = 0, firstcharacters = [ "+", "-" ] )
                                    if ( charge == 1 ) or ( charge == -1 ):
                                        multiplier = parsestring.GetInteger ( default = 1 )
                                        charge = charge * multiplier

                                    # . Check the current character which must be a "]".
                                    if parsestring.GetCharacter ( [ "]" ], default = " " ) != "]": raise SMILESReaderError ( "A closing ] is missing.", parsestring = parsestring )

                                # . Reduced atom specification.
                                else:

                                    # . Set QREDUCED and other defaults.
                                    charge   = 0
                                    hcount   = 0
                                    isotope  = 0
                                    QREDUCED = True

                                    # . Parse the individual atom fields.
                                    # . Atom symbol.
                                    smilesposition          = parsestring.position
                                    atomicNumber, isAromatic = parsestring.GetToken ( _REDUCEDELEMENTTOKENS, 2, default = ( -1, False ) )
                                    if atomicNumber == -1: raise SMILESReaderError ( "Missing or unrecognized organic-subset element symbol.", parsestring = parsestring )

                                    # . Chirality.
                                    ( chiralityclass, chiralitynumber ) = parsestring.GetToken ( _CHIRALITYREDUCED, 4, default = ( None, 0 ), firstcharacters = [ "@" ] )

                                # . Create a new atom and set its properties.
                                newatom = SMILESAtom ( atomicNumber = atomicNumber )
                                newatom.chiralityclass      = chiralityclass
                                newatom.chiralitynumber     = chiralitynumber
                                newatom.formalCharge        = charge
                                newatom.implicithydrogens   = hcount
                                newatom.isotope             = isotope
                                newatom.isAromatic           = isAromatic
                                newatom.QREDUCED            = QREDUCED
                                self.atompositions[newatom] = smilesposition

                                # . Put all the atom data in the atom list.
                                self.connectivity.AddAtom ( newatom )

                                # . Get lastatom for the first atom of a branch.
		                if lastatom is None:
		                    if len ( bracketlastatom ) > 0: lastatom = bracketlastatom[-1]

                                # . Save the current bond.
                                if ( not isinstance ( bondtype, NullBond ) ) and ( lastatom is not None ): self.connectivity.AddBond ( SMILESBond ( atom1 = newatom, atom2 = lastatom, type = bondtype ) )

                                # . Reset lastatom.
		                lastatom = newatom

                                # . Uncache any unparsed bracket pairs.
                                parsestring.UncacheBrackets ( )

	            # . Do some last checks.
                    # . Check parsestring.
                    parsestring.Finalize ( )

                    # . Check to see that the crosslink dictionary is empty.
                    if len ( crosslink ) > 0: raise SMILESReaderError ( "There are unparsed crosslink keys: " + " ".join ( [ "{:d}".format ( i ) for i in crosslink.keys ( ) ] ) )

                    # . Save the final parsed SMILES.
                    self.parsedsmiles = parsestring.ParsedString ( )

                # . Zero length string.
                else:
                    self.parsedsmiles = ""

                # . Finish up.
                self.connectivity.Finalize ( )

            # . Errors.
            except SMILESConnectivityError as error:
                if error.atoms is not None:
                    positions = []
                    for atom in error.atoms: positions.append ( self.atompositions[atom] )
                elif error.atom  is not None: positions = [ self.atompositions[error.atom] ]
                elif error.bond  is not None: positions = [ self.atompositions[error.bond.atom1], self.atompositions[error.bond.atom2] ]
                else:                         positions = None
                raise SMILESReaderError ( error.args[0], parsestring = parsestring, positions = positions )

            # . Flag the SMILES as having been parsed.
            self.QPARSED = True

    def ToSystem ( self, **keywordArguments ):
        """Return a system."""
        if self.QPARSED:
            system       = self.connectivity.ToSystem ( **keywordArguments )
            system.label = self.parsedsmiles
        else:
            system       = None
        return system

#===============================================================================
# . Helper functions.
#===============================================================================
def SMILES_ToSystem ( smiles, **keywordArguments ):
    """Helper function that reads a system from a SMILES string."""
    reader = SMILESReader ( smiles )
    reader.Parse ( )
    log = keywordArguments.get ( "log", logFile )
    reader.connectivity.Summary ( log = log )
    system = reader.ToSystem ( **keywordArguments )
    return system

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
