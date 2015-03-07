#-------------------------------------------------------------------------------
# . File      : Element.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Element classes."""

import glob, math, os, os.path

from pCore import logFile, LogFileActive, YAMLMappingFile_FromObject, YAMLMappingFile_ToObject, YAMLPickleFileExtension

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Unknown element data.
_UnknownElementNumber = -1
_UnknownElementSymbol = "*"

# . YAML tag.
_YAMLTag = "!Element"

#===================================================================================================================================
# . Element class.
#===================================================================================================================================
class Element ( object ):
    """Element defines a chemical element."""

    attributeLabelMappings = { "Atomic Number"             : "atomicNumber"             ,
                               "Bragg Radius"              : "braggRadius"              ,
                               "Coordination Angles"       : "coordinationAngles"       ,
                               "Covalent Radius"           : "covalentRadius"           ,
                               "d Quantum Number"          : "dqn"                      ,
                               "f Quantum Number"          : "fqn"                      ,
                               "Mass"                      : "mass"                     ,
                               "Name"                      : "name"                     ,
                               "Pauling Electronegativity" : "paulingElectronegativity" ,
                               "p Quantum Number"          : "pqn"                      ,
                               "RGB Color Indices"         : "rgbColorIndices"          ,
                               "Single Bond Distances"     : "singleBondDistances"      ,
                               "s Quantum Number"          : "sqn"                      ,
                               "Symbol"                    : "symbol"                   ,
                               "VdW Radius"                : "vdwRadius"                }

    defaultAttributes = { "atomicNumber"             : _UnknownElementNumber ,
                          "braggRadius"              : 0.0                   ,
                          "coordinationAngles"       : None                  ,
                          "covalentRadius"           : 0.0                   ,
                          "dqn"                      : 0                     ,
                          "fqn"                      : 0                     ,
                          "mass"                     : 0.0                   ,
                          "name"                     : ""                    ,
                          "paulingElectronegativity" : None                  ,
                          "pqn"                      : 0                     ,
                          "rgbColorIndices"          : None                  ,
                          "singleBondDistances"      : None                  ,
                          "sqn"                      : 0                     ,
                          "symbol"                   : _UnknownElementSymbol ,
                          "vdwRadius"                : 0.0                   }

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        for ( fullKey, key ) in self.__class__.attributeLabelMappings.iteritems ( ):
            attribute = getattr ( self, key, None )
            if attribute is not None: mapping[fullKey] = attribute
        return mapping

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            self.__init__ ( )
            for ( fullKey, key ) in self.__class__.attributeLabelMappings.iteritems ( ):
                attribute = mapping.pop ( fullKey, None )
                if attribute is not None: setattr ( self, key, attribute )
        except:
            raise ValueError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def GetCoordinationAngle ( self, connections ):
        """Get a coordination angle given the coordination number."""
        angle = None
        if self.coordinationAngles is not None:
            angle = self.coordinationAngles.get ( connections, None )
        return angle

    def GetSingleBondDistance ( self, atomicNumber ):
        """Get a single bond distance to another element."""
        r = None
        if self.singleBondDistances is not None:
            r = self.singleBondDistances.get ( atomicNumber, None )
        return r

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

#===================================================================================================================================
# . Element container class.
#===================================================================================================================================
class ElementContainer ( object ):
    """ElementContainer defines a container for elements."""

    def __init__ ( self, items ):
        """Constructor."""
        # . Basic indexes.
        massIndex        = []
        self.numberIndex = {}
        self.symbolIndex = {}
        for item in items:
            atomicNumber = item.atomicNumber
            symbol       = item.symbol
            if ( atomicNumber not in self.numberIndex ) and ( symbol not in self.symbolIndex ):
                massIndex.append ( ( item.mass, item ) )
                self.numberIndex[atomicNumber] = item
                self.symbolIndex[symbol]       = item
        # . Items stored in order of increasing mass.
        # . This is necessary for finding the atomic number from the mass.
        massIndex.sort ( )
        self.items = []
        for ( m, item ) in massIndex:
            self.items.append ( item )

    def AtomicNumber ( self, inString ):
        """Find the atomic number of an element from a string."""
        number = _UnknownElementNumber
        if len ( inString ) > 0:
            # . The string contains a symbol.
            if inString[0:1].isalpha ( ):
                a = inString[0:1].upper ( )
                if ( len ( inString ) > 1 ) and inString[1:2].isalpha ( ): symbol = a + inString[1:2].lower ( )
                else:                                                      symbol = a
                if symbol in self.symbolIndex: number = self.symbolIndex[symbol].atomicNumber
            # . The string is an integer.
            else:
                ndigits = 0
                for c in inString:
                    if c.isdigit ( ): ndigits += 1
                    else:             break
                if ndigits > 0: number = int ( inString[0:ndigits] )
        return number

    # . GMA: a tolerance is introduced that selects the first element whose mass differs less than this.
    # . One cannot assume that the masses increase monotonically.
    def AtomicNumberFromMass ( self, mass, tolerance = 0.0 ):
        """Find the atomic number given a mass."""
        if tolerance < 0.0: tolerance = 0.0
        useTolerance      = ( tolerance > 0.0 )
        number            = _UnknownElementNumber
        minimumDifference = math.fabs ( mass - self.items[0].mass )
        for i in range ( 1, len ( self.items ) ):
            difference = math.fabs ( mass - self.items[i].mass )
            if useTolerance:
                if difference < tolerance:
                    number = i
                    break
            else:
                if difference < minimumDifference:
                    minimumDifference = difference
                    number = i
        return number

    def AtomicNumberFromSymbol ( self, inString ):
        """Find the atomic number of an element given a symbol."""
        number = _UnknownElementNumber
        if len ( inString ) > 0:
            if inString[0:1].isalpha ( ):
                a = inString[0:1].upper ( )
                if ( len ( inString ) > 1 ) and inString[1:2].isalpha ( ): symbol = a + inString[1:2].lower ( )
                else:                                                      symbol = a
                if  symbol in self.symbolIndex: number = self.symbolIndex[symbol].atomicNumber
        return number

    def AtomicNumbers ( self ):
        """Return the atomic numbers of the elements in the container."""
        numbers = self.numberIndex.keys ( )
        numbers.sort ( )
        return numbers

    def Element ( self, atomicNumber ):
        """Return an element given an atomic number."""
        return self.numberIndex.get ( atomicNumber, None )

    @classmethod
    def FromParameterDirectory ( selfClass, path = None ):
        """Constructor from a directory of parameters."""
        if path is None: path = os.path.join ( os.getenv ( "PDYNAMO_PARAMETERS" ), "elements" )
        filePaths = glob.glob ( os.path.join ( path, "*" + YAMLPickleFileExtension ) )
        items     = []
        for filePath in filePaths:
            items.append ( YAMLMappingFile_ToObject ( filePath, Element ) )
        self = selfClass ( items )
        return self

    def Summary ( self, log = logFile ):
        """Summarizing."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "Element Container Summary" )
            summary.Entry ( "Number of Elements" , "{:d}".format ( len ( self.items ) ) )
            summary.Stop ( )

    def Symbol ( self, atomicNumber, index = None ):
        """Return a symbol given an atomic number with an optional index."""
        element = self.Element ( atomicNumber )
        if element is None: symbol = _UnknownElementSymbol
        else:               symbol = element.symbol
        if index is None: return symbol
        else:             return symbol + repr ( index )

    def Symbols ( self ):
        """Return the symbols of the elements in the container."""
        symbols = self.symbolIndex.keys ( )
        symbols.sort ( )
        return symbols

    def ToParameterDirectory ( self, path ):
        """Create a parameter directory from the container."""
        if not os.path.exists ( path ): os.mkdir ( path )
        for ( key, value ) in self.symbolIndex.iteritems ( ):
            YAMLMappingFile_FromObject ( os.path.join ( path, key + YAMLPickleFileExtension ), _YAMLTag, value, default_flow_style = False )

#===================================================================================================================================
# . Define the periodic table.
#===================================================================================================================================
PeriodicTable = ElementContainer.FromParameterDirectory ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    path = os.path.join ( os.getenv ( "PDYNAMO_SCRATCH" ), "elements" )
    PeriodicTable.ToParameterDirectory ( path )
    elements = ElementContainer.FromParameterDirectory ( path = path )
    elements.Summary ( )
