#-------------------------------------------------------------------------------
# . File      : PointGroup.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Point group classes and functions."""

import glob, math, os, os.path

from pCore import Real2DArray

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PointGroup ( object ):
    """Point group class."""

    defaultAttributes = { "label"                       : ""  , \
                          "characterSymmetryOperations" : None, \
                          "characterTable"              : None, \
                          "cInfinityRotations"          : None, \
                          "irreducibleRepresentations"  : None, \
                          "maximumDegeneracy"           :    0, \
                          "principalAxisOperation"      : None, \
                          "symmetryOperations"          : None  }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in keywordArguments.iteritems ( ):
            if key in self.__class__.defaultAttributes.keys ( ): setattr ( self, key, value )
        # . Find the maximum degeneracy.
        if self.characterTable is not None:
            e = self.characterSymmetryOperations.index ( "E" )
            m = 0.0
            for c in range ( len ( self.irreducibleRepresentations ) ):
                m = max ( m, self.characterTable[c,e] )
            self.maximumDegeneracy = int ( round ( m ) )
        # . Find the principal axis.
        if self.symmetryOperations is not None:
            if "CInfinity" in self.symmetryOperations: self.principalAxisOperation = "CInfinity"
            else:
                maximumN = 0
                for ( key, value ) in self.symmetryOperations.iteritems ( ):
                   if key.startswith ( "C" ): maximumN = max ( maximumN, int ( key[1:] ) )
                if maximumN > 0:
                    key = "C{:d}".format ( maximumN )
                    if self.symmetryOperations[key] == 1: self.principalAxisOperation = key

    def OperationKey ( self ):
        """Return the operation key."""
        if self.__dict__.get ( "_operationKey", None ) is None:
            keys = self.symmetryOperations.keys ( )
            keys.sort ( )
            items = []
            for key in keys:
                n = self.symmetryOperations[key]
                if n == 1: items.append ( key )
                else:      items.append ( repr ( n ) + "*" + key )
            self._operationKey = " ".join ( items )
        return self._operationKey

    @classmethod
    def FromText ( selfClass, path ):
        """Constructor from text file."""
        from math import cos, pi, sin
        # . Initialization.
        label                       = None
        characterSymmetryOperations = None
        cInfinityRotations          = []
        rows                        = {}
        symmetryOperations          = {}
        # . Read the file.
        inFile = open ( path, "r" )
        while True:
            try:
                tokens = next ( inFile ).split ( )
                if ( len ( tokens ) <= 0 ) or tokens[0].startswith ( "#" ): continue
                if ( tokens[0] == "Point" ) and ( len ( tokens ) >= 3 ):
                    label = tokens[2]
                elif ( tokens[0] == "Symmetry" ):
                    while True:
                        tokens = next ( inFile ).split ( )
                        if len ( tokens ) <= 0: break
                        symmetryOperations[tokens[0]] = int ( tokens[1] )
                elif ( tokens[0] == "CInfinity" ):
                    while True:
                        tokens = next ( inFile ).split ( )
                        if len ( tokens ) <= 0: break
                        cInfinityRotations.append ( tokens[0] )
                elif ( tokens[0] == "Character" ):
                    characterSymmetryOperations = next ( inFile ).split ( )
                    while True:
                        tokens = next ( inFile ).split ( )
                        if len ( tokens ) <= 0: break
                        rows[tokens[0]] = tokens[1:]
            except StopIteration:
                break
        inFile.close ( )
        # . Process character table data.
        irreducibleRepresentations = rows.keys ( )
        irreducibleRepresentations.sort ( )
        characterTable = Real2DArray.WithExtents ( len ( irreducibleRepresentations ), len ( characterSymmetryOperations ) )
        characterTable.Set ( 0.0 )
        for ( i, iR ) in enumerate ( irreducibleRepresentations ):
            row = rows[iR]
            for j in range ( len ( characterSymmetryOperations ) ):
                characterTable[i,j] = eval ( row[j] )
        # . Create the point group.
        group = selfClass ( label = label, \
                            characterSymmetryOperations = characterSymmetryOperations, \
                            characterTable              = characterTable             , \
                            cInfinityRotations          = cInfinityRotations         , \
                            irreducibleRepresentations  = irreducibleRepresentations , \
                            symmetryOperations          = symmetryOperations           )
        # . Finish up.
        return group

    def ToText ( self, path ):
        """Write to a text file."""
        # . Header.
        outFile = open ( path, "w" )
        outFile.write ( "\n" + 80 * "=" + "\n" )
        outFile.write ( "\nPoint Group:" + self.label + "\n" )
        outFile.write ( "\nSymmetry Operations:" )
        labels = self.symmetryOperations.keys ( )
        labels.sort ( )
        for label in labels:
            outFile.write ( "\n{:<20s} {:d}".format ( label, self.symmetryOperations[label] ) )
        outFile.write ( "\n" )
        if self.cInfinityRotations is not None:
            outFile.write ( "\nCInfinity Rotations:" )
            for label in self.cInfinityRotations:
                outFile.write ( "\n{:<s}".format ( label ) )
            outFile.write ( "\n" )
        outFile.write ( "\nCharacter Table:" )
        outFile.write ( "\n        " )
        for label in self.characterSymmetryOperations:
            outFile.write ( "{:s}".format ( label.center ( 16 ) ) )
        numberOperations = len ( self.characterSymmetryOperations )
        for ( i, iR ) in enumerate ( self.irreducibleRepresentations ):
            outFile.write ( "\n{:<8s}".format ( iR ) )
            for j in range ( numberOperations ): outFile.write ( "{:16.6f}".format ( self.characterTable[i,j] ) )
        outFile.write ( "\n" )
        outFile.write ( "\n" + 80 * "=" + "\n" )
        outFile.close ( )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def PointGroups_FromText ( inPath = None ):
    """Read point groups from a directory containing text files.

    Return as a dictionary with symmetry operation keys.
    """
    # . Get the path.
    if inPath is None:
        try:    inPath = os.path.join ( os.getenv ( "PDYNAMO_PARAMETERS" ), "pointGroups" )
        except: raise ( "Unable to find point group directory." )
    # . Initialization.
    pointGroups = {}
    # . Loop over all files.
    paths = glob.glob ( os.path.join ( inPath, "*.txt" ) )
    for path in paths:
        if path.find ( "00ReadMe" ) < 0:
            group = PointGroup.FromText ( path )
            key   = group.OperationKey ( )
            pointGroups[key] = group
    # . Finish up.
    return pointGroups

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
