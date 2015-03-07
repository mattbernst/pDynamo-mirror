#-------------------------------------------------------------------------------
# . File      : MMAtomTyper.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MM atom typing."""

import os

from pCore               import logFile, LogFileActive, Selection, YAMLMappingFile_ToObject, YAMLPickleFileExtension
from MMAtomTypeContainer import MMAtomTypeContainer
from MMModelError        import MMModelError
from MMPatternContainer  import MMPatternContainer
from MMSequenceLibrary   import MMSequenceLibrary

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Maximum width for table output.
_MaximumTableWidth = 100

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMAtomTyper ( object ):
    """Type atoms."""

    defaultAttributes = { "dataPath"          : None,
                          "mmAtomTypes"       : None,
                          "mmPatterns"        : None,
                          "mmSequenceLibrary" : None }

    def __init__ ( self, dataPath, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )
        self.dataPath = dataPath
        self.ReadData ( )

    def CheckUntypedAtoms ( self, connectivity, untypedAtoms, log ):
        """Check for untyped atoms."""
        # . Check for untyped atoms.
        if len ( untypedAtoms ) > 0:
            if LogFileActive ( log ):
                labels = []
                length = 0
                for i in untypedAtoms:
                    label  = connectivity.atoms[i].path
                    labels.append ( label )
                    length = max ( length, len ( label ) )
                ncolumns = _MaximumTableWidth // ( length + 2 )
                table = log.GetTable ( columns = ncolumns * [ length + 2 ] )
                table.Start ( )
                table.Title ( "Untyped Atoms" )
                for label in labels: table.Entry ( label )
                table.Stop ( )
            raise MMModelError ( "There are {:d} untyped atoms.".format ( len ( untypedAtoms ) ), untypedAtoms )

    def ReadData ( self ):
        """Read the atom typing data."""
        if self.dataPath is not None:
            # . Atom types and patterns.
            for ( attribute, tag, pathClass ) in ( ( "mmAtomTypes", "atomTypes", MMAtomTypeContainer ), ( "mmPatterns", "patterns", MMPatternContainer ) ):
                path = os.path.join ( self.dataPath, tag + YAMLPickleFileExtension )
                if os.path.exists ( path ): setattr ( self, attribute, YAMLMappingFile_ToObject ( path, pathClass ) )
            if self.mmPatterns is not None: self.mmPatterns.IndexAtomTypes ( self.mmAtomTypes )
            # . Sequence attributes.
            self.mmSequenceLibrary = MMSequenceLibrary.TestForLibrary ( self.dataPath )

    def TypeAtoms ( self, connectivity, sequence, log ):
        """Assign atom types and charges to the atoms."""
        # . Initialization.
        atomCharges  = [ 0.0  for i in range ( len ( connectivity.atoms ) ) ]
        atomTypes    = [ None for i in range ( len ( connectivity.atoms ) ) ]
        untypedAtoms = Selection.FromIterable ( range ( len ( atomTypes ) ) )
        # . Type the atoms.
        self.TypeBySequence ( sequence    , atomTypes, atomCharges, untypedAtoms )
        self.TypeByPattern  ( connectivity, atomTypes, atomCharges, untypedAtoms )
        # . Check for untyped atoms.
        self.CheckUntypedAtoms ( connectivity, untypedAtoms, log )
        # . Finish up.
        return ( atomTypes, atomCharges )

    def TypeByPattern ( self, connectivity, atomTypes, atomCharges, untypedAtoms ):
        """Type the atoms by pattern."""
        if ( self.mmAtomTypes is not None ) and ( self.mmPatterns is not None ):
            self.mmPatterns.TypeConnectivity ( connectivity, atomTypes, atomCharges, untypedAtoms )

    def TypeBySequence ( self, sequence, atomTypes, atomCharges, untypedAtoms ):
        """Type the atoms by sequence."""
        if self.mmSequenceLibrary is not None:
            self.mmSequenceLibrary.TypeSequence ( sequence, atomTypes, atomCharges, untypedAtoms )
 
#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
