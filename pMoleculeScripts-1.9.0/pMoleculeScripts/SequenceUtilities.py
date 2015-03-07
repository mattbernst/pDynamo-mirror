#-------------------------------------------------------------------------------
# . File      : SequenceUtilities.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Utilities for dealing with sequences assuming a PDB-like format."""

import string

from pCore     import Clone, logFile, LogFileActive
from pMolecule import PeriodicTable, Sequence, SequenceComponent, SequenceEntity

# . Move eventually to sequence?

#===================================================================================================================================
# . Create a sequence consisting of one atom per component named after elemental type.
#===================================================================================================================================
def CreateElementSequence ( system, entityDescription = None, entityLabel = "E" ):
    """Create a sequence of a single entity using element symbols for atom names and components and one atom per component."""
    # . Initialization.
    minorSeparator = Sequence.defaultAttributes["fieldSeparator"]
    sequence       = Sequence ( )
    entity         = SequenceEntity ( label = entityLabel )
    sequence.AddChild ( entity )
    # . Process the atoms.
    for ( index, atom ) in enumerate ( system.atoms ):
        symbol    = PeriodicTable.Symbol ( atom.atomicNumber ).upper ( )
        component = SequenceComponent ( genericLabel = symbol, label = symbol + minorSeparator + repr ( index + 1 ) )
        entity.AddChild    ( component )
        component.AddChild ( atom      )
        atom.index = index
        atom.label = symbol
    # . Finish up.
    system.__dict__["_sequence"] = sequence

#===================================================================================================================================
# . Create a homogeneous isolate sequence consisting of a single entity.
#===================================================================================================================================
def CreateHomogeneousIsolateSequence ( system, atomData = [ ( 8, "O" ), ( 1, "H1" ), ( 1, "H2" ) ], componentGenericLabel = "HOH",
                                                       entityDescription = "Water", entityLabel = "W", firstComponentNumber = 1 ):
    """Create a sequence assuming a system composed of homogeneous isolates.

    Argument values default to water.
    """
    # . Initialization.
    minorSeparator = Sequence.defaultAttributes["fieldSeparator"]
    sequence       = Sequence ( )
    entity         = SequenceEntity ( label = entityLabel )
    sequence.AddChild ( entity )
    # . Loop over the atoms per isolate.
    atomIndex = 0
    for ( index, isolate ) in enumerate ( system.connectivity.isolates ):
        if len ( isolate ) != len ( atomData ): raise ValueError ( "Incompatible atom data and isolate lengths." )
        component = SequenceComponent ( genericLabel = componentGenericLabel, label = componentGenericLabel + minorSeparator + repr ( index + 1 ) )
        entity.AddChild ( component )
        for ( ( atomicNumber, atomLabel ), i ) in zip ( atomData, isolate ):
            atom = system.atoms[i]
            if ( atomicNumber != atom.atomicNumber ): raise ValueError ( "Atomic number mismatch." )
            component.AddChild ( atom )
            atom.index = atomIndex
            atom.label = atomLabel
            atomIndex += 1
    # . Finish up.
    system.__dict__["_sequence"] = sequence

#===================================================================================================================================
# . Determine a unique entity label given a set of systems.
#===================================================================================================================================
def DetermineUniqueEntityLabel ( *arguments, **keywordArguments ):
    """Determine a unique entity label."""
    # . Get existing labels.
    labels = set ( )
    for arg in arguments:
        sequence = getattr ( arg, "sequence", None )
        if sequence is not None:
            for entity in sequence.children:
                fields = sequence.ParseLabel ( entity.label, fields = 1 )
                labels.add ( fields[0].upper ( ) )
    # . Get free labels.
    freelabels = list ( set ( string.ascii_uppercase ).difference ( labels ) )
    if len ( freelabels ) <= 0: raise ValueError ( "There are no unique entity labels remaining." )
    # . Label specified and not taken so use it.
    label = keywordArguments.get ( "label", None )
    if ( label is not None ) and ( label in freelabels ):
        unique = label
    # . Use the first free label.
    else:
        freelabels.sort ( )
        unique = freelabels.pop ( 0 )
    return unique

#===================================================================================================================================
# . Print the frequencies of individual components in the entities of a sequence.
#===================================================================================================================================
def PrintComponentFrequency ( sequence, log = logFile, title = "Sequence Component Frequency" ):
    """Print the component frequency of the entities of a sequence."""
    if LogFileActive ( log ) and isinstance ( sequence, Sequence ):
        log.Heading ( title )
        for entity in sequence.children:
            if len ( entity.children ) > 0:
                # . Frequencies.
                frequencies = {}
                for component in entity.children:
                    frequencies[component.genericLabel] = frequencies.get ( component.genericLabel, 0 ) + 1
                keys = frequencies.keys ( )
                keys.sort ( )
                # . Output.
                length = min ( 10, len ( frequencies ) )
                table  = log.GetTable ( columns = length * [ 6, 6 ] )
                table.Start ( )
                if len ( entity.label ) <= 0: title = "Component Frequency for Unnamed Entity"
                else:                         title = "Component Frequency for Entity " + entity.label
                table.Title ( title )
                for key in keys:
                    table.Entry ( key                               , alignment = "r" )
                    table.Entry ( "{:d}".format ( frequencies[key] ), alignment = "r" )
                table.Stop ( )

#===================================================================================================================================
# . Renumber the components in the entities of a sequence.
#===================================================================================================================================
def RenumberEntityComponents ( system, entityLabels = None ):
    """Renumber the components in the entities of a sequence."""
    if system.sequence is not None:
        for entity in system.sequence.children:
            if ( entityLabels is None ) or ( entity.label in entityLabels ):
                for ( i, component ) in enumerate ( entity.children ):
                    fields = system.sequence.ParseLabel ( component.label )
                    if len ( fields ) < 1: fields.append ( "" )
                    fields[1] =  repr ( i + 1 )
                    component.label = system.sequence.MakeLabel ( *fields )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
