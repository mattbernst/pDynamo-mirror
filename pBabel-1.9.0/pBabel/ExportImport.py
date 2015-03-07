#-------------------------------------------------------------------------------
# . File      : ExportImport.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for importing and exporting objects to and from external files."""

import os.path

from pCore     import Coordinates3, logFile, LogFileActive, Pickle, Unpickle, YAMLPickle, YAMLUnpickle
from pMolecule import System                                          

# . How handle object subclasses, etc.?

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
class ExportImportError ( Exception ): pass

class ExportImportHandler ( object ):
    """Class for a handler that imports or exports objects."""

    def __init__ ( self, functionIndex, identifiers, label, defaultFunction = None ):
        """Constructor."""
        self.functionIndex = functionIndex
        self.identifiers   = identifiers
        self.label         = label
        if defaultFunction is None: self.defaultFunction = functionIndex.get ( None, None )
        else:                       self.defaultFunction = defaultFunction

    def ExportObject ( self, path, object, options = None ):
        """Export an object."""
        function = self.GetFunction ( object.__class__ )
        if options is None: function ( path, object            )
        else:               function ( path, object, **options )

    def GetFunction ( self, objectClass ):
        """Get the appropriate function to handle a given object class."""
        function = self.functionIndex.get ( objectClass, self.defaultFunction )
        if function is None:
            if objectClass is None: tag = "default class"
            else:                   tag = "class \"" + objectClass.__name__ + "\""
            raise ExportImportError ( "Handler \"{:s}\" cannot treat objects of {:s}.".format ( self.label, tag ) )
        return function

    def ImportObject ( self, path, objectClass, options = None ):
        """Import an object."""
        function = self.GetFunction ( objectClass )
        if options is None: item = function ( path            )
        else:               item = function ( path, **options )
        if ( objectClass is not None ) and ( not isinstance ( item, objectClass ) ):
            raise ExportImportError ( "Invalid class for imported object: " + item.__class__.__name__ + "." )
        return item

    @property
    def identifierString ( self ):
        """Return an identifier string."""
        items = set ( self.identifiers )
        items.remove ( self.label )
        items = list ( items )
        items.sort ( )
        return ", ".join ( items )

    @property
    def objectString ( self ):
        """Return an object string."""
        items = []
        for item in self.functionIndex.keys ( ):
            if item is None: items.append ( "All"         )
            else:            items.append ( item.__name__ )
        items.sort ( )
        return ", ".join ( items )

class ExportImportHandlerContainer ( object ):
    """A container for import/export handlers."""

    def __init__ ( self ):
        """Constructor."""
        self.handlerIndex = {}
        self.handlers     = []
        self.identifiers  = set ( )

    def AddHandler ( self, functionIndex, identifiers, label, defaultFunction = None ):
        """Add a handler."""
        # . Check for duplicate identifiers.
        identifiers = set ( identifiers )
        identifiers.add ( label )
        common      = self.identifiers & identifiers
        if len ( common ) > 0:
            common = list ( common )
            common.sort ( )
            raise ExportImportError ( "Duplicate handler identifiers: {!r}.".format ( common ) )
        # . Create the object.
        handler = ExportImportHandler ( functionIndex, identifiers, label, defaultFunction = defaultFunction )
        for identifier in identifiers:
            self.handlerIndex[identifier] = handler
        self.handlers.append    ( handler     )
        self.identifiers.update ( identifiers )

    def FileFormats ( self ):
        """Get a list of file formats and their identifiers."""
        data = []
        for handler in self.handlers:
            label       = handler.label
            identifiers = list ( handler.identifiers )
            identifiers.remove ( label )
            identifiers.sort ( )
            data.append ( ( label, identifiers ) )
        data.sort ( )
        return data

    def GetHandler ( self, path, format = None ):
        """Get a handler given a path and format.

        The format is determined from the path extension if no format is given.
        """
        if format is None:
            format = os.path.splitext ( path )[-1]
            if format.startswith ( "." ): format = format[1:]
        try   : return self.handlerIndex[format]
        except: raise ExportImportError ( "Unrecognized format: " + format + "." )

    def HandlerSummary ( self, log = logFile, title = "Export/Import Handler Container" ):
        """Handler summary."""
        if LogFileActive ( log ):
            # . Gather data.
            data = []
            l0 = l1 = l2 = 18
            for handler in self.handlers:
                label       = handler.label
                identifiers = handler.identifierString
                objects     = handler.objectString
                l0 = max ( l0, len ( label       ) )
                l1 = max ( l1, len ( identifiers ) )
                l2 = max ( l2, len ( objects     ) )
                data.append ( ( label, identifiers, objects ) )
            data.sort ( )
            # . Output.
            table = log.GetTable ( columns = [ l0+2, l1+2, l2 ] )
            table.Start   ( )
            table.Title   ( title )
            table.Heading ( "Format Label"      )
            table.Heading ( "Other Identifiers" )
            table.Heading ( "Objects"           )
            for ( label, identifiers, objects ) in data:
                table.Entry ( label      , alignment = "l" )
                table.Entry ( identifiers, alignment = "l" )
                table.Entry ( objects    , alignment = "l" )
            table.Stop ( )

#===================================================================================================================================
# . Instantiate an exporter and an importer.
#===================================================================================================================================
# . Instantiation.
_Exporter = ExportImportHandlerContainer ( )
_Importer = ExportImportHandlerContainer ( )

# . Add basic formats - from outside pBabel.
# . None means that objects of all classes are handled.
_Exporter.AddHandler ( { None :       Pickle }, [ "pkl" , "PKL"  ], "Pickle"      )
_Exporter.AddHandler ( { None :   YAMLPickle }, [ "yaml", "YAML" ], "YAML Pickle" )
_Importer.AddHandler ( { None :     Unpickle }, [ "pkl" , "PKL"  ], "Pickle"      )
_Importer.AddHandler ( { None : YAMLUnpickle }, [ "yaml", "YAML" ], "YAML Pickle" )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
# . In these functions the keyword arguments include "format", "log" and the arguments
# . applicable to the export/import function found by the handler. When function-specific
# . arguments are required it is probably best to use the appropriate function directly
# . and avoid these generic export/import functions.

def ExportFileFormats ( ):
    """Return a list of the export file formats."""
    return _Exporter.FileFormats ( )

def ExportOptions ( log = logFile ):
    """List the export options."""
    _Exporter.HandlerSummary ( log = log, title = "Export Options" )

def ExportSystem ( path, system, **options ):
    """Export a system."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Exporter.GetHandler ( path, format = format )
    handler.ExportObject ( path, system, **options )
    if LogFileActive ( log ):
        log.Paragraph ( "System exported to \"{:s}\" in {:s} format.".format ( path, handler.label ) )

def ImportCoordinates3 ( path, **options ):
    """Import a set of coordinates."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Importer.GetHandler ( path, format = format )
    item    = handler.ImportObject ( path, Coordinates3, **options )
    if LogFileActive ( log ):
        log.Paragraph( "Coordinates3 imported from \"{:s}\" in {:s} format.".format ( path, handler.label ) )
    return item

def ImportFileFormats ( ):
    """Return a list of the import file formats."""
    return _Importer.FileFormats ( )

def ImportObjects ( path, **options ):
    """Import objects from the file."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Importer.GetHandler ( path, format = format )
    item    = handler.ImportObject ( path, None, **options )
    if LogFileActive ( log ):
        log.Paragraph( "Objects imported from \"{:s}\" in {:s} format.".format ( path, handler.label ) )
    return item

def ImportOptions ( log = logFile ):
    """List the import options."""
    _Importer.HandlerSummary ( log = log, title = "Import Options" )

def ImportSystem ( path, **options ):
    """Import a system."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Importer.GetHandler ( path, format = format )
    item    = handler.ImportObject ( path, System, **options )
    if LogFileActive ( log ):
        log.Paragraph( "System imported from \"{:s}\" in {:s} format.".format ( path, handler.label ) )
    return item

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
