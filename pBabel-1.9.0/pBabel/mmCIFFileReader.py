#-------------------------------------------------------------------------------
# . File      : mmCIFFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Read data from a mmCIF file."""

from pCore        import Coordinates3, logFile, LogFileActive, TextFileReader
from ExportImport import _Importer
from pMolecule    import PeriodicTable, System

from mmCIFModel import mmCIFModel, mmCIFModelAtom, mmCIFModelComponent, mmCIFModelEntity

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Continuation character.
_COMMENTCHARACTER = "#"

# . Continuation character.
_CONTINUATIONCHARACTER = ";"

# . Data token.
_DATATOKEN = "data_"

# . Default model number.
_DEFAULTMODELNUMBER = "1"

# . Key character.
_KEYCHARACTER = "_"

# . Loop token.
_LOOPTOKEN = "loop_"

# . Quotes characters.
_DOUBLEQUOTES = "\""
_SINGLEQUOTES = "'"

# . Summary page width.
_SUMMARYPAGEWIDTH = 100

# . Undefined character.
_UNDEFINEDCHARACTER = "."

# . Unknown character.
_UNKNOWNCHARACTER = "?"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFFileReader ( TextFileReader ):
    """mmCIFFileReader is the class for reading mmCIF files."""

    defaultattributes = { "nextline" : None }
    defaultattributes.update ( TextFileReader.defaultattributes )

    def GenerateModels ( self, log = logFile ):
        """Generate mmCIF models from the parsed data."""
        if self.QPARSED:
            self.mmcifmodels = {}
            for ( dataname, datablock ) in self.datablocks.iteritems ( ):
                # . Get the various tables to be processed.
                atomsite = datablock.get ( "_atom_site", None )
                entity   = datablock.get ( "_entity",    None )
                # . Do nothing if there are no atoms.
                if ( atomsite is None ) or ( atomsite.NumberOfRows ( ) <= 0 ): continue
                # . Define the basic model.
                model = mmCIFModel ( label = dataname, log = log )
                # . Add elements to the model (order is important).
                for ( elementname, attributename ) in ( ( "_entity", "EntityDefinition" ), ( "_struct_asym", "AsymmetricUnitDefinition" ), ( "_atom_site", "AtomSite" ), ( "_struct_conn", "Connection" ) ):
                    element  = datablock.get ( elementname, None )
                    function = getattr ( model, "Add" + attributename, None )
                    if ( element is not None ) and ( function is not None ):
                        for row in element.rows: function ( **row )
                # . Finish up.
                model.Finalize ( )
                self.mmcifmodels[dataname] = model
#                print self.mmcifmodels

    def GetLine ( self, QWARNING = False ):
        """Get a line of non-zero length."""
        try:
            if self.nextline is None:
                while True:
                    line = next ( self.file ).strip ( )
                    self.nlines += 1
                    if len ( line ) > 0: break
            else:
                line = self.nextline
                self.nextline = None
                self.nlines  += 1
            return line
        except:
            if QWARNING: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def GetModel ( self, datablockname = None ):
        """Get a model."""
        model = None
        if self.QPARSED:
            mmcifmodels = getattr ( self, "mmcifmodels", {} )
            if datablockname is None:
                if   len ( mmcifmodels ) == 1: model = mmcifmodels.values ( )[0]
                elif len ( mmcifmodels ) == 0: raise ValueError ( "There are no models in the data block." )
                else: raise ValueError ( "A data block name must be specified when there are multiple models." )
            else:
                model = mmcifmodels.get ( datablockname, None )
                if model is None: raise ValueError ( "Unable to find model with the name " + datablockname + "." )
        return model

    def GetTokens ( self, converters = None, separator = None, QWARNING = False ):
        """Get the tokens on a line (and continuation lines if there are any)."""
        # . Initialization.
        tokens = None
        while True:
            # . First line.
            line   = self.GetLine ( )
            tokens = self.TokenizeLine ( line )
            # . Continuation lines.
            if len ( tokens ) > 0:
                line = self.GetLine ( )
                if line.startswith ( _CONTINUATIONCHARACTER ):
                    partialtoken = line[1:].strip ( )
                    if len ( partialtoken ) > 0: continuations = [ partialtoken ]
                    else:                        continuations = [ ]
                    while True:
                        line = self.GetLine ( )
                        if line == _CONTINUATIONCHARACTER:
                            if len ( continuations ) > 0: tokens.append ( " ".join ( continuations ) )
                            break
                        else:
                            partialtoken = line.strip ( )
                            if len ( partialtoken ) > 0: continuations.append ( partialtoken )
                else:
                    self.nextline = line
                    self.nlines  -= 1
                    break
        return tokens

    def Parse ( self, log = logFile ):
        """Parse data from the file."""
        if not self.QPARSED:
            # . Initialization.
            self.datablocks = {}
            self.datablock  = {}
            if LogFileActive ( log ): self.log = log
            # . Start parsing.
            self.Open ( )
            try:
                # . Initialization.
                QLOOPHEADER = False
                QLOOPBODY   = False
                loopkeys    = None
                looptable   = None
                # . Start looping.
                while True:
                    tokens = self.GetTokens ( )
                    token  = tokens[0]
                    # . New data set.
                    if token.startswith ( _DATATOKEN ):
                        QLOOPHEADER = False
                        QLOOPBODY   = False
                        loopkeys    = None
                        looptable   = None
                        self.ProcessDataToken ( token )
                    # . Start loop.
                    elif token == _LOOPTOKEN:
                        QLOOPHEADER = True
                        QLOOPBODY   = False
                        loopkeys    = []
                        looptable   = None
                    # . In a loop header.
                    elif QLOOPHEADER:
                        if token.startswith ( _KEYCHARACTER ):
                            looptable = self.ProcessLoopHeader ( looptable, loopkeys, token )
                        else:
                            QLOOPHEADER = False
                            QLOOPBODY   = True
                            self.ProcessLoopBody ( looptable, loopkeys, tokens )
                    # . In a loop body.
                    elif QLOOPBODY:
                        if token.startswith ( _KEYCHARACTER ):
                            QLOOPHEADER = False
                            QLOOPBODY   = False
                            loopkeys    = []
                            looptable   = None
                            self.ProcessKeyValue ( token, tokens[1:] )
                        else:
                            self.ProcessLoopBody ( looptable, loopkeys, tokens )
                    # . Non-loop data.
                    elif token.startswith ( _KEYCHARACTER ):
                        self.ProcessKeyValue ( token, tokens[1:] )
                    # . Unrecognized tokens.
                    else:
                        self.Warning ( "Unrecognized data line.", False )
            except EOFError:
                pass
            # . Close the file and warning table.
            self.WarningStop ( )
            self.Close ( )
            # . Everything is now parsed.
            self.log     = None
            self.QPARSED = True

    def ParseKey ( self, key ):
        """Parse a key for table and column names."""
        items = key.split ( ".", 1 )
        if len ( items ) == 1: return ( items[0].lower ( ), ""                 )
        else:                  return ( items[0].lower ( ), items[1].lower ( ) )

    def ProcessDataToken ( self, token ):
        """Process a data token."""
        # . Get the data token key.
        key = token[len(_DATATOKEN):].upper ( )
        if key in self.datablocks: self.Warning ( "Duplicate data block keys: " + key + ".", True )
        self.datablocks[key] = self.datablock = {}

    def ProcessKeyValue ( self, key, values ):
        """Process a unique key value pair."""
        ( tablename, columnname ) = self.ParseKey ( key )
        table = self.datablock.get ( tablename, None )
        if table is None:
            table = mmCIFFileTable ( )
            self.datablock[tablename] = table
        if table.HasColumn ( columnname ): self.Warning ( "Table in single-entry mode already has column " + columnname + ".", True )
        table.SetValue ( key, values )

    def ProcessLoopBody ( self, looptable, loopkeys, values ):
        """Process loop body."""
        if looptable is None: self.Warning ( "Loop defined without data definitions.", True )
        else: looptable.AddRow ( loopkeys, values )

    def ProcessLoopHeader ( self, looptable, loopkeys, key ):
        """Process loop header."""
        # . Treat the key.
        ( tablename, columnname ) = self.ParseKey ( key )
        # . Treat the tables.
        table = self.datablock.get ( tablename, None )
        if ( table is None ) and ( looptable is None ):
            table = mmCIFFileTable ( )
            self.datablock[tablename] = table
        elif table is looptable: pass
        else: self.Warning ( "Invalid loop header definition.", True )
        # . Treat the columns.
        if table is not None:
            if table.HasColumn ( columnname ): self.Warning ( "Loop header with duplicate columns " + key + ".", True )
            table.AddColumn    ( columnname )
        loopkeys.append ( columnname )
        return table

    def Summary ( self, log = logFile ):
        """Print a summary of the stored data."""
        if self.QPARSED and LogFileActive ( log ):
            # . Heading.
            if self.name is None: log.Heading ( "mmCIF File Summary", qblankline = True )
            else:                 log.Heading ( "Summary for mmCIF File \"" + self.name + "\"", qblankline = True )
            # . Basic data.
            summary = log.GetSummary ( )
            summary.Start ( "mmCIF File Data Summary" )
            summary.Entry ( "Data Blocks" , "{:d}".format ( len ( self.datablocks ) ) )
            summary.Entry ( "Warnings"    , "{:d}".format (       self.nwarnings    ) )
            summary.Stop ( )
            # . Data blocks and models.
            mmcifmodels = getattr ( self, "mmcifmodels", {} )
            keys        = self.datablocks.keys ( )
            keys.sort ( )
            for key in keys:
                summary = log.GetSummary ( pageWidth = _SUMMARYPAGEWIDTH )
                summary.Start ( "Summary for Data Block " + key )
                tablenames = self.datablocks[key].keys ( )
                tablenames.sort ( )
                for tablename in tablenames:
                    summary.Entry ( tablename, "{:d}".format ( self.datablocks[key][tablename].NumberOfRows ( ) ) )
                summary.Stop ( )
                # . Models.
                if key in mmcifmodels: mmcifmodels[key].Summary ( log = log )


    def TokenizeLine ( self, line ):
        """Tokenize a line."""
        basictokens    = line.split ( )
        quotecharacter = None
        stringtokens   = []
        tokens         = []
        for basictoken in basictokens:
            # . Inside string.
            if len ( stringtokens ) > 0:
                # . End of string.
                if basictoken.endswith ( quotecharacter ):
                    stringtokens.append ( basictoken[0:-1].strip ( ) )
                    tokens.append ( " ".join ( stringtokens ) )
                    quotecharacter = None
                    stringtokens   = []
                # . Middle of string.
                else:
                    stringtokens.append ( basictoken )
            # . Outside string.
            # . Start (and possibly the end) of string.
            # . Make sure to check for end as well.
            elif basictoken.startswith ( _DOUBLEQUOTES ):
                if basictoken.endswith ( _DOUBLEQUOTES ):
                    token = basictoken[1:-1].strip ( )
                    if len ( token ) > 0: tokens.append ( token         )
                    else:                 tokens.append ( _UNKNOWNTOKEN )
                else:
                    quotecharacter = _DOUBLEQUOTES
                    stringtokens.append ( basictoken[1:].strip ( ) )
            elif basictoken.startswith ( _SINGLEQUOTES ):
                if basictoken.endswith ( _SINGLEQUOTES ):
                    token = basictoken[1:-1].strip ( )
                    if len ( token ) > 0: tokens.append ( token         )
                    else:                 tokens.append ( _UNKNOWNTOKEN )
                else:
                    quotecharacter = _SINGLEQUOTES
                    stringtokens.append ( basictoken[1:].strip ( ) )
            # . Comment.
            elif basictoken.find ( _COMMENTCHARACTER ) > -1:
                break
            else:
                tokens.append ( basictoken )
        # . Error if a string is not closed.
        if len ( stringtokens ) > 0:
            self.Warning ( "Unmatched quotes (" + quotecharacter + ") on line.", False )
            tokens.append ( " ".join ( stringtokens ) )
        return tokens

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFFileTable ( object ):
    """A class for storing data from a mmCIF file."""

    def __init__ ( self ):
        """Constructor."""
        self.columnnames = set ( )
        self.rows        = []

    def AddColumn ( self, columnname ):
        """Add a column to the table."""
        self.columnnames.add ( columnname )

    def AddRow ( self, keys, values ):
        """Add a row to the table."""
        row = {}
#        print "KEYS:", keys
#        print "VALUES:", values
        for ( key, value ) in zip ( keys, values ):
            if key in self.columnnames: row[key] = value
#        print "COLUMN:", self.columnnames
#        print "ROW:", row
        self.rows.append ( row )

    def GetUniqueRowValues ( self, columnname ):
        """Get the unique row values for a column."""
        values = set ( )
        for row in self.rows:
            if columnname in row: values.add ( row[columnname] )
        return values

    def HasColumn ( self, columnname ):
        """Has the table a column with this name?"""
        return ( columnname in self.columnnames )

    def NumberOfRows ( self ):
        """The number of rows in the table."""
        return len ( self.rows )

    def SetValue ( self, columnname, value ):
        """Assign a value to the table."""
        self.columnnames.add ( columnname )
        if len ( self.rows ) <= 0: self.rows.append ( {} )
        self.rows[-1][columnname] = value

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def mmCIFFile_ToSystem ( filename, altLoc = _UNDEFINEDCHARACTER, datablockname = None, log = logFile, modelNumber = _DEFAULTMODELNUMBER, embeddedHydrogens = False, useComponentLibrary = False ):
    """Helper function that returns a system defined by the atom data in a mmCIF file."""
    # . Parse the file.
    infile = mmCIFFileReader ( filename  )
    infile.Parse             ( log = log )
    infile.GenerateModels    ( log = log )
    infile.Summary           ( log = log )
    # . Get the model.
    model = infile.GetModel ( datablockname = datablockname )
    # . Make the system.
    system = model.ToSystem ( altLoc = altLoc, modelNumber = modelNumber, embeddedHydrogens = embeddedHydrogens )
    # . Print out whether there are undefined coordinates.
    if LogFileActive ( log ) and ( system.coordinates3.numberUndefined > 0 ):
        # . Determine the types of atoms with undefined coordinates.
        nheavy    = 0
        nhydrogen = 0
        for i in system.coordinates3.undefined:
            if system.atoms[i].atomicNumber == 1: nhydrogen += 1
            else:                                 nheavy    += 1
        # . Output a summary.
        summary = log.GetSummary ( )
        summary.Start ( "Undefined Coordinates" )
        summary.Entry ( "Heavy Atoms", "{:d}".format ( nheavy    ) )
        summary.Entry ( "Hydrogens",   "{:d}".format ( nhydrogen ) )
        summary.Stop ( )
    # . Finish up.
    return system

# . Importer definitions.
_Importer.AddHandler ( { System : mmCIFFile_ToSystem } , [ "mmcif", "MMCIF" ], "Macromolecular Crystallographic Information File", defaultFunction = mmCIFFile_ToSystem )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
