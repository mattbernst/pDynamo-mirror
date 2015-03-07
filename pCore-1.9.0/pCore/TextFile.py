#-------------------------------------------------------------------------------
# . File      : TextFile.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
#===============================================================================
# . Classes and functions to handle text files.
#===============================================================================

import os, string

from LogFileWriter import logFile, LogFileActive

#===============================================================================
# . Error classes.
#===============================================================================
class TextFileReaderError ( Exception ):
    """Text file reader errors."""
    pass

#===============================================================================
# . Error classes.
#===============================================================================
class TextFileWriterError ( Exception ):
    """Text file writer errors."""
    pass

#===============================================================================
# . Text file class.
#===============================================================================
class TextFile ( object ):
    """TextFile is the base class for text files."""

    defaultattributes = { "file"    : None  ,
                          "mode"    : None  ,
                          "name"    : None  ,
                          "QACTIVE" : False }

    # . Initialization.
    def __init__ ( self, name, mode ):
        """Constructor."""
        self.__dict__.update ( self.__class__.defaultattributes )
        self.mode = mode
        self.name = name

    # . Close the file.
    def Close ( self ):
        """Close the file."""
        if self.QACTIVE:
            try:    self.file.close ( )
            except: raise IOError ( "Cannot close file: " + self.name + "." )
            self.QACTIVE = False

    # . Open the file.
    def Open ( self ):
        """Open the file."""
        if not self.QACTIVE:
            try:    self.file = open ( self.name, self.mode )
            except: raise IOError ( "Cannot open file: " + self.name + "." )
            self.QACTIVE = True

#===============================================================================
# . Text file reader class.
#===============================================================================
class TextFileReader ( TextFile ):
    """TextFileReader is the base class for text files that are to be read."""

    defaultattributes = { "log"             : None,  \
                          "maximumwarnings" : 100,   \
                          "nlines"          : 0,     \
                          "nfatal"          : 0,     \
                          "nwarnings"       : 0,     \
                          "QPARSED"         : False, \
                          "warningtable"    : None   }
    defaultattributes.update ( TextFile.defaultattributes )

    def __init__ ( self, name ):
        """The constructor checks whether the file exists and can be read."""
        if not os.access ( name, os.R_OK ): raise IOError ( "Unreadable file: " + name + "." )
        super ( TextFileReader, self ).__init__ ( name, "r" )

    def GetFixedFormatArray ( self, nitems, itemsperline, itemwidth, converter = None, default = None, QWARNING = True ):
        """Parse a fixed format array that may span several lines."""
        items = None
        if ( nitems > 0 ) and ( itemsperline > 0 ) and ( itemwidth > 0 ):
            items = []
            while nitems > 0:
                itemsonline = min ( nitems, itemsperline )
                line = self.GetFixedFormatLine ( itemsonline * itemwidth, QWARNING = QWARNING )
                for i in range ( itemsonline ):
                    word = line[i*itemwidth:(i+1)*itemwidth].strip ( )
                    if converter is None:
                        items.append ( word )
                    else:
                        item = default
                        if len ( word ) > 0:
                            try:    item = converter ( word )
                            except: self.Warning ( "Unable to convert token " + repr ( i ) + ".", True )
                        items.append ( item )
                nitems -= itemsonline
        return items

    def GetFixedFormatLine ( self, linelength, QWARNING = True ):
        """Get a fixed format line with (at least) the given length."""
        try:
            line = next ( self.file ).rstrip ( ).ljust ( linelength )
            self.nlines += 1
            return line
        except:
            if QWARNING: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def GetFixedFormatTokens ( self, *arguments, **keywordArguments ):
        """Get tokens in fixed format from a line."""
        QWARNING = keywordArguments.get ( "QWARNING", True )
        tokens   = []
        length   = 0
        for arg in arguments: length = max ( length, arg[1] )
        line     = self.GetFixedFormatLine ( length, QWARNING = QWARNING )
        for ( i, ( start, stop, converter, default ) ) in enumerate ( arguments ):
            word  = line[start:stop].strip ( )
            token = default
            if len ( word ) > 0:
                if converter is None:
                    token = word
                else:
                    try:    token = converter ( word )
                    except: self.Warning ( "Unable to convert token " + repr ( i ) + ".", True )
            tokens.append ( token )
        return tokens

    def GetLine ( self, QWARNING = True ):
        """Get a line."""
        try:
            line = next ( self.file ).strip ( )
            self.nlines += 1
            return line
        except:
            if QWARNING: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def GetTokens ( self, converters = None, separator = None, QWARNING = True ):
        """Get and convert tokens on a line."""
        line = self.GetLine ( QWARNING = QWARNING )
        return self.TokenizeLine ( line, converters = converters, separator = separator )

    def Label ( self ):
        """Return a suitable label for the class."""
        return "Text File Reader"

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
	    # . Open the file.
	    self.Open ( )
	    # . Parse all the lines.
            try:
                while True:
                    line = self.GetLine ( )
                    if line is None: break
            except EOFError:
                pass
	    # . Close the file.
            self.WarningStop ( )
	    self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.QPARSED = True

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ) and self.QPARSED:
            ( head, tail ) = os.path.split ( self.name )
            summary = log.GetSummary ( )
            summary.Start ( self.Label ( ) + " Summary for " + "\"" + tail + "\"" )
            summary.Entry ( "Lines Parsed", "{:d}".format ( self.nlines ) )
            summary.Entry ( "Fatal Errors", "{:d}".format ( self.nfatal ) )
            summary.Entry ( "Warnings",     "{:d}".format ( self.nwarnings - self.nfatal ) )
            summary.Stop ( )

    def TokenizeLine ( self, line, converters = None, separator = None ):
        """Tokenize a line with optional converters and separator."""
        tokens = None
        if line is not None:
            if separator is None: tokens = line.split ( )
            else:                 tokens = line.split ( separator )
            if converters is not None:
                for ( i, ( token, converter ) ) in enumerate ( zip ( tokens, converters ) ):
                    if converter is None:
                        new = token
                    else:
                        try:
                            new = converter ( token )
                        except:
                            new = converter ( )
                            self.Warning ( "Unable to convert token " + repr ( i ) + ".", True )
                    tokens[i] = new
        return tokens

    def Warning ( self, message, QFATAL ):
        """Print a warning."""
        if ( self.log is not None ) and ( self.maximumwarnings > 0 ):
            if self.nwarnings == 0: self.WarningStart ( )
            self.nwarnings += 1
            if self.nwarnings <= self.maximumwarnings:
                self.warningtable.Entry ( "{:d}".format ( self.nlines ), alignment = "l" )
                self.warningtable.Entry ( message )
                if QFATAL:
                    self.warningtable.Entry ( "Fatal" )
                    self.nfatal += 1
                else:
                    self.warningtable.Entry ( "Warning" )

    def WarningStart ( self ):
        """Start warning printing."""
        if ( self.log is not None ):
            self.warningtable = self.log.GetTable ( columns = [ 10, 80, 10 ] )
            self.warningtable.Start ( )
            self.warningtable.Title ( "Text File Reader Warnings" )
            self.warningtable.Heading ( "Line"     )
            self.warningtable.Heading ( "Message"  )
            self.warningtable.Heading ( "Severity" )

    def WarningStop ( self ):
        """Terminate warning printing."""
        if ( self.warningtable is not None ):
            self.warningtable.Stop ( )
            self.warningtable = None
            if self.nfatal > 0:
                self.log.Paragraph ( "There have been fatal errors!" )
                raise TextFileReaderError ( "Fatal errors reading text file: " + self.name + "." )
            else:
                self.log.Paragraph ( "There have been warnings. Proceed with caution!" )
            self.log.LineBreak ( )

#===============================================================================
# . Text file writer class.
#===============================================================================
class TextFileWriter ( TextFile ):
    """TextFileWriter is the base class for text files that are to be written."""

    # . Initialization.
    def __init__ ( self, name ):
        """The constructor checks whether the file exists and can be written."""
        ( head, tail ) = os.path.split ( name )
        if len ( head ) == 0: head = "."
        if not os.access ( head, os.W_OK ): raise IOError ( "Unwriteable file: " + name + "." )
        super ( TextFileWriter, self ).__init__ ( name, "w" )

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
