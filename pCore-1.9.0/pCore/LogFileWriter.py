#-------------------------------------------------------------------------------
# . File      : LogFileWriter.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes for LogFileWriters."""

import os, os.path, sys, time

from PrintObjects import TextHeading, TextParagraph, TextSummary, TextTable, Verbatim, XHTMLHeading, XHTMLParagraph, XHTMLSummary, XHTMLTable, XHTMLVerbatim
from Time         import CPUTime

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PrintPriority ( object ):
    """A simple class to indicate print priorities.

    This class need never be instantiated by a user.
    """

    def __init__ ( self, priority ):
        """Constructor to assign a print priority."""
        self.priority = priority

"""Predefined print priorities."""
PrintPriority_None     = PrintPriority ( 0 )
PrintPriority_Verylow  = PrintPriority ( 1 )
PrintPriority_Low      = PrintPriority ( 2 )
PrintPriority_Medium   = PrintPriority ( 3 )
PrintPriority_High     = PrintPriority ( 4 )
PrintPriority_Veryhigh = PrintPriority ( 5 )
PrintPriority_Debug    = PrintPriority ( 6 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LogFileWriter ( object ):
    """Base class for all LogFileWriters."""

    defaultAttributes = { "cpuTime"          : None  ,
                          "file"             : None  ,
                          "isActive"         : True  ,
                          "isPriorityFixed"  : False ,
                          "name"             : None  ,
                          "printObjectStack" : None  ,
                          "priority"         : PrintPriority_Medium ,
                          "priorityStack"    : None  }

    # . Public methods.
    def __del__ ( self ):
        """Destructor."""
        self.Close ( )

    def __init__ ( self, fileName = None, bufferSize = 1 ):
        """Create a LogFileWriter object.

        bufferSize sets the buffering level (0 = unbuffered, 1 = line buffered, other positive value = use a buffer of that size).
        """
        # . Set default attributes.
        self.__dict__.update ( self.__class__.defaultAttributes )
        self.printObjectStack = []
        self.priorityStack    = []
        # . Assign the file object.
        if fileName is None:
            self.file = sys.stdout
            self.name = None
        else:
            self.name = fileName
            ( head, tail ) = os.path.split ( fileName )
            if len ( head ) == 0: head = "."
            if not os.access ( head, os.W_OK ): raise IOError ( "Unwriteable file: " + fileName + "." )
            try:
#                self.file = open ( fileName, "w" )
                self.file = open ( fileName, "w", bufferSize ) # . Use instead of os.fdopen.
            except:
                raise IOError ( "Cannot open file: " + self.name + "." )

    def Close ( self ):
        """Close the file."""
        self.PrintObjectsClose ( )
        if self.isActive:
            self.isActive = False
            if self.name is not None:
                try:
                    self.file.close ( )
                except:
                    raise IOError ( "Cannot close file: " + self.name + "." )

    # . Methods that return print objects.
    def GetHeading ( self, **keywordArguments ):
        pass

    def GetParagraph ( self, **keywordArguments ):
        pass

    def GetSummary ( self, **keywordArguments ):
        pass

    def GetTable ( self, **keywordArguments ):
        pass

    def GetVerbatim ( self, **keywordArguments ):
        pass

#    def Error ( self, procedure = None, message = None ):
#        """Error printing."""
#        self.ClosePrintObjects ( )
#        if self.isActive:
#            if ( self.Priority_OK ( PrintPriority_Veryhigh ) ):
#                error =
#
#                if self.qhtml:
#                    self.file.write ( "<p class = \"" + self.error_class + "\">\nGeneric Error\n<br /><br />\n<strong>Program terminating</strong>.\n</p>" )
#                else:
#                    self.file.write ( "\nGeneric Error.\n\nProgram terminating.\n" )

#    def Error_pDynamo ( self, procedure, message ):
#        """Print a pDynamo error message."""
#        if self.isActive:
#            self.Close_Environments ( )
#            if ( self.Priority_OK ( PrintPriority_Veryhigh ) ):
#                if self.qhtml:
#                    self.file.write ( "<p class = \"" + self.error_class + "\">\nError in <strong>" + procedure + "</strong>:\n<br /><br />\n<span class = \"" +
#                                        self.error_class + "\">\n" )
#                    self.Noformat_Start ( self )
#                    self.file.write ( message )
#                    self.Noformat_Stop  ( self )
#                    self.file.write ( "</span>\n<br /><br />\n<strong>Program terminating</strong>.\n</p>" )
#                else:
#                    self.file.write ( "\nError in " + procedure + ":\n\n" + message + "\n\nProgram terminating.\n" )

    # . Print object methods.
    def PrintObjectsClose ( self ):
        """Close all print objects."""
        if self.isActive:
            for item in self.printObjectStack:
                item.Stop ( )

    def PrintObjectPop ( self ):
        """Remove a print object from the stack."""
        if len ( self.printObjectStack ) > 0: self.printObjectStack.pop ( )

    def PrintObjectPush ( self, item ):
        """Add a print object to the stack."""
        self.printObjectStack.append ( item )

    # . Priority methods.
    def PriorityLock ( self ):
        """Prevent changes to the print priority.

        This is useful for debugging.
        """
        self.isPriorityFixed = True

    def PriorityOK ( self, inputpriority ):
        """Check that the print priority is reached."""
        return ( inputpriority.priority >= self.priority.priority )

    def PriorityPop ( self ):
        """Pop a priority from the stack."""
        if not self.isPriorityFixed:
            if ( len ( self.priorityStack ) > 0 ): priority = self.priorityStack.pop ( )

    def PriorityPush ( self, newPriority ):
        """Push a priority onto the stack."""
        if ( not self.isPriorityFixed ):
            self.priorityStack.append ( newPriority )
            self.priority = newPriority

    def PriorityUnlock ( self ):
        """Allow changes to the print priority."""
        self.isPriorityFixed = False

    # . Output methods.
    def Footer ( self ):
        """Write a footer and close the file."""
        if self.cpuTime is not None: self.Paragraph ( "CPU Time: " + self.cpuTime.CurrentAsString ( ) )
        self.Paragraph ( "Stop Time: " + time.asctime ( time.localtime ( None ) ) )
        self.Close ( )

    def Header ( self, title ):
        """Write a header."""
        self.cpuTime = CPUTime ( )
        if title is not None: self.Heading   ( title )
        self.Paragraph ( "Start Time: " + time.asctime ( time.localtime ( None ) ) )

    def Heading ( self, text, **keywordArguments ):
        """Heading."""
        item = self.GetHeading ( **keywordArguments )
        item.Start ( )
        item.Text  ( text )
        item.Stop  ( )

    def LineBreak ( self ):
        """Line break."""
        pass

    def Paragraph ( self, text, **keywordArguments ):
        """Paragraph."""
        item = self.GetParagraph ( **keywordArguments )
        item.Start ( )
        item.Text  ( text )
        item.Stop  ( )

    def Separator ( self ):
        """Separator."""
        pass

    def Text ( self, text ):
        """Text."""
        if self.isActive and ( text is not None ): self.file.write ( text )

    def Verbatim ( self, text, **keywordArguments ):
        """Verbatim."""
        item = self.GetVerbatim ( **keywordArguments )
        item.Start ( )
        item.Text  ( text )
        item.Stop  ( )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class TextLogFileWriter ( LogFileWriter ):
    """Text LogFileWriters."""

    def GetHeading ( self, **keywordArguments ):
        return TextHeading ( self, **keywordArguments )

    def GetParagraph ( self, **keywordArguments ):
        return TextParagraph ( self, **keywordArguments )

    def GetSummary ( self, **keywordArguments ):
        return TextSummary ( self, **keywordArguments )

    def GetTable ( self, **keywordArguments ):
        return TextTable ( self, **keywordArguments )

    def GetVerbatim ( self, **keywordArguments ):
        return Verbatim ( self, **keywordArguments )

    def LineBreak ( self ):
        """Line break."""
        self.Text ( "\n" )

    def Separator ( self ):
        """Separator."""
        self.Text ( 80 * "=" )
        self.LineBreak ( )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class XHTMLLogFileWriter ( LogFileWriter ):
    """XHTML LogFileWriters."""

    def __init__ ( self, fileName = None, styleFile = None, title = None ):
        """Constructor."""
        super ( XHTMLLogFileWriter, self ).__init__ ( fileName )
        # . Basic header.
        self.Text ( "<?xml version = \"1.0\" encoding = \"utf-8\" ?>\n" )
        self.Text ( "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n" )
        self.Text ( "<html xmlns = \"http://www.w3.org/1999/xhtml\" xml:lang = \"en\" lang = \"en\">\n" )
        self.Text ( "<head>\n" )
        if title is not None: self.Text ( "    <title>" + title + "</title>\n" )
        # . Get the style file.
        if styleFile is None: styleFile = os.getenv ( "PDYNAMO_STYLE" )
        if styleFile is None:
            self.Text ( "    <style type = \"text/css\">\n" )
            self.Text ( "    body               { background : #ffffff ; color : #000000 }\n"  )
            self.Text ( "    h1, h2, h3, h4, h5 { text-align : center }\n"                   )
            self.Text ( "    table { margin-left : auto ; margin-right : auto ; margin-top : 15px ; margin-bottom : 15px }\n" )
            self.Text ( "    </style>\n"      )
        else:
            self.Text ( "    <link href = \"" + styleFile + "\" rel = \"stylesheet\" type = \"text/css\" />\n" )
        self.Text ( "</head>\n<body>" )

    def Close ( self ):
        """Close the file."""
        self.PrintObjectsClose ( )
        if self.isActive:
            self.Text ( "</body>\n</html>\n" )
            self.isActive = False
            if self.name is not None:
                try:
                    self.file.close ( )
                except:
                    raise IOError ( "Cannot close file: " + self.name + "." )

    def GetHeading ( self, **keywordArguments ):
        return XHTMLHeading ( self, **keywordArguments )

    def GetParagraph ( self, **keywordArguments ):
        return XHTMLParagraph ( self, **keywordArguments )

    def GetSummary ( self, **keywordArguments ):
        return XHTMLSummary ( self, **keywordArguments )

    def GetTable ( self, **keywordArguments ):
        return XHTMLTable ( self, **keywordArguments )

    def GetVerbatim ( self, **keywordArguments ):
        return XHTMLVerbatim ( self, **keywordArguments )

    def LineBreak ( self ):
        """Line break."""
        self.Text ( "<br/>\n" )

    def Separator ( self ):
        """Separator."""
        self.Text ( "<hr/>\n" )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def LogFileActive ( log, pushPriority = None, referencePriority = PrintPriority_Medium ):
    """Check to see if printing is to be done on a log file.

    log is the log file.
    |pushPriority| is the new priority for the log file.
    |referencePriority| is the priority of the printing to be done.
    """
    isActive = False
    if ( log is not None ) and isinstance ( log, LogFileWriter ):
        if pushPriority is not None: log.PriorityPush ( pushPriority )
        isActive = log.PriorityOK ( referencePriority )
    return isActive

#===================================================================================================================================
# . Default log file.
#===================================================================================================================================
logFile = TextLogFileWriter ( )
#logFile = XHTMLLogFileWriter ( title = "pDynamo Output" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    pass
