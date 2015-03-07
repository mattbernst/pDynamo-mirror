#-------------------------------------------------------------------------------
# . File      : PrintObjects.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes for print objects."""

#===============================================================================
# . Various default parameters.
#===============================================================================
# . Table data.
_TABLEALIGNMENTCENTER = "c"
_TABLEALIGNMENTLEFT   = "l"
_TABLEALIGNMENTRIGHT  = "r"
_TABLECOLUMNS         = 2

# . Text.
_TEXTPAGEWIDTH         = 80
_TEXTSUMMARYVALUEWIDTH = 14
_TEXTTABLECOLUMNWIDTH  = 20

# . XHTML.
_XHTMLBLANKCHARACTER = "&nbsp;"

#===============================================================================
# . Class.
#===============================================================================
class PrintObject ( object ):
    """Base class for all print objects."""

    defaultattributes = { "owner"   : None, \
                          "QACTIVE" : None  }

    def __init__ ( self, owner, **keywordArguments ):
        """Constructor."""
        self.__dict__.update ( self.__class__.defaultattributes )
        self.owner = owner

    def Start ( self ):
        """Start a print object."""
        if self.QACTIVE is None:
            self.owner.PrintObjectPush ( self )
            self.QACTIVE = True

    def Stop ( self ):
        """Stop a print object."""
        if self.QACTIVE:
            self.owner.PrintObjectPop ()
            self.QACTIVE = False

    def LineBreak ( self ):
        """Print a line break."""
        if self.QACTIVE: self.owner.LineBreak ( )

    def Text ( self, text ):
        """Print some text."""
        if self.QACTIVE: self.owner.Text ( text )

#===============================================================================
# . Class.
#===============================================================================
class XHTMLObject ( object ):
    """Base class for XHTML objects."""

    def ClassAttributeString ( self ):
        """Return a class attribute sting."""
        if hasattr ( self, "classattribute" ) and ( self.classattribute is not None ): return " class = \"" + self.classattribute + "\""
        else:                                                                          return None

#===============================================================================
# . Heading classes.
#===============================================================================
class Heading ( PrintObject ):
    """Base class for heading writing."""

    defaultattributes = dict ( PrintObject.defaultattributes )

class TextHeading ( Heading ):
    """Text heading writing."""

    defaultattributes = { "includeBlankLine" : False,         \
                          "pageWidth"        : _TEXTPAGEWIDTH }
    defaultattributes.update ( Heading.defaultattributes )

    def __init__ ( self, owner, **keywordArguments ):
        """Constructor."""
        super ( TextHeading, self ).__init__ ( owner )
        self.includeBlankLine = keywordArguments.get ( "includeBlankLine", False          )
        self.pageWidth        = keywordArguments.get ( "pageWidth"       , _TEXTPAGEWIDTH )

    def Start ( self ):
        """Activation."""
        super ( TextHeading, self ).Start ( )
        if self.QACTIVE:
            if self.includeBlankLine: self.owner.LineBreak ( )
            self.owner.Text ( self.pageWidth * "-" + "\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.QACTIVE: self.owner.Text ( self.pageWidth * "-" + "\n" )
        super ( TextHeading, self ).Stop ( )

    def Text ( self, text ):
        """Text."""
        length = min ( self.pageWidth, len ( text ) )
        nspaces = ( self.pageWidth - length + 1 ) // 2
        self.owner.Text ( nspaces * " " + text[0:length] + "\n" )

class XHTMLHeading ( Heading, XHTMLObject ):
    """XHTML heading writing."""

    defaultattributes = { "headingTag" :"h1" }
    defaultattributes.update ( Heading.defaultattributes )

    def __init__ ( self, owner, **keywordArguments ):
        """Constructor."""
        super ( XHTMLHeading, self ).__init__ ( owner )
        self.headingTag = keywordArguments.get ( "headingTag", "h1" )

    def Start ( self ):
        """Activation."""
        super ( XHTMLHeading, self ).Start ( )
        if self.QACTIVE:
            self.owner.Text ( "<" + self.headingTag )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.QACTIVE: self.owner.Text ( "\n</" + self.headingTag + ">\n" )
        super ( XHTMLHeading, self ).Stop ( )

#===============================================================================
# . Paragraph classes.
#===============================================================================
class Paragraph ( PrintObject ):
    """Base class for paragraph writing."""
    pass

class TextParagraph ( Paragraph ):
    """Text paragraph writing."""

    def Start ( self ):
        """Activation."""
        super ( TextParagraph, self ).Start ( )
        if self.QACTIVE: self.owner.LineBreak ( )

    def Stop ( self ):
        """Deactivation."""
        if self.QACTIVE: self.owner.LineBreak ( )
        super ( TextParagraph, self ).Stop ( )

class XHTMLParagraph ( Paragraph, XHTMLObject ):
    """XHTML paragraph writing."""

    def __init__ ( self, owner ):
       super ( XHTMLParagraph, self ).__init__ ( owner )

    def Start ( self ):
        """Activation."""
        super ( XHTMLParagraph, self ).Start ( )
        if self.QACTIVE:
            self.owner.Text ( "<p" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.QACTIVE: self.owner.Text ( "\n</p>\n" )
        super ( XHTMLParagraph, self ).Stop ( )

#===============================================================================
# . Summary classes.
#===============================================================================
class Summary ( PrintObject ):
    """Base class for summary writing."""

    defaultattributes = { "nitems" : 0 }
    defaultattributes.update ( PrintObject.defaultattributes )

class TextSummary ( Summary ):
    """Text summary writing."""

    defaultattributes = { "pageWidth"  :        _TEXTPAGEWIDTH, \
                          "tagWidth"   :                     0, \
                          "valueWidth" : _TEXTSUMMARYVALUEWIDTH }
    defaultattributes.update ( Summary.defaultattributes )

    def __init__ ( self, owner, **keywordArguments ):
        """Constructor."""
        super ( TextSummary, self ).__init__ ( owner )
        self.pageWidth  = keywordArguments.get ( "pageWidth",  _TEXTPAGEWIDTH         )
        self.valueWidth = keywordArguments.get ( "valueWidth", _TEXTSUMMARYVALUEWIDTH )
        # . Adjust the sizes.
        self.pageWidth = 2 * ( ( self.pageWidth + 1 ) // 2 )
        self.tagWidth  = ( self.pageWidth - 8 ) // 2 - self.valueWidth

    def Entry ( self, tag, value ):
        """Entry."""
        if self.QACTIVE:
            QNULLENTRY = ( tag   is None ) or ( len ( tag   ) == 0 ) or \
                         ( value is None ) or ( len ( value ) == 0 )
            if not QNULLENTRY:
                if ( self.nitems == 1 ): self.owner.Text ( 2 * " " )
                self.owner.Text ( tag.ljust ( self.tagWidth ) + " = " + value.rjust ( self.valueWidth ) )
            if ( self.nitems == 1 ): self.owner.Text ( "\n" )
            self.nitems = ( self.nitems + 1 ) % 2

    def Start ( self, title ):
        """Activation."""
        super ( TextSummary, self ).Start ( )
        if self.QACTIVE:
            length = min ( self.pageWidth - 2, len ( title ) )
            n  = length + 2
            n1 = ( self.pageWidth - n + 1 ) // 2
            n2 = self.pageWidth - n - n1
            self.owner.Text ( "\n" + n1 * "-" + " " + title[0:length] + " " + n2 * "-" + "\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.QACTIVE:
            if self.nitems > 0: self.Entry ( None, None )
            self.owner.Text ( self.pageWidth * "-" + "\n" )
        super ( TextSummary, self ).Stop ( )

class XHTMLSummary ( Summary, XHTMLObject ):
    """XHTML summary writing."""

    def Entry ( self, tag, value ):
        """Entry."""
        if self.QACTIVE:
            QNULLENTRY = ( tag   == None ) or ( len ( tag   ) == 0 ) or \
                         ( value == None ) or ( len ( value ) == 0 )
            if ( self.nitems == 0 ): self.owner.Text ( "<tr>\n" )
            self.owner.Text ( "<td" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">" )
            if ( QNULLENTRY ): self.owner.Text ( _XHTMLBLANKCHARACTER )
            else:              self.owner.Text ( tag )
            self.owner.Text ( "</td>\n<td" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">" )
            if ( QNULLENTRY ): self.owner.Text ( _XHTMLBLANKCHARACTER )
            else:              self.owner.Text ( value )
            self.owner.Text ( "</td>\n" )
            if ( self.nitems == 1 ): self.owner.Text ( "</tr>\n" )
            self.nitems = ( self.nitems + 1 ) % 2

    def Start ( self, title ):
        """Activation."""
        super ( XHTMLSummary, self ).Start ( )
        if self.QACTIVE:
            self.owner.Text ( "<div" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n<table border = \"1\" cellpadding = \"10\" cellspacing = \"0\">\n<tr><th" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( " colspan = \"4\">" + title + "</th></tr>\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.QACTIVE:
            if self.nitems > 0: self.Entry ( None, None )
            self.owner.Text ( "</table></div>\n" )
        super ( XHTMLSummary, self ).Stop ( )

#===============================================================================
# . Table classes.
#===============================================================================
class Table ( PrintObject ):
    """Base class for table writing."""

    defaultattributes = { "ncolumns" : 0, \
                          "nitems"   : 0  }
    defaultattributes.update ( PrintObject.defaultattributes )

    def Heading ( self, text, columnSpan = 1 ):
        """Heading."""
        self.Entry ( text, alignment = "c", columnSpan = columnSpan, isHeader = True )

    def Title ( self, text ):
        """Title spanning the table."""
        self.EndRow ( )
        self.Entry  ( text, alignment = "c", columnSpan = self.ncolumns, isHeader = True )

class TextTable ( Table ):
    """Text table writing."""

    defaultattributes = { "columnWidths" : 0, \
                          "pageWidth"    : 0  }
    defaultattributes.update ( Table.defaultattributes )

    def __init__ ( self, owner, **keywordArguments ):
        """Constructor."""
        super ( TextTable, self ).__init__ ( owner )
        if "columns" in keywordArguments: self.columnWidths = keywordArguments["columns"]
        else:                   self.columnWidths = [ _TEXTTABLECOLUMNWIDTH for i in range ( _TABLECOLUMNS ) ]
        self.ncolumns  = len ( self.columnWidths )
        self.nitems    = 0
        self.pageWidth = sum ( self.columnWidths )

    def EndRow ( self ):
        """Terminate a row."""
        if self.QACTIVE:
            if self.nitems > 0:
#                self.owner.Text ( "\n" )
                self.nitems = 0

    def Entry ( self, text, alignment = _TABLEALIGNMENTRIGHT, columnSpan = 1, isHeader = False ):
        """Entry."""
        if self.QACTIVE and ( columnSpan > 0 ):
            align  = alignment[0:1].lower ( )
            nstart = self.nitems + 1
            nstop  = min ( self.ncolumns, self.nitems + columnSpan )
            n = sum ( self.columnWidths[nstart-1:nstop] )
            if ( nstart == 1 ): self.owner.Text ( "\n" )
            if ( text is None ) or ( len ( text ) == 0 ): self.owner.Text ( n * " " )
            else:
                nafter  = 0
                nbefore = 0
                nblanks = 0
                if len ( text ) <= n: nblanks = n - len ( text )
                if   align == _TABLEALIGNMENTLEFT : nafter  = nblanks
                elif align == _TABLEALIGNMENTRIGHT: nbefore = nblanks
                else:
                    nbefore = ( nblanks + 1 ) // 2
                    nafter  = nblanks - nbefore
                self.owner.Text ( nbefore * " " + text )
                if nstop < self.ncolumns: self.owner.Text ( nafter * " " )
            if ( nstop == self.ncolumns ) and isHeader: self.owner.Text ( "\n" + self.pageWidth * "-" )
            self.nitems = nstop
            if self.nitems == self.ncolumns: self.nitems = 0

    def Start ( self ):
        """Activation."""
        super ( TextTable, self ).Start ( )
        self.owner.Text ( "\n" + self.pageWidth * "-" )

    def Stop ( self ):
        """Deactivation."""
        if self.QACTIVE:
            self.owner.Text ( "\n" + self.pageWidth * "-" + "\n" )
        super ( TextTable, self ).Stop ( )

class XHTMLTable ( Table, XHTMLObject ):
    """XHTML table writing."""

    def __init__ ( self, owner, **keywordArguments ):
        """Constructor."""
        super ( XHTMLTable, self ).__init__ ( owner )
        if "columns" in keywordArguments: self.ncolumns = len ( keywordArguments["columns"] )
        else:                   self.ncolumns = _TABLECOLUMNS
        self.nitems = 0

    def EndRow ( self ):
        """Terminate a row."""
        if self.QACTIVE:
            if self.nitems > 0:
                for i in range ( self.nitems, self.ncolumns ): self.Entry ( None )
                self.nitems = 0

    def Entry ( self, text, alignment = _TABLEALIGNMENTRIGHT, columnSpan = 1, isHeader = False ):
        """Entry."""
        if self.QACTIVE and ( columnSpan > 0 ):
            align  = alignment[0:1].lower ( )
            nstart = self.nitems + 1
            nstop  = min ( self.ncolumns, self.nitems + columnSpan )
            if nstart == 1: self.owner.Text ( "\n<tr>" )
            if isHeader: self.owner.Text ( "<th" )
            else:       self.owner.Text ( "<td" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( " colspan = \"{:d}\">".format ( columnSpan ) )
            if ( text is None ) or ( len ( text ) == 0 ): self.owner.Text ( _XHTMLBLANKCHARACTER )
            else:                                         self.owner.Text ( text )
            if ( isHeader ): self.owner.Text ( "</th>" )
            else:           self.owner.Text ( "</td>" )
            self.owner.Text ( "\n" )
            if nstop == self.ncolumns: self.owner.Text ( "</tr>\n" )
            self.nitems = nstop
            if self.nitems == self.ncolumns: self.nitems = 0

    def Start ( self ):
        """Activation."""
        super ( XHTMLTable, self ).Start ( )
        if self.QACTIVE:
            self.owner.Text ( "<div" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n<table border = \"1\" cellpadding = \"10\" cellspacing = \"0\">\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.QACTIVE:
            self.EndRow ( )
            self.owner.Text ( "</table></div>\n" )
        super ( XHTMLTable, self ).Stop ( )

#===============================================================================
# . Verbatim classes.
#===============================================================================
class Verbatim ( PrintObject ):
    """Base class for verbatim writing.

    This can be used directly for text output.
    """

    pass

class XHTMLVerbatim ( Verbatim, XHTMLObject ):
    """XHTML verbatim writing."""

    def __init__ ( self, owner ):
       super ( XHTMLVerbatim, self ).__init__ ( owner )

    def Start ( self ):
        """Activation."""
        super ( XHTMLVerbatim, self ).Start ( )
        if self.QACTIVE:
            self.owner.Text ( "<pre" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.QACTIVE: self.owner.Text ( "</pre>\n" )
        super ( XHTMLVerbatim, self ).Stop ( )

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__":

    pass
