#-------------------------------------------------------------------------------
# . File      : Testing.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Testing classes and functions."""

# . Extends unittests although for the moment only a few of its capabilities are used. This may change if needs get more complicated.

import math, os, unittest

from pCore import logFile, LogFileActive, Pickle, TextLogFileWriter, Unpickle

#===================================================================================================================================
# . Error class.
#===================================================================================================================================
class TestError ( Exception ):
    """Test errors."""
    pass

#===================================================================================================================================
# . Subclasses of unittest.
#===================================================================================================================================
class TestCase ( unittest.TestCase ):
    """A test case."""

    def __init__ ( self, *arguments ):
        """Constructor."""
        super ( TestCase, self ).__init__ ( *arguments )
        self.log                     = logFile
        self.doLong                  = False
        self.fullVerificationSummary = False
        self.generateReferenceData   = False
        self.label                   = self.__module__
        self.outputPath              = None
        self.referenceDataPath       = None
        self.resultPath              = None

    def GenerateReferenceData ( self ): return False

    def GetLog ( self ):
        """Get the log file."""
        if self.outputPath is not None:
            self.log = TextLogFileWriter ( fileName = self.outputPath )
        log = self.log
        if not LogFileActive ( log ): log = None
        return log

    def MakeShort ( self ): return True

#-----------------------------------------------------------------------------------------------------------------------------------
class TestResult ( unittest.TestResult ):
    """Test case results."""
    pass

#===================================================================================================================================
# . Test data.
#===================================================================================================================================
class TestDatum ( object ):
    """Base class for an item of test data."""

    def __init__ ( self, label, value, parent, **keywordArguments ):
        """Constructor."""
	self.label  = label
        self.parent = parent
	self.value  = value

    def __len__ ( self ): return 1

    def CheckInputValues ( self ):
        """Check input values."""
	pass

    def DeviationAsString ( self, value ): return ""

    def ToleranceAsString ( self ): return ""

    def ValueAsString     ( self ): return ""

    def VerifyAgainst ( self, value ):
        """Verify the datum against an input value."""
	return False

#-----------------------------------------------------------------------------------------------------------------------------------
class TestReal ( TestDatum ):
    """A real item of test data."""

    def __init__ ( self, label, value, parent, **keywordArguments ):
        """Constructor."""
	super ( TestReal, self ).__init__ ( label, value, parent )
	for attribute in ( "absoluteErrorTolerance", "percentErrorTolerance", "toleranceFormat", "units", "valueFormat" ):
	    setattr ( self, attribute, keywordArguments.get ( attribute, None ) )
	self.CheckInputValues ( )

    def CheckInputValues ( self ):
        """Check input values."""
        # . Initialization.
        self.usePercentError = False
        # . Get the absoluteErrorTolerance.
        if self.percentErrorTolerance is not None:
            if self.value != 0.0:
                absoluteErrorTolerance = ( self.percentErrorTolerance * self.value ) / 100.0
                if ( self.absoluteErrorTolerance is None ) or ( absoluteErrorTolerance < self.absoluteErrorTolerance ):
                    self.absoluteErrorTolerance = absoluteErrorTolerance
                    self.usePercentError        = True
        # . Overall checks.
	isOK = isinstance ( self.value, float ) and ( self.absoluteErrorTolerance is not None )
        if not isOK: raise TestError ( "Invalid real datum input values." )
        # . Formats.
        if self.toleranceFormat is None: self.toleranceFormat = "{:.4e}"
        if self.valueFormat     is None: self.valueFormat     = "{:.4e}"

    def Deviation ( self, value ):
        """Get the absolute deviation."""
        return ( value - self.value )

    def DeviationAsString ( self, value ):
        """Get the deviation as a string."""
        deviation = self.Deviation ( value )
        if self.usePercentError:
            deviation *= ( 100.0 / self.value )
            string = "{:.1%}".format ( deviation )
        else:
            string = self.toleranceFormat.format ( deviation )
        return string

    def ToleranceAsString ( self ):
        """Get the tolerance as a string."""
        if self.usePercentError: string = "{:.1%}".format ( self.percentErrorTolerance )
        else:                    string = self.toleranceFormat.format ( self.absoluteErrorTolerance )
        return string

    def ValueAsString ( self, value = None ):
        """Get the value as a string."""
        if value is None: value = self.value
        return self.valueFormat.format ( value )

    def VerifyAgainst ( self, value, results ):
        """Verify the datum against an input value."""
	if math.fabs ( self.Deviation ( value ) ) <= self.absoluteErrorTolerance:
            results.AddSuccess ( self.label, self.parent, value )
        else:
            results.AddFailure ( self.label, self.parent, value )

#===================================================================================================================================
# . Test data results.
#===================================================================================================================================
class TestDataResult ( object ):
    """An object to hold the results of data set verifications."""

    def __init__ ( self, dataSet ):
        """Constructor."""
        self.dataSet    = dataSet
	self.failures   = {}
        self.missing    = {}
        self.successes  = {}
        self.unverified = len ( dataSet )

    def AddFailure ( self, label, parent, value ):
        """Add a failure."""
        items = self.failures.get ( parent, {} )
        items[label] = value
        self.failures[parent] = items
        self.unverified -= 1

    def AddMissing ( self, label, parent ):
        """Add a missing result."""
        items = self.missing.get ( parent, [] )
        items.append ( label )
        self.missing[parent] = items

    def AddSuccess ( self, label, parent, value ):
        """Add a success."""
        items = self.successes.get ( parent, {} )
        items[label] = value
        self.successes[parent] = items
        self.unverified -= 1

    @staticmethod
    def Counter ( dictionary ):
        """Count up the number of items in the dictionary."""
        length = 0
        for items in dictionary.values ( ): length += len ( items )
        return length

    def NumberFailures  ( self ): return TestDataResult.Counter ( self.failures  )
    def NumberMissing   ( self ): return TestDataResult.Counter ( self.missing   )
    def NumberSuccesses ( self ): return TestDataResult.Counter ( self.successes )

    def Summary ( self, fullSummary = False, log = logFile ):
        """Write a summary of the results."""
        if LogFileActive ( log ):
            if fullSummary:
                summary = log.GetSummary ( )
                summary.Start ( "Verification Summary for " + self.dataSet.label )
                summary.Entry ( "Failures"  , "{:d}".format ( self.NumberFailures  ( ) ) )
                summary.Entry ( "Missing"   , "{:d}".format ( self.NumberMissing   ( ) ) )
                summary.Entry ( "Successes" , "{:d}".format ( self.NumberSuccesses ( ) ) )
                summary.Entry ( "Unverified", "{:d}".format ( self.unverified          ) )
                summary.Stop ( )
                self.dataSet.ResultsSummary ( log = log, failures = self.failures, missing = self.missing, successes = self.successes )
            else:
                if self.WasSuccessful ( ):
                    log.Paragraph ( "The data set " + self.dataSet.label + " was successfully verified." )
                else:
                    self.dataSet.ResultsSummary ( log = log, failures = self.failures, missing = self.missing )

    def WasSuccessful ( self ):
        """Was the verification successful?"""
        return ( ( len ( self.failures ) == 0 ) and ( len ( self.missing ) == 0 ) )

#===================================================================================================================================
# . Test data sets.
#===================================================================================================================================
class TestDataSet ( object ):
    """A collection of test data items."""

    def __init__ ( self, label, parent = None ):
        """Constructor."""
        self.children = {}
	self.data     = {}
	self.label    = label
        self.parent   = parent

    def __len__ ( self ):
        """Length."""
        length = 0
        for item in self.children.values ( ): length += len ( item )
        for item in self.data.values     ( ): length += len ( item )
        return length

    def AddDatum ( self, datum ):
        """Add a datum."""
        uniqueLabel = ( datum.label not in self.children ) and ( datum.label not in self.data )
        if   isinstance ( datum, TestDataSet ) and uniqueLabel: self.children[datum.label] = datum
        elif isinstance ( datum, TestDatum   ) and uniqueLabel: self.data    [datum.label] = datum
	else: raise TestError ( "Invalid datum added to data set." )

    def Path ( self ):
        """Return the path for the data set."""
        items  = [ self.label ]
        parent = self.parent
        while ( parent is not None ):
            items.append ( parent.label )
            parent = parent.parent
        items.reverse ( )
        return ":".join ( items )

    def ResultsSummary ( self, failures = {}, log = logFile, missing = {}, successes = None ):
        """Write a results summary of the data set."""
        if LogFileActive ( log ):
            if len ( self ) > 0:
                # . Initialization.
                addFailureColumn = False
                # . Gather data.
                localFailures = failures.get ( self, {} )
                items         = dict ( localFailures )
                if successes is not None:
                    localSuccesses = successes.get ( self, {} )
                    items.update ( localSuccesses )
                    addFailureColumn = ( len ( localSuccesses ) < len ( items ) )
                # . Do the summary of successes and failures.
                if len ( items ) > 0:
                    columns = [ 20, 20, 20, 20 ]
                    if addFailureColumn: columns.append ( 8 )
                    table = log.GetTable ( columns = columns )
                    table.Start   ( )
                    table.Title   ( "Results for " + self.Path ( ) )
                    table.Heading ( "Label"     )
                    table.Heading ( "Observed"  )
                    table.Heading ( "Reference" )
                    table.Heading ( "Deviation" )
                    if addFailureColumn: table.Heading ( "Fail" )
                    keys = items.keys ( )
                    keys.sort ( )
                    for key in keys:
                        datum = self.data[key]
                        value = items[key]
                        table.Entry ( datum.label, alignment = "left" )
                        table.Entry ( datum.ValueAsString     ( value ) )
                        table.Entry ( datum.ValueAsString     ( ) )
                        table.Entry ( datum.DeviationAsString ( value ) )
                        if addFailureColumn:
                            if key in localFailures: table.Entry ( "*", alignment = "center" )
                            else:                    table.Entry ( ""  )
                    table.Stop ( )
                # . Missing.
                localMissing = missing.get ( self, [] )
                if len ( localMissing ) > 0:
                    localMissing.sort ( )
                    maximumLabelLength = 0
                    for label in localMissing: maximumLabelLength = max ( maximumLabelLength, len ( label ) )
                    columnWidth   = max ( 20, maximumLabelLength + 2 )
                    numberColumns = min ( len ( localMissing ), 4 )
                    table = log.GetTable ( columns = numberColumns * [ columnWidth ] )
                    table.Start   ( )
                    table.Title   ( "Missing Items in " + self.Path ( ) )
                    for label in localMissing: table.Entry ( label )
                    table.Stop ( )
                # . Header.
#                else:
#                    log.Heading ( self.Path ( ), includeBlankLine = True )
                # . Children.
                if len ( self.children ) > 0:
                    keys = self.children.keys ( )
                    keys.sort ( )
                    for key in keys:
                        self.children[key].ResultsSummary ( failures = failures, log = log, missing = missing, successes = successes )

    def Summary ( self, log = logFile ):
        """Write a summary of the data set."""
        if LogFileActive ( log ):
            if len ( self ) == 0:
                log.Paragraph ( "Data set " + self.label + " is empty." )
            else:
                # . Data.
                if len ( self.data ) > 0:
                    table = log.GetTable ( columns = [ 20, 20, 20 ] )
                    table.Start   ( )
                    table.Title   ( self.Path ( ) )
                    table.Heading ( "Label"     )
                    table.Heading ( "Reference" )
                    table.Heading ( "Tolerance" )
                    keys = self.data.keys ( )
                    keys.sort ( )
                    for key in keys:
                        datum = self.data[key]
                        table.Entry ( datum.label, alignment = "left" )
                        table.Entry ( datum.ValueAsString     ( ) )
                        table.Entry ( datum.ToleranceAsString ( ) )
                    table.Stop ( )
                # . Header.
#                else:
#                    log.Heading ( self.Path ( ), QBLANKLINE = True )
                # . Children.
                if len ( self.children ) > 0:
                    keys = self.children.keys ( )
                    keys.sort ( )
                    for key in keys:
                        self.children[key].Summary ( log = log )

    def VerifyAgainst ( self, dataSet, results = None ):
        """Verify the data set against a data set input as a dictionary."""
        if results is None: results = TestDataResult ( self )
        for ( label, value ) in dataSet.iteritems ( ):
            if   label in self.children: results = self.children[label].VerifyAgainst ( value, results = results )
            elif label in self.data    : self.data[label].VerifyAgainst ( value, results )
	    else: results.AddMissing ( label, self )
        return results

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
