#-------------------------------------------------------------------------------
# . File      : EnergyTerms.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes for manipulating potential energy terms."""

from pCore import logFile, LogFileActive

#===============================================================================
# . Class.
#===============================================================================
class EnergyTerms ( object ):
    """EnergyTerms defines a container for potential energy terms."""

    def __init__ ( self ):
        """Constructor."""
        self.__grms  = None
        self.__terms = []
        self.__total = 0.0

    def Append ( self, toappend ):
        """Append energy terms to those already defined."""
        if isinstance ( toappend, tuple ) and  ( len ( toappend ) == 2 ) and isinstance ( toappend[0], basestring ) and isinstance ( toappend[1], float ):
            self.__terms.append ( toappend )

    def Extend ( self, toextend ):
        """Add energy terms to those already defined."""
        for term in toextend: self.Append ( term )

    def RMSGradient ( self, value ):
        """Assign a value to the RMS gradient."""
        if isinstance ( value, float ): self.__grms = value
        else:                           self.__grms = None

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( pageWidth = 90, valueWidth = 16 )
            summary.Start ( "Summary of Energy Terms" )
            if self.__total is None: summary.Entry ( "Potential Energy" , "None" )
            else:                    summary.Entry ( "Potential Energy" , "{:16.4f}".format ( self.__total ) )
            if self.__grms  is None: summary.Entry ( "RMS Gradient"     , "None" )
            else:                    summary.Entry ( "RMS Gradient"     , "{:16.4f}".format ( self.__grms  ) )
            for ( name, value ) in self.__terms: summary.Entry ( name, "{:16.4f}".format ( value ) )
            summary.Stop ( )

    def Terms ( self, asDictionary = False ):
        """Return the terms."""
        if asDictionary:
            terms = {}
            for ( name, value ) in self.__terms: terms[name] = value
            return terms
        else:
            return self.__terms

    def Total ( self ):
        """Sum the energy terms to find the total energy."""
        self.__total = 0.0
        for ( name, value ) in self.__terms: self.__total += value

    # . Properties.
    def __GetPotentialEnergy ( self ): return self.__total
    potentialEnergy = property ( __GetPotentialEnergy, None, None, "Potential Energy." )
    def __GetRMSGradient     ( self ): return self.__grms
    rmsGradient     = property ( __GetRMSGradient,     None, None, "RMS Gradient."     )

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
