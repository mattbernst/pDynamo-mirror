#-------------------------------------------------------------------------------
# . File      : pCore.Correlation.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Functions for handling correlations."""

from LogFileWriter import logFile, LogFileActive

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def Correlation_AutoCorrelation ( Real1DArray x, Integer correlationLength = -1, normalize = True, removeMean = True, useFFT = False ):
    """Calculate an auto-correlation function."""
    cdef Boolean     cNormalize, cRemoveMean, cUseFFT
    cdef Real1DArray correlation
    # . Get all options.
    if correlationLength <= 0: correlationLength = len ( x ) - 1
    if normalize:  cNormalize  = CTrue
    else:          cNormalize  = CFalse
    if removeMean: cRemoveMean = CTrue
    else:          cRemoveMean = CFalse
    if useFFT:     cUseFFT     = CTrue
    else:          cUseFFT     = CFalse
    # . Calculate the function.
    correlation         = Real1DArray.Raw ( )
    correlation.cObject = Correlation_MakeSimple ( x.cObject, NULL, NULL, cUseFFT, cNormalize, cRemoveMean, correlationLength, NULL, NULL )
    correlation.isOwner = True
    # . Finish up.
    return correlation

def Correlation_CrossCorrelation ( Real1DArray x, Real1DArray y, Integer correlationLength = -1, normalize = True, removeMean = True, useFFT = False ):
    """Calculate a cross-correlation function."""
    cdef Boolean     cNormalize, cRemoveMean, cUseFFT
    cdef Real1DArray correlation
    # . Get all options.
    if correlationLength <= 0: correlationLength = min ( len ( x ) - 1, len ( y ) - 1 )
    if normalize:  cNormalize  = CTrue
    else:          cNormalize  = CFalse
    if removeMean: cRemoveMean = CTrue
    else:          cRemoveMean = CFalse
    if useFFT:     cUseFFT     = CTrue
    else:          cUseFFT     = CFalse
    # . Calculate the function.
    correlation         = Real1DArray.Raw ( )
    correlation.cObject = Correlation_MakeSimple ( x.cObject, y.cObject, NULL, cUseFFT, cNormalize, cRemoveMean, correlationLength, NULL, NULL )
    correlation.isOwner = True
    # . Finish up.
    return correlation

def Correlation_DotProductAutoCorrelation ( Real2DArray x, Integer correlationLength = -1, normalize = True, removeMean = True, useFFT = False ):
    """Calculate a dot-product auto-correlation function."""
    cdef Boolean     cNormalize, cRemoveMean, cUseFFT
    cdef Real1DArray correlation
    # . Get all options.
    if correlationLength <= 0: correlationLength = x.Length ( 0 ) - 1
    if normalize:  cNormalize  = CTrue
    else:          cNormalize  = CFalse
    if removeMean: cRemoveMean = CTrue
    else:          cRemoveMean = CFalse
    if useFFT:     cUseFFT     = CTrue
    else:          cUseFFT     = CFalse
    # . Calculate the function.
    correlation         = Real1DArray.Raw ( )
    correlation.cObject = Correlation_MakeDotProduct ( x.cObject, NULL, NULL, cUseFFT, cNormalize, cRemoveMean, correlationLength, NULL, NULL )
    correlation.isOwner = True
    # . Finish up.
    return correlation

def Correlation_DotProductCrossCorrelation ( Real2DArray x, Real2DArray y, Integer correlationLength = -1, normalize = True, removeMean = True, useFFT = False ):
    """Calculate a dot-product cross-correlation function."""
    cdef Boolean     cNormalize, cRemoveMean, cUseFFT
    cdef Real1DArray correlation
    # . Get all options.
    if correlationLength <= 0: correlationLength = min ( x.Length ( 0 ) - 1, y.Length ( 0 ) - 1 )
    if normalize:  cNormalize  = CTrue
    else:          cNormalize  = CFalse
    if removeMean: cRemoveMean = CTrue
    else:          cRemoveMean = CFalse
    if useFFT:     cUseFFT     = CTrue
    else:          cUseFFT     = CFalse
    # . Calculate the function.
    correlation         = Real1DArray.Raw ( )
    correlation.cObject = Correlation_MakeDotProduct ( x.cObject, y.cObject, NULL, cUseFFT, cNormalize, cRemoveMean, correlationLength, NULL, NULL )
    correlation.isOwner = True
    # . Finish up.
    return correlation
