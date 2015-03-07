#-------------------------------------------------------------------------------
# . File      : pMolecule.PairwiseInteraction.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Pairwise interaction class."""

import math

from pCore import AdaptiveSimpsonsRule, logFile, LogFileActive, RawObjectConstructor, \
                  UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE, UNITS_LENGTH_ANGSTROMS_TO_BOHRS

# . Have some sort of check or test on spline construction?

#===================================================================================================================================
# . Functions for creation of splines.
#===================================================================================================================================
# . Exponential distributions could be added but the formulae are more complicated and there are many more cases
# . (i.e. width1 = 0.0, width1 = width2, width1 != width2, the limiting case of the latter with width1 ~ width2,
# . and the r->0 versions of all of them).

# . Delta/Delta Coulomb.
def FDelta ( r ): return (   1.0 /   r       )
def GDelta ( r ): return ( - 1.0 / ( r * r ) )

# . Delta/Gaussian and Gaussian/Gaussian Coulomb.
# . f(r) = (a/pi)**(3/2) * exp ( - a * r**2 ).
# . <r> = 2 / sqrt ( a * pi ) so a = 4 / ( pi * <r>**2 ).
# . The p arguments in these functions are "sqrt ( a )" for a Delta/Gaussian interaction
# . or "sqrt ( a * b / ( a + b ) )" for a Gaussian/Gaussian interaction.
# . Small x expansions up to 10th order. At switch over point differences are of the order of 10^-16 (F) and 10^-14/15 (G).
def FGaussian ( p ):
    def F ( r ):
        x  = p * r
        x2 = x * x
        if x > 1.0e-2: f = math.erf ( x ) / x
        else:          f = ( 2.0 - x2 * ( 2.0 / 3.0 - x2 * ( 0.2 - x2 * ( 1.0 / 21.0 - x2 * ( 1.0 / 108.0 - x2 / 660.0 ) ) ) ) ) / math.sqrt ( math.pi )
        return ( f * p )
    return F

def GGaussian ( p ):
    def G ( r ):
        sPi = math.sqrt ( math.pi )
        x   = p * r
        x2  = x * x
        if x > 1.0e-2: g = 2.0 * math.exp ( - x2 ) / ( sPi * x ) - math.erf ( x ) / x2
        else:          g = - x * ( 4.0 / 3.0 - x2 * ( 0.8 - x2 * ( 2.0 / 7.0 - x2 * ( 2.0 / 27.0 - x2 / 66.0 ) ) ) ) / sPi
        return ( g * p * p )
    return G

# . Switching function.
def Switch ( r, rOn, rOff ): return ( ( rOff**2 - r**2 )**2 * ( rOff**2 + 2.0 * r**2 - 3.0 * rOn**2 ) / ( rOff**2 - rOn**2 )**3 )

# . Integration.
# . Seems to be OK for most cases with precisions of at most ~ 1.0e-15.
def Integrate ( F, G, rDamp, rOn, rOff, epsilon, maximumIterations, pointDensity ):
    """Integration."""
    # . Function definition.
    def GS ( r ): return ( G ( r ) * Switch ( r, rOn, rOff ) )
    # . Initialization.
    numberOfPoints = int ( math.ceil ( float ( pointDensity ) * rOff + 1.0 ) )
    x = Real1DArray.WithExtent ( numberOfPoints )
    y = Real1DArray.WithExtent ( numberOfPoints )
    # . Calculate x.
    increment = rOff / float ( numberOfPoints - 1 )
    x[ 0] = 0.0
    for i in range ( 1, numberOfPoints - 1 ): x[i] = increment * float ( i )
    x[-1] = rOff
    # . Calculate y.
    # . Switching region - rOn to rOff.
    integral = 0.0
    iUpper   = numberOfPoints - 2
    y[-1]    = 0.0
    for i in range ( numberOfPoints - 2, -1, -1 ):
        rLower = x[i  ]
        rUpper = x[i+1]
        if rLower < rOn:
            iUpper = i
            if rUpper > rOn:
                ( value, state ) = AdaptiveSimpsonsRule ( GS, rOn, rUpper, epsilon, maximumIterations )
                integral -= value
            break
        ( value, state ) = AdaptiveSimpsonsRule ( GS, rLower, rUpper, epsilon, maximumIterations )
        integral -= value
        y[i]      = integral
    # . Normal and damping regions - 0 to rOn.
    # . Determine some constants.
    offSet = integral - F ( rOn )
    if rDamp > 0:
        fDamp = F ( rDamp ) + offSet
        gDamp = G ( rDamp )
        alpha =   0.5 * gDamp / rDamp
        beta  = - 0.5 * gDamp * rDamp + fDamp
    else:
        alpha = 0.0
        beta  = 0.0
    # . Remaining points.
    for i in range ( iUpper + 1 ):
        r = x[i]
        if r < rDamp: y[i] = alpha * r**2 + beta
        else:         y[i] = F ( r ) + offSet
    # . Square x.
    for i in range ( numberOfPoints ): x[i] = x[i]**2
    # . Finish up.
    return ( x, y )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionABFS ( object ):
    """An ABFS pairwise interaction."""

    def __copy__ ( self ):
        """Copying."""
        options = self.__getstate__ ( )
        new     = self.__class__ ( **options )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner: PairwiseInteractionABFS_Deallocate ( &self.cObject )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.PairwiseInteractionABFS"

    def __getstate__ ( self ):
        """Return the state."""
        return { "dampingCutoff"      : self.cObject.dampingCutoff            ,
                 "electrostaticModel" : self.electrostaticModel               ,
                 "innerCutoff"        : self.cObject.innerCutoff              ,
                 "outerCutoff"        : self.cObject.outerCutoff              ,
                 "splinePointDensity" : self.cObject.splinePointDensity       ,
                 "useAnalyticForm"    : self.cObject.useAnalyticForm == CTrue ,
                 "width1"             : self.width1                           ,
                 "width2"             : self.width2                           }
                 
    def __init__ ( self, **options ):
        """Constructor with options."""
        self._Initialize ( )
        self._Allocate   ( )
        self.SetOptions ( **options )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        self.SetOptions ( **state )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = PairwiseInteractionABFS_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject            = NULL
        self.electrostaticModel = "Delta/Delta"
        self.isOwner            = False
        self.width1             = 0.0
        self.width2             = 0.0

    def CheckOptions ( self ):
        """Check the options."""
        if self.electrostaticModel is None:
            self.electrostaticModel = "Delta/Delta"
        elif self.electrostaticModel not in ( "Delta/Delta", "Delta/Gaussian", "Gaussian/Gaussian" ):
            raise ValueError ( "Invalid pairwise interaction electrostatic model: " + self.electrostaticModel + "." )

    @classmethod
    def FromOptions ( selfClass, **options ):
        """Constructor from options."""
        return selfClass ( **options )

    def MakeElectrostaticSplines ( self, useAtomicUnits = False ):
        """Make the non-delta/delta electrostatic splines using numerical integration."""
        cdef CReal1DArray *cX, *cY
        cdef  Real1DArray   x,   y
        factor = 2.0 / math.sqrt ( math.pi )
        if   self.electrostaticModel == "Delta/Gaussian"    : p = factor / self.width2
        elif self.electrostaticModel == "Gaussian/Gaussian" : p = factor / math.sqrt ( self.width1**2 + self.width2**2 )
        ( x, y ) = Integrate ( FGaussian ( p ), GGaussian ( p ), self.cObject.dampingCutoff, self.cObject.innerCutoff, self.cObject.outerCutoff, 1.0e-15, 1000, self.cObject.splinePointDensity )
        if useAtomicUnits: y.Scale ( 1.0e+00 / UNITS_LENGTH_ANGSTROMS_TO_BOHRS       )
        else:              y.Scale ( UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE )
        cX = x.cObject ; x.isOwner = False ; x = None
        cY = y.cObject ; y.isOwner = False ; y = None
        CubicSpline_Deallocate           ( &self.cObject.electrostaticSpline )
        CubicSpline_MakeFromReal1DArrays ( &self.cObject.electrostaticSpline, &cX, &cY, 1, 0.0, 1, 0.0 )

    def MakeSplines ( self, electrostatic = True, lennardJones = True, useAtomicUnits = False ):
        """Make the splines for the interaction."""
        cdef Boolean cUseAtomicUnits
        if self.cObject != NULL:
            if electrostatic:
                if self.electrostaticModel == "Delta/Delta":
                    if useAtomicUnits: cUseAtomicUnits = CTrue
                    else:              cUseAtomicUnits = CFalse
                    CubicSpline_Deallocate ( &self.cObject.electrostaticSpline )
                    self.cObject.electrostaticSpline = PairwiseInteractionABFS_MakeElectrostaticSpline ( self.cObject, cUseAtomicUnits, NULL )
                else:
                    self.MakeElectrostaticSplines ( useAtomicUnits = useAtomicUnits )
            if lennardJones:
                CubicSpline_Deallocate ( &self.cObject.lennardJonesASpline )
                CubicSpline_Deallocate ( &self.cObject.lennardJonesBSpline )
                self.cObject.lennardJonesASpline = PairwiseInteractionABFS_MakeLennardJonesASpline ( self.cObject, NULL )
                self.cObject.lennardJonesBSpline = PairwiseInteractionABFS_MakeLennardJonesBSpline ( self.cObject, NULL )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        if "dampingCutoff"      in keywordArguments: self.cObject.dampingCutoff      = keywordArguments.pop ( "dampingCutoff"      )
        if "electrostaticModel" in keywordArguments: self.electrostaticModel         = keywordArguments.pop ( "electrostaticModel" )
        if "innerCutoff"        in keywordArguments: self.cObject.innerCutoff        = keywordArguments.pop ( "innerCutoff"        )
        if "outerCutoff"        in keywordArguments: self.cObject.outerCutoff        = keywordArguments.pop ( "outerCutoff"        )
        if "splinePointDensity" in keywordArguments: self.cObject.splinePointDensity = keywordArguments.pop ( "splinePointDensity" )
        if "width1"             in keywordArguments: self.width1                     = keywordArguments.pop ( "width1"             )
        if "width2"             in keywordArguments: self.width2                     = keywordArguments.pop ( "width2"             )
        if "useAnalyticForm"    in keywordArguments:
            value = keywordArguments.pop ( "useAnalyticForm" )
            if value: self.cObject.useAnalyticForm = CTrue
            else:     self.cObject.useAnalyticForm = CFalse
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )
        # . Check models.
        self.CheckOptions ( )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "Pairwise Interaction Summary" )
            summary.Entry ( "Damping Cutoff"       , "{:.3f}".format ( self.cObject.dampingCutoff            ) )
            summary.Entry ( "Electrostatic Model"  , self.electrostaticModel )
            summary.Entry ( "Inner Cutoff"         , "{:.3f}".format ( self.cObject.innerCutoff              ) )
            summary.Entry ( "Outer Cutoff"         , "{:.3f}".format ( self.cObject.outerCutoff              ) )
            summary.Entry ( "Spline Point Density" , "{:d}"  .format ( self.cObject.splinePointDensity       ) )
            summary.Entry ( "Use Analytic Form"    , "{!r}"  .format ( self.cObject.useAnalyticForm == CTrue ) )
            summary.Entry ( "Width 1"              , "{:.3f}".format ( self.width1                           ) )
            summary.Entry ( "Width 2"              , "{:.3f}".format ( self.width2                           ) )
            summary.Stop ( )
