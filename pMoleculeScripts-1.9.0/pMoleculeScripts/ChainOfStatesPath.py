#-------------------------------------------------------------------------------
# . File      : ChainOfStatesPath.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Chain-of-states path class and utilities."""

import math

from pCore import CubicSpline, logFile, LogFileActive, MultiCubicSpline, Real1DArray

from ChainOfStatesObjectiveFunction import ChainOfStatesObjectiveFunction

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ChainOfStatesPath ( object ):
    """A chain-of-states path."""

    # . Objectfunction attributes.
    defaultAttributes = { "angles"                 : None ,
                          "arcLengths"             : None ,
                          "curvatures"             : None ,
                          "dA"                     : None ,
                          "dB"                     : None ,
                          "distances1"             : None ,
                          "distances2"             : None ,
                          "functionProfile"        : None ,
                          "functionSplineA"        : None ,
                          "functionSplineI"        : None ,
                          "functionValues"         : None ,
                          "g"                      : None ,
                          "gI"                     : None ,
                          "gIp"                    : None ,
                          "gP"                     : None ,
                          "imageSpline"            : None ,
                          "maximaA"                : None ,
                          "maximaI"                : None ,
                          "maximaIStructures"      : None ,
                          "numberOfImages"         : 0    ,
                          "numberOfImageVariables" : 0    ,
                          "rmsGradients"           : None ,
                          "rmsGradientSpline"      : None ,
                          "x"                      : None ,
                          "xI"                     : None }

    def _Allocate ( self ):
        """Allocation."""
        # . Image arrays.
        self.angles         = Real1DArray ( self.numberOfImages ) ; self.angles.Set         ( 0.0 )
        self.arcLengths     = Real1DArray ( self.numberOfImages ) ; self.arcLengths.Set     ( 0.0 )
        self.curvatures     = Real1DArray ( self.numberOfImages ) ; self.curvatures.Set     ( 0.0 )
        self.distances1     = Real1DArray ( self.numberOfImages ) ; self.distances1.Set     ( 0.0 )
        self.distances2     = Real1DArray ( self.numberOfImages ) ; self.distances2.Set     ( 0.0 )
        self.functionValues = Real1DArray ( self.numberOfImages ) ; self.functionValues.Set ( 0.0 )
        self.rmsGradients   = Real1DArray ( self.numberOfImages ) ; self.rmsGradients.Set   ( 0.0 )
        # . Work arrays.
        self.dA             = Real1DArray ( self.numberOfImageVariables ) ; self.dA.Set ( 0.0 )
        self.dB             = Real1DArray ( self.numberOfImageVariables ) ; self.dB.Set ( 0.0 )
        # . Variable arrays - full.
        self.g              = Real1DArray ( self.numberOfImages * self.numberOfImageVariables ) ; self.g.Set  ( 0.0 )
        self.gP             = Real1DArray ( self.numberOfImages * self.numberOfImageVariables ) ; self.gP.Set ( 0.0 )
        self.x              = Real1DArray ( self.numberOfImages * self.numberOfImageVariables ) ; self.x.Set  ( 0.0 )
        # . Variable arrays - slices.
        self.gI = [] ; self.gIp = [] ; self.xI = []
        for image in range ( self.numberOfImages ):
            start = image * self.numberOfImageVariables
            stop  = start + self.numberOfImageVariables
            self.gI.append  ( self.g [start:stop] )
            self.gIp.append ( self.gP[start:stop] )
            self.xI.append  ( self.x [start:stop] )

    def _Initialize ( self ):
        """Initialization."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def __init__ ( self, numberOfImages, numberOfImageVariables ):
        """Constructor."""
        self._Initialize ( )
        self.numberOfImages         = max ( 0, numberOfImages         )
        self.numberOfImageVariables = max ( 0, numberOfImageVariables )
        self._Allocate   ( )

    def BuildFunctionProfile ( self, numberOfPoints = None ):
        """Construct a function profile as a function of arc length."""
        if self.imageSpline is not None:
            if self.functionSplineI is None: self.BuildFunctionSplines ( )
            # . Automatically include optimized images.
            profile = []
            for i in range ( self.numberOfImages ):
                s = float ( i )
                profile.append ( ( self.imageSpline.ArcLength ( 0.0, s ), self.functionSplineI.Evaluate ( s )[0], "*" ) )
            # . Interpolated points.
            if ( numberOfPoints is None ) or ( numberOfPoints <= 0 ): numberOfPoints = 2 * self.numberOfImages
            arcStart  = self.functionSplineA.lower
            arcStop   = self.functionSplineA.upper
            step      = ( arcStop - arcStart ) / float ( numberOfPoints + 1 )
            for i in range ( numberOfPoints ):
                a = arcStart + float ( i + 1 ) * step
                profile.append ( ( a, self.functionSplineA.Evaluate ( a )[0], "" ) )
            # . Sort the profile by increasing arc length.
            profile.sort ( )
            # . Finish up.
            self.functionProfile = profile

    def BuildFunctionSplines ( self ):
        """Build function splines - with respect to arc length and image number."""
        if self.imageSpline is not None:
            # . Get cumulative arclength.
            totalLength = 0.0
            arcLengths  = []
            for a in self.arcLengths:
                arcLengths.append ( totalLength )
                totalLength += a
            # . Match end-point first derivatives.
            self.imageSpline.Tangent ( 0.0                              , tangent = self.dA ) ; dF0 = self.dA.Dot ( self.gIp[0] )
            self.imageSpline.Tangent ( float ( self.numberOfImages - 1 ), tangent = self.dB ) ; dFn = self.dB.Dot ( self.gIp[self.numberOfImages-1] )
            self.functionSplineA   = CubicSpline.FromFunctionValues ( self.functionValues, lowerDerivative = 1, lowerValue = dF0, upperDerivative = 1, upperValue = dFn, x = arcLengths )
            self.functionSplineI   = CubicSpline.FromFunctionValues ( self.functionValues, lowerDerivative = 1, lowerValue = dF0, upperDerivative = 1, upperValue = dFn )
            self.rmsGradientSpline = CubicSpline.FromFunctionValues ( self.rmsGradients )

    def BuildImageSpline ( self ):
        """Build the image spline using image number as the x-values and zero second derivatives at the end points."""
        self.imageSpline = MultiCubicSpline ( self.xI )

    def ComputeAngles ( self ):
        """Compute the angles for the images."""
        self.angles[ 0] = 0.0
        self.angles[-1] = 0.0
        for i in range ( 1, self.numberOfImages - 1 ):
            self.xI[i+1].CopyTo ( self.dA ) ; self.dA.AddScaledArray ( -1.0, self.xI[i  ] )
            self.xI[i  ].CopyTo ( self.dB ) ; self.dB.AddScaledArray ( -1.0, self.xI[i-1] )
            self.dA.Normalize ( ) ; self.dB.Normalize ( )
            dotAB = - self.dA.Dot ( self.dB )
            if   dotAB < -1.0: dotAB = -1.0
            elif dotAB >  1.0: dotAB =  1.0
            self.angles[i] = math.degrees ( math.acos ( dotAB ) )

    def ComputeArcLengths ( self ):
        """Compute arc lengths between images."""
        if self.imageSpline is not None:
            for i in range ( self.numberOfImages - 1 ):
                self.arcLengths[i] = self.imageSpline.ArcLength ( float ( i ), float ( i + 1 ) )
            self.arcLengths[-1] = 0.0

    def ComputeCurvatures ( self ):
        """Compute the curvatures for the images."""
        if self.imageSpline is not None:
            for i in range ( self.numberOfImages ):
                self.curvatures[i] = self.imageSpline.Curvature ( float ( i ) )

    def ComputeDistances ( self ):
        """Compute nearest and next-nearest neighbor distances between images."""
        self.distances1[-1] = 0.0
        self.distances2[-1] = 0.0
        for i in range ( self.numberOfImages - 1 ):
            # . Nearest neighbor (i,i+1).
            self.xI[i].CopyTo ( self.dA )
            self.dA.AddScaledArray ( -1.0, self.xI[i+1] )
            self.distances1[i] = self.dA.Norm2 ( )
            # . Next-nearest neighbor (i,i+2).
            if i < self.numberOfImages - 2:
                self.xI[i].CopyTo ( self.dA )
                self.dA.AddScaledArray ( -1.0, self.xI[i+2] )
                self.distances2[i] = self.dA.Norm2 ( )
            else:
                self.distances2[i] = 0.0

    def Dump ( self, objectiveFunction, images = None ):
        """Dump the path to an objective function."""
        if images is None: images = range ( self.numberOfImages )
        for image in images:
            objectiveFunction.VariablesPut ( self.xI[image] )
            objectiveFunction.DumpImage ( image )

    def Finalize ( self ):
        """Finalization."""
        pass

    def FindMaxima ( self, computeStructures = False ):
        """Find maxima along the path."""
        if self.imageSpline is not None:
            if self.functionSplineA is None: self.BuildFunctionSplines ( )
            self.maximaA = self.functionSplineA.FindMaxima ( )
            self.maximaI = self.functionSplineI.FindMaxima ( )
            if computeStructures:
                structures = []
                for ( s, fs ) in self.maximaI:
                    v = Real1DArray.WithExtent ( self.imageSpline.NumberOfSplines ( ) )
                    v.Set ( 0.0 )
                    self.imageSpline.Evaluate ( s, f = v )
                    structures.append ( ( s, fs, v ) )
                self.maximaIStructures = structures

    def Irregularity ( self ):
        """Return the maximum distance ratio between images."""
        if self.imageSpline is None: return ( max ( self.distances1 ) ) / ( min ( self.distances1[:-1] ) )
        else:                        return ( max ( self.arcLengths ) ) / ( min ( self.arcLengths[:-1] ) )

    def Load ( self, objectiveFunction, images = None ):
        """Load the path from an objective function."""
        if images is None: images = range ( self.numberOfImages )
        for image in images:
            objectiveFunction.LoadImage ( image )
            objectiveFunction.VariablesGet ( self.xI[image] )

    # . Need to save new structures on trajectory.
    def RedistributeImages ( self ):
        """Redistribute images along the spline."""
        if self.imageSpline is not None:
            # . Get the new image locations.
            totalLength = sum ( self.arcLengths )
            # . Calculate the new interpolated structures (except at the end points).
            for i in range ( 1, self.numberOfImages - 1 ):
                arcLength = float ( i ) * totalLength / float ( self.numberOfImages - 1 )
                try:
                    s = self.imageSpline.FindAbscissaFromArcLength ( arcLength )
                except:
                    print ( "{:25.15f} {:25.15f} {:25.15f}".format ( totalLength, arcLength, sum ( self.imageSpline.arclengths ) ) )
                    self.arcLengths.Print ( )
                    self.imageSpline.arclengths.Print ( )
                    raise
                self.imageSpline.Evaluate ( s, f = self.xI[i] )
            # . Update the spline data.
            self.UpdateSplineData ( )

    def Summary ( self, log = logFile ):
        """Summary of a log of path data."""
        if LogFileActive ( log ):
            self.BuildFunctionSplines ( )
            self.BuildFunctionProfile ( )
            self.FindMaxima           ( )
            self.SummaryPath            ( log = log )
            self.SummaryFunctionProfile ( log = log )
            self.SummaryMaxima          ( log = log )

    def SummaryFunctionProfile ( self, log = logFile ):
        """Output the function profile."""
        # . Output.
        if ( self.functionProfile is not None ) and LogFileActive ( log ):
            table  = log.GetTable ( columns = [ 10, 20, 10 ] )
            table.Start   ( )
            table.Title   ( "Arc Length Interpolation" )
            table.Heading ( "Arc Length" )
            table.Heading ( "Function"   )
            table.Heading ( "Optimized"  )
            for ( a, f, tag ) in self.functionProfile:
                table.Entry ( "{:.3f}".format ( a ) )
                table.Entry ( "{:.3f}".format ( f ) )
                table.Entry ( tag )
            table.Stop ( )

    def SummaryMaxima ( self, log = logFile ):
        """Output the maxima."""
        # . Output.
        if LogFileActive ( log ) and ( self.maximaA is not None ) and ( self.maximaI is not None ):
            # . Consistency between two profiles.
            if len  ( self.maximaA ) == len ( self.maximaI ):
                log.Paragraph ( "Number of maxima found = {:d}".format ( len ( self.maximaI ) ) )
                table  = log.GetTable ( columns = [ 20, 20, 20 ] )
                table.Start   ( )
                table.Title   ( "Maxima Summary" )
                table.Heading ( "Structure"      )
                table.Heading ( "Arc Length"     )
                table.Heading ( "Function"       )
                for ( ( a, fa ), ( s, fs ) ) in zip ( self.maximaA, self.maximaI ):
                    table.Entry ( "{:.3f}".format ( s  ) )
                    table.Entry ( "{:.3f}".format ( a  ) )
                    table.Entry ( "{:.3f}".format ( fs ) )
                table.Stop ( )
            # . Inconsistency between profiles.
            else:
                for ( tag, maxima ) in ( ( "Arc Length", self.maximaA ), ( "Structure", self.maximaI ) ):
                    log.Paragraph ( "Number of {:s} maxima found = {:d}".format ( tag.lower ( ), len ( maxima ) ) )
                    table  = log.GetTable ( columns = [ 20, 20 ] )
                    table.Start   ( )
                    table.Title   ( tag + " Maxima Summary" )
                    table.Heading ( "Abscissa" )
                    table.Heading ( "Function" )
                    for ( x, f ) in maxima:
                        table.Entry ( "{:.3f}".format ( x ) )
                        table.Entry ( "{:.3f}".format ( f ) )
                    table.Stop ( )

    def SummaryPath ( self, log = logFile ):
        """Output a path summary."""
        if LogFileActive ( log ):
            # . Initialization.
            fMinimum       = min ( self.functionValues )
            hasSpline      = ( self.imageSpline is not None )
            totalArcLength = 0.0
            totalDistance  = 0.0
            # . Output.
            columns = [ 10 ] + 8 * [ 18 ]
            if hasSpline: columns.extend ( [ 18, 18 ] )
            table  = log.GetTable ( columns = columns )
            table.Start   ( )
            table.Title   ( "Path Summary"      )
            table.Heading ( "Image"             )
            table.Heading ( "Abs. Function"     )
            table.Heading ( "Rel. Function"     )
            table.Heading ( "RMS Gradient"      )
            table.Heading ( "Angle"             )
            if hasSpline: table.Heading ( "Arc Length (i,i+1)" )
            table.Heading ( "Curvature"         )
            table.Heading ( "Distance  (i,i+1)" )
            table.Heading ( "Distance  (i,i+2)" )
            if hasSpline: table.Heading ( "Total Arc Length"   )
            table.Heading ( "Total Distance"    )
            for i in range ( self.numberOfImages ):
                table.Entry ( "{:d}".format ( i ) )
                table.Entry ( "{:.3f}".format ( self.functionValues[i]            ) )
                table.Entry ( "{:.3f}".format ( self.functionValues[i] - fMinimum ) )
                table.Entry ( "{:.3f}".format ( self.rmsGradients  [i]            ) )
                table.Entry ( "{:.1f}".format ( self.angles        [i]            ) )
                if hasSpline: table.Entry ( "{:.3f}".format ( self.arcLengths[i]  ) )
                table.Entry ( "{:.3f}".format ( self.curvatures    [i]            ) )
                table.Entry ( "{:.3f}".format ( self.distances1    [i]            ) )
                table.Entry ( "{:.3f}".format ( self.distances2    [i]            ) )
                if hasSpline:
                    table.Entry ( "{:.3f}".format ( totalArcLength ) )
                    totalArcLength += self.arcLengths[i]
                table.Entry ( "{:.3f}".format ( totalDistance ) )
                totalDistance += self.distances1[i]
            table.Stop ( )

    # . Need more tangent options.
    def Tangent ( self, image, objectiveFunction, tangent ):
        """Calculate the tangent at a point."""
        # . End points.
        if image == 0:
            self.xI[1].CopyTo ( tangent )
            tangent.AddScaledArray ( -1.0, self.xI[0] )
        elif image == self.numberOfImages - 1:
            self.xI[image].CopyTo ( tangent )
            tangent.AddScaledArray ( -1.0, self.xI[image-1] )
        # . Intermediate points.
        else:
            # . Function values.
            v1 = self.functionValues[image-1]
            v2 = self.functionValues[image  ]
            v3 = self.functionValues[image+1]
            # . Get displacments dA = R(i+1) - R(i), dB = R(i) - R(i-1).
            self.xI[image+1].CopyTo ( self.dA ) ; self.dA.AddScaledArray ( -1.0, self.xI[image  ] )
            self.xI[image  ].CopyTo ( self.dB ) ; self.dB.AddScaledArray ( -1.0, self.xI[image-1] )
            # . Henkelman-Jonsson method.
            v23  = abs ( v3 - v2  ) ; v12  = abs ( v2 - v1  )
            vMax = max ( v12, v23 ) ; vMin = min ( v12, v23 )
            if ( v3 >= v2 ) and ( v2 >= v1 ):
                self.dA.CopyTo ( tangent )
            elif ( v1 >= v2 ) and ( v2 >= v3 ):
                self.dB.CopyTo ( tangent )
            elif ( v3 > v1 ):
                self.dA.CopyTo ( tangent )
                tangent.Scale ( vMax )
                tangent.AddScaledArray ( vMin, self.dB )
            elif ( v3 < v1 ):
                self.dA.CopyTo ( tangent )
                tangent.Scale ( vMin )
                tangent.AddScaledArray ( vMax, self.dB )
            else: raise ValueError ( "Unable to calculate tangent - probably due to numerical errors." )
        # . Project out the linear constraints and normalize.
        objectiveFunction.ApplyLinearConstraints ( tangent )
        tangent.Normalize ( )

    def UpdateNonSplineData ( self, images = None ):
        """Update non-spline path data."""
        self.ComputeAngles    ( )
        self.ComputeDistances ( )

    def UpdateSplineData ( self, images = None ):
        """Update spline path data."""
        self.BuildImageSpline  ( )
        self.ComputeArcLengths ( )
        self.ComputeCurvatures ( )

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def ChainOfStatesPath_RedistributeImages ( target, imageTrajectory, log = logFile ):
    """Redistribute the images along a path in a trajectory."""
    # . Set up the objective function.
    if isinstance ( target, ChainOfStatesObjectiveFunction ): of = target
    else: of = ChainOfStatesObjectiveFunction.FromSystem ( target )
    of.InitializeImages ( imageTrajectory )
    # . Set up the path.
    path = ChainOfStatesPath ( of.NumberOfImages ( ), of.NumberOfVariables ( ) )
    path.Load ( of )
    # . Redistribute.
    path.RedistributeImages ( )
    # . Save.
    path.Dump ( of )
    # . Finish up.
    if LogFileActive ( log ): log.Paragraph ( "Images redistributed along chain-of-states path." )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
