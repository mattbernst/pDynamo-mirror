#-------------------------------------------------------------------------------
# . File      : pCore.CubicSpline.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle cubic splines."""

#from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CubicSpline:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef CubicSpline new
        new = self.__class__.Raw ( )
        CubicSpline_Clone ( self.cObject, &new.cObject )
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            CubicSpline_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.CubicSpline"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )

#    def __reduce_ex__ ( self, protocol ):
#        """Pickling protocol."""
#        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def AbscissaValue ( self, Integer i ):
        """Return an abscissa value."""
        if ( i < 0 ) or ( i >= self.NumberOfPoints ( ) ): raise IndexError ( "Point index out of range." )
        return self.cObject.x.data[i]

    # . Better handling of errors here if x is out of range.
    def Evaluate ( self, Real x ):
        """Evaluate the spline given an abscissa value."""
        cdef Real   f, g, h
        cdef Status status
        status = CubicSpline_Evaluate ( self.cObject, x, &f, &g, &h )
#        if status != Status_Null: print "Status>", status, x
        return ( f, g, h )

    def FindMaxima ( self ):
        """Find the maxima."""
        cdef Integer i
        cdef CReal1DArray *cmaxima
        maxima  = []
        cmaxima = NULL
        CubicSpline_FindExtrema ( self.cObject, &cmaxima, NULL )
        if cmaxima != NULL:
            for i from 0 <= i < cmaxima.length:
                x = cmaxima.data[i]
                f = self.Evaluate ( x )[0]
                maxima.append ( ( x, f ) )
            Real1DArray_Deallocate ( &cmaxima )
        return maxima

    @classmethod
    def FromFunctionValues ( selfClass, y, lowerDerivative = 2, lowerValue = 0.0, upperDerivative = 2, upperValue = 0.0, x = None ):
        """Constructor."""
        cdef CReal1DArray *xArray, *yArray
        cdef CubicSpline   self
        cdef Integer       i, lD, size, uD
        cdef Real          lV, uV
        cdef Status        status
        # . Initialization.
        self = selfClass.Raw ( )
        # . Length of spline.
        size = len ( y )
        # . Get X.
        xArray = NULL
        if x is not None:
            if len ( x ) != size: raise ValueError ( "Incompatible x and y dimensions." )
            xArray = Real1DArray_Allocate ( size, NULL )
            if xArray == NULL: raise ValueError ( "Out of memory - x." )
            for ( i, v ) in enumerate ( x ): xArray.data[i] = v
        # . Get Y.
        yArray = Real1DArray_Allocate ( size, NULL )
        if yArray == NULL: raise ValueError ( "Out of memory - y." )
        for ( i, v ) in enumerate ( y ): yArray.data[i] = v
        # . Options.
        lD = lowerDerivative
        lV = lowerValue
        uD = upperDerivative
        uV = upperValue
        # . Make the spline.
        status = CubicSpline_MakeFromReal1DArrays ( &self.cObject, &xArray, &yArray, lD, lV, uD, uV )
#        if status != Status_Null: print "Status>", status
        if self.cObject != NULL: self.isOwner = True
        # . Finish up.
        return self

    def Integrate ( self, Real a, Real b ):
        """Integrate the spline between the limits [a,b]."""
        return CubicSpline_Integrate ( self.cObject, a, b, NULL )

    def IntegrateFull ( self ):
        """Integrate the spline over its full range."""
        return CubicSpline_IntegrateFull ( self.cObject, NULL )

    def NumberOfPoints ( self ):
        """The number of points."""
        n = 0
        if self.cObject != NULL: n = self.cObject.length
        return n

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    # . Properties.
    property lower:
        def __get__ ( self ): return self.AbscissaValue ( 0 )
    property upper:
        def __get__ ( self ): return self.AbscissaValue ( self.NumberOfPoints ( ) - 1 )
