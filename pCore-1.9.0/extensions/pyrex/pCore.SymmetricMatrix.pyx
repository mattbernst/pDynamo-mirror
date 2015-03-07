#-------------------------------------------------------------------------------
# . File      : pCore.SymmetricMatrix.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle symmetric matrices."""

from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

# . Symmetric matrices currently do not support views although they could, in principle,
# . be implemented (e.g. with iterators).

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . These data must correspond to the enums in the C source.
# . Updating options.
_UpdatingOptions = { "BFGS"   : SYMMETRICMATRIXUPDATING_BFGS,
                     "BOFILL" : SYMMETRICMATRIXUPDATING_BOFILL,
                     "MS"     : SYMMETRICMATRIXUPDATING_MS,
                     "POWELL" : SYMMETRICMATRIXUPDATING_POWELL }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetricMatrix:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef SymmetricMatrix new
        new = self.__class__.WithExtent ( self.rows )
        self.CopyTo ( new )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            SymmetricMatrix_Deallocate ( &self.cObject )
            self.isOwner = False
        self.owner   = None

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getitem__ ( self, indices ):
        """Get an item."""
        cdef Integer i, j
        cdef Real    value
        cdef Status  status
        try:    ( i, j ) = indices
        except: raise TypeError ( "Expecting two-tuple of integers as indices." )
        status = Status_Continue
        value  = SymmetricMatrix_GetItem ( self.cObject, i, j, &status )
        if status != Status_Continue: raise IndexError ( "Indices ({:d},{:d}) out of range.".format ( i, j ) )
        return value

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.SymmetricMatrix"

    def __getstate__ ( self ):
        """Return the state."""
        items = []
        for i from 0 <= i < self.size:
            items.append ( self.cObject.data[i] )
        return { "items" : items, "shape" : [ self.cObject.dimension, self.cObject.dimension ], "storage" : "LowerTriangleRowMajor" }

    def __init__ ( self, extent ):
        """Constructor with extent."""
        self._Initialize ( )
        self._Allocate ( extent )

    def __len__ ( SymmetricMatrix self ):
        """Return the size of the symmetricmatrix."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, indices, Real value ):
        """Set an item."""
        cdef Integer i, j
        cdef Status  status
        try:    ( i, j ) = indices
        except: raise TypeError ( "Expecting two-tuple of integers as indices." )
        status = Status_Continue
        SymmetricMatrix_SetItem ( self.cObject, i, j, value, &status )
        if status != Status_Continue: raise IndexError ( "Indices ({:d},{:d}) out of range.".format ( i, j ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        items                = state["items"]
        ( extent0, extent1 ) = state["shape"]
        self._Allocate ( extent0 )
        for i from 0 <= i < self.size:
            self.cObject.data[i] = items[i]

    def _Allocate ( self, Integer extent ):
        """Allocation."""
        self.cObject = SymmetricMatrix_Allocate ( extent )
        self.isOwner = True
        self.owner   = None

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.owner   = None

    def AddScaledMatrix ( self, Real alpha, SymmetricMatrix other ):
        """Add a scaled matrix."""
        SymmetricMatrix_AddScaledMatrix ( self.cObject, alpha, other.cObject )

    def Copy_Column_From_Vector ( SymmetricMatrix self, Integer column, Integer mstart, Real1DArray vector, Integer vstart, Integer n ):
        """Set the elements of a column from a vector."""
        SymmetricMatrix_Set_Column_From_Vector ( self.cObject, column, mstart, vector.cObject, vstart, n )

    def Copy_DB_From_DB ( SymmetricMatrix self, Integer istart, SymmetricMatrix other, Integer jstart, Integer n ):
        """Copy the elements of the DB of one matrix to another."""
        SymmetricMatrix_Set_DB_From_SM_DB ( self.cObject, istart, other.cObject, jstart, n )

    def CopyTo ( self, SymmetricMatrix other ):
        """Copying."""
        SymmetricMatrix_CopyTo ( self.cObject, other.cObject )

    def Diagonalize ( SymmetricMatrix self, Real1DArray eigenValues, Real2DArray eigenVectors = None ):
        """Find the eigenvalues and, optionally, the eigenvectors.

        Note that the input matrix is destroyed by this operation.
        """
        if eigenVectors is None: SymmetricMatrix_Diagonalize ( self.cObject, eigenValues.cObject,                 NULL, NULL )
        else:                    SymmetricMatrix_Diagonalize ( self.cObject, eigenValues.cObject, eigenVectors.cObject, NULL )

    def Extent ( self ):
        """Return the extent."""
        return SymmetricMatrix_Dimension ( self.cObject )

    def MaximumAbsoluteValue ( self ):
        """Return the maximum absolute value in the matrix."""
        return SymmetricMatrix_MaximumAbsoluteValue ( self.cObject )

#    def PostMultiply ( self, other, result ):
#        """Post-multiply a matrix by something else."""
#        cdef Matrix          r
#        cdef SymmetricMatrix a
#        if isinstance ( other, SymmetricMatrix ) and isinstance ( result, Matrix ):
#            a = other
#            r = result
#            SymmetricMatrix_Multiply2 ( self.cObject, a.cObject, r.cObject )

    def Print ( self, indexWidth = 6    , itemFormat = "{:18.8f}" , itemsPerRow =  6 , itemWidth  = 19 ,
                      labels     = None , log        = logFile    , title       = None ):
        """Printing."""
        if LogFileActive ( log ):
            dimension = self.cObject.dimension
            if ( labels is None ) or ( len ( labels ) != dimension ):
                labels = []
                for i in range ( dimension ): labels.append ( "{:d}".format ( i ) )
                indexWidth = len ( "{:d}".format ( dimension ) )
            else:
                indexWidth = 0
                for label in labels: indexWidth = max ( len ( label ), indexWidth )
            indexWidth = max ( itemWidth, indexWidth + 2 )
            ( ntables, nlast ) = divmod ( dimension, itemsPerRow )
            if nlast > 0:
                ntables = ntables + 1
            for itable in range ( ntables ):
                if ( itable == ntables - 1 ) and ( nlast > 0 ): nelements = nlast
                else:                                           nelements = itemsPerRow
                columns = [ indexWidth ]
                for i in range ( nelements ): columns.append ( itemWidth )
                table = log.GetTable ( columns = columns )
                table.Start ( )
                if title is not None: table.Title ( title )
                table.Entry ( None )
                for i in range ( nelements ): table.Entry ( labels[i + itable * itemsPerRow], alignment = "center" )
                for irow in range ( dimension ):
                    table.Entry ( labels[irow], alignment = "left" )
                    for i in range ( nelements ): table.Entry (  itemFormat.format ( SymmetricMatrix_GetItem ( self.cObject, irow, i + itable * itemsPerRow, NULL ) ) )
                table.Stop ( )

    def PrintWithCondition ( self, conditionFunction, indexWidth = 6    , itemFormat = "{:18.8f}" , itemWidth  = 19   ,
                                                      labels     = None , log        = logFile    , title      = None ):
        """Printing of items which satisfy a particular condition."""
        if LogFileActive ( log ):
            # . Label widths.
            dimension = self.cObject.dimension
            if ( labels is None ) or ( len ( labels ) != dimension ):
                labels = []
                for i in range ( dimension ): labels.append ( "{:d}".format ( i ) )
                indexWidth = len ( "{:d}".format ( dimension ) )
            else:
                indexWidth = 0
                for label in labels: indexWidth = max ( len ( label ), indexWidth )
            indexWidth = max ( itemWidth, indexWidth + 2 )
            # . Table.
            columns = [ indexWidth, indexWidth, itemWidth ]
            table = log.GetTable ( columns = columns )
            table.Start ( )
            if title is not None: table.Title ( title )
            for i in range ( dimension ):
                for j in range ( i + 1 ):
                    x = SymmetricMatrix_GetItem ( self.cObject, i, j, NULL )
                    if conditionFunction ( x, i, j ):
                        table.Entry ( labels[i], alignment = "left" )
                        table.Entry ( labels[j], alignment = "left" )
                        table.Entry ( itemFormat.format ( x ) )
            table.Stop ( )

    def ProjectOutVectors ( SymmetricMatrix self, Real2DArray vectors ):
        """Project a set of vectors from the matrix."""
        SymmetricMatrix_ProjectOut ( &self.cObject, vectors.cObject, NULL )

    def Raise ( SymmetricMatrix self, Real2DArray vectors, Real value ):
        """Make 'vectors' eigenvectors of the matrix with eigenvalues 'value'."""
        SymmetricMatrix_Raise ( self.cObject, vectors.cObject, value, NULL )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Scale ( SymmetricMatrix self, Real value ):
        """Scale all the elements of a symmetricmatrix."""
        SymmetricMatrix_Scale ( self.cObject, value )

    def Scale_DB ( SymmetricMatrix self, Integer start, Real scale, Integer n ):
        """Scale the elements of a diagonal block."""
        SymmetricMatrix_Scale_DB ( self.cObject, start, n, scale )

    def Scale_OB ( SymmetricMatrix self, Integer istart, Integer jstart, Real scale, Integer ni, Integer nj ):
        """Scale the elements of an off-diagonal block."""
        SymmetricMatrix_Scale_OB ( self.cObject, istart, ni, jstart, nj, scale )

    def Set ( SymmetricMatrix self, Real value ):
        """Set all the elements of a symmetricmatrix."""
        SymmetricMatrix_Set ( self.cObject, value )

    def Update ( SymmetricMatrix self, Real1DArray dx, Real1DArray dg, option = "BFGS", tolerance = None ):
        """Update the matrix."""
        cdef Real                         tol
        cdef SYMMETRICMATRIXUPDATING_OPTION formula
        formula = _UpdatingOptions.get ( option.upper ( ), SYMMETRICMATRIXUPDATING_BFGS )
        if tolerance is None:
            SymmetricMatrix_Update ( self.cObject, dx.cObject, dg.cObject, formula, NULL )
        else:
            tol = tolerance
            SymmetricMatrix_Update ( self.cObject, dx.cObject, dg.cObject, formula, &tol )

    def VectorMultiply ( SymmetricMatrix self, Real1DArray other, Real1DArray result ):
        """Multiply by a vector."""
        SymmetricMatrix_VectorMultiply ( self.cObject, other.cObject, result.cObject, NULL )

    @classmethod
    def WithExtent ( selfClass, extent ):
        """Constructor with extent."""
        return selfClass ( extent )

    # . Properties.
    property columns:
        def __get__ ( self ): return self.cObject.dimension

    property isSquare:
        def __get__ ( self ): return ( self.rows == self.columns )

    property rank:
        def __get__ ( self ): return 2

    property rows:
        def __get__ ( self ): return self.cObject.dimension

    property shape:
        def __get__ ( self ):
            return [ self.rows, self.columns ]

    property size:
        def __get__ ( self ):
            return SymmetricMatrix_Size ( self.cObject )
