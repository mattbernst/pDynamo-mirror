"""Linear algebra (matrix) benchmarks."""

from pCore import CPUTime, logFile, NormalDeviateGenerator, RandomNumberGenerator, Real1DArray, Real2DArray, SymmetricMatrix

# . Matrix sizes.
eExtents = ( 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000 )
mExtents = ( 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000 )

# . Functions.
def Eigenvalues ( extent, cpuTimer, ndg ):
    s = SymmetricMatrix ( extent         ) ; s.Set ( 0.0 )
    e = Real1DArray     ( extent         ) ; e.Set ( 0.0 )
    v = Real2DArray     ( extent, extent ) ; v.Set ( 0.0 )
    for i in range ( extent ):
        for j in range ( i+1 ): s[i,j] = ndg.NextDeviate ( )
        s[i,i] *= 2.0
    tStart = cpuTimer.Current ( )
    s. Diagonalize ( e, eigenVectors = v )
    return ( cpuTimer.Current ( ) - tStart )

def MatrixMultiply ( extent, cpuTimer, ndg ):
    a = Real2DArray ( extent, extent ) ; a.Set ( 0.0 )
    b = Real2DArray ( extent, extent ) ; b.Set ( 0.0 )
    c = Real2DArray ( extent, extent ) ; c.Set ( 0.0 )
    for i in range ( extent ):
        for j in range ( extent ):
            a[i,j] = ndg.NextDeviate ( )
            b[i,j] = ndg.NextDeviate ( )
    tStart = cpuTimer.Current ( )
    c.MatrixMultiply ( a, b )
    return ( cpuTimer.Current ( ) - tStart )

# . Initialization.
cpuTimer = CPUTime ( )
rng      = RandomNumberGenerator.WithSeed ( 314159 )
ndg      = NormalDeviateGenerator.WithRandomNumberGenerator ( rng, mu = 0.0, sigma = 5.0 )

# . Calculation.
for ( extents, function, tag ) in ( ( eExtents, Eigenvalues   , "Eigenvalue"      ) ,
                                    ( mExtents, MatrixMultiply, "Matrix Multiply" ) ):
    times = [ function ( extent, cpuTimer, ndg ) for extent in extents ]
    table = logFile.GetTable ( columns = [ 10, 20, 20 ] )
    table.Start   ( )
    table.Title   ( tag + " Timings" )
    table.Heading ( "Extent" )
    table.Heading ( "Times", columnSpan = 2 )
    for ( extent, time ) in zip ( extents, times ):
        table.Entry ( "{:d}"  .format ( extent ) )
        table.Entry ( "{:.3f}".format ( time   ) )
        table.Entry ( CPUTime.TimeToString ( time ) )
    table.Stop ( )
