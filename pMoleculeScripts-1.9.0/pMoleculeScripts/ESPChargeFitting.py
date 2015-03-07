#-------------------------------------------------------------------------------
# . File      : ESPChargeFitting.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""ESP charge fitting."""

import math

from pCore     import Coordinates3, logFile, LogFileActive, Real1DArray, Real2DArray, UNITS_LENGTH_ANGSTROMS_TO_BOHRS, Vector3
from pMolecule import LebedevLaikovGrid_GetGridPoints, PeriodicTable

#===================================================================================================================================
# . Parameter definitions.
#===================================================================================================================================
# . Default values - RESP parameters.
_DEFAULT_ARESP = 0.0005 # . Weak constraints.
_DEFAULT_BRESP = 0.1

# . Default values - RESP iterator.
_DEFAULT_FTOLERANCE        = 1.0e-10
_DEFAULT_LOGFREQUENCY      = 1
_DEFAULT_MAXIMUMITERATIONS = 250
_DEFAULT_QTOLERANCE        = 1.0e-8

#===================================================================================================================================
# . Perform the charge fitting.
#===================================================================================================================================
def ESPChargeFitting ( system, aresp = _DEFAULT_ARESP, bresp = _DEFAULT_BRESP, log = logFile, QRESP = False ):
    """Perform an ESP charge fit."""

    # . Check for a system with a qcModel.
    if not ( hasattr ( system, "energyModel" ) and hasattr ( system.energyModel, "qcModel" ) ): raise ValueError ( "System does not have a QC energy model." )

    # . Initialization.
    natoms = len ( system.energyModel.qcAtoms )
    ndim   = natoms + 1
    qtotal = 0.0
    if hasattr ( system, "electronicState" ): qtotal = float ( system.electronicState.charge )

#    atomicNumbers = system.atoms.GetItemAttributes ( "atomicNumber" )
#    qtotal = float ( sum ( atomicNumbers ) )

    # . Get the grid points for the QC atoms only and convert to bohrs.
    gridPoints = GenerateVanDerWaalsSurface ( system, log = log, qcAtomsOnly = True, scalingfactors = [ 1.4, 1.6, 1.8, 2.0 ] )
    gridPoints.Scale ( UNITS_LENGTH_ANGSTROMS_TO_BOHRS )

    # . Allocate space - one larger than necessary for fInteraction.
    fInteraction = Real2DArray.WithExtents ( ndim, gridPoints.rows )
    fInteraction.Set ( 0.0 )

    # . Get the observed potentials at the grid points.
    phi = system.energyModel.qcModel.GridPointPotentials ( system.configuration, gridPoints )

    # . Get the interaction terms for each atom at the grid points.
    coordinates3 = system.energyModel.qcAtoms.GetCoordinates3 ( system.coordinates3, toBohrs = True )
    GetInteractionTerms ( coordinates3, gridPoints, fInteraction )

    # . Get the A matrix and the B vector.
    A = Real2DArray.WithExtents ( ndim, ndim ) ; A.Set ( 0.0 )
    B = Real1DArray.WithExtent  ( ndim )       ; B.Set ( 0.0 )
    A.MatrixMultiply ( fInteraction, fInteraction, yTranspose = True )
    fInteraction.VectorMultiply ( phi, B )

    # . Add the total charge constraint terms.
    for i in range ( natoms ):
        A[i,ndim-1] = 1.0
        A[ndim-1,i] = 1.0
    B[ndim-1] = qtotal

#    A.Print ( )
#    B.Print ( )
#    phi.Print ( )

    # . Get |phi^2|.
    phi2 = phi.Dot ( phi )

    # . Iterative solution.
    if QRESP:
        ( isConverged, sos, solution, condition, rank ) = RESPIterator ( natoms, A, B, fInteraction, phi, aresp, bresp, log = log )

    # . Non-iterative solution.
    else:

        # . Solve the equations.
        ( condition, rank ) = A.SolveLinearEquationsBySVD ( B )
        solution = B

        # . Determine the sum of squares.
        fInteraction.VectorMultiply ( solution, phi, alpha = -1.0, beta = 1.0, transpose = True )
        sos = phi.Dot ( phi )

    # . Determine the error measure.
    error = math.sqrt ( sos / phi2 )

    # . Do some printing.
    if LogFileActive ( log ):
        summary = log.GetSummary ( )
        summary.Start ( "ESP Charge Fitting Summary" )
        summary.Entry ( "Charges"   , "{:d}"  .format ( natoms    ) )
        summary.Entry ( "Rank"      , "{:d}"  .format ( rank      ) )
        summary.Entry ( "Error"     , "{:.5g}".format ( error     ) )
        summary.Entry ( "Condition" , "{:.5g}".format ( condition ) )
        if QRESP:
            if isConverged: summary.Entry ( "Converged" , "True"  )
            else:           summary.Entry ( "Converged" , "False" )
        summary.Stop ( )

    # . Return the charges.
    charges = Real1DArray.WithExtent ( natoms )
    for i in range ( natoms ):
        charges[i] = solution[i]
    return charges 

#===================================================================================================================================
# . Generate grid points on the van der Waals surface for a system.
# . Maybe could simplify this by scaling a single complete surface?
#===================================================================================================================================
def GenerateVanDerWaalsSurface ( system, gridangularmomentum = 21, log = logFile, qcAtomsOnly = False, scalingfactors = [ 1.0 ] ):
    """Generate a superposition of van der Waals surfaces represented by grid points."""

    # . QC atoms only (including link atoms).
    if qcAtomsOnly:
        qcAtoms       = system.energyModel.qcAtoms
        atomicNumbers = qcAtoms.GetAtomicNumbers ( )
        coordinates3  = qcAtoms.GetCoordinates3 ( system.coordinates3 )
    # . All atoms.
    else:
        atomicNumbers = system.atoms.GetItemAttributes ( "atomicNumber" )
        coordinates3  = system.coordinates3

    # . Find the number of atoms.
    natoms = len ( atomicNumbers )

    # . Set radii.
    radii = Real1DArray.WithExtent ( natoms )
    for ( i, n ) in enumerate ( atomicNumbers ):
        radii[i] = PeriodicTable.Element ( n ).vdwRadius

    # . Get the grid points for a single center.
    ( basicgrid, weights ) = LebedevLaikovGrid_GetGridPoints ( gridangularmomentum )

    # . Allocate space for the possible number of grid points.
    npossible  = natoms * basicgrid.rows * len ( scalingfactors )
    gridPoints = Coordinates3.WithExtent ( npossible )
    gridPoints.Set ( 0.0 )

    # . Initialization.
    nfound      = 0
    atomgrid    = Coordinates3.WithExtent ( basicgrid.rows )
    translation = Vector3.Null ( )
    atomgrid.Set ( 0.0 )

    # . Loop over scaling factors.
    for factor in scalingfactors:

        # . Loop over points.
        for i in range ( natoms ):

            # . Get the radius.
            iradius = factor * radii[i]

            # . Get the translation.
            translation[0] = coordinates3[i,0]
            translation[1] = coordinates3[i,1]
            translation[2] = coordinates3[i,2]

            # . Get the scaled grid centered at the point.
            basicgrid.CopyTo   ( atomgrid    )
            atomgrid.Scale     ( iradius     )
            atomgrid.Translate ( translation )

            # . Remove points that are within the current scaled radii of other points.
            for p in range ( atomgrid.rows ):
                QOK = True
                x   = atomgrid[p,0]
                y   = atomgrid[p,1]
                z   = atomgrid[p,2]
                for j in range ( natoms ):
                    if j != i:
                        dx       = coordinates3[j,0] - x
                        dy       = coordinates3[j,1] - y
                        dz       = coordinates3[j,2] - z
                        jradius2 = ( factor * radii[j] )**2
                        if ( dx**2 + dy**2 + dz**2 ) <= jradius2:
                            QOK = False
                            break
                if QOK:
                    gridPoints[nfound,0] = x
                    gridPoints[nfound,1] = y
                    gridPoints[nfound,2] = z
                    nfound += 1

    # . Reduce the array size if necessary.
    if nfound < npossible:
        newpoints = Coordinates3.WithExtent ( nfound )
        for p in range ( nfound ):
            newpoints[p,0] = gridPoints[p,0]
            newpoints[p,1] = gridPoints[p,1]
            newpoints[p,2] = gridPoints[p,2]
        gridPoints = newpoints

    # . Create a system.
#    from pBabel  import XYZFile_FromSystem
#    from pMolecule import System
#    ngrid = gridPoints.rows
#    junk = System ( ngrid * [ 1 ] )
#    junk.coordinates3 = gridPoints
#    XYZFile_FromSystem ( "junk.xyz", junk )

    # . Do some printing.
    if LogFileActive ( log ):
        summary = log.GetSummary ( )
        summary.Start ( "van der Waals Surface Generation Summary" )
        summary.Entry ( "Atoms"           , "{:d}".format ( natoms                 ) )
        summary.Entry ( "Surfaces"        , "{:d}".format ( len ( scalingfactors ) ) )
        summary.Entry ( "Found Points"    , "{:d}".format ( nfound                 ) )
        summary.Entry ( "Possible Points" , "{:d}".format ( npossible              ) )
        summary.Stop ( )

    # . Finish up.
    return gridPoints

#===================================================================================================================================
# . Get the interaction terms between the model charges and grid points.
#===================================================================================================================================
def GetInteractionTerms ( coordinates3, gridPoints, fInteraction ):
    """Get interaction terms."""
    # . To start with use simple 1/r terms. Should be generalized to allow non-delta function (e.g. Gaussian) densities for atoms.
    # . In which case would require radii or Gaussian widths as well.
    fInteraction.Set ( 0.0 )
    for i in range ( coordinates3.rows ):
        x = coordinates3[i,0]
        y = coordinates3[i,1]
        z = coordinates3[i,2]
        for g in range ( gridPoints.rows ):
            dx = gridPoints[g,0] - x
            dy = gridPoints[g,1] - y
            dz = gridPoints[g,2] - z
            r  = math.sqrt ( dx**2 + dy**2 + dz**2 )
            if r != 0.0:
                fInteraction[i,g] = 1.0 / r

#===================================================================================================================================
# . Get the RESP constraints.
#===================================================================================================================================
def GetRESPConstraints ( natoms, aresp, bresp, charges, A ):
    """Add in RESP constraints to the diagonal elements of the A matrix."""
    for i in range ( natoms ):
        q = charges[i]
        A[i,i] += aresp * q / math.sqrt ( q**2 + bresp**2 )

#===================================================================================================================================
# . RESP iterator function.
#===================================================================================================================================
def RESPIterator ( natoms, A, rhs, fInteraction, phi, aresp, bresp, ftolerance = _DEFAULT_FTOLERANCE, log = logFile, logFrequency = _DEFAULT_LOGFREQUENCY, maximumIterations = _DEFAULT_MAXIMUMITERATIONS, qtolerance = _DEFAULT_QTOLERANCE ):
    """Solve the RESP equations by simple iteration."""

    # . Allocate space.
    n       = len ( rhs )
    Atemp   = Real2DArray.WithExtents ( n, n )
    Btemp   = Real1DArray.WithExtent  ( n )
    phitemp = Real1DArray.WithExtent  ( len ( phi ) )
    q       = Real1DArray.WithExtent  ( n, initializer = 0.0 ) # . With initial values.
    qold    = Real1DArray.WithExtent  ( n )

    # . Check for printing.
    QPRINT = ( logFrequency > 0 ) and ( logFrequency < maximumIterations ) and LogFileActive ( log )
    if QPRINT:
        table = log.GetTable ( columns = [ 10, 20, 20, 20, 20, 10 ] )
        table.Start ( )
        table.Title ( "RESP Solver" )
        table.Heading ( "Iteration"   )
        table.Heading ( "Function"    )
        table.Heading ( "Change in F" )
        table.Heading ( "Change in Q" )
        table.Heading ( "Condition"   )
        table.Heading ( "Rank"        )

    # . Determine the sum of squares with initial zero charges.
    f = phi.Dot ( phi )

    # . Initialization.
    niterations = 0
    isConverged  = False

    # . Perform the iterations.
    while ( niterations < maximumIterations ) and ( not isConverged ):

        # . Save old data.
        fold = f
        q.CopyTo ( qold )

        # . Fill new A and RHS.
        A.CopyTo   ( Atemp )
        rhs.CopyTo ( Btemp )

        # . Add in the constraints.
        GetRESPConstraints ( natoms, aresp, bresp, q, Atemp )

        # . Solve the equations.
        ( condition, rank ) = Atemp.SolveLinearEquationsBySVD ( Btemp )
        Btemp.CopyTo ( q )

        # . Determine the sum of squares.
        phi.CopyTo ( phitemp )
        fInteraction.VectorMultiply ( q, phitemp, alpha = -1.0, beta = 1.0, transpose = True )
        f = phitemp.Dot ( phitemp )

        # . Check for convergence.
        qold.AddScaledArray ( -1.0, q )
        fdifference = math.fabs ( fold - f )
        qdifference = qold.AbsoluteMaximum ( )
        isConverged  = ( fdifference < ftolerance ) and ( qdifference < qtolerance )

        # . Printing.
        if QPRINT and ( niterations % logFrequency == 0 ):
            table.Entry ( "{:d}"  .format ( niterations ) )
            table.Entry ( "{:.6g}".format ( f           ) )
            table.Entry ( "{:.6g}".format ( fdifference ) )
            table.Entry ( "{:.6g}".format ( qdifference ) )
            table.Entry ( "{:.6g}".format ( condition   ) )
            table.Entry ( "{:d}"  .format ( rank        ) )

        # . End of loop.
        niterations += 1

    # . Stop printing.
    if QPRINT:
        table.Stop ( )
        if isConverged: log.Paragraph ( "RESP procedure converged." )
        else:           log.Paragraph ( "Warning: RESP procedure not converged." )

    # . Finish up.
    return ( isConverged, f, q, condition, rank )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
