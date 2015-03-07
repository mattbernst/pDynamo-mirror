#-------------------------------------------------------------------------------
# . File      : JaguarScripts.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Scripts for manipulating files from the Jaguar program."""

import math

from pBabel    import JaguarInputFileReader, JaguarOutputFileReader
from pCore     import Clone, logFile, LogFileActive, Real1DArray, Real2DArray, SymmetricMatrix
from pMolecule import PeriodicTable

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The tolerance for bond order values.
_DEFAULTBONDORDERTOLERANCE = 0.1

# . The tolerance for charge deviations.
_DEFAULTCHARGETOLERANCE = 1.0e-3

# . The occupancy tolerance for including orbitals in the density matrix calculation.
_OCCUPANCYTOLERANCE = 1.0e-10

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def JaguarBondOrders ( infile, outfile, bondOrdertolerance = _DEFAULTBONDORDERTOLERANCE, chargetolerance = _DEFAULTCHARGETOLERANCE, log = logFile, QCROSS = False ):
    """Calculate bond orders given Jaguar input and output files."""

    # . Parse the input file.
    infile = JaguarInputFileReader ( infile )
    infile.Parse ( )
    infile.Summary ( log = logFile )

    # . Parse the output file.
    outfile = JaguarOutputFileReader ( outfile )
    outfile.Parse   ( )
    outfile.Summary ( log = logFile )

    # . Get data from the files.
    # . Input.
    atomicNumbersi = getattr ( infile,  "atomicNumbers", None )
    coordinates3   = getattr ( infile,  "coordinates3",  None )
    orbitalsets    = getattr ( infile,  "orbitalsets",   None )

    # . Output.
    atomicNumberso = getattr ( outfile, "atomicNumbers", None )
    nfunctions     = getattr ( outfile, "nfunctions",    None )
    overlap        = getattr ( outfile, "overlap",       None )

    # . Check for coordinates.
    if coordinates3 is None: raise ValueError ( "The coordinate data is missing from the input file." )

    # . Check the orbital data.
    QOK = ( orbitalsets is not None ) and ( len ( orbitalsets ) > 0 )
    if QOK:
        if   ( len ( orbitalsets ) == 1 ) and ( ""      in orbitalsets )                               : spinDensitiesRESTRICTED = True
        elif ( len ( orbitalsets ) == 2 ) and ( "alpha" in orbitalsets ) and ( "beta" in orbitalsets ) : spinDensitiesRESTRICTED = False
        else: QOK = False
    if not QOK: raise ValueError ( "Invalid orbital data on input file." )
    if spinDensitiesRESTRICTED: nbasisi = infile.orbitalsets[""]     [1]
    else:               nbasisi = infile.orbitalsets["alpha"][1]

    # . Check the overlap.
    if ( overlap is None ): raise ValueError ( "The overlap matrix is missing from the output file." )
    nbasiso = overlap.Dimension ( )

    # . Check the array giving the number of functions per atom.
#    print nfunctions, len ( nfunctions ), len ( atomicNumberso ), sum ( nfunctions ), nbasiso
    if ( nfunctions is None ) or ( len ( nfunctions ) != len ( atomicNumberso ) or ( sum ( nfunctions ) != nbasiso ) ): raise ValueError ( "Basis function data on the output file is missing or invalid." )

    # . Create the function index array.
    findices = []
    first = 0
    last  = 0
    for f in nfunctions:
        last  += f
        findices.append ( ( first, last ) )
        first += f

    # . Check for compatibility between the data.
    QOK = ( atomicNumbersi is not None ) and ( atomicNumberso is not None ) and ( len ( atomicNumbersi ) == len ( atomicNumberso ) ) and ( nbasisi == nbasiso )
    if QOK:
        for ( i, j ) in zip ( atomicNumbersi, atomicNumberso ):
            if i != j:
                QOK = False
                break
    if not QOK: raise ValueError ( "The systems on the input and output files are incompatible." )

    # . Set the keys for the calculation.
    if spinDensitiesRESTRICTED: keys = [ "" ]
    else:               keys = [ "alpha", "beta" ]

    # . Get the densities multiplied by the overlap.
    ps = {}
    for key in keys:
        p      = _MakeDensity ( orbitalsets[key] )
#        p.Print ( title = "**1**" )
        result = Real2DArray.WithExtents ( nbasisi, nbasisi )
        result.Set ( 999.0 )
        p.PostMultiply ( overlap, result )
#        overlap.Print ( title = "**2**" )
#        result.Print ( title = "**3**" )
        ps[key] = result
#    f = STOP
    # . Scale ps correctly for the spin-restricted case.
    if spinDensitiesRESTRICTED: ps[""].Scale ( 2.0 )

    # . If cross terms are not required condense the ps matrices.
    if ( not QCROSS ) and ( not spinDensitiesRESTRICTED ):
        tps  = ps.pop ( keys[0] )
        for key in keys[1:]: tps.AddScaledMatrix ( 1.0, ps[key] )
        ps   = { "": tps }
        keys = [ "" ]

    # . Get the bond-orders.
    bondOrders = {}
    for key1 in keys:
        for key2 in keys:
            bondOrders[ ( key1, key2 ) ] = _MakeBondOrders ( ps[key1], ps[key2], findices )

    # . Make the total bond-order if necessary.
    bokeys = bondOrders.keys ( )
    if len ( bokeys ) > 1:
        bokeys.sort ( )
        tbo    = Clone ( bondOrders[bokeys[0]] )
        for key in bokeys[1:]: tbo.AddScaledMatrix ( 1.0, bondOrders[key] )
        tkey   = ( "", "" )
        bokeys.append ( tkey )
        bondOrders[tkey] = tbo

    # . Compute the electronic contribution to the Mulliken charges.
    qmulliken = Real1DArray.WithExtent ( len ( atomicNumbersi ) )
    qmulliken.Set ( 0.0 )
    for key in keys: _MakeElectronicMullikenCharges ( ps[key], findices, qmulliken )

    # . Determine atom valencies.
    free      = Real1DArray.WithExtent ( len ( atomicNumbersi ) ) ; free.Set      ( 0.0 )
    valencies = Real1DArray.WithExtent ( len ( atomicNumbersi ) ) ; valencies.Set ( 0.0 )
    tbo       = bondOrders[ ( "", "" ) ]
    for i in range ( len ( atomicNumbersi ) ):
        valencies[i] = ( - 2.0 * qmulliken[i] ) - tbo[i,i]
        totalbo = 0.0
        for j in range ( len ( atomicNumbersi ) ): totalbo += tbo[i,j]
        free[i] = ( - 2.0 * qmulliken[i] ) - totalbo

    # . Add in the core contributions to the Mulliken charges.
    for ( i, q ) in enumerate ( atomicNumbersi ): qmulliken[i] += float ( q )
    if outfile.QECP:
        for ( i, q ) in enumerate ( outfile.ecpelectrons ): qmulliken[i] -= float ( q )

    # . Output the results.
    if LogFileActive ( log ):

        # . Get the spin label.
        if spinDensitiesRESTRICTED: spinlabel = "Spin Restricted"
        else:               spinlabel = "Spin Unrestricted"

        # . Create the atom names.
        atomnames = []
        for ( i, n ) in enumerate ( atomicNumbersi ):
            atomnames.append ( PeriodicTable.Symbol ( n, index = i + 1 ) )

        # . Atom data.
        columns = [ 10, 20, 20, 20 ]
        if not spinDensitiesRESTRICTED: columns.append ( 20 )
        table   = log.GetTable ( columns = columns )
        table.Start ( )
        table.Title ( "Atom Data (" + spinlabel + ")" )
        table.Heading ( "Atom"            )
        table.Heading ( "Charge"          )
        table.Heading ( "Self Bond Order" )
        table.Heading ( "Valence"         )
        if not spinDensitiesRESTRICTED: table.Heading ( "Free Valence" )
        for ( i, ni ) in enumerate ( atomnames ):
            table.Entry ( ni )
            table.Entry ( "{:.3f}".format ( qmulliken[i] ) )
            table.Entry ( "{:.3f}".format ( tbo[i,i]     ) )
            table.Entry ( "{:.3f}".format ( valencies[i] ) )
            if not spinDensitiesRESTRICTED: table.Entry ( "{:.3f}".format ( free[i] ) )
        table.Stop ( )

        # . Bond orders.
        for key in bokeys:
            orders = bondOrders[key]
            table  = log.GetTable ( columns = [ 10, 10, 20, 20 ] )
            table.Start ( )
            if key == ( "", "" ): table.Title ( "Total Bond Orders" )
            else:                 table.Title ( key[0].title ( ) + "/" + key[1].title ( ) + " Bond Orders" )
            table.Heading ( "Atom 1"   )
            table.Heading ( "Atom 2"   )
            table.Heading ( "Order"    )
            table.Heading ( "Distance" )
            for ( i, ni ) in enumerate ( atomnames ):
                for ( j, nj ) in enumerate ( atomnames[0:i] ):
                    b = orders[i,j]
                    if math.fabs ( b ) > bondOrdertolerance:
                        table.Entry ( ni )
                        table.Entry ( nj )
                        table.Entry ( "{:.3f}".format ( b ) )
                        table.Entry ( "{:.3f}".format ( coordinates3.Distance ( i, j ) ) )
            table.Stop ( )

        # . Checks on the calculation.
        # . Free valence.
        if spinDensitiesRESTRICTED:
            deviation = free.AbsoluteMaximum ( )
            if deviation > chargetolerance: log.Paragraph ( "Warning: the largest deviation between the free valence values and zero is {:.3f}.".format ( deviation ) )

        # . Total charge.
        deviation = math.fabs ( qmulliken.Sum ( ) - float ( infile.charge ) )
        if deviation > chargetolerance: log.Paragraph ( "Warning: the total charge deviates from the formal charge by {:.3f}.".format ( deviation ) )

        # . Check for a valid reference set of Mulliken charges.
        qreference = getattr ( outfile, "qmulliken", None )
        if ( qreference is not None ) and ( len ( qreference ) == len ( atomicNumbersi ) ):
            qmulliken.AddScaledArray ( -1.0, qreference )
            deviation = qmulliken.AbsoluteMaximum ( )
            if deviation > chargetolerance: log.Paragraph ( "Warning: the largest deviation between calculated and read Mulliken charges is {:.3f}.".format ( deviation ) )

    # . Finish up.
    return bondOrders

#===================================================================================================================================
# . Private functions.
# . These functions could be made more efficient.
#===================================================================================================================================
def _MakeBondOrders ( ps1, ps2, findices ):
    """Make bond orders."""
    natoms = len ( findices )
    orders = SymmetricMatrix.WithExtent ( natoms )
    orders.Set ( 0.0 )
    for ( i, ( ifirst, ilast ) ) in enumerate ( findices ):
        for ( j, ( jfirst, jlast ) ) in enumerate ( findices[0:i+1] ):
            b = 0.0
            for m in range ( ifirst, ilast ):
                for n in range ( jfirst, jlast ):
                    b += ps1[m,n] * ps2[n,m]
            orders[i,j] = b
    return orders

def _MakeDensity ( orbitalset ):
    """Make a density."""
    ( norbitals, nbasis, energies, occupancies, vectors ) = orbitalset
    p = SymmetricMatrix.WithExtent ( nbasis )
    p.Set ( 0.0 )
    for i in range ( norbitals ):
        o = occupancies[i]
        if math.fabs ( o ) > _OCCUPANCYTOLERANCE:
            for m in range ( nbasis ):
                for n in range ( m + 1 ):
                    p[m,n] += o * vectors[m,i] * vectors[n,i]
    return p

def _MakeElectronicMullikenCharges ( ps, findices, q ):
    """Determine an electronic contribution to the Mulliken charges."""
    for ( i, ( first, last ) ) in enumerate ( findices ):
        for j in range ( first, last ): q[i] -= ps[j,j]

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
