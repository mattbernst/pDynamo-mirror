#-------------------------------------------------------------------------------
# . File      : NormalModes.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Helper functions for performing normal mode and associated thermodynamics calculations."""

import math

from pCore              import Clone,  CONSTANT_AVOGADRO_NUMBER, CONSTANT_BOLTZMANN, CONSTANT_MOLAR_GAS, CONSTANT_PLANCK, CONSTANT_SPEED_OF_LIGHT, \
                               Coordinates3, logFile, LogFileActive, Real1DArray, Real2DArray, UNITS_MASS_AMU_TO_KG, UNITS_PRESSURE_ATMOSPHERES_TO_PASCALS
from pMolecule          import SystemGeometryObjectiveFunction
from TrajectoryAnalysis import CovarianceMatrix

#===================================================================================================================================
# . Parameter definitions.
#===================================================================================================================================
# . Frequency output options.
_FREQUENCYFORMAT   = "{:10.3f}"
_FREQUENCYWIDTHS   = 12
_NFREQUENCYCOLUMNS =  8

# . Hessian modify options.
_MODIFYOPTIONS   = ( "PROJECT", "RAISE" )
_RAISEEIGENVALUE = 50000.0

# . Mode output options.
_NMODECOLUMNS       = 6
_MODEFORMAT         = "{:11.5f}"
_MODEWIDTHCOMPONENT =  3
_MODEWIDTHELEMENT   = 12
_MODEWIDTHINDEX     =  6
_MODEWIDTHNAME      =  4

# . Mode trajectory options.
_LOWFREQUENCY       = 0.1

# . The conversion factor from cm^-1 to ps^-1.
_TO_HZ = CONSTANT_SPEED_OF_LIGHT / 1.0e+10

# . The conversion factor from K to kJ mol^-1.
_TO_KJMOL = CONSTANT_AVOGADRO_NUMBER * CONSTANT_BOLTZMANN * 1.0e-03

# . The conversion factor from internal units to cm^-1.
_TO_WAVENUMBERS = 1.0e+11 / ( 2.0 * math.pi * CONSTANT_SPEED_OF_LIGHT )

# . The conversion factor from wavenumbers to Joules.
_WAVENUMBERS_TO_JOULES = 1.0e+12 * CONSTANT_PLANCK * _TO_HZ

# . Thermodynamic options.
_EXPONENTIALUNDERFLOW    = 75.0
_ROTATIONLINEARTOLERANCE = 1.0e-3
_THERMODYNAMICPROPERTIES = ( "Constant Pressure Heat Capacity", \
                             "Constant Volume Heat Capacity",   \
                             "Enthalpy",                        \
                             "Entropy",                         \
                             "Gibbs Free Energy",               \
                             "Helmholtz Free Energy",           \
                             "Internal Energy",                 \
                             "Log Partition Function"           )

# . Tolerances.
_QHTolerance = 1.0e-05

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NormalModeState ( object ):
    """A class to contain the results of a normal mode calculation."""

    # . Attributes.
    attributes = { "covariance"  : None, \
                   "dimension"   :    0, \
                   "freeAtoms"   : None, \
                   "frequencies" : None, \
                   "hessian"     : None, \
                   "modes"       : None, \
                   "nrtmodes"    :    0, \
                   "weights"     : None  }

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        # . Set defaults for all attributes.
        for ( key, value ) in self.__class__.attributes.iteritems ( ): setattr ( self, key, value )
        # . Process keywords.
        for ( key, value ) in          keywordArguments.iteritems ( ): setattr ( self, key, value )

#===================================================================================================================================
# . Normal modes.
#===================================================================================================================================
def NormalModes_SystemGeometry ( system, hessian = None, log = logFile, modify = None, title = "Harmonic Frequencies (cm^(-1))" ):
    """Determine the normal modes for a system."""

    # . Get some options.
    if modify is None: modopt = None
    else:              modopt = modify.upper ( )

    # . Get the Hessian with mass-weighting.
    of = SystemGeometryObjectiveFunction.FromSystem ( system )
    of.DefineWeights ( )
    n  = of.NumberOfVariables ( )
    if hessian is None:
        x  = Real1DArray.WithExtent ( n )
        x.Set ( 0.0 )
        of.VariablesGet ( x )
        hessian = of.NumericalHessian ( x )

    # . Get the mass-weighted rotation-translation vectors and count their number.
    of.RemoveRotationTranslation ( )
    if of.linearScalars is None: nrtmodes = 0
    else:                        nrtmodes = len ( of.linearScalars )

    # . Modify the Hessian.
    if modopt in _MODIFYOPTIONS:
        if   modopt == "PROJECT" : hessian.ProjectOutVectors ( of.linearVectors                   )
        elif modopt == "RAISE"   : hessian.Raise             ( of.linearVectors, _RAISEEIGENVALUE )

    # . Diagonalization.
    # . Maybe should save hessian here as it is destroyed by the diagonalization.
    eigenvalues  = Real1DArray.WithExtent  ( n )    ; eigenvalues.Set  ( 0.0 )
    eigenvectors = Real2DArray.WithExtents ( n, n ) ; eigenvectors.Set ( 0.0 )
    hessian.Diagonalize ( eigenvalues, eigenvectors )

    # . Convert eigenvalues to frequencies.
    for ( i, e ) in enumerate ( eigenvalues ):
        f = math.sqrt ( math.fabs ( e ) ) * _TO_WAVENUMBERS
        if e < 0.0: f *= -1.0
        eigenvalues[i] = f

    # . Un-mass-weight the modes.
    for r in range ( n ):
        w = 1.0 / of.variableWeights[r]
        for c in range ( n ): eigenvectors[r,c] *= w

    # . Do some printing.
    if LogFileActive ( log ):
        table = log.GetTable ( columns = _NFREQUENCYCOLUMNS * [ _FREQUENCYWIDTHS ] )
        table.Start ( )
        table.Title ( title )
        for f in eigenvalues: table.Entry ( _FREQUENCYFORMAT.format ( f ) )
        table.Stop ( )

    # . Save all data.
    state = NormalModeState ( dimension = n, freeAtoms = of.freeAtoms, frequencies = eigenvalues, modes = eigenvectors, nrtmodes = nrtmodes, weights = of.variableWeights )
    system.configuration.SetTemporaryAttribute ( "nmState", state )

    # . Finish up.
    return state

#===================================================================================================================================
# . Normal mode printing.
#===================================================================================================================================
def NormalModesPrint_SystemGeometry ( system, log = logFile, modes = None, selection = None, state = None, title = "Normal Mode Eigenvectors" ):
    """Print the normal modes."""
    # . Check for printing.
    if LogFileActive ( log ):

        # . Get the state.
        if state is None: state = system.configuration.nmState

        # . Get state-related information.
        n            = state.dimension
        eigenvectors = state.modes
        frequencies  = state.frequencies

        # . Get the atom selection as the intersection of the input selection and the free atoms.
        # . Input selection.
        if selection == None: atoms = range ( len ( system.atoms ) )
        else:                 atoms = selection

        # . Free atoms.
        if state.freeAtoms == None: freeAtoms = range ( len ( system.atoms ) )
        else:                       freeAtoms = state.freeAtoms

        # . The intersection.
        atomselection = list ( set ( atoms ).intersection ( set ( freeAtoms ) ) )
        atomselection.sort ( )

        # . Get the mode selection.
        if modes == None: modeselection = range ( n )
        else:             modeselection = modes
        nmodes = len ( modeselection )

        # . There are atoms and modes.
        if ( len ( atomselection ) > 0 ) and ( nmodes > 0 ):

            # . Get the free atom positions.
            freeindices = {}
            for ( i, f ) in enumerate ( freeAtoms ): freeindices[f] = i

            # . Get atom columns and names.
            atomcolumns = []
            atomnames   = []
            for atom in atomselection:
                atomcolumns.append ( 3 * freeindices[atom] )
                atomnames.append ( system.atoms[atom].symbol )

            # . Do the printing.
            for start in range ( 0, nmodes, _NMODECOLUMNS ):
                stop    = min ( start + _NMODECOLUMNS, nmodes )
                columns = [ _MODEWIDTHINDEX, _MODEWIDTHNAME, _MODEWIDTHCOMPONENT ]
                for row in range ( start, stop ): columns.append ( _MODEWIDTHELEMENT )
                table = log.GetTable ( columns = columns )
                table.Start ( )
                table.Title ( title )
                table.Entry ( "Mode",          alignment = "l", columnSpan = 3 )
                for i in range ( start, stop ): table.Entry ( "{:d}".format ( modeselection[i] ) )
                table.Entry ( "Freq. (cm^-1)", alignment = "l", columnSpan = 3 )
                for i in range ( start, stop ): table.Entry ( _MODEFORMAT.format ( frequencies[modeselection[i]]          ) )
                table.Entry ( "Freq. (ps^-1)", alignment = "l", columnSpan = 3 )
                for i in range ( start, stop ): table.Entry ( _MODEFORMAT.format ( frequencies[modeselection[i]] * _TO_HZ ) )
                for ( atom, column, name ) in zip ( atomselection, atomcolumns, atomnames ):
                    table.Entry ( "{:d}".format ( atom ) )
                    table.Entry ( name   )
                    table.Entry ( "x", alignment = "c" )
                    for i in range ( start, stop ): table.Entry ( _MODEFORMAT.format ( eigenvectors[modeselection[i],column  ] ) )
                    table.Entry ( None )
                    table.Entry ( None )
                    table.Entry ( "y", alignment = "c" )
                    for i in range ( start, stop ): table.Entry ( _MODEFORMAT.format ( eigenvectors[modeselection[i],column+1] ) )
                    table.Entry ( None )
                    table.Entry ( None )
                    table.Entry ( "z", alignment = "c" )
                    for i in range ( start, stop ): table.Entry ( _MODEFORMAT.format ( eigenvectors[modeselection[i],column+2] ) )
                table.Stop ( )

#===================================================================================================================================
# . Generate a normal mode trajectory.
#===================================================================================================================================
def NormalModesTrajectory_SystemGeometry ( system, trajectory, cycles = 10, frames = 21, mode = 0, state = None, temperature = 300.0 ):
    """Generate a normal mode trajectory."""

    # . Get the state.
    if state is None: state = system.configuration.nmState

    # . Get state-related information.
    if state.freeAtoms == None: freeAtoms = range ( len ( system.atoms ) )
    else:                       freeAtoms = state.freeAtoms
    frequencies = state.frequencies
    modes       = state.modes

    # . Get the mode frequency.
    omega = math.fabs ( frequencies[mode] )

    # . Calculate the number of frames.
    total = cycles * frames

    # . Check for a calculation.
    if ( mode >= 0 ) and ( mode < state.dimension ) and ( omega > _LOWFREQUENCY ) and ( total > 0 ):

        # . Calculate the amplitude (in Angstroms).
        amplitude = math.sqrt ( 2.0 * 1.0e-3 * CONSTANT_AVOGADRO_NUMBER * CONSTANT_BOLTZMANN * temperature ) * ( _TO_WAVENUMBERS / omega )

        # . Allocate space for the coordinates and mode.
        coordinates3 = Clone ( system.coordinates3 )
        displacement = Coordinates3.WithExtent ( len ( system.atoms ) )
        displacement.Set ( 0.0 )

        # . Get the displacement.
        for f in freeAtoms:
            displacement[f,0] = modes[mode,3*f  ]
            displacement[f,1] = modes[mode,3*f+1]
            displacement[f,2] = modes[mode,3*f+2]

        # . Loop over the cycles and frames.
        # . Calculate the displacement prefactor using sine instead of cosine.
        trajectory.WriteHeader ( temperature = temperature )
        for c in range ( cycles ):
            for f in range ( frames ):
                factor = amplitude * math.sin ( 2.0 * math.pi * float ( f ) / float ( frames ) )
                coordinates3.CopyTo ( system.coordinates3  )
                system.coordinates3.AddScaledMatrix ( factor, displacement )
                trajectory.WriteOwnerData ( )

        # . Finish up.
        trajectory.WriteFooter ( )
        trajectory.Close       ( )

#===================================================================================================================================
# . Quasi-harmonic modes.
#===================================================================================================================================
# . Here H = T * inverse ( C ) so equations are T * inverse ( C ) * A = w^2 * M * A.
# . Rearranging gives: C * M * A = ( T / w^2 ) * A or C' * A' = e * A' where
# . C' = (M^1/2) * C * (M^1/2), A' = (M^1/2) * A and e = ( T / w^2 ). Therefore,
# . w = Sqrt ( T / e ) and A = A'/(M^1/2).
def QuasiHarmonic_SystemGeometry ( system, log = logFile, modify = None, temperature = 300.0, title = "Quasi-Harmonic Frequencies (cm^(-1))", trajectories = None ):
    """Determine the quasi-harmonic modes for a system."""

    # . Initialization.
    state = None

    # . Determine if any atoms are fixed.
    hc = system.hardConstraints
    if ( hc is not None ) and ( hc.fixedAtoms is not None ) and ( len ( hc.fixedAtoms ) > 0 ): fixedAtoms = hc.fixedAtoms
    else: fixedAtoms = None

    # . Get the covariance matrix.
    covariance = CovarianceMatrix ( trajectories, selection = fixedAtoms )

    # . Proceed with the analysis.
    if covariance is not None:

        # . Get some options.
        if modify is None: modopt = None
        else:              modopt = modify.upper ( )

        # . Mass-weight the covariance matrix.
        # . Weights are square roots of masses.
        of = SystemGeometryObjectiveFunction.FromSystem ( system )
        of.DefineWeights ( )
        n  = of.NumberOfVariables ( )
        for i in range ( n ):
            wI = of.variableWeights[i]
            for j in range ( i + 1 ):
                wJ = of.variableWeights[j]
                covariance[i,j] *= ( wI * wJ )

        # . Get the mass-weighted rotation-translation vectors and count their number.
        of.RemoveRotationTranslation ( )
        if of.linearScalars is None: nrtmodes = 0
        else:                        nrtmodes = len ( of.linearScalars )

        # . Modify the Hessian.
        if modopt in _MODIFYOPTIONS:
            if   modopt == "PROJECT" : covariance.ProjectOutVectors ( of.linearVectors      )
            elif modopt == "RAISE"   : covariance.Raise             ( of.linearVectors, 0.0 )

        # . Diagonalization.
        eigenValues  = Real1DArray.WithExtent  ( n )    ; eigenValues.Set  ( 0.0 )
        eigenVectors = Real2DArray.WithExtents ( n, n ) ; eigenVectors.Set ( 0.0 )
        covariance.Diagonalize ( eigenValues, eigenVectors )

        # . Convert eigenvalues to frequencies.
        conversionFactor = math.sqrt ( _TO_KJMOL * temperature ) * _TO_WAVENUMBERS
        numberZero       = 0
        for ( i, e ) in enumerate ( eigenValues ):
            eAbs = math.fabs ( e )
            if eAbs <= _QHTolerance:
                f = 0.0
                numberZero += 1
            else:
                f = math.sqrt ( 1.0 / eAbs ) * conversionFactor
                if e < 0.0: f *= -1.0
            eigenValues[i] = f

        # . Un-mass-weight the modes.
        for r in range ( n ):
            w = 1.0 / of.variableWeights[r]
            for c in range ( n ): eigenVectors[r,c] *= w

        # . Reverse in place (excluding zero modes).
        temporary = Real1DArray.WithExtent ( n )
        for i in range ( ( n - numberZero ) // 2 ):
            # . Indices.
            lower = i + numberZero
            upper = n - i - 1
            # . Eigenvalues.
            e = eigenValues[upper]
            eigenValues[upper] = eigenValues[lower]
            eigenValues[lower] = e
            # . Eigenvectors.
            for j in range ( n ): temporary[j]          = eigenVectors[j,upper]
            for j in range ( n ): eigenVectors[j,upper] = eigenVectors[j,lower]
            for j in range ( n ): eigenVectors[j,lower] = temporary[j]

        # . Do some printing.
        if LogFileActive ( log ):
            table = log.GetTable ( columns = _NFREQUENCYCOLUMNS * [ _FREQUENCYWIDTHS ] )
            table.Start ( )
            table.Title ( title )
            for f in eigenValues: table.Entry ( _FREQUENCYFORMAT.format ( f ) )
            table.Stop ( )

        # . Save all data.
        state = NormalModeState ( dimension = n, freeAtoms = of.freeAtoms, frequencies = eigenValues, modes = eigenVectors, nrtmodes = nrtmodes, weights = of.variableWeights )
        system.configuration.SetTemporaryAttribute ( "qhState", state )

    # . Finish up.
    return state

#===================================================================================================================================
# . Thermodynamical quantities within the RRHO approximation.
#===================================================================================================================================
def ThermodynamicsRRHO_SystemGeometry ( system, pressure = 1.0, state = None, symmetryNumber = 1, temperature = 300.0 ):
    """Determine thermodynamical quantities within the RRHO approximation."""

    # . Allowing multiple P and T values would be more efficient.

    # . Get the state.
    if state is None: state = system.configuration.nmState

    # . Calculate pressure, R, RT and volume.
    p  = UNITS_PRESSURE_ATMOSPHERES_TO_PASCALS * pressure
    R  = CONSTANT_MOLAR_GAS / 1.0e+3
    RT = R * temperature
    v  = CONSTANT_BOLTZMANN * temperature / p

    # . Get some counters and the masses.
    natoms = len ( system.atoms )
    masses = system.atoms.GetItemAttributes ( "mass" )
    if state.freeAtoms == None:
        nfree = natoms
    else:
        nfree = len ( state.freeAtoms )
        for i in system.hardConstraints.fixedatoms: masses[i] = 0.0
    totmas = UNITS_MASS_AMU_TO_KG * sum ( masses )

    # . Electronic contributions.
    # . Assume there is one state with a multiplicity of one.
    electronic = {}

    # . Rotational contributions.
    # . Initialization.
    rotation = {}
    # . Polyatomic system.
    if nfree > 1:

        # . Get the moments of inertia factor (in amu angstroms**2).
        coordinates3  = Clone ( system.coordinates3 )
        coordinates3.ToPrincipalAxes ( weights = masses )
        inertiamatrix = coordinates3.InertiaMatrix ( weights = masses )
        mproduct      = 1.0
        nzeromoments  = 0
        for i in range ( 3 ):
            m = inertiamatrix[i,i]
            if math.fabs ( m ) <= _ROTATIONLINEARTOLERANCE: nzeromoments += 1
            else:                                           mproduct     *= m

        # . Calculate some factors.
        factor = ( 8.0 * math.pi * math.pi * CONSTANT_BOLTZMANN * temperature * UNITS_MASS_AMU_TO_KG ) / ( CONSTANT_PLANCK**2 * 1.0e+20 )

        # . Linear molecule.
        if nzeromoments >= 2:
            z = factor * mproduct / float ( symmetryNumber )
            rotation["Log Partition Function"       ] = math.log ( z )
            rotation["Constant Volume Heat Capacity"] = R
            rotation["Entropy"                      ] = R * ( rotation["Log Partition Function"] + 1.0 )
            rotation["Internal Energy"              ] = RT
        # . Non-linear molecule.
        else:
            z = math.sqrt ( math.pi * factor**3 * mproduct ) / float ( symmetryNumber )
            rotation["Log Partition Function"]        = math.log ( z )
            rotation["Constant Volume Heat Capacity"] = 1.5 * R
            rotation["Entropy"                      ] = R * ( rotation["Log Partition Function"] + 1.5 )
            rotation["Internal Energy"              ] = 1.5 * RT

        # . Calculate the common terms.
        rotation["Helmholtz Free Energy"          ] = - RT * rotation["Log Partition Function"]
        rotation["Constant Pressure Heat Capacity"] = rotation["Constant Volume Heat Capacity"]
        rotation["Helmholtz Free Energy"          ] = rotation["Internal Energy"]
        rotation["Gibbs Free Energy"              ] = rotation["Helmholtz Free Energy"] - temperature * rotation["Entropy"]

    # . Translational contributions.
    z = ( math.sqrt ( 2.0 * math.pi * totmas * CONSTANT_BOLTZMANN * temperature ) / CONSTANT_PLANCK )**3 * v
    translation                                    = {}
    translation["Log Partition Function"         ] = math.log ( z )
    translation["Helmholtz Free Energy"          ] = - RT * translation["Log Partition Function"]
    translation["Constant Pressure Heat Capacity"] = 2.5 * R
    translation["Constant Volume Heat Capacity"  ] = 1.5 * R
    translation["Helmholtz Free Energy"          ] = 2.5 * RT
    translation["Entropy"                        ] = R * ( translation["Log Partition Function"] + 1.5 )
    translation["Internal Energy"                ] = 1.5 * RT
    translation["Gibbs Free Energy"              ] = translation["Helmholtz Free Energy"] - temperature * translation["Entropy"]

    # . Vibrational contributions.
    # . Initialization.
    vibration = {}
    # . Polyatomic system.
    if natoms > 1:

        # . Get some state-related information.
        frequencies = state.frequencies
        nrtmodes    = state.nrtmodes

        # . Remove imaginary and rotation-translation modes.
        start = len ( frequencies )
        zero  = []
        for ( i, f ) in enumerate ( frequencies ):
            zero.append ( f )
            if len ( zero ) > nrtmodes:
                old = math.fabs ( zero.pop ( 0 ) )
                if f > old:
                    start = i
                    break

        # . Loop over the frequencies.
        cv  = 0.0
        lnz = 0.0
        u   = 0.0
        zpe = 0.0
        for i in range ( start, len ( frequencies ) ):
            omega = _WAVENUMBERS_TO_JOULES * frequencies[i]
            hvkt  = omega / ( CONSTANT_BOLTZMANN * temperature )
            lnz  -= 0.5 * hvkt
            zpe  += 0.5 * 1.0e-3 * CONSTANT_AVOGADRO_NUMBER * omega
            if hvkt <= _EXPONENTIALUNDERFLOW:
                expm  = math.exp ( - hvkt )
                expp  = 1.0 - expm
                lnz  -= math.log ( expp )
                expp  = expm / expp
                cv   += R * hvkt * hvkt * expp * ( 1.0 + expp )
                u    += RT * hvkt * expp

        # . Assign the terms.
        vibration["Log Partition Function"         ] = lnz
        vibration["Helmholtz Free Energy"          ] = - RT * lnz
        vibration["Constant Pressure Heat Capacity"] = cv
        vibration["Constant Volume Heat Capacity"  ] = cv
        vibration["Internal Energy"                ] = u + zpe
        vibration["Entropy"                        ] = R * lnz + vibration["Internal Energy"] / temperature
        vibration["Gibbs Free Energy"              ] = vibration["Internal Energy"] - temperature * vibration["Entropy"]
        vibration["Helmholtz Free Energy"          ] = vibration["Internal Energy"]

    # . Get the totals.
    tdics = {}
    for key in _THERMODYNAMICPROPERTIES:
        tdics[key] = electronic.get ( key, 0.0 ) + rotation.get ( key, 0.0 ) + translation.get ( key, 0.0 ) + vibration.get ( key, 0.0 )

    # . Add in some extra terms pertaining to the ideal gas.
    factor = math.log ( CONSTANT_AVOGADRO_NUMBER ) - 1.0
    tdics["Gibbs Free Energy"]     -= RT * factor
    tdics["Entropy"]               -= R  * factor
    tdics["Helmholtz Free Energy"] += RT * factor

    # . Return.
    return tdics

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
