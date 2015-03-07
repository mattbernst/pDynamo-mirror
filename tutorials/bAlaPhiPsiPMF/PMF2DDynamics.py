"""2-D PMFs - data generation."""

from Definitions import *

# . Dynamics parameters.
_LogFrequency     = 1000
_NSteps0          = 10000
_NSteps1          = 10000
_TimeStep         = 0.001

# . Window parameters.
_ForceConstant    = 0.1
_NumberOfWindows  = 36
_WindowIncrement  = 360.0 / float ( _NumberOfWindows )
seed0             = 215689

# . Get the system.
system              = Unpickle ( os.path.join ( outPath, "bAla.pkl" ) )
system.coordinates3 = XYZFile_ToCoordinates3 ( os.path.join ( outPath, "bAla_md.xyz" ) )
system.Summary ( )
system.Energy  ( )

# . Define constraints.
constraints = SoftConstraintContainer ( )
system.DefineSoftConstraints ( constraints )

# . Initial dihedral angles.
phi0 = system.coordinates3.Dihedral ( *phiAtomIndices )
psi0 = system.coordinates3.Dihedral ( *psiAtomIndices )

# . Generator options.
randomNumberGenerator  = RandomNumberGenerator.WithSeed ( seed0 )
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )
seed0 += 1

# . Loop over phi windows.
for i in range ( _NumberOfWindows ):

    # . Phi constraint.
    phi        = phi0 + float ( i ) * _WindowIncrement
    phiSCModel = SoftConstraintEnergyModelHarmonic ( phi, _ForceConstant )
    constraint = SoftConstraintDihedral ( * ( list ( phiAtomIndices ) + [ phiSCModel ] ) )
    constraints["Phi"] = constraint

    # . Loop over psi windows.
    for j in range ( _NumberOfWindows ):

        # . Psi constraint.
        psi        = psi0 + float ( j ) * _WindowIncrement
        psiSCModel = SoftConstraintEnergyModelHarmonic ( psi, _ForceConstant )
        constraint = SoftConstraintDihedral ( * ( list ( psiAtomIndices ) + [ psiSCModel ] ) )
        constraints["Psi"] = constraint

        # . Equilibration.
        LangevinDynamics_SystemGeometry ( system                                 ,
                                          collisionFrequency     =          25.0 ,
                                          logFrequency           = _LogFrequency ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  =      _NSteps0 ,
                                          temperature            =         300.0 ,
                                          timeStep               =     _TimeStep )

        # . Data-collection.
        trajectoryPath = os.path.join ( outPath, "bAla_phi_{:d}_psi_{:d}.trj".format ( i, j ) )
        trajectory     = SystemSoftConstraintTrajectory ( trajectoryPath, system, mode = "w" )
        LangevinDynamics_SystemGeometry ( system                                 ,
                                          collisionFrequency     =          25.0 ,
                                          logFrequency           = _LogFrequency ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  =      _NSteps1 ,
                                          temperature            =         300.0 ,
                                          timeStep               =     _TimeStep ,
                                          trajectories = [ ( trajectory, 1 ) ] )
