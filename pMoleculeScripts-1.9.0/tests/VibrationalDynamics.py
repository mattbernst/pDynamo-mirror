"""Test for vibrational dynamics."""

import glob, math, os

from pBabel           import AmberTrajectoryFileReader, AmberTrajectoryFileWriter, MOLFile_ToSystem
from pCore            import Clone, NormalDeviateGenerator, RandomNumberGenerator, Selection, TestCase
from pMolecule        import DIISSCFConverger, MMModelOPLS, NBModelFull, QCModelMNDO
from pMoleculeScripts import LangevinDynamics_SystemGeometry, LBFGSMinimize_SystemGeometry, NormalModes_SystemGeometry, QuasiHarmonic_SystemGeometry

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Dynamics options.
_CollisionFrequency =     1.0
_LogFrequency       =   10000
_NSteps0            =    1000
_NSteps1            =  100000
_Seed               =  156451
_SaveFrequency      =     100
_Temperature        =  1000.0

# . Molecules.
_MoleculeLabels = ( "hydrogenFluoride", "water", "benzene", "cyclohexane" )
_QCModels       = ( "hydrogenFluoride", )

# . Optimization options.
_Iterations   = 1000
_Tolerance    = 1.0e-05

# . Paths.
_Destination = "vibrationalDynamics"

# . Tolerances.
_RMSAbsoluteErrorTolerance = 1.0e-03

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class VibrationalDynamicsTest ( TestCase ):
    """Test vibrational dynamics."""

    def runTest ( self ):
        """The test."""

        # . Paths.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "mol" )
        if self.resultPath is None: outPath  = os.path.join ( os.getenv ( "PDYNAMO_SCRATCH" ), _Destination )
        else:                       outPath  = os.path.join ( self.resultPath                , _Destination )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )
        log = self.GetLog ( )

        # . Energy models.
        mmModel = MMModelOPLS ( "protein" )
        nbModel = NBModelFull ( )
        qcModel = QCModelMNDO ( converger = DIISSCFConverger ( densityTolerance = 1.0e-10 ) )

        # . Initialization.
        numberFailures = 0

        # . Loop over molecules.
        for moleculeLabel in _MoleculeLabels:

            # . Get the system.
            system = MOLFile_ToSystem ( os.path.join ( dataPath, moleculeLabel + ".mol" ) )
            if moleculeLabel in _QCModels:
                system.DefineQCModel ( qcModel )
            else:
                system.DefineMMModel ( mmModel, log = log )
                system.DefineNBModel ( nbModel )
            system.Summary ( log = log )
            system.Energy  ( log = log )

            # . Minimize well.
            LBFGSMinimize_SystemGeometry ( system,
                                           log                  =           log ,
                                           logFrequency         = _LogFrequency ,
                                           maximumIterations    =   _Iterations ,
                                           rmsGradientTolerance =    _Tolerance )

            # . Normal mode analysis.
            nmState = NormalModes_SystemGeometry ( system, log = log )

            # . Do a dynamics simulation: equilibration and then data collection.
            normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( _Seed ) )
            LangevinDynamics_SystemGeometry ( system                                          ,
                                              collisionFrequency     = _CollisionFrequency    ,
                                              log                    = log                    ,
                                              logFrequency           = _LogFrequency          ,
                                              normalDeviateGenerator = normalDeviateGenerator ,
                                              steps                  = _NSteps0               ,
                                              temperature            = _Temperature           ,
                                              timeStep               = 0.001                  )
            reference3 = Clone ( system.coordinates3 )
            trajectory = AmberTrajectoryFileWriter ( os.path.join ( outPath, moleculeLabel + ".crd" ), system )
            LangevinDynamics_SystemGeometry ( system                                          ,
                                              collisionFrequency     = _CollisionFrequency    ,
                                              log                    = log                    ,
                                              logFrequency           = _LogFrequency          ,
                                              normalDeviateGenerator = normalDeviateGenerator ,
                                              steps                  = _NSteps1               ,
                                              temperature            = _Temperature           ,
                                              timeStep               = 0.001                  ,
                                              trajectories           = [ ( trajectory, _SaveFrequency ) ] )

            # . Check RMSs.
            masses = system.atoms.GetItemAttributes ( "mass" )
            rms0   = system.coordinates3.RMSDeviation ( reference3, weights = masses )
            system.coordinates3.Superimpose ( reference3, weights = masses )
            rms1 = system.coordinates3.RMSDeviation ( reference3, weights = masses )
            if ( math.fabs ( rms1 - rms0 ) >= _RMSAbsoluteErrorTolerance ): numberFailures += 1

            # . Do a quasi-harmonic analysis.
            trajectory = AmberTrajectoryFileReader ( os.path.join ( outPath, moleculeLabel + ".crd" ), system )
            qhState    = QuasiHarmonic_SystemGeometry ( system, log = log, temperature = _Temperature, trajectories = [ trajectory ] )

        # . Success/failure.
        self.assertTrue ( ( numberFailures == 0 ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = VibrationalDynamicsTest ( )
    test.runTest ( )
