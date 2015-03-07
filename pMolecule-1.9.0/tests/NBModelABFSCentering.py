"""Test centering for the ABFS NB model."""

import glob, math, os, os.path

from pBabel           import DCDTrajectoryFileReader, DCDTrajectoryFileWriter, ImportSystem
from pCore            import NormalDeviateGenerator, RandomNumberGenerator, TestCase
from pMolecule        import MMModelOPLS, NBModelABFS
from pMoleculeScripts import LangevinDynamics_SystemGeometry

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_Destination = "NBModelABFSCentering"
_NLog        =  100
_NSave       =  100
_NSteps      = 1000
_Tolerance   = 0.001

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NBABFSCenteringTest ( TestCase ):
    """Test centering for the ABFS NB model."""

    def runTest ( self ):
        """The test."""

        # . Paths.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data", "mol2" )
        if self.resultPath is None: outPath  = os.path.join ( os.getenv ( "PDYNAMO_SCRATCH" ), _Destination )
        else:                       outPath  = os.path.join ( self.resultPath, _Destination )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )
        log = self.GetLog ( )

        # . NB models.
        nbModelNoC = NBModelABFS ( useCentering = False )
        nbModelC   = NBModelABFS ( useCentering = True  )

        # . Set up the system.
        system = ImportSystem ( os.path.join ( dataPath, "waterBox.mol2" ), log = log )
        system.DefineMMModel ( MMModelOPLS ( "bookSmallExamples" ) )
        system.DefineNBModel ( nbModelC )
        system.Summary ( log = log )
        system.Energy  ( log = log )

        # . Do a short dynamics.
        # . Define a normal deviate generator in a given state.
        normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 614108 ) )

        # . Dynamics.
        trajectoryPath = os.path.join ( outPath, "waterBox_C_1ps.dcd" )
        trajectory     = DCDTrajectoryFileWriter ( trajectoryPath, system )
        LangevinDynamics_SystemGeometry ( system                            ,
                                          collisionFrequency     =     25.0 ,
                                          log                    =      log ,
                                          logFrequency           =    _NLog ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  =  _NSteps ,
                                          temperature            =    300.0 ,
                                          timeStep               =    0.001 ,
                                          trajectories = [ ( trajectory, _NSave ) ] )

        # . Calculate trajectory energies with different NB models.
        energies = []
        for ( i, nbModel ) in enumerate ( ( nbModelC, nbModelNoC ) ):
            system.DefineNBModel ( nbModel )
            trajectory = DCDTrajectoryFileReader ( trajectoryPath, system )
            trajectory.ReadHeader ( )
            e = []
            while trajectory.RestoreOwnerData ( ):
                e.append ( system.Energy ( log = None ) )
            trajectory.Close ( )
            energies.append ( e )

        # . Check deviations.
        ( e0, e1 ) = energies
        maximumDeviation = 0.0
        for i in range ( len ( e0 ) ):
            maximumDeviation = max ( maximumDeviation, math.fabs ( e0[i] - e1[i] ) )

        # . Summary of results.
        if log is not None:
            log.Paragraph ( "Maximum deviation = {:.5f}".format ( maximumDeviation ) )

        # . Success/failure.
        self.assertTrue ( ( maximumDeviation <= _Tolerance ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = NBABFSCenteringTest ( )
    test.run ( )
