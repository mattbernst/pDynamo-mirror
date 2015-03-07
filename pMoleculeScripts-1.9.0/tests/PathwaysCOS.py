"""Pathways tests."""

import math, os, os.path

from pBabel           import AmberTrajectory_FromSystemGeometryTrajectory, SystemGeometryTrajectory
from pCore            import CPUTime, TestCase
from pMoleculeScripts import ChainOfStatesOptimizePath_SystemGeometry, SGOFProcessPoolFactory

from PathTestSystems  import PathTestReportsSummary, PathTestSystems

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Pathways directory.
_pathways = "pathways"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PathwaysCOSTest ( TestCase ):
    """Calculate pathways between two conformations of a variety of systems."""

    def __init__ ( self, *args ):
        """Constructor."""
        super ( PathwaysCOSTest, self ).__init__ ( *args )
        # . Options.
        poolFactory = None # SGOFProcessPoolFactory ( maximumProcesses = 6, poolType = "Multiprocessing" )
        # . Attributes.
        self.doLong            = False
        self.doShort           = True
        self.doSuperimposition = False
        self.numberOfImages    = 11
        self.optimizerOptions  = { "fixedTerminalImages"                        : True        ,
                                   "forceOneSingleImageOptimization"            : False       ,
                                   "forceSingleImageOptimizations"              : False       ,
                                   "forceSplineRedistributionCheckPerIteration" : False       ,
                                   "freezeRMSGradientTolerance"                 :    0.0      ,
                                   "logFrequency"                               :    1        ,
                                   "maximumIterations"                          :  5000       ,
                                   "optimizer"                                  : None        ,
                                   "poolFactory"                                : poolFactory ,
                                   "rmsGradientTolerance"                       :    0.1      ,
                                   "splineRedistributionTolerance"              :    1.5      ,
                                   "springForceConstant"                        :    0.0      ,
                                   "useSplineRedistribution"                    : True        }

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK = True

        # . Paths.
        outPath = None
        if self.resultPath is None: outPath = os.path.join ( os.getenv ( "PDYNAMO_SCRATCH" ), _pathways )
        else:                       outPath = os.path.join ( self.resultPath, _pathways )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )
        log = self.GetLog ( )

        # . Loop over the tests.
        reports = {}
        for testSystem in PathTestSystems:

            # . Check whether to skip this system.
            if ( testSystem.isLong and ( not self.doLong ) ) or ( ( not testSystem.isLong ) and ( not self.doShort ) ): continue

            # . Get the system.
            system = testSystem.GetSystem ( log = log )
            system.Summary ( log = log )

            # . Create a starting trajectory.
            trajectoryPath = os.path.join ( outPath, "COS_" + testSystem.tag + ".trj" )
            testSystem.GetGrowingStringInitialPath ( self.numberOfImages, trajectoryPath, log = log )

            # . Superimpose the trajectory structures.
            if self.doSuperimposition:
                trajectory  = SystemGeometryTrajectory ( trajectoryPath, system, mode = "a+" )
                rmsResults0 = trajectory.Superimpose ( )

            # . Pathway.
            cpuTime        = CPUTime ( )
            options        = dict ( self.optimizerOptions )
            options["log"] = log
            trajectory     = SystemGeometryTrajectory ( trajectoryPath, system, mode = "a+" )
            report         = ChainOfStatesOptimizePath_SystemGeometry ( system, trajectory, **options )
            report["Time"] = cpuTime.Current ( )
            reports[testSystem.tag] = report

            # . Superimpose the trajectory structures.
            if self.doSuperimposition:
                trajectory  = SystemGeometryTrajectory ( trajectoryPath, system, mode = "a+" )
                rmsResults1 = trajectory.Superimpose ( )
                for ( ( a, b ), ( c, d ) ) in zip ( rmsResults0, rmsResults1 ):
                    print ( "{:10.3f} {:10.3f} {:10.3f}   | {:10.3f} {:10.3f} {:10.3f}".format ( a, b, math.fabs ( a - b ), c, d, math.fabs ( c - d ) ) )

            # . Convert the trajectory.
            amberPath = trajectoryPath[0:-4] + ".crd"
            AmberTrajectory_FromSystemGeometryTrajectory ( amberPath, trajectoryPath, system )

        # . Finish up.
        numberConverged = PathTestReportsSummary ( reports, log = log )
        self.assertTrue ( numberConverged == len ( reports ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = PathwaysCOSTest ( )
    test.run ( )
