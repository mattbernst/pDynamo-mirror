"""pDynamo system benchmarks."""

import glob, os, os.path

from pBabel           import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, XYZFile_ToCoordinates3
from pCore            import logFile, CPUTime, Selection, YAMLUnpickle
from pMolecule        import CrystalClassCubic, ElectronicState, NBModelABFS, QCModelMNDO, SystemWithTimings
from pMoleculeScripts import LangevinDynamics_SystemGeometry, VelocityVerletDynamics_SystemGeometry

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . File and path names.
_TestDataFileName = "systemData.yaml"
_TestRootPath     = "data"

# . Other options.
_CrystalClasses       = { "Cubic" : CrystalClassCubic ( ) }
_DoDynamics           = True
_ForceNoQC            = False
_NBModel              = NBModelABFS ( )
_QCModel              = QCModelMNDO ( )
_QCRegionKey          = "Large QC Region" # "Small QC Region"
_Steps                = 1000
_UseLangevin          = True
_UseSystemWithTimings = True

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def FindTests ( ):
    """Find the tests."""
    paths = glob.glob ( os.path.join ( _TestRootPath, "*" ) )
    tests = []
    for path in paths:
        if os.path.exists ( os.path.join ( path, _TestDataFileName ) ): tests.append ( path )
    tests.sort ( )
    return tests

def SetUpSystem ( path, forceNoQC = False, useSystemWithTimings = True ):
    """Set up the system."""
    # . Get system data.
    systemData = YAMLUnpickle ( os.path.join ( path, _TestDataFileName ) )
    # . Get the parameters.
    parameters = CHARMMParameterFiles_ToParameters ( glob.glob ( os.path.join ( path, "*.prm" ) ) )
    # . Get the test name.
    name       = os.path.split ( path )[-1]
    # . Generate the system.
    system              = CHARMMPSFFile_ToSystem ( os.path.join ( path, name + ".psfx" ), isXPLOR = True, parameters = parameters )
    if useSystemWithTimings: system = SystemWithTimings.FromSystem ( system )
    system.coordinates3 = XYZFile_ToCoordinates3 ( os.path.join ( path, name + ".xyz" ) )
    system.label        = systemData["Label"]
    # . Symmetry.
    if systemData.get ( "Crystal Class", None ) is not None:
        symmetryOptions                 = systemData["Symmetry Parameters"]
        symmetryOptions["crystalClass"] = _CrystalClasses[systemData["Crystal Class"]]
        system.DefineSymmetry ( **symmetryOptions )
    # . QC data.
    if not forceNoQC:
        qcData = systemData.get ( _QCRegionKey, None )
        if qcData is not None:
            # . Electronic state.
            system.electronicState = ElectronicState ( charge       = qcData.get ( "Charge"      , 0 ) ,
                                                       multiplicity = qcData.get ( "Multiplicity", 1 ) )
            # . QC atoms.
            qcAtoms = set ( )
            for path in qcData["Atoms"]:
                index = system.sequence.AtomIndex ( path )
                qcAtoms.add ( index )
            system.DefineQCModel ( _QCModel, qcSelection = Selection ( qcAtoms ) )
    # . Finish set up.
    system.DefineNBModel ( _NBModel )
    return system

#===================================================================================================================================
# . Execution.
#===================================================================================================================================
# . Find the tests to run.
tests = FindTests ( )

# . Loop over tests.
results = []
for test in tests:

    # . Get the system.
    system = SetUpSystem ( test, forceNoQC = _ForceNoQC, useSystemWithTimings = _UseSystemWithTimings )

    # . Output.
    print ( "\n" + 80 * "=" + "\n" + "Test for \"" + system.label + "\":\n" + 80 * "=" )
    system.Summary ( )

    # . Do the test.
    cpu = CPUTime ( )
    if _UseSystemWithTimings: system.TimingStart ( )
    system.Energy ( doGradients = True )
    if _DoDynamics:
        if _UseLangevin:
            LangevinDynamics_SystemGeometry       ( system                                 ,
                                                    collisionFrequency        =       25.0 ,
                                                    logFrequency              =        100 ,
                                                    steps                     =     _Steps ,
                                                    temperature               =      300.0 ,
                                                    timeStep                  =      0.001 )
        else:
            VelocityVerletDynamics_SystemGeometry ( system                                 ,
                                                    logFrequency              =        100 ,
                                                    steps                     =     _Steps ,
                                                    temperatureScaleFrequency =        100 ,
                                                    temperatureScaleOption    = "constant" ,
                                                    temperatureStart          =      300.0 ,
                                                    timeStep                  =      0.001 )
    results.append ( ( system.label, cpu.CurrentAsString ( ) ) )
    if _UseSystemWithTimings: system.TimingSummary ( orderByMagnitude = True )
    system.configuration.nbState.StatisticsSummary ( )

# . Finish up.
table = logFile.GetTable ( columns = [ 20, 30 ] )
table.Start   ( )
table.Title   ( "Test Results" )
table.Heading ( "Test" )
table.Heading ( "Time" )
for ( label, time ) in results:
    table.Entry ( label, alignment = "left" )
    table.Entry ( time )
table.Stop ( )
