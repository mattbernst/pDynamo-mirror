"""Definition of simple QC systems for testing."""

import glob, math, os.path

from pBabel    import MOLFile_ToSystem, XYZFile_ToSystem
from pCore     import logFile, LogFileActive, Selection
from pMolecule import DIISSCFConverger, ElectronicState

#===================================================================================================================================
# . Basic parameters.
#===================================================================================================================================
# . Paths.
_dataPath       = os.path.join ( os.getenv ( "PDYNAMO_PMOLECULE" ), "data" )
_structuresPath = os.path.join ( os.getenv ( "PDYNAMO_ROOT" ), "molecularStructures" )

#===================================================================================================================================
# . Basic closed shell tests.
#===================================================================================================================================
# . Charges.
_closedShellCharges = { "cadmiumComplex"  :  2  ,
                        "fch3cl"          : -1  ,
                        "histamine"       :  1  ,
                        "methylPhosphate" : -2  ,
                        "rhodiumComplex"  : -2  }

# . Sources.
_closedShellSources = ( os.path.join ( _dataPath                                      , "xyz" ), \
                        os.path.join ( _structuresPath, "difficultSCFCases"           , "xyz" ), \
                        os.path.join ( _structuresPath, "gaussianGeometryOptimization", "xyz" )  )

#===================================================================================================================================
# . Basic radical tests (includes closed and open shell cases).
#===================================================================================================================================
# . Systems.
_radicalKeywordLabels = ( "label", "fileName", "fileFormat", "dataPath", "charge", "multiplicity" )
_radicalSystems       = ( ( "Tyrosine Dipeptide Singlet", "tyrosineDipeptide", "xyz", "xyz"     , 0, 1 ), \
                          ( "Tyrosine Dipeptide Triplet", "tyrosineDipeptide", "xyz", "xyz"     , 0, 3 ), \
                          ( "Allyl Radical"             , "allyl"            , "xyz", "radicals", 0, 2 ), \
                          ( "Methylene Radical"         , "methylene"        , "xyz", "radicals", 0, 3 )  )

# . QC model options.
_radicalQCModelKeywords = { "Tyrosine Dipeptide Singlet" : ( { "isSpinRestricted" : True , "occupancyType" : "Cardinal"            }, \
                                                             { "isSpinRestricted" : False, "occupancyType" : "Cardinal"            }, \
                                                             { "isSpinRestricted" : True , "occupancyType" : "Fixed Fractional"    }, \
                                                             { "isSpinRestricted" : False, "occupancyType" : "Fixed Fractional"    }, \
                                                             { "isSpinRestricted" : True , "occupancyType" : "Variable Fractional" }, \
                                                             { "isSpinRestricted" : False, "occupancyType" : "Variable Fractional" } ), \
                            "Tyrosine Dipeptide Triplet" : ( { "isSpinRestricted" : False, "occupancyType" : "Cardinal"            }, \
                                                             { "isSpinRestricted" : False, "occupancyType" : "Fixed Fractional"    }, \
                                                             { "isSpinRestricted" : False, "occupancyType" : "Variable Fractional" } ), \
                            "Allyl Radical"              : ( { "isSpinRestricted" : False, "occupancyType" : "Cardinal"            },
                                                             { "isSpinRestricted" : False, "occupancyType" : "Fixed Fractional"    },
                                                             { "isSpinRestricted" : False, "occupancyType" : "Variable Fractional" } ), \
                            "Methylene Radical"          : ( { "isSpinRestricted" : False, "occupancyType" : "Cardinal"            }, \
                                                             { "isSpinRestricted" : False, "occupancyType" : "Fixed Fractional"    }, \
                                                             { "isSpinRestricted" : False, "occupancyType" : "Variable Fractional" } ) }

#===================================================================================================================================
# . Functions to return test sets.
#===================================================================================================================================
# . Closed shell systems.
def GetClosedShellMoleculeData ( ):
    """Get data for closed shell systems."""
    # . Get paths.
    paths = []
    for source in _closedShellSources:
        if os.path.exists ( source ):
            paths.extend ( glob.glob ( os.path.join ( source, "*.xyz" ) ) )
    # . Gather data.
    data = {}
    for path in paths:
        ( head, tail ) = os.path.split ( path )
        label    = tail[0:-4]
        fileName = label
        data[label] = { "label" : label, "fileName" : label, "fileFormat" : "xyz", "dataPath" : head, "charge" : _closedShellCharges.get ( label, 0 ), "multiplicity" : 1 }
    # . Finish up.
    labels = data.keys ( )
    labels.sort ( )
    closedShellMoleculeData = []
    for label in labels: closedShellMoleculeData.append ( data[label] )
    return closedShellMoleculeData

# . Radical systems.
def GetRadicalMoleculeData ( ):
    """Get data for radical systems."""
    radicalMoleculeData = []
    for values in _radicalSystems:
        label = values[0]
        for qcModelKeywords in _radicalQCModelKeywords[label]:
            items = { key : value for ( key, value ) in zip ( _radicalKeywordLabels, values ) }
            items["dataPath"] = os.path.join ( _dataPath, items["dataPath"] )
            items["qcModelKeywords"] = qcModelKeywords
            radicalMoleculeData.append ( items )
    return radicalMoleculeData

#===================================================================================================================================
# . Class for a QC test system.
#===================================================================================================================================
class QCTestSystem ( object ):
    """QC test system."""

    def __init__ ( self, **kwargs ):
        """Constructor."""
        for ( attribute, value ) in kwargs.iteritems ( ):
            setattr ( self, attribute, value )

    def GetSystem ( self, log = logFile, maximumAtoms = None ):
        """Get the system with the energy model defined."""
        # . Get the QC model options.
        convergerKeywords = getattr ( self, "convergerKeywords",   {} )
        qcModelArguments  = getattr ( self, "qcModelArguments" ,   [] )
        qcModelClass      = getattr ( self, "qcModelClass"     , None )
        qcModelKeywords   = getattr ( self, "qcModelKeywords"  ,   {} )
        # . Basic setup.
        if   self.fileFormat == "mol": molecule = MOLFile_ToSystem  ( os.path.join ( self.dataPath, self.fileName + ".mol" ) )
        elif self.fileFormat == "xyz": molecule = XYZFile_ToSystem  ( os.path.join ( self.dataPath, self.fileName + ".xyz" ) )
        molecule.electronicState = ElectronicState ( charge = getattr ( self, "charge", 0 ), multiplicity = getattr ( self, "multiplicity", 1 ) )
        molecule.label           = self.label
        # . Only keep the molecule if it is not too large.
        if ( maximumAtoms is None ) or ( ( maximumAtoms is not None ) and ( len ( molecule.atoms ) <= maximumAtoms ) ):
            # . Define the QC model.
            if qcModelClass is not None:
                converger           = DIISSCFConverger ( **convergerKeywords )
                kwargs              = dict ( qcModelKeywords )
                kwargs["converger"] = converger
                qcModel = qcModelClass ( *qcModelArguments, **kwargs )
                molecule.DefineQCModel ( qcModel, log = log )
            # . Summary.
            if LogFileActive ( log ):
                molecule.Summary ( log = log )
                log.Paragraph ( "\nFormula = " + molecule.atoms.FormulaString ( ) + "." )
        # . Molecule rejected.
        else: molecule = None
        # . Finish up.
        return molecule

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
