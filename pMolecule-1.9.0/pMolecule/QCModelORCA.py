#-------------------------------------------------------------------------------
# . File      : QCModelORCA.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""The ORCA QC model."""

import glob, math, os, os.path, subprocess, re

from pCore            import Coordinates3, logFile, LogFileActive, Real1DArray, Vector3, \
                             UNITS_DIPOLE_ATOMIC_UNITS_TO_DEBYES, UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE, UNITS_LENGTH_ANGSTROMS_TO_BOHRS
from Element          import PeriodicTable
from NBModelORCAState import NBModelORCAState
from QCModel          import QCModel

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Default cleanup option.
_DefaultCleanUp = True

# . Default command.
# . The directory where the programs are found must be on the bin path.
_DefaultCommand = "orca"

# . Default error suffix.
_DefaultErrorPrefix = "error_"

# . Default job name.
_DefaultJobName = "job"

# . Default keywords.
_DefaultKeywords = [ "BP86", "DefBas4", "FinalGrid4", "Grid3", "RI", "SCFConv8", "TightSCF" ]

# . Default label.
_DefaultLabel       = "BP86:DefBas4"
_DefaultLabelPrefix = "ORCA:"

# . Default scratch.
_DefaultScratch = os.getenv ( "PDYNAMO_SCRATCH" )

# . Default use random job.
_DefaultUseRandomJob = False

# . Default use random scratch.
_DefaultUseRandomScratch = False

#===================================================================================================================================
# . Utility function.
#===================================================================================================================================
import random, string
def RandomString ( characters = ( string.ascii_lowercase + string.digits ), size = 12, startingCharacters = None ):
    """Generate a random string."""
    if startingCharacters is None: startingCharacters = string.ascii_lowercase
    return "".join ( [ random.choice ( startingCharacters ) ] + [ random.choice ( characters ) for x in range ( size - 1 ) ] )

#===================================================================================================================================
# . Error class.
#===================================================================================================================================
class ORCAError ( Exception ):
    pass

#===================================================================================================================================
# . Class for the model state.
#===================================================================================================================================
class QCModelORCAState ( object ):

    def __del__ ( self ):
        """Deallocation."""
        self.DeleteJobFiles ( )
#        super ( QCModelORCAState, self ).__del__ ( )

    def __init__ ( self ):
        """Constructor."""
        pass

    def DeleteJobFiles ( self ):
        """Delete job files."""
        if self.deleteJobFiles:
            try:
                jobFiles = glob.glob ( os.path.join ( self.jobPath + ".*" ) )
                for jobFile in jobFiles: os.remove ( jobFile )
                if self.scratchPath is not None: os.rmdir ( self.scratchPath )
            except:
                pass

    def Finalize ( self, qcAtoms ):
        """Finish up a calculation."""
        # . Gradients.
        if self.gradients3 is not None:
            qcAtoms.SetGradients3 ( self.coordinates3, self.qcGradients3, self.gradients3, toInternalUnits = True )

    def Initialize ( self, qcAtoms, coordinates3, gradients3 ):
        """Initialize a calculation."""
        # . Initialize results.
        self.cycles          = 0
        self.energy          = 0.0
        self.isConverged     = False
        self.spinSquared     = 0.0
        self.mayerBondOrders = []
        self.chelpgCharges.Set   ( 0.0 )
        self.dipole.Set          ( 0.0 )
        self.lowdinCharges.Set   ( 0.0 )
        self.lowdinSpins.Set     ( 0.0 )
        self.mullikenCharges.Set ( 0.0 )
        self.mullikenSpins.Set   ( 0.0 )
        # . Coordinates.
        self.coordinates3 = coordinates3
        qcAtoms.FillCoordinates3 ( coordinates3, self.qcCoordinates3 )
        # . Gradients.
        self.gradients3 = gradients3
        self.qcGradients3.Set ( 0.0 )

    def SaveErrorFiles ( self, message  ):
        """Save the input and output files for inspection if there is an error."""
        for path in ( self.engradPath, self.inputPath, self.outputPath, self.pcgradPath, self.pcPath ):
            if os.path.exists ( path ):
                ( head, tail ) = os.path.split ( path )
                os.rename ( path, os.path.join ( head, _DefaultErrorPrefix + tail ) )
        raise ORCAError ( message + "\nCheck the files \"{:s}/{:s}*\".".format ( self.scratchPath, _DefaultErrorPrefix ) )

    @classmethod
    def Setup ( selfClass, qcAtoms, electronicState, nbState, job, scratch, deleteJobFiles = False, useRandomJob = False, useRandomScratch = False ):
        """Set up a state given the relevant data."""
        self = selfClass ( )
        # . QC atoms.
        # . Basic attributes.
        self.charge          = electronicState.charge
        self.multiplicity    = electronicState.multiplicity
        self.numberOfQCAtoms = len ( qcAtoms )
        # . File paths.
        if useRandomJob: job = RandomString ( )
        if useRandomScratch:
            self.scratchPath = os.path.join ( scratch, RandomString ( ) )
            if not os.path.exists ( self.scratchPath ): os.mkdir ( self.scratchPath )
            self.jobPath     = os.path.join ( self.scratchPath, job )
        else:
            self.scratchPath = None
            self.jobPath     = os.path.join ( scratch, job )
        self.engradPath = self.jobPath + ".engrad"
        self.inputPath  = self.jobPath + ".inp"
        self.outputPath = self.jobPath + ".log"
	self.pcgradPath = self.jobPath + ".pcgrad"
	self.pcPath     = self.jobPath + ".pc"
        # . Options.
        self.deleteJobFiles = deleteJobFiles
        # . Symbols.
        atomicNumbers = qcAtoms.GetAtomicNumbers ( )
        self.symbols  = []
        for i in atomicNumbers:
            self.symbols.append ( PeriodicTable.Symbol ( i ) )
        # . Array attributes.
        self.qcCoordinates3  = Coordinates3.WithExtent ( self.numberOfQCAtoms ) ; self.qcCoordinates3.Set  ( 0.0 )
        self.dipole          = Vector3.Null ( )
        self.qcGradients3    = Coordinates3.WithExtent ( self.numberOfQCAtoms ) ; self.qcGradients3.Set    ( 0.0 )
        self.chelpgCharges   = Real1DArray.WithExtent  ( self.numberOfQCAtoms ) ; self.chelpgCharges.Set   ( 0.0 )
        self.lowdinCharges   = Real1DArray.WithExtent  ( self.numberOfQCAtoms ) ; self.lowdinCharges.Set   ( 0.0 )
        self.lowdinSpins     = Real1DArray.WithExtent  ( self.numberOfQCAtoms ) ; self.lowdinSpins.Set     ( 0.0 )
        self.mullikenCharges = Real1DArray.WithExtent  ( self.numberOfQCAtoms ) ; self.mullikenCharges.Set ( 0.0 )
        self.mullikenSpins   = Real1DArray.WithExtent  ( self.numberOfQCAtoms ) ; self.mullikenSpins.Set   ( 0.0 )
        # . Orbital data.
        self.HOMO            = -1
        self.LUMO            = -1
        self.orbitalEnergies = []
        # . NB state.
        self.nbState = nbState
        return self

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelORCA ( QCModel ):
    """The ORCA QC model class."""

    defaultAttributes = { "command"          : _DefaultCommand          ,
                          "deleteJobFiles"   : _DefaultCleanUp          ,
                          "job"              : _DefaultJobName          ,
                          "keywords"         : _DefaultKeywords         ,
                          "label"            : _DefaultLabelPrefix + _DefaultLabel ,
                          "scratch"          : _DefaultScratch          ,
                          "useRandomJob"     : _DefaultUseRandomJob     ,
                          "useRandomScratch" : _DefaultUseRandomScratch }

    defaultAttributeNames = { "Command"             : "command"          ,
                              "Delete Job Files"    : "deleteJobFiles"   ,
                              "Job Name"            : "job"              ,
                              "Keywords"            : "keywords"         ,
                              "Label"               : "label"            ,
                              "Scratch"             : "scratch"          ,
                              "Use Random Job Name" : "useRandomJob"     ,
                              "Use Random Scratch"  : "useRandomScratch" }

    def __getstate__ ( self ): return self.__dict__

    def _Initialize ( self ):
        """Initialization."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def AtomicCharges ( self, configuration, chargeModel = "Mulliken", spinDensities = False ):
        """Atomic charges."""
        qcState = getattr ( configuration, "qcState", None )
        if qcState is None:
            return None
        else:
            if spinDensities:
                if   chargeModel == "Lowdin": return getattr ( qcState, "lowdinSpins"    , None )
                else:                         return getattr ( qcState, "mullikenSpins"  , None )
            else:
                if   chargeModel == "Chelpg": return getattr ( qcState, "chelpgCharges"  , None )
                elif chargeModel == "Lowdin": return getattr ( qcState, "lowdinCharges"  , None )
                else:                         return getattr ( qcState, "mullikenCharges", None )

    def Clear ( self, configuration ):
        """Clear up temporary data."""
        if configuration is not None:
            qcState = getattr ( configuration, "qcState", None )
            if qcState is not None:
                qcState.DeleteJobFiles ( )
                delattr ( configuration, "qcState" )

    def DipoleMoment ( self, configuration, center = None ):
        """Dipole Moment."""
        qcState = getattr ( configuration, "qcState", None )
        if qcState is None: return None
        else:               return getattr ( qcState, "dipole", None )

    def Energy ( self, qcAtoms, qcParameters, electronicState, configuration, log = logFile ):
        """Calculate the quantum chemical energy."""
        # . Check for a valid NB state.
        nbState = getattr ( configuration, "nbState", None )
        if ( nbState is not None ) and not isinstance ( nbState, NBModelORCAState ): raise ORCAError ( "Invalid non-bonding model state." )
        # . Check for a valid state.
        qcState = getattr ( configuration, "qcState", None )
        if qcState is None:
            qcState = QCModelORCAState.Setup ( qcAtoms, electronicState, nbState, self.job, self.scratch, deleteJobFiles = self.deleteJobFiles, \
                                                                     useRandomJob = self.useRandomJob, useRandomScratch = self.useRandomScratch )
            setattr ( configuration, "qcState", qcState )
        # . Get the coordinates.
        coordinates3 = getattr ( configuration, "coordinates3", None )
        gradients3   = getattr ( configuration, "gradients3",   None )
        qcState.Initialize ( qcAtoms, coordinates3, gradients3 )
        # . Do the calculation.
        self.WriteInputFile ( qcState )
        isOK = self.Execute ( qcState )
        if not isOK: qcState.SaveErrorFiles ( "Error executing program." )
        # . Read the output files.
        if qcState.gradients3 is not None:
            isOK = self.ReadEngradFile ( qcState )
            if not isOK: qcState.SaveErrorFiles ( "Error reading engrad file." )
        isOK = self.ReadOutputFile ( qcState )
        if not isOK: qcState.SaveErrorFiles ( "Error reading output file." )
        # . Deal with the gradients.
        qcState.Finalize ( qcAtoms )
        # . Return.
        return [ ( "ORCA QC", qcState.energy * UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE ) ]

    # . File executing/reading/writing methods.
    def Execute ( self, qcState ):
        """Execute the job."""
        try:
#            currentdirectory = os.getcwd ( )
#            os.chdir ( _DefaultDIRECTORY )
            outFile = open ( qcState.outputPath, "w" )
            subprocess.check_call ( [ self.command, qcState.inputPath ], stderr = outFile, stdout = outFile )
            outFile.close ( )
#            os.chdir ( currentdirectory )
            return True
        except:
            return False

    def IsSpinRestricted ( self ):
        """Is the model spin restricted?"""
        return False

    def MakeLabel ( self ):
        """Construct a model label."""
        pass

    def MayerBondOrders ( self, qcAtoms, qcParameters, configuration ):
        """Mayer Bond Orders."""
        qcState = getattr ( configuration, "qcState", None )
        if qcState is None: return None
        else:               return getattr ( qcState, "mayerBondOrders", None )

    def OrbitalEnergies ( self, configuration ):
        """Orbital energies and HOMO and LUMO indices."""
        qcState = getattr ( configuration, "qcState", None )
        if qcState is None: return ( None, -1, -1 )
        else:
            return ( getattr ( qcState, "orbitalEnergies", None ), \
                     getattr ( qcState, "HOMO", -1 ), \
                     getattr ( qcState, "LUMO", -1 ) )

    def ProcessArguments ( self, *arguments ):
        """Process constructor arguments."""
        if len ( arguments ) > 0:
            self.keywords = arguments[0].split ( ":" )
            self.label    = _DefaultLabelPrefix + arguments[0]
            if len ( arguments ) > 1: self.keywords.extend ( arguments[1:] )

    def ReadEngradFile ( self, qcState ):
        """Read an engrad file."""
        try:
            egFile = open ( qcState.engradPath, "r" )
            # . Skip the number of atoms section and the energy header.
            for i in range ( 7 ): next ( egFile )
            # . Get the energy.
            qcState.energy = float ( next ( egFile ) )
            # . Skip the gradients header.
            for i in range ( 3 ): ( next ( egFile ) )
            # . Get the gradients.
            for i in range ( qcState.numberOfQCAtoms ):
                for j in range ( 3 ):
                    qcState.qcGradients3[i,j] = float ( ( next ( egFile ) ) )
            egFile.close ( )
            if qcState.nbState is not None: qcState.nbState.ReadPCgradFile ( qcState.pcgradPath )
            return True
        except:
            return False

    def ReadOutputFile ( self, qcState ):
        """Read an output file."""
        try:
            qcState.isConverged = False
            outFile = open ( qcState.outputPath, "r" )
            while True:
                try:
                    line = next ( outFile ).strip ( )
                    # . CHELPG charges.
                    if line == "CHELPG Charges":
                        line = next ( outFile )
                        for i in range ( qcState.numberOfQCAtoms ):
                            words = next ( outFile ).split ( ":", 1 )
                            qcState.chelpgCharges[i] = float ( words[-1] )
                    # . Convergence OK.
                    elif line.find ( "SCF CONVERGED AFTER" ) >= 0:
                        words              = line.split ( )
                        qcState.cycles     = int ( words[-3] )
                        qcState.isConverged = True
                    # . Convergence not OK.
                    elif line.find ( "SCF NOT CONVERGED AFTER" ) >= 0:
                        words          = line.split ( )
                        qcState.cycles = int ( words[-3] )
                    # . Dipole.
                    elif line.startswith ( "Total Dipole Moment" ):
                        words = line.split ( )
                        for ( i, word ) in enumerate ( words[-3:] ):
                            qcState.dipole[i] = UNITS_DIPOLE_ATOMIC_UNITS_TO_DEBYES * float ( word )
                    # . Energy.
                    elif line.startswith ( "FINAL SINGLE POINT ENERGY" ):
                        qcState.energy = float ( line.split ( )[-1] )
                    # . Lowdin charges.
                    elif line == "LOEWDIN ATOMIC CHARGES":
                        next ( outFile )
                        for i in range ( qcState.numberOfQCAtoms ):
                            words = next ( outFile ).split ( ":", 1 )
                            qcState.lowdinCharges[i] = float ( words[-1] )
                    # . Lowdin charges and spin densities.
                    elif line.startswith ( "LOEWDIN ATOMIC CHARGES AND SPIN " ):
                        next ( outFile )
                        for i in range ( qcState.numberOfQCAtoms ):
                            words = next ( outFile ).split ( )
                            qcState.lowdinCharges[i] = float ( words[-2] )
                            qcState.lowdinSpins[i]   = float ( words[-1] )
                    # . Mayer bond orders.
                    elif line.startswith ( "Mayer bond orders larger than" ):
                        qcState.mayerBondOrders = []
                        while True:
			    line = next ( outFile )
			    atompairs=re.findall(r'(\d+)-\w+\s*,\s*(\d+)-\w+', line)
			    bondOrders=p=re.findall(r'(\d\.\d+)', line)
			    if atompairs:
			    	for (pair, bond) in zip(atompairs, bondOrders):
			    	    i=int(pair[0])
				    j=int(pair[1])
				    qcState.mayerBondOrders.append ( ( i, j, float ( bond ) ) )
                            else:
                                break
                    # . Mulliken charges.
                    elif line == "MULLIKEN ATOMIC CHARGES":
                        line = next ( outFile )
                        for i in range ( qcState.numberOfQCAtoms ):
                            words = next ( outFile ).split ( ":", 1 )
                            qcState.mullikenCharges[i] = float ( words[-1] )
                    # . Mulliken charges and spin densities.
                    elif line.startswith ( "MULLIKEN ATOMIC CHARGES AND SPIN " ):
                        next ( outFile )
                        for i in range ( qcState.numberOfQCAtoms ):
                            words = next ( outFile ).split ( )
                            qcState.mullikenCharges[i] = float ( words[-2] )
                            qcState.mullikenSpins[i]   = float ( words[-1] )
                    # . Orbital energies.
                    elif line == "ORBITAL ENERGIES":
                        HOMO = -1
                        LUMO = -1
                        orbitalEnergies = []
                        for i in range ( 3 ): next ( outFile )
                        index = 0
                        while True:
                            tokens = next ( outFile ).split ( )
                            if len ( tokens ) < 3: break
                            else:
                                occupancy = float ( tokens[1] )
                                if ( HOMO == -1 ) and ( LUMO == -1 ) and ( occupancy <= 1.0e-6 ):
                                    HOMO = index - 1
                                    LUMO = index
                                orbitalEnergies.append ( float ( tokens[2] ) )
                                index += 1
                        energies = Real1DArray.WithExtent ( len ( orbitalEnergies ) )
                        energies.Set ( 0.0 )
                        for ( i, e ) in enumerate ( orbitalEnergies ): energies[i] = e
                        qcState.HOMO = HOMO
                        qcState.LUMO = LUMO
                        qcState.orbitalEnergies = energies
                    # . <S**2>.
                    elif line.startswith ( "Expectation value of <S**2>" ):
                        qcState.spinSquared = float ( line.split ( ":", 1 )[-1] )
                except StopIteration:
                    break
            outFile.close ( )
            return qcState.isConverged
        except Exception as e:
            print ( e[0] )
            return False

    def SetOptions ( self, **options ):
        """Set options for the model."""
        if "command"          in options: self.command          = options.pop ( "command"          )
        if "deleteJobFiles"   in options: self.deleteJobFiles   = options.pop ( "deleteJobFiles"   )
        if "job"              in options: self.job              = options.pop ( "job"              )
        if "keywords"         in options: self.keywords         = options.pop ( "keywords"         )
        if "label"            in options: self.label            = options.pop ( "label"            )
        if "scratch"          in options: self.scratch          = options.pop ( "scratch"          )
        if "useRandomJob"     in options: self.useRandomJob     = options.pop ( "useRandomJob"     )
        if "useRandomScratch" in options: self.useRandomScratch = options.pop ( "useRandomScratch" )
        if len ( options ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( options.keys ( ) ) ) + "." )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "QC Model ORCA" )
            keys = self.__class__.defaultAttributeNames.keys ( )
            keys.sort ( )
            for key in keys:
                value = getattr ( self, self.__class__.defaultAttributeNames[key] )
                if   isinstance ( value, bool  ):
                    if value: valuestring = "True"
                    else:     valuestring = "False"
                elif isinstance ( value, float      ): valuestring = "{:g}".format ( value )
                elif isinstance ( value, basestring ): valuestring =  value
                else:                                  valuestring = str ( value )
                summary.Entry ( key, valuestring )
            summary.Stop ( )

    def WriteInputFile ( self, qcState, mode = None ):
        """Write an input file."""
        inFile = open ( qcState.inputPath, "w" )
        inFile.write ( "#\n" )
        inFile.write ( "# ORCA Job.\n" )
        inFile.write ( "#\n" )
        if qcState.nbState is not None: inFile.write ( '%pointcharges "' + qcState.pcPath + '"\n' )
        if mode is None:
            if qcState.gradients3 is not None: mode = "ENGRAD"
            else:                              mode = "ENERGY"
        inFile.write ( "! " + mode + " " + " ".join ( self.keywords ) + "\n" )
        inFile.write ( "* xyz {:d} {:d}\n".format ( qcState.charge, qcState.multiplicity ) )
        for i in range ( qcState.numberOfQCAtoms ):
            inFile.write ( "{:<12s}{:20.10f}{:20.10f}{:20.10f}\n".format ( qcState.symbols[i], qcState.qcCoordinates3[i,0], qcState.qcCoordinates3[i,1], qcState.qcCoordinates3[i,2] ) )
        if qcState.nbState is not None: qcState.nbState.WriteInputFile ( qcState.pcPath )
        inFile.write ( "*\n" )
        inFile.close ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

