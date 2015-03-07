#-------------------------------------------------------------------------------
# . File      : JaguarOutputFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for reading Jaguar output files."""

import os.path

from pCore        import Coordinates3, logFile, LogFileActive, Real1DArray, SymmetricMatrix, TextFileReader, Vector3
from ExportImport import _Importer
from pMolecule    import PeriodicTable, System

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class JaguarOutputFileReader ( TextFileReader ):
    """JaguarOutputFileReader is the class for Jaguar output files that are to be read."""

    defaultattributes = { "atomicNumbers"       : None     ,
                          "bestenergy"          : 1.0e+300 ,
                          "bestcoordinates3"    : None     ,
                          "charge"              : 0        ,
                          "coordinates3"        : None     ,
                          "dipole"              : None     ,
                          "displacementmaximum" : None     ,
                          "displacementrms"     : None     ,
                          "ecpelectrons"        : None     ,
                          "energy"              : 0.0      ,
                          "energychange"        : None     ,
                          "gradientmaximum"     : None     ,
                          "gradientrms"         : None     ,
                          "mp2energy"           : 0.0      ,
                          "multiplicity"        : 1        ,
                          "natoms"              : 0        ,
                          "nbasis"              : 0        ,
                          "nfunctions"          : None     ,
                          "ngeometries"         : 0        ,
                          "overlap"             : None     ,
                          "QCOMPLETED"          : False    ,
                          "QECP"                : False    ,
                          "qesp"                : None     ,
                          "QFATALERROR"         : False    ,
                          "qmulliken"           : None     ,
                          "QOPTIMIZATION"       : False    ,
                          "qspin"               : None     ,
                          "scfenergy"           : None     }
    defaultattributes.update ( TextFileReader.defaultattributes )

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Loop over the lines.
                while True:
                    line = self.GetLine ( QWARNING = False )
                    if   len ( line ) == 0: continue
                    elif line.find ( "Atomic charges from electrostatic potential:"      ) >= 0: self.qesp                = self.ParseCharges ( )
                    elif line.find ( "Atomic charges from Mulliken population analysis:" ) >= 0: self.qmulliken           = self.ParseCharges ( )
                    elif line.find ( "Atomic Spin Densities from Mulliken analysis:"     ) >= 0: self.qspin               = self.ParseCharges ( )
                    elif line.find ( "  energy change:"                                  ) >= 0: self.energychange        = self.ParseConvergenceCriterium ( line )
                    elif line.find ( "  gradient maximum:"                               ) >= 0: self.gradientmaximum     = self.ParseConvergenceCriterium ( line )
                    elif line.find ( "  gradient rms:"                                   ) >= 0: self.gradientrms         = self.ParseConvergenceCriterium ( line )
                    elif line.find ( "  displacement maximum:"                           ) >= 0: self.displacementmaximum = self.ParseConvergenceCriterium ( line )
                    elif line.find ( "  displacement rms:"                               ) >= 0: self.displacementrms     = self.ParseConvergenceCriterium ( line )
                    elif line.find ( "Effective Core Potential"                          ) >= 0: self.ParseECPs      ( )
                    elif line.find ( "Gaussian Functions - Normalized coefficients"      ) >= 0: self.ParseFunctions ( )
                    elif line.find ( "geometry:"                                         ) >= 0: self.ParseGeometry  ( )
                    elif line.find ( "Moments from quantum mechanical wavefunction:"     ) >= 0: self.ParseDipole    ( )
                    elif line.find ( "olap"                                              ) >= 0: self.ParseOverlap   ( )
                    elif line.find ( "ERROR: fatal error -- debug information follows"   ) >= 0: self.QFATALERROR   = True
                    elif line.find ( "  energy:"                                         ) >= 0: self.energy        = float ( line.split ( )[1]  )
                    elif line.find ( "number of alpha electrons...."                     ) >= 0: self.nalpha        = int   ( line.split ( )[-1] )
                    elif line.find ( "number of beta electrons....."                     ) >= 0: self.nbeta         = int   ( line.split ( )[-1] )
                    elif line.find ( "functions...."                                     ) >= 0: self.nbasis        = int   ( line.split ( )[-1] )
                    elif line.find ( "Geometry optimization complete"                    ) >= 0: self.QCOMPLETED    = True
                    elif line.find ( "geometry optimization step"                        ) >= 0: self.QOPTIMIZATION = True
                    elif line.find ( "multiplicity:"                                     ) >= 0: self.multiplicity  = int   ( line.split ( )[-1] )
                    elif line.find ( "net molecular charge:"                             ) >= 0: self.charge        = int   ( line.split ( )[-1] )
                    elif line.find ( "SCFE: SCF energy: DFT"                             ) >= 0: self.scfenergy     = float ( line.split ( )[4]  )
                    elif line.find ( "SCFE: SCF energy: UDFT"                            ) >= 0: self.scfenergy     = float ( line.split ( )[4]  )
                    elif line.find ( "SCFE:  Solution phase energy: DFT"                 ) >= 0: self.scfenergy     = float ( line.split ( )[5]  )
                    elif line.find ( "SCFE:  Solution phase energy: UDFT"                ) >= 0: self.scfenergy     = float ( line.split ( )[5]  )
                    elif line.find ( "stopping now - optimization seems to be stuck"     ) >= 0: self.QFATALERROR   = True
                    elif line.find ( "Total LMP2 Energy......."                          ) >= 0: self.mp2energy     = float ( line.split ( )[3]  )
                # . Shuffle data if there is a geometry optimization.
                if self.QOPTIMIZATION:
                    if self.energy < self.bestenergy: # . Note not <=.
                        self.bestcoordinates3 = self.coordinates3
                        self.bestenergy       = self.energy
                # . Check the ecp data.
                if self.QECP:
                    if len ( self.ecpelectrons ) != self.natoms: self.Warning ( "Mismatch between number of atoms and ECP atom data.", False )
            except EOFError:
                pass
            except Exception as e:
                print ( e )
	    # . Close the file.
            self.WarningStop ( )
	    self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.QPARSED = True

    def ParseCharges ( self ):
        """Parse charges from a charge analysis."""
        if ( self.natoms > 0 ):
            nall = self.natoms + len ( self.dummies )
            q    = Real1DArray.WithExtent ( self.natoms )
            q.Set ( 0.0 )
            i    = -1
            n    = -1
            for iblock in range ( ( nall + 4 ) // 5 ):
                self.GetLine ( ) # Blank line.
                self.GetLine ( ) # Atom heading.
                line   = self.GetLine ( )
                tokens = line[7:].split ( )
                for token in tokens:
                    n += 1
                    if n not in self.dummies:
                        i += 1
                        q[i] = float ( token )
            return q
        else: return None

    def ParseConvergenceCriterium ( self, line ):
        """Parse a convergence criterium - a pair of numbers giving the actual value and the required value."""
        words = line.split ( ) ;
        return [ float ( words[2] ), float ( words[5] ) ]

    def ParseDipole ( self ):
        """Parse a dipole."""
        line = self.GetLine ( )
        if line.find ( "Dipole Moments" ) >= 0:
            tokens = self.GetTokens ( converters = [ None, float, None, float, None, float ] )
            self.dipole = Vector3.Null ( )
            self.dipole[0] = tokens[1]
            self.dipole[1] = tokens[3]
            self.dipole[2] = tokens[5]

    def ParseECPs ( self ):
        """Parse the ECP section."""
        # . Initialization.
        ecpelectrons = []
        # . Header.
        self.GetLine ( ) # . Blank line.
        self.GetLine ( ) # . Table header.
        # . Loop over atom lines.
        while True:
            tokens = self.GetTokens ( )
            # . Blank line.
            if len ( tokens )   == 0: break
            # . Atom line.
            elif len ( tokens ) == 2:
                atomicNumber = PeriodicTable.AtomicNumber ( tokens[0] )
                if ( atomicNumber > 0 ): ecpelectrons.append ( int ( tokens[1] ) )
        # . Finish up.
        self.ecpelectrons = ecpelectrons
        self.QECP         = True

    def ParseFunctions ( self ):
        """Parse a function list."""
        for i in range ( 7 ): self.GetLine ( ) # . All header lines.
        nfunctions = []
        tokens     = self.GetTokens ( converters = [ None, None, None, int ] )
        oldname    = tokens[0]
        nstart     = tokens[3]
        nstop      = nstart
        while True:
            tokens = self.GetTokens ( )
            if ( len ( tokens ) == 8 ):
                if tokens[0] != oldname:
                    nfunctions.append ( nstop - nstart + 1 )
                    nstart  = int ( tokens[3] )
                    nstop   = nstart
                    oldname = tokens[0]
                else:
                    nstop = int ( tokens[3] )
            elif ( len ( tokens ) == 4 ):
                nstop = int ( tokens[1] )
            else:
                nfunctions.append ( nstop - nstart + 1 )
                break
        self.nfunctions = nfunctions

    def ParseGeometry ( self ):
        """Parse a geometry."""
        self.GetLine ( ) # .The Angstroms line.
        self.GetLine ( ) # . The atoms heading.
        atomicNumbers = []
        xyz           = []
        dummies       = set ( )
        n             = 0
        while True:
            tokens = self.GetTokens ( )
            if ( len ( tokens ) == 4 ):
                atomicNumber = PeriodicTable.AtomicNumber ( tokens[0] )
                if ( atomicNumber > 0 ):
                    atomicNumbers.append ( atomicNumber )
                    x = float ( tokens[1] )
                    y = float ( tokens[2] )
                    z = float ( tokens[3] )
                    xyz.append ( ( x, y, z ) )
                else:
                    dummies.add ( n )
                n += 1
            else:
                break
        if ( len ( atomicNumbers ) > 0 ):
            self.atomicNumbers = atomicNumbers
            if self.coordinates3 is None: self.coordinates3 = Coordinates3.WithExtent ( len ( xyz ) )
            for ( i, ( x, y, z ) ) in enumerate ( xyz ):
                self.coordinates3[i,0] = x
                self.coordinates3[i,1] = y
                self.coordinates3[i,2] = z
            self.natoms = len ( self.atomicNumbers )
            self.ngeometries += 1
            self.dummies      = dummies

    def ParseOverlap ( self ):
        """Parse the overlap matrix."""
        # . The number of basis functions is unknown.
        data = [ ]
        # . Loop over the first block.
        self.GetLine ( ) # . Blank line.
        self.GetLine ( ) # . Integer heading.
        ncolumns = 0
        while True:
            tokens = self.GetTokens ( )
            if len ( tokens ) <= 0: break
            i = int ( tokens[0] )
            l = len ( tokens )
            for j in range ( 1, l ): data.append ( ( i, j, float ( tokens[j] ) ) )
        ncolumns = ncolumns + l - 1
        # . Loop over subsequent blocks.
        nbasis  = i
        nblocks = ( ( nbasis + 4 ) // 5 ) - 1
        if nblocks > 0:
            for iblock in range ( nblocks ):
                self.GetLine ( ) # . Integer heading.
                while True:
                    tokens = self.GetTokens ( )
                    if len ( tokens ) <= 0: break
                    i = int ( tokens[0] )
                    l = len ( tokens )
                    for j in range ( 1, l ): data.append ( ( i, ncolumns + j, float ( tokens[j] ) ) )
                ncolumns = ncolumns + l - 1
        nbastr = ( nbasis * ( nbasis + 1 ) ) // 2
        if len ( data ) == nbastr:
            self.nbasis  = nbasis
            self.overlap = SymmetricMatrix.WithExtent ( nbasis )
            self.overlap.Set ( 0.0 )
            for ( i, j, o ) in data:
                self.overlap[i-1,j-1] = o
        else:
            self.Warning ( "Invalid number of overlap matrix elements: {:d} and {:d}.".format ( len ( data ), nbastr ), False )

    def Summary ( self, log = logFile ):
        """Print a summary of the stored data."""
        if self.QPARSED and LogFileActive ( log ):
            ( head, tail ) = os.path.split ( self.name )
            summary = log.GetSummary ( )
            summary.Start ( "Jaguar Output File Summary for " + "\"" + tail + "\"" )
            summary.Entry ( "Number of Atoms"     , "{:d}".format ( self.natoms       ) )
            summary.Entry ( "Number of Geometries", "{:d}".format ( self.ngeometries  ) )
            summary.Entry ( "Net Charge"          , "{:d}".format ( self.charge       ) )
            summary.Entry ( "Multiplicity"        , "{:d}".format ( self.multiplicity ) )
            if ( self.displacementmaximum is not None ) and ( self.displacementrms is not None ) and ( self.energychange is not None ) and ( self.gradientmaximum is not None ) and ( self.gradientrms is not None ):
                summary.Entry ( "Optimization Finished", "{!r}"  .format ( self.QCOMPLETED             ) )
                summary.Entry ( "Optimization Stuck"   , "{!r}"  .format ( self.QFATALERROR            ) )
                summary.Entry ( "Maximum Displacement" , "{:.6e}".format ( self.displacementmaximum[0] ) )
                summary.Entry ( "Target"               , "{:.6e}".format ( self.displacementmaximum[1] ) )
                summary.Entry ( "RMS Displacement"     , "{:.6e}".format ( self.displacementrms[0]     ) )
                summary.Entry ( "Target"               , "{:.6e}".format ( self.displacementrms[1]     ) )
                summary.Entry ( "Energy Change"        , "{:.6e}".format ( self.energychange[0]        ) )
                summary.Entry ( "Target"               , "{:.6e}".format ( self.energychange[1]        ) )
                summary.Entry ( "Maximum Gradient"     , "{:.6e}".format ( self.gradientmaximum[0]     ) )
                summary.Entry ( "Target"               , "{:.6e}".format ( self.gradientmaximum[1]     ) )
                summary.Entry ( "RMS Gradient"         , "{:.6e}".format ( self.gradientrms[0]         ) )
                summary.Entry ( "Target"               , "{:.6e}".format ( self.gradientrms[1]         ) )
                summary.Entry ( "Best Energy"          , "{:.6e}".format ( self.bestenergy             ) )
                summary.Entry ( "Last Energy"          , "{:.6e}".format ( self.energy                 ) )
            elif self.scfenergy is not None:
                summary.Entry ( "SCF Energy"           , "{:.6e}".format ( self.scfenergy              ) )
            summary.Entry ( "Lines Parsed", "{:d}".format ( self.nlines ) )
            summary.Entry ( "Fatal Errors", "{:d}".format ( self.nfatal ) )
            summary.Entry ( "Warnings",     "{:d}".format ( self.nwarnings - self.nfatal ) )
            summary.Stop ( )

    def ToCoordinates ( self, QBEST = True ):
        """Return a coordinate set from the file."""
        coordinates3 = None
        if self.QPARSED:
            if QBEST and ( self.bestcoordinates3 is not None ): coordinates3 = self.bestcoordinates3
            else:                                               coordinates3 = self.coordinates3
        return coordinates3

    def ToSystem ( self, QBEST = True ):
        """Return a system constructed from the file data."""
        system = None
        if self.QPARSED and len ( self.atomicNumbers > 0 ):
            system              = System.FromAtoms ( self.atomicNumbers )
            system.coordinates3 = self.ToCoordinates3 ( QBEST = QBEST )
            system.label        = "System from Jaguar Output File"
        return system

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def JaguarOutputFile_ToCoordinates3 ( filename ):
    """Helper function that reads the coordinates from a Jaguar output file."""
    outfile = JaguarOutputFileReader ( filename )
    outfile.Parse ( )
    coordinates3 = outfile.ToCoordinates3 ( )
    return coordinates3

def JaguarOutputFile_ToSystem ( filename ):
    """Helper function that reads a system from a Jaguar outputfile."""
    outfile = JaguarOutputFileReader ( filename )
    outfile.Parse ( )
    system = outfile.ToSystem ( )
    return system

# . Importer definitions.
_Importer.AddHandler ( { Coordinates3 : JaguarOutputFile_ToCoordinates3 ,
                         System       : JaguarOutputFile_ToSystem       } , [ "jagout", "jout", "JAGOUT", "JOUT" ], "Jaguar Output", defaultFunction = JaguarOutputFile_ToSystem )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
