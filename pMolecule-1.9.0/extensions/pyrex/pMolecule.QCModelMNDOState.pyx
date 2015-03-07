#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModelMNDOState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines a state for a MNDO quantum chemical model."""

from pCore import CLibraryError, logFile, LogFileActive, UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE

import math

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModelMNDOState:
    """Define a state for a semi-empirical MNDO-type quantum chemical model."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            QCModelMNDOState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.QCModelMNDOState"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = QCModelMNDOState_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def __len__ ( self ):
        """Return the length as the number of QC atoms."""
        length = 0
        if self.cObject != NULL:
            length = QCAtomContainer_Size ( self.cObject.qcAtoms )
        return length

    def BackupFock ( self ):
        """Backup the Fock matrices."""
        cdef Status status
        status = Status_Continue
        QCModelMNDOState_BackupFock ( self.cObject, &status )
        if status != Status_Continue: raise CLibraryError ( "Error backing up Fock matrices." )

# #ifdef MNDOCI
    def CIVectorsTable ( self, log = logFile, nvectors = 20, vectorsPerRow = 8 ):
        """Output a table of CI vectors."""
        cdef Integer i, index, j
        cdef MNDOCIState *cistate
        if LogFileActive ( log ):
            # . Get the CI state.
            cistate = self.cObject.cistate
            if cistate == NULL:
                log.Paragraph ( "No MNDO CI information found." )
            else:
                # . Get the total number of configurations and the number to write.
                if nvectors is None: n = cistate.numberOfStates
                else:                n = min ( nvectors, cistate.numberOfStates )
                # . Write the vectors.
                if n > 0:
                    # . Find the number of tables to write.
                    ntables = n / vectorsPerRow
                    if ( n % vectorsPerRow > 0 ): ntables += 1
                    # . Write the tables.
                    for itable in range ( ntables ):
                        # . Number of vectors to write.
                        ncolumns = vectorsPerRow
                        if itable == ntables - 1: ncolumns = n - itable * vectorsPerRow
                        # . Find width of first column.
                        nwidth = max ( 10, 2 * ( cistate.nactive + 1 ) )
                        # . Start tables.
                        columns = [ 6, nwidth ] + ncolumns * [ 15 ]
                        table   = log.GetTable ( columns = columns )
                        table.Start ( )
                        table.Title ( "MOPAC CI Vectors" )
                        table.Heading ( "" )
                        table.Heading ( "" )
                        for i in range ( ncolumns ): table.Heading ( "{:d}".format ( i+1+itable*vectorsPerRow ) )
                        table.Entry ( "Energy", alignment = "l" )
                        table.Entry ( "" )
                        for i in range ( ncolumns ):
                            index = i+itable*vectorsPerRow
                            e     = Real1DArray_GetItem ( cistate.ciEnergies, index, NULL ) + cistate.baseline
                            table.Entry ( "{:.6f}".format ( e ) )
                        table.Entry ( "<S^2>", alignment = "l" )
                        table.Entry ( "" )
                        for i in range ( ncolumns ):
                            index = i+itable*vectorsPerRow
                            x     = Real1DArray_GetItem ( cistate.spins, index, NULL )
                            table.Entry ( "{:.6f}".format ( x ) )
                        table.Entry ( "S", alignment = "l" )
                        table.Entry ( "" )
                        for i in range ( ncolumns ):
                            index = i+itable*vectorsPerRow
                            x     = Real1DArray_GetItem ( cistate.spins, index, NULL )
                            y     = 0.5 * ( -1.0 + math.sqrt ( 1.0 + 4.0 * x ) )
                            table.Entry ( "{:.6f}".format ( y ) )
                        for i in range ( cistate.nconfigurations ):
                            alphas = []
                            betas  = []
                            for j from 0 <= j < cistate.nactive: alphas.append ( "{:d}".format ( Integer1DArray_GetItem ( cistate.configurations[i].alphas, j, NULL ) ) )
                            for j from 0 <= j < cistate.nactive: betas.append  ( "{:d}".format ( Integer1DArray_GetItem ( cistate.configurations[i].betas , j, NULL ) ) )
                            table.Entry ( "{:d}".format ( i+1 ) )
                            table.Entry ( "".join ( alphas ) + " " + "".join ( betas ) )
                            for j in range ( ncolumns ):
                                index = j+itable*vectorsPerRow
                                table.Entry ( "{:.6f}".format ( Real2DArray_GetItem ( cistate.ciVectors, index, i, NULL ) ) )
                        table.Stop ( )

    def CIWavefunctionSummary ( self, log = logFile, ncoefficients = None ):
        """Write a summary of the CI wavefunction."""
        cdef Integer i
        cdef MNDOCIState *cistate
        if LogFileActive ( log ):
            # . Get the CI state.
            cistate = self.cObject.cistate
            if cistate == NULL:
                log.Paragraph ( "No MNDO CI information found." )
            else:
                # . Get the total number of configurations and the number to write.
                if ncoefficients is None: n = cistate.nconfigurations
                else:                     n = min ( ncoefficients, cistate.nconfigurations )
                # . Sort the configurations.
                sorted = []
                for i from 0 <= i < cistate.nconfigurations: sorted.append ( ( math.fabs ( Real1DArray_GetItem ( cistate.ciVector, i, NULL ) ), i ) )
                sorted.sort ( reverse = True )
                # . Write the table.
                doneFundamental = False
                columns = [ 7, 7, 14 ] + 2 * cistate.nactive * [ 3 ]
                table   = log.GetTable ( columns = columns )
                table.Start ( )
                table.Title ( "MOPAC CI Wavefunction Coefficients" )
                table.Heading ( "Order"       )
                table.Heading ( "Conf."       )
                table.Heading ( "Coefficient" )
                table.Heading ( "Alphas", columnSpan = cistate.nactive )
                table.Heading ( "Betas" , columnSpan = cistate.nactive )
                for ( i, ( c, index ) ) in enumerate ( sorted ):
                    table.Entry   ( "{:d}".format ( i )     )
                    table.Entry   ( "{:d}".format ( index ) )
                    table.Entry   ( "{:12.6f}".format ( Real1DArray_GetItem ( cistate.ciVector, index, NULL ) ) )
                    for i from 0 <= i < cistate.nactive: table.Entry ( "{:d}".format ( Integer1DArray_GetItem ( cistate.configurations[index].alphas, i, NULL ) ) )
                    for i from 0 <= i < cistate.nactive: table.Entry ( "{:d}".format ( Integer1DArray_GetItem ( cistate.configurations[index].betas , i, NULL ) ) )
                    if index == 0: doneFundamental = True
                if not doneFundamental:
                    table.Heading ( "HF Ground State", columnSpan = len ( columns ) )
                    table.Entry   ( "{:d}".format ( n ) )
                    table.Entry   ( "{:d}".format ( 0 ) )
                    table.Entry   ( "{:12.6f}".format ( Real1DArray_GetItem ( cistate.ciVector, 0, NULL ) ) )
                    for i from 0 <= i < cistate.nactive: table.Entry ( "{:d}".format ( Integer1DArray_GetItem ( cistate.configurations[0].alphas, i, NULL ) ) )
                    for i from 0 <= i < cistate.nactive: table.Entry ( "{:d}".format ( Integer1DArray_GetItem ( cistate.configurations[0].betas , i, NULL ) ) )
                table.Stop ( )
# #endif /*MNDOCI*/

    def GetEnergies ( self ):
        """Return the final energies."""
        energies = []
        if self.cObject != NULL:
# #ifdef MNDOCI
            if self.cObject.cistate != NULL:
                energies.append ( ( "Pure QC", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * ( self.cObject.cistate.ciEnergy                    + self.cObject.qcAtoms.energybaseline ) ) )
            else:
                energies.append ( ( "Pure QC", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * ( self.cObject.eelectronic + self.cObject.enuclear + self.cObject.qcAtoms.energybaseline ) ) )
                if self.cObject.qcmmstate != NULL:
                    energies.append ( ( "QC/MM Elect.", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * self.cObject.eqcmm ) )
                    energies.append ( ( "QC/QC Elect.", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * self.cObject.eqcqc ) )
# #else /*MNDOCI*/
#            energies.append ( ( "Pure QC", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * ( self.cObject.eelectronic + self.cObject.enuclear + self.cObject.qcAtoms.energybaseline ) ) )
#            if self.cObject.qcmmstate != NULL:
#                energies.append ( ( "QC/MM Elect.", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * self.cObject.eqcmm ) )
#                energies.append ( ( "QC/QC Elect.", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * self.cObject.eqcqc ) )
# #endif /*MNDOCI*/
        return energies

    def GetItem ( self, name ):
        """Return a scalar value."""
        if   name == "Electronic Energy" : return self.cObject.eelectronic
        elif name == "Nuclear Energy"    : return self.cObject.enuclear
        elif name == "SCF Cycles"        : return self.cObject.ncycles
        else                             : return None

    def MakeDensitiesFromFock ( self ):
        """Make the densities from the Fock matrices."""
        QCModelMNDOState_MakeDensities ( self.cObject )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def RestoreFock ( self ):
        """Restore the Fock matrices."""
        QCModelMNDOState_RestoreFock ( self.cObject )

    def SetConvergenceData ( self, isConverged, ncycles ):
        """Set convergence data for the state."""
        if self.cObject != NULL:
            if isConverged: self.cObject.isConverged = CTrue
            else:           self.cObject.isConverged = CFalse
            self.cObject.ncycles = ncycles

    def SpinExpectationValues ( self ):
        """Calculate <Sz> and <S2>."""
        cdef Real Sz, S2
        Sz = 0.0
        S2 = 0.0
        if self.cObject.densityq != NULL:
            QCOnePDM_SpinExpectationValues ( self.cObject.densityp, self.cObject.densityq, NULL, &Sz, &S2 )
        return ( Sz, S2 )

    def Summary ( self, log = logFile, nabove = 10 ):
        """Summary."""
        cdef Integer i
# #ifdef MNDOCI
        cdef MNDOCIState *cistate
# #endif /*MNDOCI*/
        if LogFileActive ( log ):
            # . SCF.
            summary = log.GetSummary ( )
            summary = log.GetSummary ( valueWidth = 14 )
            summary.Start ( "MNDO QC Model State Summary" )
            summary.Entry ( "Electronic Energy",   "{:14.8g}".format ( self.cObject.eelectronic + self.cObject.eqcmm + self.cObject.eqcqc ) )
            summary.Entry ( "Nuclear Energy",      "{:14.8g}".format ( self.cObject.enuclear ) )
            summary.Entry ( "One-Electron Energy", "{:14.8g}".format ( self.cObject.eoei     ) )
            summary.Entry ( "Two-Electron Energy", "{:14.8g}".format ( self.cObject.etei     ) )
            if self.cObject.qcmmstate != NULL:
                summary.Entry ( "QC/MM Elect.", "{:14.8g}".format ( self.cObject.eqcmm ) )
                summary.Entry ( "QC/QC Elect.", "{:14.8g}".format ( self.cObject.eqcqc ) )
            if self.cObject.densityp.occupancyType == QCOnePDMOccupancyType_FractionalVariable:
                summary.Entry ( "Occupancy Energy", "{:14.8g}".format ( self.cObject.eocc ) )
                summary.Entry ( None, None )
                if self.cObject.densityq == NULL:
                    summary.Entry ( "Occupied Orbitals",  "{:14.8g}".format ( self.cObject.densityp.numberOccupied ) )
                    summary.Entry ( "Fermi Energy",       "{:14.8g}".format ( self.cObject.densityp.fermiEnergy    ) )
                else:
                    summary.Entry ( "Alpha Occupied",     "{:14.8g}".format ( self.cObject.densityp.numberOccupied ) )
                    summary.Entry ( "Alpha Fermi Energy", "{:14.8g}".format ( self.cObject.densityp.fermiEnergy    ) )
                    summary.Entry ( "Beta  Occupied",     "{:14.8g}".format ( self.cObject.densityq.numberOccupied ) )
                    summary.Entry ( "Beta  Fermi Energy", "{:14.8g}".format ( self.cObject.densityq.fermiEnergy    ) )
            if self.cObject.densityq != NULL:
                ( Sz, S2 ) = self.SpinExpectationValues ( )
                summary.Entry ( "<Sz>", "{:14.3f}".format ( Sz ) )
                summary.Entry ( "<S2>", "{:14.3f}".format ( S2 ) )
            if self.cObject.isConverged == CTrue: summary.Entry ( "SCF Convergence", "True"  )
            else:                                 summary.Entry ( "SCF Convergence", "False" )
            summary.Entry ( "SCF Cycles", "{:d}".format ( self.cObject.ncycles ) )
            summary.Stop ( )
# #ifdef MNDOCI
            # . CI.
            cistate = self.cObject.cistate
            if cistate != NULL:
                # . Basic summary.
                summary = log.GetSummary ( )
                summary.Start ( "MOPAC CI Calculation Summary" )
                summary.Entry ( "Configurations"     , "{:d}"   .format ( cistate.nconfigurations          ) )
                summary.Entry ( "Active Electrons"   , "{:d}"   .format ( cistate.nelectrons               ) )
                summary.Entry ( "Core Orbitals"      , "{:d}"   .format ( cistate.ncore                    ) )
                summary.Entry ( "Active Orbitals"    , "{:d}"   .format ( cistate.nactive                  ) )
                summary.Entry ( "CI Matrix Sparsity" , "{:5.1%}".format ( cistate.ciMatrixSparsity / 100.0 ) )
                summary.Entry ( "Number Of States"   , "{:d}"   .format ( cistate.numberOfStates           ) )
                # . Solver data.
                if cistate.doSparse == CTrue:
                    summary.Entry ( "Solver Converged" , "{:d}".format ( cistate.eigenvalueSolverReport.isConverged == CTrue              ) )
                    summary.Entry ( "Converged Pairs"  , "{:d}".format ( cistate.eigenvalueSolverReport.convergedPairs                    ) )
                    summary.Entry ( "Solver Iterations", "{:d}".format ( cistate.eigenvalueSolverReport.numberMatrixVectorMultiplications ) )
                    summary.Entry ( "Return Code"      , "{:d}".format ( cistate.eigenvalueSolverReport.returnCode                        ) )
                    if cistate.eigenvalueSolverReport.solutionChecked == CTrue:
                        summary.Entry ( "Solution Checked"         , "{!r}"    .format ( cistate.eigenvalueSolverReport.solutionChecked == CTrue ) )
                        summary.Entry ( "Eigenvalue Error"         , "{:14.6g}".format ( cistate.eigenvalueSolverReport.eigenvalueError          ) )
                        summary.Entry ( "Eigenvector Error"        , "{:14.6g}".format ( cistate.eigenvalueSolverReport.eigenvectorError         ) )
                        summary.Entry ( "Normalization Error"      , "{:14.6g}".format ( cistate.eigenvalueSolverReport.normalizationError       ) )
                # . Localization options.
                if cistate.localizeStart > -1:
                    summary.Entry ( "Localize Start", "{:d}".format ( cistate.localizeStart ) )
                    summary.Entry ( "Localize Stop" , "{:d}".format ( cistate.localizeStop  ) )
                # . Gradient data.
                if cistate.doGradients == CTrue:
                    summary.Entry ( "Non-Redundant"       , "{:d}".format ( cistate.numberNonRedundant        ) )
                    summary.Entry ( "Redundant"           , "{:d}".format ( cistate.numberRedundant           ) )
                    summary.Entry ( "Degenerate Redundant", "{:d}".format ( cistate.numberDegenerateRedundant ) )
                    summary.Entry ( "CPHF Iterations"     , "{:d}".format ( cistate.cphfIterations            ) )
                    summary.Entry ( "CPHF Status Flag"    , "{:d}".format ( cistate.cphfErrorFlag             ) )
                summary.Stop ( )
	        # . Details of configurations.
                spinlabels = { 1 : "Singlet", 2 : "Doublet", 3 : "Triplet", 4 : "Quartet", 5 : "Quintet", 6 : "Sextet", 7 : "Septet", 8 : "Octet", 9 : "Nonet", 10 : "Decuplet", 11 : "Undecuplet" }
                table      = log.GetTable ( columns = [ 8 ] + 4 * [ 18 ] )
                table.Start ( )
                table.Title ( "CI Configurations" )
                table.Heading ( "Conf."  )
                table.Heading ( "Energy" )
                table.Heading ( "Type"   )
                table.Heading ( "<S^2>"  )
                table.Heading ( "S"      )
                for i in range ( min ( cistate.numberOfStates, cistate.root + nabove ) ):
                    e     = Real1DArray_GetItem ( cistate.ciEnergies, i, NULL ) + cistate.baseline
                    x     = Real1DArray_GetItem ( cistate.spins, i, NULL )
                    y     = 0.5 * ( -1.0 + math.sqrt ( 1.0 + 4.0 * x ) )
                    j     = int ( round ( 2.0 * y + 1.0 ) )
                    label = spinlabels.get ( j, "{:d}".format ( j ) )
                    table.Entry ( "{:d}".format ( i ) )
                    table.Entry ( "{:.6f}".format ( e ) )
                    table.Entry ( label )
                    table.Entry ( "{:.5f}".format ( x ) )
                    table.Entry ( "{:.5f}".format ( y ) )
                table.Stop ( )
                # . Degenerate eigenvalues?
                if cistate.orbitalDegeneracies == CTrue: log.Paragraph ( "** There were \"quasi\"-degenerate orbitals in the CI calculation **" )
                # . Fractionally-occupied inactive orbitals?
                if cistate.fractionallyOccupiedInactive == CTrue: log.Paragraph ( "** There were fractionally-occupied inactive orbitals (either core or virtual) in the CI calculation **" )
                # . Possible gradient error.
                if ( cistate.ciGradientError == CTrue ) or ( cistate.cphfErrorFlag != 0 ): log.Paragraph ( "** Possible gradient error (inappropriate model or CPHF problem). **" )
                # . Root not found.
                if cistate.rootNotFound == CTrue: log.Paragraph ( "** It was not possible to locate the appropriate root. **" )
                # . Final energy.
                else: log.Paragraph ( "Selected State CI Energy = {:.8g} Hartrees.".format ( cistate.baseline + cistate.rootenergy ) )
# #endif /*MNDOCI*/

    # . Properties.
# #ifdef MNDOCI
    property ciStateEnergies:
        def __get__ ( self ):
            cdef MNDOCIState *cistate
            cdef Real1DArray  item
            item = None
            cistate = self.cObject.cistate
            if ( cistate != NULL ) and ( cistate.ciEnergies != NULL ):
                item = Real1DArray.Raw ( )
                item.cObject = cistate.ciEnergies
                item.isOwner = False
                item.owner   = self
            return item
    property ciStateSpins:
        def __get__ ( self ):
            cdef MNDOCIState *cistate
            cdef Real1DArray  item
            item = None
            cistate = self.cObject.cistate
            if ( cistate != NULL ) and ( cistate.spins != NULL ):
                item = Real1DArray.Raw ( )
                item.cObject = cistate.spins
                item.isOwner = False
                item.owner   = self
            return item
# #endif /*MNDOCI*/
    property densityp:
        def __get__ ( self ):
            cdef QCOnePDM item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.densityp != NULL ):
                item = QCOnePDM.Raw ( )
                item.cObject = self.cObject.densityp
                item.isOwner = False
                item.owner   = self
            return item
    property densityq:
        def __get__ ( self ):
            cdef QCOnePDM item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.densityq != NULL ):
                item = QCOnePDM.Raw ( )
                item.cObject = self.cObject.densityq
                item.isOwner = False
                item.owner   = self
            return item
# #ifdef MNDOCI
    property numberCIConfigurations:
        def __get__ ( self ):
            cdef MNDOCIState *cistate
            item = 0
            cistate = self.cObject.cistate
            if cistate != NULL:
                item = cistate.nconfigurations
            return item
# #endif /*MNDOCI*/
    property oneElectronMatrix:
        def __get__ ( self ):
            cdef SymmetricMatrix item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.oneelectronmatrix != NULL ):
                item = SymmetricMatrix.Raw ( )
                item.cObject = self.cObject.oneelectronmatrix
                item.isOwner = False
            return item
    property overlap:
        def __get__ ( self ): return None
    property orthogonalizer:
        def __get__ ( self ): return None
