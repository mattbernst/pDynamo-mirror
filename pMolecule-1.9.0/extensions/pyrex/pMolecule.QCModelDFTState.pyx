#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModelDFTState.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines a state for a DFT quantum chemical model."""

from pCore import logFile, LogFileActive, UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModelDFTState:
    """Define a state for a DFT quantum chemical model."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            QCModelDFTState_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.QCModelDFTState"

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = QCModelDFTState_Allocate ( )
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

    def GetEnergies ( self ):
        """Return the final energies."""
        energies = []
        if self.cObject != NULL:
            energies.append ( ( "Pure QC", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * ( ( self.cObject.eelectronic + self.cObject.enuclear ) + self.cObject.qcAtoms.energybaseline ) ) )
            if self.cObject.qcmmstate != NULL:
                energies.append ( ( "QC/MM Elect.", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * self.cObject.eqcmm ) )
                energies.append ( ( "QC/QC Elect.", UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE * self.cObject.eqcqc ) )
        return energies

    def GetItem ( self, name ):
        """Return a scalar value."""
        if name == "SCF Cycles":
            return self.cObject.ncycles
        else:
            return None

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetConvergenceData ( self, isConverged, ncycles ):
        """Set convergence data for the state."""
        if self.cObject != NULL:
            if isConverged: self.cObject.isConverged = CTrue
            else:          self.cObject.isConverged = CFalse
            self.cObject.ncycles = ncycles

    def Summary ( self, log = logFile ):
        """Summary."""
        cdef Real Sz, S2
        if LogFileActive ( log ):
            # . Intermediate quantities.
            nf = float ( self.cObject.qcAtoms.nfbasisw )
            no = float ( self.cObject.qcAtoms.nobasisw )
            # . Summary.
            summary = log.GetSummary ( )
            summary = log.GetSummary ( valueWidth = 14 )
            summary.Start ( "DFT QC Model State Summary" )
            summary.Entry ( "Electronic Energy"   , "{:14.8g}".format ( self.cObject.eelectronic + self.cObject.eqcmm + self.cObject.eqcqc ) )
            summary.Entry ( "Nuclear Energy"      , "{:14.8g}".format ( self.cObject.enuclear ) )
            summary.Entry ( "One-Electron Energy" , "{:14.8g}".format ( self.cObject.eoei     ) )
            summary.Entry ( "Two-Electron Energy" , "{:14.8g}".format ( self.cObject.etei     ) )
            summary.Entry ( "Quadrature Energy"   , "{:14.8g}".format ( self.cObject.equad    ) )
            summary.Entry ( "Quadrature Density"  , "{:14.8g}".format ( self.cObject.rhoquad  ) )
            summary.Entry ( "Quadrature Points"   , "{:d}".format ( self.cObject.dftgrid.numberOfPoints ) )
            if DFTGrid_HasFunctionData ( self.cObject.dftgrid ):
                n = DFTGrid_NumberOfFunctionValues ( self.cObject.dftgrid )
                summary.Entry ( "Quadrature BF Values"   , "{:d}"       .format ( n ) )
                summary.Entry ( "Quadrature BF Sparsity" , "{:10.2%}"   .format ( ( 1.0 - float ( n ) / ( float ( self.cObject.dftgrid.numberOfPoints ) * no ) ) ) )
                summary.Entry ( "Quadrature BF Storage"  , "{:10.4f} GB".format ( DFTGrid_FunctionDataSize ( self.cObject.dftgrid ) ) )
            summary.Entry ( "Fit Integrals"         , "{:d}"       .format ( self.cObject.fitintegrals.ndata ) )
            summary.Entry ( "Fit Integral Sparsity" , "{:10.2%}"   .format ( ( 1.0 - float ( self.cObject.fitintegrals.ndata ) / ( nf * ( no * ( no + 1.0 ) ) / 2.0 ) ) ) )
            summary.Entry ( "Fit Integral Storage"  , "{:10.4f} GB".format ( BlockStorage_Size ( self.cObject.fitintegrals ) ) )
            if self.cObject.qcmmstate != NULL:
                summary.Entry ( "QC/MM Elect." , "{:14.8g}".format ( self.cObject.eqcmm ) )
                summary.Entry ( "QC/QC Elect." , "{:14.8g}".format ( self.cObject.eqcqc ) )
            if self.cObject.densityp.occupancyType == QCOnePDMOccupancyType_FractionalVariable:
                summary.Entry ( "Occupancy Energy", "{:14.8g}".format ( self.cObject.eocc, ) )
                summary.Entry ( None, None )
                if self.cObject.densityq == NULL:
                    summary.Entry ( "Occupied Orbitals"  , "{:14.8g}".format ( self.cObject.densityp.numberOccupied ) )
                    summary.Entry ( "Fermi Energy"       , "{:14.8g}".format ( self.cObject.densityp.fermiEnergy    ) )
                else:
                    summary.Entry ( "Alpha Occupied"     , "{:14.8g}".format ( self.cObject.densityp.numberOccupied ) )
                    summary.Entry ( "Alpha Fermi Energy" , "{:14.8g}".format ( self.cObject.densityp.fermiEnergy    ) )
                    summary.Entry ( "Beta  Occupied"     , "{:14.8g}".format ( self.cObject.densityq.numberOccupied ) )
                    summary.Entry ( "Beta  Fermi Energy" , "{:14.8g}".format ( self.cObject.densityq.fermiEnergy    ) )
            if self.cObject.densityq != NULL:
                QCOnePDM_SpinExpectationValues ( self.cObject.densityp, self.cObject.densityq, self.cObject.overlap, &Sz, &S2 )
                summary.Entry ( "<Sz>" , "{:14.3f}".format ( Sz ) )
                summary.Entry ( "<S2>" , "{:14.3f}".format ( S2 ) )
            if self.cObject.isConverged == CTrue: summary.Entry ( "SCF Convergence" , "True"  )
            else:                                 summary.Entry ( "SCF Convergence" , "False" )
            summary.Entry ( "SCF Cycles", "{:d}".format ( self.cObject.ncycles ) )
            n = self.cObject.qcAtoms.nobasis - Real2DArray_Length ( self.cObject.orthogonalizer, 1 )
            if n != 0: summary.Entry ( "Basis Linear Dependence" , "{:d}".format ( n ) )
            summary.Stop ( )

    # . Properties.
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
    property overlap:
        def __get__ ( self ):
            cdef SymmetricMatrix item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.overlap != NULL ):
                item = SymmetricMatrix.Raw ( )
                item.cObject = self.cObject.overlap
                item.isOwner = False
            return item
    property orthogonalizer:
        def __get__ ( self ):
            cdef Real2DArray item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.orthogonalizer != NULL ):
                item = Real2DArray.Raw ( )
                item.cObject = self.cObject.orthogonalizer
                item.isOwner = False
                item.owner   = self
            return item
