#-------------------------------------------------------------------------------
# . File      : pMolecule.QCParameters.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle the parameter data necessary for a quantum chemical calculation."""

import os, os.path

from pCore   import logFile, LogFileActive, RawObjectConstructor, YAMLMappingFile_ToObject, YAMLPickleFileExtension
from Element import PeriodicTable

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCParameters:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef QCParameters new
        new         = self.__class__.Raw ( )
        new.cObject = QCParameters_Clone ( self.cObject )
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        if self.isOwner:
            QCParameters_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.QCParameters"

    def __getstate__ ( self ):
        """Return the state."""
        cdef GaussianBasis  basis
        cdef MNDOParameters parameters
        centers = []
        for i from 0 <= i < self.cObject.ncenters:
            center = { "atomicNumber" : self.cObject.centers[i].atomicNumber }
            if self.cObject.centers[i].densitybasis != NULL:
                basis = GaussianBasis.Raw ( )
                basis.cObject = self.cObject.centers[i].densitybasis
                basis.isOwner = False
                center["densityBasis"] = basis
            if self.cObject.centers[i].orbitalbasis != NULL:
                basis = GaussianBasis.Raw ( )
                basis.cObject = self.cObject.centers[i].orbitalbasis
                basis.isOwner = False
                center["orbitalBasis"] = basis
            if self.cObject.centers[i].poissonbasis != NULL:
                basis = GaussianBasis.Raw ( )
                basis.cObject = self.cObject.centers[i].poissonbasis
                basis.isOwner = False
                center["poissonBasis"] = basis
            if self.cObject.centers[i].mndoparameters != NULL:
                parameters = MNDOParameters.Raw ( )
                parameters.cObject = self.cObject.centers[i].mndoparameters
                parameters.isOwner = False
                center["mndoParameters"] = parameters
            centers.append ( center )
        return { "centers" : centers }

    def __init__ ( self, numberOfCenters ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfCenters )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef GaussianBasis  basis
        cdef MNDOParameters parameters
        # . Allocate the object.
        centers = state["centers"]
        self._Allocate ( len ( centers ) )
        # . Fill the object.
        for ( i, center ) in enumerate ( centers ):
            self.cObject.centers[i].atomicNumber = center["atomicNumber"]
            if "densityBasis" in center:
                basis = center["densityBasis"]
                self.cObject.centers[i].densitybasis = basis.cObject
                basis.isOwner                        = False
            if "orbitalBasis" in center:
                basis = center["orbitalBasis"]
                self.cObject.centers[i].orbitalbasis = basis.cObject
                basis.isOwner                        = False
            if "poissonBasis" in center:
                basis = center["poissonBasis"]
                self.cObject.centers[i].poissonbasis = basis.cObject
                basis.isOwner                        = False
            if "mndoParameters" in center:
                parameters = center["mndoParameters"]
                self.cObject.centers[i].mndoparameters = parameters.cObject
                parameters.isOwner                     = False
            self.ScaleMNDOOrbitalBasisExponents ( i )

    def _Allocate ( self, numberOfCenters ):
        """Constructor."""
        self.cObject = QCParameters_Allocate ( numberOfCenters )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SummaryEntry ( self, object summary ):
        """Summary entry."""
        if summary is not None:
            summary.Entry ( "QC Parameter Centers", "{:d}".format ( self.cObject.ncenters ) )

    def ScaleMNDOOrbitalBasisExponents ( self, Integer icenter ):
        """Scale the exponents of a MNDO orbital basis."""
        if ( self.cObject.centers[icenter].mndoparameters != NULL ) and ( self.cObject.centers[icenter].orbitalbasis != NULL ):
            if ( self.cObject.centers[icenter].mndoparameters.norbitals > 0 ):
                GaussianBasis_ScaleExponents ( self.cObject.centers[icenter].orbitalbasis         , \
                                               self.cObject.centers[icenter].mndoparameters.zetas , \
                                               self.cObject.centers[icenter].mndoparameters.zetap , \
                                               self.cObject.centers[icenter].mndoparameters.zetad   )
            else:
                GaussianBasis_Deallocate ( &(self.cObject.centers[icenter].orbitalbasis) )

#===================================================================================================================================
# . Class methods.
#===================================================================================================================================
def QCParameters_Define ( atomicNumbers, densityPath = None, log = logFile, mndoPath = None, orbitalPath = None ):
    """Constructor given a list of atomic numbers and parameter directories."""
    # . Declarations.
    cdef GaussianBasis  basis
    cdef MNDOParameters mndo
    cdef QCParameters   self
    # . Allocate the object.
    self = QCParameters.Raw ( )
    self._Allocate ( len ( atomicNumbers ) )
    # . Extract data for each atomic number.
    missingpairs = set ( )
    warnings     = []
    for ( icenter, atomicNumber ) in enumerate ( atomicNumbers ):
        self.cObject.centers[icenter].atomicNumber = atomicNumber
        symbol            = PeriodicTable.Symbol ( atomicNumber )
        orbitalBasisLabel = symbol
        warning           = []
        # . MNDO parameters.
        if mndoPath is not None:
            filename = os.path.join ( mndoPath, symbol + YAMLPickleFileExtension )
            try:
                mndo = YAMLMappingFile_ToObject ( filename, MNDOParameters )
                self.cObject.centers[icenter].mndoparameters = mndo.cObject
                mndo.isOwner = False
                # . Create path for orbital basis.
                labels = []
                if mndo.cObject.qns != 0: labels.append ( "{:d}s".format ( mndo.cObject.qns ) )
                if mndo.cObject.qnp != 0: labels.append ( "{:d}p".format ( mndo.cObject.qnp ) )
                if mndo.cObject.qnd != 0: labels.append ( "{:d}d".format ( mndo.cObject.qnd ) )
                orbitalBasisLabel = "".join ( labels ).lower ( )
                # . Check for missing diatomic core-core pairs.
                if mndo.cObject.QDIATOMIC == CTrue:
                    for ( jcenter, other ) in enumerate ( atomicNumbers ):
                        if ( icenter == jcenter ) or ( ( other < mndo.cObject.ndiatomic ) and ( mndo.cObject.QDIATOMICFLAGS[other] == CTrue ) ):
                            pass
                        else:
                            if atomicNumber >= other: key = ( atomicNumber, other )
                            else:                     key = ( other, atomicNumber )
                            missingpairs.add ( key )
            except:
                warnings.append ( ( mndoPath, symbol ) )
        # . Gaussian basis sets.
        if densityPath is not None:
            filename = os.path.join ( densityPath, symbol + YAMLPickleFileExtension )
            try:
                basis = YAMLMappingFile_ToObject ( filename, GaussianBasis )
                self.cObject.centers[icenter].densitybasis = basis.cObject
                basis.cObject.type = GaussianBasisType_Coulomb
                basis.isOwner      = False
                basis.Normalize ( ) # . Fudge for renormalization.
            except:
                warnings.append ( ( densityPath, symbol ) )
        if orbitalPath is not None:
            filename = os.path.join ( orbitalPath, orbitalBasisLabel + YAMLPickleFileExtension )
            try:
                basis = YAMLMappingFile_ToObject ( filename, GaussianBasis )
                self.cObject.centers[icenter].orbitalbasis = basis.cObject
                basis.cObject.type = GaussianBasisType_Orbital
                if mndoPath is not None: basis.cObject.atomicNumber = atomicNumber
                basis.isOwner      = False
            except:
                warnings.append ( ( orbitalPath, orbitalBasisLabel ) )
        # . Scaling.
        if mndoPath is not None: self.ScaleMNDOOrbitalBasisExponents ( icenter )
    # . MNDO warnings.
    if ( mndoPath is not None ) and ( len ( missingpairs ) > 0 ):
         if LogFileActive ( log ):
             log.Paragraph ( "Possible MNDO missing diatomic core-core interaction parameters: {!r}.".format ( missingpairs ) )
    # . Finish up.
    if len ( warnings ) > 0: raise ValueError ( "Error processing some QC parameter files.", warnings )
    return self
