#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModelDFT.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines a DFT quantum chemical model."""

import os

from pCore                       import CLibraryError, Clone, logFile, LogFileActive, UNITS_LENGTH_ANGSTROMS_TO_BOHRS
from DFTFunctional               import DFTFunctionalIDsFromKeyword
from QCChargeModelOptions        import QCChargeModel_ToEnum, QCChargeModel_ToString
from QCOnePDM                    import QCOnePDMOccupancyType_ToEnum, QCOnePDMOccupancyType_ToString
from QCParameters                import QCParameters_Define
from QCMMLinkAtomCouplingOptions import QCMMLinkAtomCoupling_ToEnum, QCMMLinkAtomCoupling_ToString

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Default values.
_DefaultDensityBasis  = "demon"
_DefaultFunctional    = "lda"
_DefaultOrbitalBasis  = "321g"

# . Accuracies - these data must correspond to the enums in the C source.
DFTGridAccuracy_ToEnum   = { "Very Low" : DFTGridAccuracy_VeryLow, "Low" : DFTGridAccuracy_Low, "Medium" : DFTGridAccuracy_Medium, "High" : DFTGridAccuracy_High, "Very High" : DFTGridAccuracy_VeryHigh }
DFTGridAccuracy_ToString = { DFTGridAccuracy_VeryLow : "Very Low", DFTGridAccuracy_Low : "Low", DFTGridAccuracy_Medium : "Medium", DFTGridAccuracy_High : "High", DFTGridAccuracy_VeryHigh : "Very High" }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModelDFT ( QCModel ):
    """Define a DFT quantum chemical model."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            QCModelDFT_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.QCModelDFT"

    def __getstate__ ( self ):
        """Return the state."""
        return { "converger"                 : self.converger                         ,
                 "label"                     : self.label                             ,
                 "densityBasis"              : self.densityBasis                      ,
                 "orbitalBasis"              : self.orbitalBasis                      ,
                 "inCore"                    : self.cObject.inCore                    ,
                 "accuracy"                  : DFTGridAccuracy_ToString[self.cObject.accuracy],
                 "functional"                : self.functional                        ,
                 "fermiBroadening"           : self.cObject.fermiBroadening           ,
                 "keepOrbitalData"           : self.cObject.keepOrbitalData           ,
                 "isSpinRestricted"          : self.cObject.isSpinRestricted          ,
                 "linkAtomRatio"             : self.cObject.linkAtomRatio             ,
                 "numberFractionalAlphaHOOs" : self.cObject.numberFractionalAlphaHOOs ,
                 "numberFractionalAlphaLUOs" : self.cObject.numberFractionalAlphaLUOs ,
                 "numberFractionalBetaHOOs"  : self.cObject.numberFractionalBetaHOOs  ,
                 "numberFractionalBetaLUOs"  : self.cObject.numberFractionalBetaLUOs  ,
                 "numberFractionalHOOs"      : self.cObject.numberFractionalHOOs      ,
                 "numberFractionalLUOs"      : self.cObject.numberFractionalLUOs      ,
                 "occupancyType"             : QCOnePDMOccupancyType_ToString[self.cObject.occupancyType] ,
                 "qcChargeModel"             : QCChargeModel_ToString[self.cObject.qcChargeModel]         }

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = QCModelDFT_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        super ( QCModelDFT, self )._Initialize ( )
        self.cObject      = NULL
        self.isOwner      = False
        self.densityBasis = _DefaultDensityBasis
        self.functional   = _DefaultFunctional
        self.orbitalBasis = _DefaultOrbitalBasis

    def AtomicCharges ( self, configuration, chargeModel = None, spinDensities = False ):
        """Atomic charges."""
        cdef Boolean            cSpinDensities
        cdef QCChargeModel   model
        cdef QCModelDFTState qcState
        cdef Real1DArray     qcCharges
        # . Initialization.
        qcCharges = None
        # . Check the model.
        if chargeModel is None:
            model = self.cObject.qcChargeModel
        else:
            try:    model = QCChargeModel_ToEnum[chargeModel]
            except: raise TypeError ( "Unknown charge model - {:s}.".format ( chargeModel ) )
        # . Get the state object.
        if hasattr ( configuration, "qcState" ):
            qcState = configuration.qcState
            if spinDensities: cSpinDensities = CTrue
            else:             cSpinDensities = CFalse
            qcCharges = Real1DArray.WithExtent ( len ( qcState ) )
            qcCharges.Set ( 0.0 )
            QCModelDFT_AtomicCharges ( self.cObject, qcState.cObject, &model, cSpinDensities, qcCharges.cObject )
        return qcCharges

    def BondOrders ( self, configuration, chargeModel = None ):
        """Bond orders."""
        cdef QCChargeModel   model
        cdef QCModelDFTState qcState
        cdef Real1DArray     charges, freevalence, totalvalence
        cdef SymmetricMatrix bondOrders
        # . Initialization.
        bondOrders   = None
        charges      = None
        freevalence  = None
        totalvalence = None
        # . Check the model.
        if chargeModel is None:
            model = self.cObject.qcChargeModel
        else:
            try:    model = QCChargeModel_ToEnum[chargeModel]
            except: raise TypeError ( "Unknown charge model - {:s}.".format ( chargeModel ) )
        # . Get the state object.
        if hasattr ( configuration, "qcState" ):
            qcState      = configuration.qcState
            bondOrders   = SymmetricMatrix.WithExtent ( len ( qcState ) ) ; bondOrders.Set   ( 0.0 )
            charges      = Real1DArray.WithExtent     ( len ( qcState ) ) ; charges.Set      ( 0.0 )
            freevalence  = Real1DArray.WithExtent     ( len ( qcState ) ) ; freevalence.Set  ( 0.0 )
            totalvalence = Real1DArray.WithExtent     ( len ( qcState ) ) ; totalvalence.Set ( 0.0 )
            QCModelDFT_MayerBondOrders ( self.cObject, qcState.cObject, &model, bondOrders.cObject, charges.cObject, freevalence.cObject, totalvalence.cObject )
        return ( bondOrders, charges, freevalence, totalvalence )

    def CheckOptions ( self ):
        """Check the options."""
        # . Ensure that a functional model exists.
        if ( self.cObject != NULL ) and ( self.cObject.functionalModel == NULL ):
            self.DefineFunctionalModel ( self.functional )

    def DefineFunctionalModel ( self, keyword ):
        """Define the functional model given a keyword."""
        cdef Integer1DArray ids
        cdef Status         status
        ids     = DFTFunctionalIDsFromKeyword ( keyword )
        DFTFunctionalModel_Deallocate ( &(self.cObject.functionalModel) )
        status  = Status_Continue
        self.cObject.functionalModel = DFTFunctionalModel_MakeFromIDs ( ids.cObject, self.cObject.isSpinRestricted, &status )
        if status != Status_Continue: raise CLibraryError ( "Error creating DFT functional model." )
        self.functional = keyword
        self.label      = None

    def DefineParameters ( self, anlist, log = logFile ):
        """Define parameters."""
        path      = os.getenv ( "PDYNAMO_PARAMETERS" )
        if path is None: raise ValueError ( "Parameter path environmental variable missing." )
        basisPath = os.path.join ( path, "gaussianBasisSets" )
        return QCParameters_Define ( anlist, densityPath = os.path.join ( basisPath, self.densityBasis ),
                                             orbitalPath = os.path.join ( basisPath, self.orbitalBasis ),
                                             log = log )

    def DipoleMoment ( self, configuration, Vector3 center = None ):
        """Dipole Moment."""
        cdef Coordinates3    coordinates3
        cdef QCModelDFTState qcState
        cdef Vector3         dipole, pcenter
        cdef CVector3       *ccenter
        # . Initialization.
        dipole = None
        # . Get the coordinate and state objects.
        if hasattr ( configuration, "coordinates3" ): coordinates3 = configuration.coordinates3
        else:                                         coordinates3 = None
        if hasattr ( configuration, "qcState"      ): qcState      = configuration.qcState
        else:                                         qcState      = None
        if ( coordinates3 is not None ) and ( qcState is not None ):
            if center is None:
                ccenter = NULL
            else:
                pcenter = Clone ( center )
                pcenter.Scale ( UNITS_LENGTH_ANGSTROMS_TO_BOHRS )
                ccenter = pcenter.cObject
            dipole         = Vector3.Raw ( )
            dipole.cObject = QCModelDFT_DipoleMoment ( self.cObject, qcState.cObject, coordinates3.cObject, ccenter )
            dipole.isOwner = True
        return dipole

    def EnergyFinalize ( self, configuration ):
        """Finalization."""
        cdef QCModelDFTState qcState
        if hasattr ( configuration, "qcState" ):
            qcState = configuration.qcState
            QCModelDFTState_Finalize ( qcState.cObject, self.cObject.keepOrbitalData )

    def EnergyFock ( self, configuration ):
        """Fock matrix construction."""
        cdef Real          eelectronic
        cdef QCModelDFTState qcState
        eelectronic = 0.0
        if hasattr ( configuration, "qcState" ):
            qcState = configuration.qcState
            QCModelDFT_Fock     ( self.cObject, qcState.cObject, &eelectronic )
            QCModelDFT_QCMMFock ( self.cObject, qcState.cObject, &eelectronic )
        return eelectronic

    def EnergyGradients ( self, configuration ):
        """Gradients."""
        cdef QCModelDFTState qcState
        if hasattr ( configuration, "qcState" ):
            qcState = configuration.qcState
            QCModelDFT_Gradients ( self.cObject, qcState.cObject )

    def EnergyInitialize ( self, configuration ):
        """Integrals and space allocation."""
        cdef Coordinates3           coordinates3, pgradients3
        cdef QCMMInteractionState   pqcmmstate
        cdef QCModelDFTState        qcState
        cdef CCoordinates3         *cgradients3
        cdef CQCMMInteractionState *cqcmmstate
        # . Get the coordinate and state objects.
        if hasattr ( configuration, "coordinates3" ): coordinates3 = configuration.coordinates3
        else:                                         coordinates3 = None
        if hasattr ( configuration, "qcState"      ): qcState      = configuration.qcState
        else:                                         qcState      = None
        if ( coordinates3 is not None ) and ( qcState is not None ):
            # . Check for gradients.
            if hasattr ( configuration, "gradients3" ): pgradients3 = configuration.gradients3
            else:                                       pgradients3 = None
            if pgradients3 is None: cgradients3 = NULL
            else:                   cgradients3 = pgradients3.cObject
            # . Check for QC/MM interactions.
            if hasattr ( configuration, "qcmmstate" ): pqcmmstate = configuration.qcmmstate
            else:                                      pqcmmstate = None
            if pqcmmstate is None: cqcmmstate = NULL
            else:                  cqcmmstate = pqcmmstate.cObject
            # . Initialize the state.
            QCModelDFTState_Initialize ( qcState.cObject, coordinates3.cObject, cqcmmstate, cgradients3 )
            # . DFT grid.
            qcState.cObject.dftgrid = DFTGrid_Construct ( self.cObject.accuracy, qcState.cObject.qcAtoms, qcState.cObject.qccoordinates3 )
            # . Calculate the integrals and nuclear energy.
            QCModelDFT_Integrals ( self.cObject, qcState.cObject )
            # . Orthogonalizer.
#            qcState.cObject.orthogonalizer = SymmetricMatrix_Orthogonalizing_Transformation ( qcState.cObject.overlap )
# . Reorthogonalize?
#        QCOnePDM_Reorthogonalize       ( pdensityp.cObject, porthogonalizer.cObject )
#        QCOnePDM_Reorthogonalize       ( cdensityq,         porthogonalizer.cObject )

    def GetOption ( self, optionName ):
        """Get an option for the model."""
        value = None
        if   optionName == "isSpinRestricted"         : value = ( self.cObject.isSpinRestricted == CTrue )
        elif optionName == "numberFractionalAlphaHOOs": value = self.cObject.numberFractionalAlphaHOOs
        elif optionName == "numberFractionalAlphaLUOs": value = self.cObject.numberFractionalAlphaLUOs
        elif optionName == "numberFractionalBetaHOOs" : value = self.cObject.numberFractionalBetaHOOs
        elif optionName == "numberFractionalBetaLUOs" : value = self.cObject.numberFractionalBetaLUOs
        elif optionName == "numberFractionalHOOs"     : value = self.cObject.numberFractionalHOOs
        elif optionName == "numberFractionalLUOs"     : value = self.cObject.numberFractionalLUOs
        elif optionName == "occupancyType"            : value = QCOnePDMOccupancyType_ToString[self.cObject.occupancyType]
        else: raise AttributeError ( "Unknown option: " + optionName + "." )
        return value

    def GridPointDensities ( self, configuration, Coordinates3 gridPoints, spinDensities = False ):
        """Densities at grid points - grid points are already in atomic units."""
        cdef Boolean            cSpinDensities
        cdef Coordinates3    coordinates3
        cdef QCModelDFTState qcState
        cdef Real1DArray     data
        # . Initialization.
        data = None
        # . Options.
        if spinDensities: cSpinDensities = CTrue
        else:             cSpinDensities = CFalse
        # . Get the coordinate and state objects.
        if hasattr ( configuration, "coordinates3" ): coordinates3 = configuration.coordinates3
        else:                                         coordinates3 = None
        if hasattr ( configuration, "qcState"      ): qcState      = configuration.qcState
        else:                                         qcState      = None
        if ( coordinates3 is not None ) and ( qcState is not None ):
            data = Real1DArray.WithExtent ( gridPoints.rows )
            QCModelDFT_GridPointDensities ( self.cObject, qcState.cObject, coordinates3.cObject, gridPoints.cObject, data.cObject, cSpinDensities )
        return data

    def GridPointOrbitals ( self, configuration, Coordinates3 gridPoints, orbitals = None, useDensityP = True ):
        """Orbitals at grid points - grid points are already in atomic units."""
        cdef Boolean            useDensityPOrbitals
        cdef Integer             i, n, *orbitalindices
        cdef Coordinates3    coordinates3
        cdef QCModelDFTState qcState
        cdef Real1DArray     data
        # . Initialization.
        data           = None
        orbitalindices = NULL
        # . Options.
        if orbitals is not None:
            n = len ( orbitals )
            orbitalindices = Memory_Allocate_Array_Integer ( n )
            for i from 0 <= i < n: orbitalindices[i] = orbitals[i]
        else:
            n = 1
        if useDensityP: useDensityPOrbitals = CTrue
        else:      useDensityPOrbitals = CFalse
        # . Get the coordinate and state objects.
        if hasattr ( configuration, "coordinates3" ): coordinates3 = configuration.coordinates3
        else:                                         coordinates3 = None
        if hasattr ( configuration, "qcState"      ): qcState      = configuration.qcState
        else:                                         qcState      = None
        if ( coordinates3 is not None ) and ( qcState is not None ):
            data = Real1DArray.WithExtent ( gridPoints.rows * n )
            QCModelDFT_GridPointOrbitals ( self.cObject, qcState.cObject, coordinates3.cObject, gridPoints.cObject, data.cObject, n, orbitalindices, useDensityPOrbitals )
            Memory_Deallocate_Integer ( &orbitalindices )
        return data

    def GridPointPotentials ( self, configuration, Coordinates3 gridPoints ):
        """Electrostatic data at grid points - grid points are already in atomic units."""
        cdef Coordinates3    coordinates3
        cdef QCModelDFTState qcState
        cdef Real1DArray     data
        # . Initialization.
        data = None
        # . Get the coordinate and state objects.
        if hasattr ( configuration, "coordinates3" ): coordinates3 = configuration.coordinates3
        else:                                         coordinates3 = None
        if hasattr ( configuration, "qcState"      ): qcState      = configuration.qcState
        else:                                         qcState      = None
        if ( coordinates3 is not None ) and ( qcState is not None ):
            data = Real1DArray.WithExtent ( gridPoints.rows )
            QCModelDFT_GridPointPotentials ( self.cObject, qcState.cObject, coordinates3.cObject, gridPoints.cObject, data.cObject )
        return data

    def IsSpinRestricted ( self ):
        """Is the model spin restricted?"""
        isSpinRestricted = True
        if self.cObject != NULL: isSpinRestricted = ( self.cObject.isSpinRestricted == CTrue )
        return isSpinRestricted

    def MakeLabel ( self ):
        """Construct a model label."""
        if self.label is None:
            self.label = self.functional.upper ( ) + "/" + self.orbitalBasis.upper ( ) + "/" + self.densityBasis.upper ( )

    def OrbitalEnergies ( self, configuration, useDensityP = True ):
        """Get the orbital energies and HOMO and LUMO indices."""
        cdef Boolean         useDensityPOrbitals
        cdef Integer         homo, lumo
        cdef QCModelDFTState qcState
        cdef Real1DArray     data
        # . Initialization.
        data = None
        homo = -1
        lumo = -1
        # . Options.
        if useDensityP: useDensityPOrbitals = CTrue
        else:           useDensityPOrbitals = CFalse
        # . Get the coordinate and state objects.
        if hasattr ( configuration, "qcState"): qcState = configuration.qcState
        else:                                   qcState = None
        if qcState is not None:
            data         = Real1DArray.Raw ( )
            data.cObject = QCModelDFT_OrbitalEnergies ( self.cObject, qcState.cObject, useDensityPOrbitals, &homo, &lumo )
            data.isOwner = True
            if data.cObject == NULL: data = None
        return ( data, homo, lumo )

    def Orbitals ( self, configuration, useDensityP = True ):
        """Get the orbitals."""
        cdef Boolean         useDensityPOrbitals
        cdef Real2DArray     data
        cdef QCModelDFTState qcState
        # . Initialization.
        data = None
        # . Options.
        if useDensityP: useDensityPOrbitals = CTrue
        else:           useDensityPOrbitals = CFalse
        # . Get the state object.
        if hasattr ( configuration, "qcState"): qcState = configuration.qcState
        else:                                   qcState = None
        if qcState is not None:
            data         = Real2DArray.Raw ( )
            data.cObject = QCModelDFT_Orbitals ( self.cObject, qcState.cObject, useDensityPOrbitals )
            data.isOwner = True
            if data.cObject == NULL: data = None
        return data

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        if "converger"                 in keywordArguments: self.converger                         = keywordArguments.pop ( "converger"                 )
        if "fermiBroadening"           in keywordArguments: self.cObject.fermiBroadening           = keywordArguments.pop ( "fermiBroadening"           )
        if "label"                     in keywordArguments: self.label                             = keywordArguments.pop ( "label"                     )
        if "densityBasis"              in keywordArguments: self.densityBasis                      = keywordArguments.pop ( "densityBasis"              )
        if "orbitalBasis"              in keywordArguments: self.orbitalBasis                      = keywordArguments.pop ( "orbitalBasis"              )
        if "numberFractionalAlphaHOOs" in keywordArguments: self.cObject.numberFractionalAlphaHOOs = keywordArguments.pop ( "numberFractionalAlphaHOOs" )
        if "numberFractionalAlphaLUOs" in keywordArguments: self.cObject.numberFractionalAlphaLUOs = keywordArguments.pop ( "numberFractionalAlphaLUOs" )
        if "numberFractionalBetaHOOs"  in keywordArguments: self.cObject.numberFractionalBetaHOOs  = keywordArguments.pop ( "numberFractionalBetaHOOs"  )
        if "numberFractionalBetaLUOs"  in keywordArguments: self.cObject.numberFractionalBetaLUOs  = keywordArguments.pop ( "numberFractionalBetaLUOs"  )
        if "numberFractionalHOOs"      in keywordArguments: self.cObject.numberFractionalHOOs      = keywordArguments.pop ( "numberFractionalHOOs"      )
        if "numberFractionalLUOs"      in keywordArguments: self.cObject.numberFractionalLUOs      = keywordArguments.pop ( "numberFractionalLUOs"      )
        if "inCore"           in keywordArguments:
            if keywordArguments.pop ( "inCore"           ): self.cObject.inCore           = CTrue
            else:                                           self.cObject.inCore           = CFalse
        if "isSpinRestricted" in keywordArguments:
            if keywordArguments.pop ( "isSpinRestricted" ): self.cObject.isSpinRestricted = CTrue
            else:                                           self.cObject.isSpinRestricted = CFalse
        if "keepOrbitalData" in keywordArguments:
            if keywordArguments.pop ( "keepOrbitalData"  ): self.cObject.keepOrbitalData  = CTrue
            else:                                           self.cObject.keepOrbitalData  = CFalse
        if "linkAtomRatio"   in keywordArguments:
            if keywordArguments.pop ( "linkAtomRatio"    ): self.cObject.linkAtomRatio    = CTrue
            else:                                           self.cObject.linkAtomRatio    = CFalse
        if "accuracy"        in keywordArguments:
            keyword = keywordArguments.pop ( "accuracy" )
            option  = DFTGridAccuracy_ToEnum.get ( keyword.lower ( ).title  ( ), None )
            if option is None: raise TypeError ( "Unrecognized DFT grid accuracy option: " + keyword + "." )
            else:              self.cObject.accuracy = option
        if "functional"      in keywordArguments:
            keyword = keywordArguments.pop ( "functional" )
            self.DefineFunctionalModel ( keyword )
        if "occupancyType"   in keywordArguments:
            keyword = keywordArguments.pop ( "occupancyType" )
            option  = QCOnePDMOccupancyType_ToEnum.get ( keyword, None )
            if option is None: raise TypeError ( "Unrecognized occupancy type option: " + keyword + "." )
            else:              self.cObject.occupancyType = option
        if "qcChargeModel"   in keywordArguments:
            keyword = keywordArguments.pop ( "qcChargeModel" )
            option  = QCChargeModel_ToEnum.get ( keyword, None )
            if option is None: raise TypeError ( "Unrecognized QC charge model option: " + keyword + "." )
            else:              self.cObject.qcChargeModel = option
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )
        self.CheckOptions ( )

    def SetUp ( self, QCAtomContainer qcAtoms, QCParameters qcParameters, electronicState, configuration ):
        """Set up the QC model state."""
        cdef Real                alphaCharge, betaCharge
        cdef QCModelDFTState       qcState
        cdef CQCAtomContainer *cqcAtoms
        cdef CQCParameter     *cqcParameters
        # . There must be a configuration.
        if ( configuration is not None ):
            # . Get data for the new state.
            if qcAtoms      is None: cqcAtoms      = NULL
            else:                    cqcAtoms      = qcAtoms.cObject
            if qcParameters is None: cqcParameters = NULL
            else:                    cqcParameters = qcParameters.cObject
            # . Check the consistency of the model with the system.
            ( alphaCharge, betaCharge ) = self.CrossCheckModelAndState ( electronicState.charge, electronicState.multiplicity, cqcAtoms.nuclearCharge, cqcAtoms.nobasis )
            # . Create the new state.
            qcState         = QCModelDFTState.Raw ( )
            qcState.cObject = QCModelDFTState_Setup ( cqcAtoms, cqcParameters, alphaCharge, betaCharge, self.cObject.occupancyType            , self.cObject.isSpinRestricted         , \
                                                                                                        self.cObject.numberFractionalHOOs     , self.cObject.numberFractionalLUOs     , \
                                                                                                        self.cObject.numberFractionalAlphaHOOs, self.cObject.numberFractionalAlphaLUOs, \
                                                                                                        self.cObject.numberFractionalBetaHOOs , self.cObject.numberFractionalBetaLUOs , \
                                                                                                                                                           self.cObject.fermiBroadening )
            qcState.isOwner = True
            if qcState.cObject == NULL: raise ValueError ( "Unable to create QCModel state either due to a memory problem or because the model is inapplicable." )
            setattr ( configuration, "qcState", qcState )
