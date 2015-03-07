#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModelMNDO.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines a MNDO quantum chemical model."""

import math, os

from pCore                       import logFile, LogFileActive

from QCChargeModelOptions        import QCChargeModel_ToEnum, QCChargeModel_ToString
from QCOnePDM                    import QCOnePDMOccupancyType_ToEnum, QCOnePDMOccupancyType_ToString
from QCParameters                import QCParameters_Define
from QCMMLinkAtomCouplingOptions import QCMMLinkAtomCoupling_ToEnum, QCMMLinkAtomCoupling_ToString

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Default values.
_DEFAULT_MNDOMODEL    = "am1"
#_DEFAULT_ORBITALBASIS = "mndosto6g"
_DEFAULT_ORBITALBASIS = "mndostong"

# #ifdef MNDOCI
# . Dictionaries for interconverting between C and Python representations.
MNDOCIAlgorithm_ToEnum   = { "Direct" : MNDOCIAlgorithm_Direct, "Full" : MNDOCIAlgorithm_Full, "Sparse" : MNDOCIAlgorithm_Sparse }
MNDOCIAlgorithm_ToString = { MNDOCIAlgorithm_Direct : "Direct", MNDOCIAlgorithm_Full : "Full", MNDOCIAlgorithm_Sparse : "Sparse" }

MNDOCIMethod_ToEnum   = { "Doubles" : MNDOCIMethod_Doubles, "Full" : MNDOCIMethod_Full, "Singles" : MNDOCIMethod_Singles, "Singles/Doubles" : MNDOCIMethod_SinglesDoubles, "User" : MNDOCIMethod_UserSpecified }
MNDOCIMethod_ToString = { MNDOCIMethod_Doubles : "Doubles", MNDOCIMethod_Full : "Full", MNDOCIMethod_Singles : "Singles", MNDOCIMethod_SinglesDoubles : "Singles/Doubles", MNDOCIMethod_UserSpecified : "User" }
# #endif /*MNDOCI*/

# . Tolerances.
_CardinalOccupancyTolerance = 1.0e-12

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModelMNDO ( QCModel ):
    """Define a semi-empirical MNDO-type quantum chemical model."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            QCModelMNDO_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.QCModelMNDO"

    def __getstate__ ( self ):
        """Return the state."""
# #ifdef MNDOCI
        cdef Integer i, j
        cdef MNDOCIModel *ciModel
        ciModel = self.cObject.ciModel
        hasCI   = ( ciModel != NULL )
# #endif /*MNDOCI*/
        state = { "converger"                 : self.converger                         ,
                  "label"                     : self.label                             ,
                  "mndoModel"                 : self.mndoModel                         ,
                  "orbitalBasis"              : self.orbitalBasis                      ,
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
                  "occupancyType"             : QCOnePDMOccupancyType_ToString[self.cObject.occupancyType],
# #ifdef MNDOCI
                  "qcChargeModel"             : QCChargeModel_ToString[self.cObject.qcChargeModel],
                  "hasCI"                     : hasCI }
        if hasCI:
            state.update ( { "checkAlgorithm"      : ciModel.checkAlgorithm      ,
                             "doAllStates"         : ciModel.doAllStates         ,
                             "identifyRootSpin"    : ciModel.identifyRootSpin    ,
                             "localizeOrbitals"    : ciModel.localizeOrbitals    ,
                             "activeElectrons"     : ciModel.activeElectrons     ,
                             "activeOrbitals"      : ciModel.activeOrbitals      ,
                             "localizeStart"       : ciModel.localizeStart       ,
                             "localizeStop"        : ciModel.localizeStop        ,
                             "minimalMultiplicity" : ciModel.minimalMultiplicity ,
                             "numberOfStates"      : ciModel.numberOfStates      ,
                             "requiredRoot"        : ciModel.requiredRoot        ,
                             "rootMultiplicity"    : ciModel.rootMultiplicity    ,
                             "CIAlgorithm"         : MNDOCIAlgorithm_ToString[ciModel.algorithm],
                             "CIMethod"            : MNDOCIMethod_ToString   [ciModel.method   ],
                             "eigenvalueSolverIterations"      : ciModel.eigenvalueSolver.maximumMatrixVectorMultiplications,
                             "eigenvalueSolverPreconditioning" : ciModel.eigenvalueSolver.usePreconditioning,
                             "eigenvalueSolverTolerance"       : ciModel.eigenvalueSolver.errorTolerance } )
            # . User-specified microstates.
            if ciModel.method == MNDOCIMethod_UserSpecified:
                microstates = []
                for i from 0 <= i < ciModel.microstates.length0:
                    microstate = []
                    for j from 0 <= j < ciModel.microstates.length1: microstate.append ( "{:d}".format ( Integer2DArray_GetItem ( ciModel.microstates, i, j, NULL ) ) )
                    microstates.append ( "".join ( microstate ) )
                state["microStates"] = microstates
# #else /*MNDOCI*/
#                  "qcChargeModel"             : QCChargeModel_ToString[self.cObject.qcChargeModel] }
# #endif /*MNDOCI*/
        return state

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = QCModelMNDO_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        super ( QCModelMNDO, self )._Initialize ( )
        self.cObject      = NULL
        self.isOwner      = False
        self.mndoModel    = _DEFAULT_MNDOMODEL
        self.orbitalBasis = _DEFAULT_ORBITALBASIS

    def AtomicCharges ( self, configuration, chargeModel = None, spinDensities = False ):
        """Atomic charges."""
        cdef Boolean          cSpinDensities
        cdef QCModelMNDOState qcState
        cdef Real1DArray      qcCharges
        # . Initialization.
        qcCharges = None
        # . Get the state object.
        if hasattr ( configuration, "qcState" ):
            qcState = configuration.qcState
            if spinDensities: cSpinDensities = CTrue
            else:             cSpinDensities = CFalse
            qcCharges = Real1DArray.WithExtent ( len ( qcState ) )
            qcCharges.Set ( 0.0 )
            QCModelMNDO_LowdinCharges ( self.cObject, qcState.cObject, cSpinDensities, qcCharges.cObject )
        return qcCharges

    def BondOrders ( self, configuration, chargeModel = None ):
        """Bond orders."""
        cdef QCModelMNDOState qcState
        cdef SymmetricMatrix  bondOrders
        cdef Real1DArray      charges, freevalence, totalvalence
        # . Initialization.
        bondOrders   = None
        charges      = None
        freevalence  = None
        totalvalence = None
        # . Get the state object.
        if hasattr ( configuration, "qcState" ):
            qcState      = configuration.qcState
            bondOrders   = SymmetricMatrix.WithExtent ( len ( qcState ) ) ; bondOrders.Set   ( 0.0 )
            charges      = Real1DArray.WithExtent     ( len ( qcState ) ) ; charges.Set      ( 0.0 )
            freevalence  = Real1DArray.WithExtent     ( len ( qcState ) ) ; freevalence.Set  ( 0.0 )
            totalvalence = Real1DArray.WithExtent     ( len ( qcState ) ) ; totalvalence.Set ( 0.0 )
            QCModelMNDO_MayerBondOrders ( self.cObject, qcState.cObject, bondOrders.cObject, charges.cObject, freevalence.cObject, totalvalence.cObject )
        return ( bondOrders, charges, freevalence, totalvalence )

    def ChargeConstraintFunction ( self, QCModelMNDOState qcState, ChargeConstraintContainer chargeConstraints, Real1DArray lambdas, Real1DArray gradients, SymmetricMatrix hessian, buildFock = True ):
        """Calculate the charge constraint function and, optionally, the gradient and hessian."""
        cdef Boolean           cBuildFock
        cdef CReal1DArray     *cGradients
        cdef CSymmetricMatrix *cHessian
        if buildFock: cBuildFock = CTrue
        else:         cBuildFock = CFalse
        if gradients is None: cGradients = NULL
        else:                 cGradients = gradients.cObject
        if hessian   is None: cHessian   = NULL
        else:                 cHessian   = hessian.cObject
        return QCModelMNDO_CCFunction ( self.cObject, qcState.cObject, chargeConstraints.cObject, lambdas.cObject, cBuildFock, cGradients, cHessian, NULL )

# #ifdef MNDOCI
    def CIStateCharacters ( self, configuration, Matrix33 rotation, Integer1DArray mapping, Selection stateIndices, includeCoreOrbitals = False ):
        """Determine the characters of CI states under a rotation."""
        cdef Boolean          cIncludeCoreOrbitals
        cdef Real1DArray      characters
        cdef QCModelMNDOState state
        # . Initialization.
        characters = None
        # . Options.
        if includeCoreOrbitals: cIncludeCoreOrbitals = CTrue
        else:                   cIncludeCoreOrbitals = CFalse
        # . Get the state object.
        if hasattr ( configuration, "qcState"): state = configuration.qcState
        else:                                   state = None
        # . Determine the characters.
        if state is not None:
            characters         = Real1DArray.Raw ( )
            characters.cObject = QCModelMNDO_CIStateCharacters ( self.cObject, state.cObject, rotation.cObject, mapping.cObject, stateIndices.cObject, cIncludeCoreOrbitals, NULL )
            characters.isOwner = True
            if characters.cObject == NULL: characters = None
        return characters

    def CrossCheckCIModelAndState ( self, multiplicity, alphaCharge, betaCharge ):
        """Extra checks for a CI model."""
        cdef Integer rootMultiplicity
        # . Basic checks.
        if not self.GetOption ( "isSpinRestricted" ): raise ValueError ( "A CI calculation requires a spin-restricted QC model." )
        numberElectrons = int ( round ( alphaCharge + betaCharge ) )
        if ( math.fabs ( float ( numberElectrons ) - ( alphaCharge + betaCharge ) ) > _CardinalOccupancyTolerance ): raise ValueError ( "A CI calculation requires an integer number of electrons." )
        # . Ground-state multiplicity.
        numberAlpha = ( numberElectrons + ( multiplicity - 1 ) ) // 2
        numberBeta  = ( numberElectrons - ( multiplicity - 1 ) ) // 2
        if ( ( ( numberAlpha + numberBeta ) != numberElectrons ) or ( numberBeta < 0 ) ): raise ValueError ( "Impossible spin multiplicity for CI ground state." )
        # . Required root multiplicity.
        rootMultiplicity = self.cObject.ciModel.rootMultiplicity
        if rootMultiplicity > 0:
            if ( ( ( numberElectrons + rootMultiplicity - 1 ) % 2 != 0 ) or ( numberElectrons < ( rootMultiplicity - 1 ) ) ): raise ValueError ( "Impossible spin multiplicity for the required CI root." )
        return ( numberAlpha, numberBeta )
# #endif /*MNDOCI*/

    def DefineParameters ( self, anlist, log = logFile ):
        """Define parameters."""
        path           = os.getenv ( "PDYNAMO_PARAMETERS" )
        if path is None: raise ValueError ( "Parameter path environmental variable missing." )
        parametersPath = os.path.join ( path, "mndoParameters" )
        return QCParameters_Define ( anlist, mndoPath    = os.path.join ( parametersPath, self.mndoModel    ),
                                             orbitalPath = os.path.join ( parametersPath, self.orbitalBasis ),
                                             log = log )

    def DipoleMoment ( self, configuration, Vector3 center = None ):
        """Dipole Moment."""
        cdef Coordinates3     coordinates3
        cdef QCModelMNDOState qcState
        cdef Vector3          dipole
        cdef CVector3        *ccenter
        # . Initialization.
        dipole = None
        # . Get the coordinate and state objects.
        if hasattr ( configuration, "coordinates3" ): coordinates3 = configuration.coordinates3
        else:                                         coordinates3 = None
        if hasattr ( configuration, "qcState"      ): qcState      = configuration.qcState
        else:                                         qcState      = None
        if ( coordinates3 is not None ) and ( qcState is not None ):
            if center is None: ccenter = NULL
            else:              ccenter = center.cObject
            dipole         = Vector3.Raw ( )
            dipole.cObject = QCModelMNDO_DipoleMoment ( self.cObject, qcState.cObject, coordinates3.cObject, ccenter )
            dipole.isOwner = True
        return dipole

    def EnergyFinalize ( self, configuration ):
        """Finalization."""
        cdef QCModelMNDOState qcState
        if hasattr ( configuration, "qcState" ):
            qcState = configuration.qcState
            QCModelMNDOState_Finalize ( qcState.cObject, self.cObject.keepOrbitalData )

    def EnergyFock ( self, configuration ):
        """Fock matrix construction."""
        cdef Real           eelectronic
        cdef QCModelMNDOState qcState
        eelectronic = 0.0
        if hasattr ( configuration, "qcState" ):
            qcState = configuration.qcState
            QCModelMNDO_Fock           ( self.cObject, qcState.cObject, &eelectronic )
            QCModelMNDO_QCMMFockLowdin ( self.cObject, qcState.cObject, &eelectronic )
        return eelectronic

    def EnergyGradients ( self, configuration ):
        """Gradients."""
        cdef QCModelMNDOState qcState
        if hasattr ( configuration, "qcState" ):
            qcState = configuration.qcState
            QCModelMNDO_Gradients ( self.cObject, qcState.cObject )

    def EnergyInitialize ( self, configuration ):
        """Integrals and space allocation."""
        cdef Coordinates3           coordinates3, pgradients3
        cdef QCMMInteractionState   pqcmmstate
        cdef QCModelMNDOState       qcState
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
            QCModelMNDOState_Initialize ( qcState.cObject, coordinates3.cObject, cqcmmstate, cgradients3 )
            # . Calculate the integrals and nuclear energy.
            QCModelMNDO_Integrals ( self.cObject, qcState.cObject )

# #ifdef MNDOCI
    def EnergyPostProcessing ( self, configuration ):
        """Do a CI calculation if required."""
        cdef QCModelMNDOState qcState
        if hasattr ( configuration, "qcState" ): qcState = configuration.qcState
        else:                                    qcState = None
        if qcState is not None: QCModelMNDO_CIEnergy ( self.cObject, qcState.cObject )
# #endif /*MNDOCI*/

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
        cdef Boolean             cSpinDensities
        cdef Coordinates3     coordinates3
        cdef QCModelMNDOState qcState
        cdef Real1DArray      data
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
            QCModelMNDO_GridPointDensities ( self.cObject, qcState.cObject, coordinates3.cObject, gridPoints.cObject, data.cObject, cSpinDensities )
        return data

    def GridPointOrbitals ( self, configuration, Coordinates3 gridPoints, orbitals = None, useDensityP = True ):
        """Orbitals at grid points - grid points are already in atomic units."""
        cdef Boolean             useDensityPOrbitals
        cdef Integer              i, n, *orbitalindices
        cdef Coordinates3     coordinates3
        cdef QCModelMNDOState qcState
        cdef Real1DArray      data
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
            QCModelMNDO_GridPointOrbitals ( self.cObject, qcState.cObject, coordinates3.cObject, gridPoints.cObject, data.cObject, n, orbitalindices, useDensityPOrbitals )
            Memory_Deallocate_Integer ( &orbitalindices )
        return data

    def GridPointPotentials ( self, configuration, Coordinates3 gridPoints ):
        """Electrostatic data at grid points - grid points are already in atomic units."""
        cdef Coordinates3     coordinates3
        cdef QCModelMNDOState qcState
        cdef Real1DArray      data
        # . Initialization.
        data = None
        # . Get the coordinate and state objects.
        if hasattr ( configuration, "coordinates3" ): coordinates3 = configuration.coordinates3
        else:                                         coordinates3 = None
        if hasattr ( configuration, "qcState"      ): qcState      = configuration.qcState
        else:                                         qcState      = None
        if ( coordinates3 is not None ) and ( qcState is not None ):
            data = Real1DArray.WithExtent ( gridPoints.rows )
            QCModelMNDO_GridPointPotentials ( self.cObject, qcState.cObject, coordinates3.cObject, gridPoints.cObject, data.cObject )
        return data

    def IsSpinRestricted ( self ):
        """Is the model spin restricted?"""
        isSpinRestricted = True
        if self.cObject != NULL: isSpinRestricted = ( self.cObject.isSpinRestricted == CTrue )
        return isSpinRestricted

    def LocalizeOrbitals ( self, configuration, doBeta = False, startOrbital = -1, stopOrbital = -1 ):
        """Localize a set of orbitals."""
        cdef Boolean          doQ
        cdef QCModelMNDOState qcState
        # . Options.
        if doBeta: doQ = CTrue
        else:      doQ = CFalse
        # . Get the coordinate and state objects.
        if hasattr ( configuration, "qcState"): qcState = configuration.qcState
        else:                                   qcState = None
        if qcState is not None:
            QCModelMNDO_LocalizeOrbitals ( self.cObject, qcState.cObject, doQ, startOrbital, stopOrbital, NULL )

    def MakeLabel ( self ):
        """Construct a model label."""
        if self.label is None: self.label = self.mndoModel.upper ( )

    def OrbitalCharacters ( self, configuration, Matrix33 rotation, Integer1DArray mapping, Selection orbitalIndices, useDensityP = True ):
        """Determine the characters of orbitals under a rotation."""
        cdef Boolean          cUseDensityP
        cdef Real1DArray      characters
        cdef QCModelMNDOState state
        # . Initialization.
        characters = None
        # . Options.
        if useDensityP: cUseDensityP = CTrue
        else:           cUseDensityP = CFalse
        # . Get the state object.
        if hasattr ( configuration, "qcState"): state = configuration.qcState
        else:                                   state = None
        # . Determine the characters.
        if state is not None:
            characters         = Real1DArray.Raw ( )
            characters.cObject = QCModelMNDO_OrbitalCharacters ( self.cObject, state.cObject, rotation.cObject, mapping.cObject, cUseDensityP, orbitalIndices.cObject, NULL )
            characters.isOwner = True
        return characters

    def OrbitalEnergies ( self, configuration, useDensityP = True ):
        """Get the orbital energies and HOMO and LUMO indices."""
        cdef Boolean          useDensityPOrbitals
        cdef Integer          homo, lumo
        cdef QCModelMNDOState qcState
        cdef Real1DArray      data
        # . Initialization.
        data = None
        homo = -1
        lumo = -1
        # . Options.
        if useDensityP: useDensityPOrbitals = CTrue
        else:           useDensityPOrbitals = CFalse
        # . Get the state object.
        if hasattr ( configuration, "qcState"): qcState = configuration.qcState
        else:                                   qcState = None
        if qcState is not None:
            data         = Real1DArray.Raw ( )
            data.cObject = QCModelMNDO_OrbitalEnergies ( self.cObject, qcState.cObject, useDensityPOrbitals, &homo, &lumo )
            data.isOwner = True
            if data.cObject == NULL: data = None
        return ( data, homo, lumo )

    def Orbitals ( self, configuration, useDensityP = True ):
        """Get the orbitals."""
        cdef Boolean          useDensityPOrbitals
        cdef Real2DArray      data
        cdef QCModelMNDOState qcState
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
            data.cObject = QCModelMNDO_Orbitals ( self.cObject, qcState.cObject, useDensityPOrbitals )
            data.isOwner = True
            if data.cObject == NULL: data = None
        return data

    def PrintOnePDMArrays ( self, configuration, log = logFile ):
        """Print the arrays from the one-particle density matrix data structures."""
        cdef QCAtomContainer  qcAtoms
        cdef QCModelMNDOState qcState
        cdef QCParameters     qcParameters
        if LogFileActive ( logFile ):
            if hasattr ( configuration, "qcState"): qcState = configuration.qcState
            else:                                   qcState = None
            if qcState is not None:
                # . Basis function labels.
                qcAtoms              = QCAtomContainer.Raw ( )
                qcAtoms.cObject      = qcState.cObject.qcAtoms
                qcAtoms.isOwner      = False
                qcParameters         = QCParameters.Raw ( )
                qcParameters.cObject = qcState.cObject.qcParameters
                qcParameters.isOwner = False
                labels = qcAtoms.GetOrbitalBasisFunctionLabels ( qcParameters )
                # . Densities.
                densityp = qcState.densityp
                densityq = qcState.densityq
                if densityq is None: onepdms = [ ( "Total", densityp ) ]
                else:                onepdms = [ ( "Alpha", densityp ), ( "Beta", densityq ) ]
                for ( tag, onepdm ) in onepdms:
                    log.Heading   ( "Data for " + tag + " One-Particle Density Matrix", includeBlankLine = True )
                    density = onepdm.density
                    if density is not None: density.Print ( labels = labels, log = log, title = tag + " Density Matrix" )
                    fock = onepdm.fock
                    if fock is not None: fock.Print ( labels = labels, log = log, title = tag + " Fock Matrix" )
                    orbitals = onepdm.orbitals
                    if orbitals is not None:
                        columnLabels = []
                        occupancies  = onepdm.occupancies
                        if occupancies is not None: columnLabels.append ( ( "Occupancies", occupancies ) )
                        energies     = onepdm.energies
                        if energies    is not None: columnLabels.append ( ( "Energies"   , energies    ) )
                        orbitals.Print ( columnLabels = columnLabels, rowLabels = labels, log = log, title = "Orbitals from " + tag + " Density Matrix" )

    def ProcessArguments ( self, *arguments ):
        """Process constructor arguments."""
        if len ( arguments ) > 0: setattr ( self, "mndoModel", arguments[0] )

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
# #ifdef MNDOCI
        cdef Integer          i, j
        cdef CInteger2DArray *cmicrostates
        cdef MNDOCIModel     *ciModel
# #endif /*MNDOCI*/
        # . Basic options.
        if "converger"                 in keywordArguments: self.converger                         = keywordArguments.pop ( "converger"                 )
        if "fermiBroadening"           in keywordArguments: self.cObject.fermiBroadening           = keywordArguments.pop ( "fermiBroadening"           )
        if "label"                     in keywordArguments: self.label                             = keywordArguments.pop ( "label"                     )
        if "mndoModel"                 in keywordArguments: self.mndoModel                         = keywordArguments.pop ( "mndoModel"                 )
        if "orbitalBasis"              in keywordArguments: self.orbitalBasis                      = keywordArguments.pop ( "orbitalBasis"              )
        if "numberFractionalAlphaHOOs" in keywordArguments: self.cObject.numberFractionalAlphaHOOs = keywordArguments.pop ( "numberFractionalAlphaHOOs" )
        if "numberFractionalAlphaLUOs" in keywordArguments: self.cObject.numberFractionalAlphaLUOs = keywordArguments.pop ( "numberFractionalAlphaLUOs" )
        if "numberFractionalBetaHOOs"  in keywordArguments: self.cObject.numberFractionalBetaHOOs  = keywordArguments.pop ( "numberFractionalBetaHOOs"  )
        if "numberFractionalBetaLUOs"  in keywordArguments: self.cObject.numberFractionalBetaLUOs  = keywordArguments.pop ( "numberFractionalBetaLUOs"  )
        if "numberFractionalHOOs"      in keywordArguments: self.cObject.numberFractionalHOOs      = keywordArguments.pop ( "numberFractionalHOOs"      )
        if "numberFractionalLUOs"      in keywordArguments: self.cObject.numberFractionalLUOs      = keywordArguments.pop ( "numberFractionalLUOs"      )
        if "isSpinRestricted" in keywordArguments:
            if keywordArguments.pop ( "isSpinRestricted" ): self.cObject.isSpinRestricted = CTrue
            else:                                           self.cObject.isSpinRestricted = CFalse
        if "keepOrbitalData" in keywordArguments:
            if keywordArguments.pop ( "keepOrbitalData"  ): self.cObject.keepOrbitalData  = CTrue
            else:                                           self.cObject.keepOrbitalData  = CFalse
        if "linkAtomRatio"   in keywordArguments:
            if keywordArguments.pop ( "linkAtomRatio"    ): self.cObject.linkAtomRatio    = CTrue
            else:                                           self.cObject.linkAtomRatio    = CFalse
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
# #ifdef MNDOCI
        # . CI options.
        hasCI = keywordArguments.pop ( "hasCI", False ) or ( keywordArguments.get ( "CIMethod", None ) is not None )
        if hasCI:
            ciModel = self.cObject.ciModel
            if ciModel == NULL:
                ciModel = MNDOCIModel_Allocate ( NULL )
                self.cObject.ciModel = ciModel
            if "activeElectrons"            in keywordArguments: ciModel.activeElectrons     = keywordArguments.pop ( "activeElectrons"     )
            if "activeOrbitals"             in keywordArguments: ciModel.activeOrbitals      = keywordArguments.pop ( "activeOrbitals"      )
            if "localizeStart"              in keywordArguments: ciModel.localizeStart       = keywordArguments.pop ( "localizeStart"       )
            if "localizeStop"               in keywordArguments: ciModel.localizeStop        = keywordArguments.pop ( "localizeStop"        )
            if "minimalMultiplicity"        in keywordArguments: ciModel.minimalMultiplicity = keywordArguments.pop ( "minimalMultiplicity" )
            if "numberOfStates"             in keywordArguments: ciModel.numberOfStates      = keywordArguments.pop ( "numberOfStates"      )
            if "requiredRoot"               in keywordArguments: ciModel.requiredRoot        = keywordArguments.pop ( "requiredRoot"        )
            if "rootMultiplicity"           in keywordArguments: ciModel.rootMultiplicity    = keywordArguments.pop ( "rootMultiplicity"    )
            if "doAllStates"         in keywordArguments:
                if keywordArguments.pop ( "doAllStates"      ): ciModel.doAllStates      = CTrue
                else:                                           ciModel.doAllStates      = CFalse
            if "identifyRootSpin" in keywordArguments:
                if keywordArguments.pop ( "identifyRootSpin" ): ciModel.identifyRootSpin = CTrue
                else:                                           ciModel.identifyRootSpin = CFalse
            if "localizeOrbitals" in keywordArguments:
                if keywordArguments.pop ( "localizeOrbitals" ): ciModel.localizeOrbitals = CTrue
                else:                                           ciModel.localizeOrbitals = CFalse
            # . The CI algorithm.
            if "CIAlgorithm" in keywordArguments:
                keyword = keywordArguments.pop ( "CIAlgorithm" )
                option  = MNDOCIAlgorithm_ToEnum.get ( keyword, None )
                if option is None: raise TypeError ( "Unrecognized CI Algorithm: " + keyword + "." )
                else:              ciModel.algorithm = option
            if "checkAlgorithm" in keywordArguments:
                if keywordArguments.pop ( "checkAlgorithm" ): ciModel.checkAlgorithm = CTrue
                else:                                         ciModel.checkAlgorithm = CFalse
            # . Solver.
            if "eigenvalueSolverIterations" in keywordArguments: ciModel.eigenvalueSolver.maximumMatrixVectorMultiplications = keywordArguments.pop ( "eigenvalueSolverIterations" )
            if "eigenvalueSolverTolerance"  in keywordArguments: ciModel.eigenvalueSolver.errorTolerance                     = keywordArguments.pop ( "eigenvalueSolverTolerance"  )
            if "eigenvalueSolverPreconditioning" in keywordArguments:
                if keywordArguments.pop ( "eigenvalueSolverPreconditioning" ): ciModel.eigenvalueSolver.usePreconditioning = CTrue
                else:                                                          ciModel.eigenvalueSolver.usePreconditioning = CFalse
            # . Basic checks.
            if ciModel.rootMultiplicity > 0: ciModel.identifyRootSpin = CTrue
            else:                            ciModel.identifyRootSpin = CFalse
            # . The CI method.
            if "CIMethod" in keywordArguments:
                keyword = keywordArguments.pop ( "CIMethod" )
                option  = MNDOCIMethod_ToEnum.get ( keyword, None )
                if option is None: raise TypeError ( "Unrecognized CI Method: " + keyword + "." )
                else:              ciModel.method = option
                # . Generate microstates.
                microstates = keywordArguments.pop ( "microStates", None )
                if option == MNDOCIMethod_UserSpecified:
                    if microstates is None: raise TypeError ( "User specified CI method requested but microstates were not found." )
                    cmicrostates = Integer2DArray_Allocate ( len ( microstates ), 2 * ciModel.activeOrbitals, NULL )
                    isOK         = False
                    if cmicrostates != NULL:
                        isOK = True
                        Integer2DArray_Set ( cmicrostates, 0 )
                        for ( i, microstate ) in enumerate ( microstates ):
                            indices = []
                            for m in microstate: indices.append ( int ( m ) )
                            if ( len ( indices ) != ( 2 * ciModel.activeOrbitals ) ) or ( sum ( indices ) != ciModel.activeElectrons ):
                                isOK = False
                                break
                            for j from 0 <= j < ( 2 * ciModel.activeOrbitals ): Integer2DArray_SetItem ( cmicrostates, i, j, indices[j], NULL )
                        ciModel.microstates = cmicrostates
                    if not isOK: raise TypeError ( "Invalid user microstate specification." )
# #endif /*MNDOCI*/
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )

    def SetUp ( self, QCAtomContainer qcAtoms, QCParameters qcParameters, electronicState, configuration ):
        """Set up the QC model state."""
        cdef Real                alphaCharge, betaCharge
        cdef Integer               numberAlpha, numberBeta
# #ifdef MNDOCI
        cdef MNDOCIModel          *ciModel
        cdef MNDOCIState          *cistate
# #endif /*MNDOCI*/
        cdef QCModelMNDOState      qcState
        cdef CQCAtomContainer *cqcAtoms
        cdef CQCParameter     *cqcParameters
# #ifdef MNDOCI
        # . Check for a CI model.
        hasCI = ( self.cObject.ciModel != NULL )
# #endif /*MNDOCI*/
        # . There must be a configuration.
        if ( configuration is not None ):
            # . Get data for the new state.
            if qcAtoms      is None: cqcAtoms      = NULL
            else:                    cqcAtoms      = qcAtoms.cObject
            if qcParameters is None: cqcParameters = NULL
            else:                    cqcParameters = qcParameters.cObject
            # . Check the consistency of the model with the system.
# #ifdef MNDOCI
            ( alphaCharge, betaCharge ) = self.CrossCheckModelAndState   ( electronicState.charge, electronicState.multiplicity, cqcAtoms.nuclearCharge, cqcAtoms.nobasis, permitFractionalRestrictedNonSinglet = hasCI )
            if hasCI: ( numberAlpha, numberBeta ) = self.CrossCheckCIModelAndState ( electronicState.multiplicity, alphaCharge, betaCharge )
# #else /*MNDOCI*/
#            ( alphaCharge, betaCharge ) = self.CrossCheckModelAndState   ( electronicState.charge, electronicState.multiplicity, cqcAtoms.nuclearCharge, cqcAtoms.nobasis )
# #endif /*MNDOCI*/
            # . Create the new state.
            qcState         = QCModelMNDOState.Raw ( )
            qcState.cObject = QCModelMNDOState_Setup ( cqcAtoms, cqcParameters, alphaCharge, betaCharge, self.cObject.occupancyType            , self.cObject.isSpinRestricted         , \
                                                                                                         self.cObject.numberFractionalHOOs     , self.cObject.numberFractionalLUOs     , \
                                                                                                         self.cObject.numberFractionalAlphaHOOs, self.cObject.numberFractionalAlphaLUOs, \
                                                                                                         self.cObject.numberFractionalBetaHOOs , self.cObject.numberFractionalBetaLUOs , \
                                                                                                                                                            self.cObject.fermiBroadening )
            qcState.isOwner = True
            isOK = ( qcState.cObject != NULL )
            if isOK:
                setattr ( configuration, "qcState", qcState )
# #ifdef MNDOCI
                # . CI state.
                if hasCI:
                    ciModel = self.cObject.ciModel
                    cistate = MNDOCIModel_MakeState ( ciModel, electronicState.multiplicity, numberAlpha, numberBeta, cqcAtoms.nobasis, NULL )
                    if cistate == NULL: isOK = False
                    else: qcState.cObject.cistate = cistate
# #endif /*MNDOCI*/
            # . Check for an error.
            if not isOK: raise ValueError ( "Unable to create QCModel state either due to a memory problem or because the model is inapplicable." )
