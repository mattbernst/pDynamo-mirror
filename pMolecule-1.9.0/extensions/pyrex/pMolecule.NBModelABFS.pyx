#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelABFS.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines an atom-based force-switching NB model."""

from pCore                       import CLibraryError, logFile, LogFileActive, UNITS_LENGTH_BOHRS_TO_ANGSTROMS
from QCMMLinkAtomCouplingOptions import QCMMLinkAtomCoupling_ToEnum, QCMMLinkAtomCoupling_ToString

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelABFS ( NBModel ):
    """Define an ABFS NB model."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner: NBModelABFS_Deallocate ( &self.cObject )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModelABFS"

    def __getstate__ ( self ):
        """Return the state."""
        state =  { "checkForInverses"     : ( self.cObject.checkForInverses == CTrue ) ,
                   "imageExpandFactor"    :   self.cObject.imageExpandFactor           , 
                   "dampingCutoff"        :   self.cObject.dampingCutoff               , 
                   "dielectric"           :   self.cObject.dielectric                  , 
                   "electrostaticScale14" :   self.cObject.electrostaticScale14        , 
                   "innerCutoff"          :   self.cObject.innerCutoff                 , 
                   "listCutoff"           :   self.cObject.listCutoff                  , 
                   "outerCutoff"          :   self.cObject.outerCutoff                 , 
                   "qcmmCoupling"         : QCMMLinkAtomCoupling_ToString[self.cObject.qcmmCoupling] ,
                   "useCentering"         : ( self.cObject.useCentering == CTrue )     }
        if self.generator               is not None: state["generator"              ] = self.generator
        if self.mmmmPairwiseInteraction is not None: state["mmmmPairwiseInteraction"] = self.mmmmPairwiseInteraction
        if self.qcmmPairwiseInteraction is not None: state["qcmmPairwiseInteraction"] = self.qcmmPairwiseInteraction
        if self.qcqcPairwiseInteraction is not None: state["qcqcPairwiseInteraction"] = self.qcqcPairwiseInteraction
        return state

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = NBModelABFS_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject                 = NULL
        self.generator               = None
        self.isOwner                 = False
        self.label                   = "ABFS"
        self.mmmmPairwiseInteraction = None
        self.qcmmPairwiseInteraction = None
        self.qcqcPairwiseInteraction = None

    def CheckGenerator ( self ):
        """Check the pairlist generator."""
        # . Create a suitable default generator.
        if self.generator is None: self.generator = PairListGenerator.FromOptions ( cutoff               =  self.cObject.listCutoff ,
                                                                                    cutoffCellSizeFactor =  0.5   ,
                                                                                    minimumCellExtent    =  2     ,
                                                                                    minimumCellSize      =  3.0   ,
                                                                                    minimumExtentFactor  =  1.5   ,
                                                                                    minimumPoints        =  500   ,
                                                                                    sortIndices          =  False ,
                                                                                    useGridByCell        =  True  )
        # . Ensure compatibility between the generator and the list cutoff.
        else: self.generator.SetOptions ( cutoff = self.cObject.listCutoff )

    def CheckPairwiseInteractions ( self ):
        """Check the pairwise interactions."""
        # . Create suitable default interactions - these are the same as in the old GABFS module.
        # . Reasonable equivalents using only damping are achieved with values of 1.25 (qcmm) and 1.75 (qcqc) for default 8/12 inner/outer cutoffs.
        # . MM/MM.
        if self.mmmmPairwiseInteraction is None:
            self.mmmmPairwiseInteraction = PairwiseInteractionABFS.FromOptions ( dampingCutoff = self.cObject.dampingCutoff ,
                                                                                 innerCutoff   = self.cObject.innerCutoff   ,
                                                                                 outerCutoff   = self.cObject.outerCutoff   )
        self.mmmmPairwiseInteraction.MakeSplines ( )
        # . QC/MM.
        if self.qcmmPairwiseInteraction is None:
            self.qcmmPairwiseInteraction = PairwiseInteractionABFS.FromOptions ( dampingCutoff      = 0.0                        ,
                                                                                 electrostaticModel = "Delta/Gaussian"           ,
                                                                                 innerCutoff        = self.cObject.innerCutoff   ,
                                                                                 outerCutoff        = self.cObject.outerCutoff   ,
                                                                                 width2             = 2.0 * UNITS_LENGTH_BOHRS_TO_ANGSTROMS )
        self.qcmmPairwiseInteraction.MakeSplines ( lennardJones = False, useAtomicUnits = True )
        # . QC/QC.
        if self.qcqcPairwiseInteraction is None:
            self.qcqcPairwiseInteraction = PairwiseInteractionABFS.FromOptions ( dampingCutoff      = 0.0                        ,
                                                                                 electrostaticModel = "Gaussian/Gaussian"        ,
                                                                                 innerCutoff        = self.cObject.innerCutoff   ,
                                                                                 outerCutoff        = self.cObject.outerCutoff   ,
                                                                                 width1             = 2.0 * UNITS_LENGTH_BOHRS_TO_ANGSTROMS ,
                                                                                 width2             = 2.0 * UNITS_LENGTH_BOHRS_TO_ANGSTROMS )
        self.qcqcPairwiseInteraction.MakeSplines ( lennardJones = False, useAtomicUnits = True )

    def ClearPairwiseInteractions ( self ):
        """Clear the pairwise interactions."""
        for label in ( "mmmmPairwiseInteraction" ,
                       "qcmmPairwiseInteraction" ,
                       "qcqcPairwiseInteraction" ):
            setattr ( self, label, None )

    def Energy ( self, configuration ):
        """Energy and gradients."""
        cdef CNBModelABFSState *cnbState
        cdef  NBModelABFSState   nbState
        # . Initialization.
        energies = []
        # . Get the state object.
        if hasattr ( configuration, "nbState" ):
            nbState  = configuration.nbState
            cnbState = nbState.cObject
            # . Calculate the energies and gradients.
            NBModelABFS_MMMMEnergy   ( self.cObject, self.mmmmPairwiseInteraction.cObject, cnbState )
            NBModelABFS_QCMMEnergyLJ ( self.cObject, self.mmmmPairwiseInteraction.cObject, cnbState )
            nbState.GetEnergies ( energies )
        return energies

    def QCMMGradients ( self, configuration ):
        """Calculate the QC/MM electrostatic gradients."""
        cdef NBModelABFSState nbState
        # . Get the state object.
        if hasattr ( configuration, "nbState" ):
            nbState = configuration.nbState
            NBModelABFS_QCMMGradients ( self.cObject, self.qcmmPairwiseInteraction.cObject, self.qcqcPairwiseInteraction.cObject, nbState.cObject )

    def QCMMPotentials ( self, configuration ):
        """Calculate the QC/MM electrostatic potentials."""
        cdef NBModelABFSState nbState
        # . Get the state object.
        if hasattr ( configuration, "nbState" ):
            nbState = configuration.nbState
            NBModelABFS_QCMMPotentials ( self.cObject, self.qcmmPairwiseInteraction.cObject, self.qcqcPairwiseInteraction.cObject, nbState.cObject )

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        # . Simple options.
        if "dampingCutoff"           in keywordArguments: self.cObject.dampingCutoff        = keywordArguments.pop ( "dampingCutoff"           )
        if "dielectric"              in keywordArguments: self.cObject.dielectric           = keywordArguments.pop ( "dielectric"              )
        if "electrostaticScale14"    in keywordArguments: self.cObject.electrostaticScale14 = keywordArguments.pop ( "electrostaticScale14"    )
        if "generator"               in keywordArguments: self.generator                    = keywordArguments.pop ( "generator"               )
        if "imageExpandFactor"       in keywordArguments: self.cObject.imageExpandFactor    = keywordArguments.pop ( "imageExpandFactor"       )
        if "innerCutoff"             in keywordArguments: self.cObject.innerCutoff          = keywordArguments.pop ( "innerCutoff"             )
        if "listCutoff"              in keywordArguments: self.cObject.listCutoff           = keywordArguments.pop ( "listCutoff"              )
        if "outerCutoff"             in keywordArguments: self.cObject.outerCutoff          = keywordArguments.pop ( "outerCutoff"             )
        if "mmmmPairwiseInteraction" in keywordArguments: self.mmmmPairwiseInteraction      = keywordArguments.pop ( "mmmmPairwiseInteraction" )
        if "qcmmPairwiseInteraction" in keywordArguments: self.qcmmPairwiseInteraction      = keywordArguments.pop ( "qcmmPairwiseInteraction" )
        if "qcqcPairwiseInteraction" in keywordArguments: self.qcqcPairwiseInteraction      = keywordArguments.pop ( "qcqcPairwiseInteraction" )
        # . Booleans.
        if "checkForInverses"        in keywordArguments:
            if keywordArguments.pop ( "checkForInverses" ): self.cObject.checkForInverses = CTrue
            else:                                           self.cObject.checkForInverses = CFalse
        if "useCentering"            in keywordArguments:
            if keywordArguments.pop ( "useCentering"     ): self.cObject.useCentering     = CTrue
            else:                                           self.cObject.useCentering     = CFalse
        # . QC/MM coupling.
        if "qcmmCoupling" in keywordArguments:
            coupling = keywordArguments.pop ( "qcmmCoupling" )
            option   = QCMMLinkAtomCoupling_ToEnum.get ( coupling, None )
            if option is None: raise TypeError ( "Unrecognized QC/MM coupling option: " + coupling + "." )
            else:              self.cObject.qcmmCoupling = option
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )
        # . Do minimum checks on cutoff values.
        if ( self.cObject.dampingCutoff < 0.0                        ) or \
           ( self.cObject.innerCutoff   < self.cObject.dampingCutoff ) or \
           ( self.cObject.outerCutoff   < self.cObject.innerCutoff   ) or \
           ( self.cObject.listCutoff    < self.cObject.outerCutoff   ):
            raise TypeError ( "Invalid cutoff values: damping - {:.3f}; inner - {:.3f}; outer - {:.3f}; list - {:.3f}.".format ( self.cObject.dampingCutoff , \
                                                                                                                                 self.cObject.innerCutoff   , \
                                                                                                                                 self.cObject.outerCutoff   , \
                                                                                                                                 self.cObject.listCutoff    ) )
        # . Check the generator and pairwise interactions.
        self.CheckGenerator            ( )
        self.CheckPairwiseInteractions ( )

    def SetUp ( self, MMAtomContainer mmAtoms, QCAtomContainer qcAtoms, LJParameterContainer ljParameters, LJParameterContainer ljParameters14, Selection fixedAtoms,
                                                            SelfPairList interactions14, SelfPairList exclusions, symmetry, isolates, configuration, log = logFile ):
        """Set up the energy calculation."""
        cdef Coordinates3                 coordinates3
        cdef Coordinates3                 gradients3
        cdef NBModelABFSState             nbState
        cdef QCMMInteractionState         qcmmstate
        cdef Status                       status
        cdef SymmetryParameters           symmetryParameters
        cdef SymmetryParameterGradients   symmetryParameterGradients
        cdef Transformation3Container     transformations
        cdef CCoordinates3               *ccoordinates3, *cgradients3
        cdef CLJParameterContainer       *cljParameters, *cljParameters14
        cdef CPairList                   *cexclusions, *cinteractions14
        cdef CQCAtomContainer            *cqcAtoms
        cdef CReal1DArray                *cQCCharges, *cQCMMPotentials
        cdef CSelection                  *cfixedAtoms
        cdef CSymmetricMatrix            *cQCQCPotentials
        cdef CSymmetryParameters         *cSymmetryParameters
        cdef CSymmetryParameterGradients *cSymmetryParameterGradients
        cdef CTransformation3Container   *ctransformations
        # . There must be a configuration.
        if ( configuration is not None ):
            # . If there is an existing state it must (or should be) correct.
            # . qcmmstate should also be present in the case of QC/MM interactions.
            # . Create a new state.
            if not hasattr ( configuration, "nbState" ):
                # . Symmetry information.
                ctransformations = NULL
                if symmetry is not None:
                    if hasattr ( symmetry, "transformations" ):
                        transformations = symmetry.transformations
                        if transformations is not None: ctransformations = transformations.cObject
                # . Check for QC atoms.
                QQC = ( qcAtoms is not None ) and ( qcAtoms.size > 0 )
                # . Create a new QC/MM state if there are QC atoms.
                if QQC:
                    qcmmstate = QCMMInteractionState.WithExtent ( qcAtoms.size, includeQCQC = ( ctransformations != NULL ) )
                    setattr ( configuration, "qcmmstate", qcmmstate )
                    cQCCharges      = qcmmstate.cObject.qcCharges
                    cQCMMPotentials = qcmmstate.cObject.qcmmPotentials
                    cQCQCPotentials = qcmmstate.cObject.qcqcPotentials
                else:
                    cQCCharges      = NULL
                    cQCMMPotentials = NULL
                    cQCQCPotentials = NULL
                # . Get all data for the new state.
                if exclusions is None:     cexclusions     = NULL
                else:                      cexclusions     = exclusions.cObject
                if fixedAtoms is None:     cfixedAtoms     = NULL
                else:                      cfixedAtoms     = fixedAtoms.cObject
                if interactions14 is None: cinteractions14 = NULL
                else:                      cinteractions14 = interactions14.cObject
                if ljParameters   is None: cljParameters   = NULL
                else:                      cljParameters   = ljParameters.cObject
                if ljParameters14 is None: cljParameters14 = NULL
                else:                      cljParameters14 = ljParameters14.cObject
                if qcAtoms is None:        cqcAtoms        = NULL
                else:                      cqcAtoms        = qcAtoms.cObject
                # . Create the new state.
                nbState         = NBModelABFSState.Raw ( )
                nbState.cObject = NBModelABFSState_SetUp ( mmAtoms.cObject, cqcAtoms, cfixedAtoms, cexclusions, cinteractions14, cljParameters, cljParameters14,
                                                                     cQCCharges, cQCMMPotentials, cQCQCPotentials, ctransformations, self.cObject.qcmmCoupling )
                status          = Status_Continue
                NBModelABFSState_SetUpCentering ( nbState.cObject, self.cObject.useCentering, &status )
                nbState.isOwner = True
                setattr ( configuration, "nbState", nbState )
                if ( nbState.cObject == NULL ) or ( status != Status_Continue ): raise CLibraryError ( "Unable to create NB state." )
            # . Get the state.
            nbState = configuration.nbState
            # . Initialize the NB state (including gradients).
            if hasattr ( configuration, "coordinates3"               ): coordinates3               = configuration.coordinates3
            else:                                                       coordinates3               = None
            if hasattr ( configuration, "gradients3"                 ): gradients3                 = configuration.gradients3
            else:                                                       gradients3                 = None
            if hasattr ( configuration, "symmetryParameterGradients" ): symmetryParameterGradients = configuration.symmetryParameterGradients
            else:                                                       symmetryParameterGradients = None
            if hasattr ( configuration, "symmetryParameters"         ): symmetryParameters         = configuration.symmetryParameters
            else:                                                       symmetryParameters         = None
            if coordinates3               is None: ccoordinates3               = NULL
            else:                                  ccoordinates3               = coordinates3.cObject
            if gradients3                 is None: cgradients3                 = NULL
            else:                                  cgradients3                 = gradients3.cObject
            if symmetryParameterGradients is None: cSymmetryParameterGradients = NULL
            else:                                  cSymmetryParameterGradients = symmetryParameterGradients.cObject
            if symmetryParameters         is None: cSymmetryParameters         = NULL
            else:                                  cSymmetryParameters         = symmetryParameters.cObject
            NBModelABFSState_Initialize ( nbState.cObject, ccoordinates3, cSymmetryParameters, cgradients3, cSymmetryParameterGradients )
            # . Check for an update.
            status     = Status_Continue
            updateDone = ( NBModelABFS_Update ( self.cObject, self.generator.cObject, nbState.cObject, &status ) == CTrue )
            if status != Status_Continue: raise CLibraryError ( "Unable to create NB lists." )
            if updateDone: nbState.Summary ( log = log )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "ABFS NB Model Summary" )
            summary.Entry ( "Dielectric"      , "{:.6f}".format ( self.cObject.dielectric           ) )
            summary.Entry ( "El. 1-4 Scaling" , "{:.6f}".format ( self.cObject.electrostaticScale14 ) )
            summary.Entry ( "Damping Cutoff"  , "{:.6f}".format ( self.cObject.dampingCutoff        ) )
            summary.Entry ( "Inner Cutoff"    , "{:.6f}".format ( self.cObject.innerCutoff          ) )
            summary.Entry ( "List Cutoff"     , "{:.6f}".format ( self.cObject.listCutoff           ) )
            summary.Entry ( "Outer Cutoff"    , "{:.6f}".format ( self.cObject.outerCutoff          ) )
            summary.Entry ( "QC/MM Coupling"  ,  QCMMLinkAtomCoupling_ToString[self.cObject.qcmmCoupling] )
            summary.Entry ( "Use Centering"   , "{!r}".format ( self.cObject.useCentering == CTrue ) )
            summary.Stop ( )
            if self.generator is not None: self.generator.Summary ( log = log )
            # . Interactions?

    # . Properties - readonly.
    property dielectric:
        def __get__ ( self ): return self.cObject.dielectric
    property electrostaticScale14:
        def __get__ ( self ): return self.cObject.electrostaticScale14
    property innerCutoff:
        def __get__ ( self ): return self.cObject.innerCutoff
    property listCutoff:
        def __get__ ( self ): return self.cObject.listCutoff
    property outerCutoff:
        def __get__ ( self ): return self.cObject.outerCutoff
