#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelSSBP.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the SSBP NB model."""

from pCore import logFile, LogFileActive

#from QCMMLinkAtomCouplingOptions import QCMMLinkAtomCoupling_ToEnum, QCMMLinkAtomCoupling_ToString

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelSSBP ( NBModelFull ):
    """Define a SSBP NB model."""

    def __dealloc__ ( self ):
        """Finalization."""
# . Why doesn't this work?
#        super ( NBModelSSBP, self ).__dealloc__ ( )
        if self.isOwner:
            NBModelFull_Deallocate ( &self.cObject )
            self.isOwner = False
        if self.isOwner2:
            SSBPModel_Deallocate ( &self.cObject2 )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModelSSBP"

    def __getstate__ ( self ):
        """Return the state."""
        state = super ( NBModelSSBP, self ).__getstate__ ( )
        state.update ( { "doAngularPotential"      : ( self.cObject2.doAngularPotential      == CTrue ) ,
                         "doCavityPotential"       : ( self.cObject2.doCavityPotential       == CTrue ) ,
                         "doEmpiricalCorrection"   : ( self.cObject2.doEmpiricalCorrection   == CTrue ) ,
                         "doHardSphereRestriction" : ( self.cObject2.doHardSphereRestriction == CTrue ) ,
                         "doKirkwood"              : ( self.cObject2.doKirkwood              == CTrue ) ,
                         "doNonBondingTerms"       :   self.doNonBondingTerms                           ,
                         "fixCavityRadius"         : ( self.cObject2.fixCavityRadius         == CTrue ) ,
                         "maximumL"                :   self.cObject2.maximumL                           ,
                         "cavityRadius"            :   self.cObject2.cavityRadius                       ,
                         "cavityRadiusIncrement"   :   self.cObject2.cavityRadiusIncrement              ,
                         "dielectricInside"        :   self.cObject2.dielectricInside                   ,
                         "dielectricOutside"       :   self.cObject2.dielectricOutside                  ,
                         "empirical1"              :   self.cObject2.empirical1                         ,
                         "empirical2"              :   self.cObject2.empirical2                         ,
                         "kirkwoodRadiusIncrement" :   self.cObject2.kirkwoodRadiusIncrement            ,
                         "pressure"                :   self.cObject2.pressure                           ,
                         "surfaceTension"          :   self.cObject2.surfaceTension                   } )
        return state

    def _Allocate ( self ):
        """Allocation."""
        super ( NBModelSSBP, self )._Allocate ( )
        self.cObject2 = SSBPModel_Allocate ( )
        self.isOwner2 = True

    def _Initialize ( self ):
        """Initialization."""
        super ( NBModelSSBP, self )._Initialize ( )
        self.cObject2          = NULL
        self.isOwner2          = False
        self.label             = "SSBP"
        self.doNonBondingTerms = True

    def Energy ( self, configuration ):
        """Energy and gradients.

        This method also does the work of QCMMPotentials when there are QC atoms.
        """
        cdef CSSBPModelState *cState
# # ifdef SSBPTESTQC
#        cdef Real             energy
# # endif
        cdef SSBPModelState   state
        # . Initialization.
        if self.doNonBondingTerms: energies = super ( NBModelSSBP, self ).Energy ( configuration )
        else:                      energies = []
        # . Get the state object.
        if hasattr ( configuration, "ssbpState" ):
            state  = configuration.ssbpState
            cState = state.cObject
            # . Determine gradient flags.
            doGradients = ( cState.gradients3 != NULL )
            doQC        = ( cState.qcAtoms    != NULL )
            # . Calculate the energies and gradients (MM case only).
            if doQC or ( not doGradients ): SSBPModel_Energy ( self.cObject2, CFalse, cState )
            else:                           SSBPModel_Energy ( self.cObject2, CTrue , cState )
            state.GetEnergies ( energies )
# # ifdef SSBPTESTQC
#            # . Gradients for QC case.
#            if doQC and doGradients:
#                energy = SSBPModel_Kirkwood ( self.cObject2, CTrue, cState )
#                print ( "\nKirkwood Energy = {:20.5f}\n".format ( energy ) )
# # endif
        return energies

# . The following is not needed if SSBPTESTQC

    def QCMMGradients ( self, configuration ):
        """Calculate the QC/MM electrostatic gradients.

        This method also calculates the MM gradient terms when there are QC atoms.
        """
        cdef CSSBPModelState *cState
        cdef SSBPModelState   state
        # . Initialization.
        if self.doNonBondingTerms:
            super ( NBModelSSBP, self ).QCMMGradients ( configuration )
        # . Get the state object.
        if hasattr ( configuration, "ssbpState" ):
            state  = configuration.ssbpState
            cState = state.cObject
            SSBPModel_Kirkwood ( self.cObject2, CTrue, cState )

    def QCMMPotentials ( self, configuration ):
        """Calculate the QC/MM electrostatic potentials."""
        if self.doNonBondingTerms:
            super ( NBModelSSBP, self ).QCMMPotentials ( configuration )

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        keywordArguments["skipInvalidOptionCheck"] = True
        keywordArguments = super ( NBModelSSBP, self ).SetOptions ( **keywordArguments )
        if "doAngularPotential"      in keywordArguments:
            if keywordArguments.pop ( "doAngularPotential"      ): self.cObject2.doAngularPotential      = CTrue
            else:                                                  self.cObject2.doAngularPotential      = CFalse
        if "doCavityPotential"       in keywordArguments:
            if keywordArguments.pop ( "doCavityPotential"       ): self.cObject2.doCavityPotential       = CTrue
            else:                                                  self.cObject2.doCavityPotential       = CFalse
        if "doEmpiricalCorrection"   in keywordArguments:
            if keywordArguments.pop ( "doEmpiricalCorrection"   ): self.cObject2.doEmpiricalCorrection   = CTrue
            else:                                                  self.cObject2.doEmpiricalCorrection   = CFalse
        if "doHardSphereRestriction" in keywordArguments:
            if keywordArguments.pop ( "doHardSphereRestriction" ): self.cObject2.doHardSphereRestriction = CTrue
            else:                                                  self.cObject2.doHardSphereRestriction = CFalse
        if "doKirkwood"              in keywordArguments:
            if keywordArguments.pop ( "doKirkwood"              ): self.cObject2.doKirkwood              = CTrue
            else:                                                  self.cObject2.doKirkwood              = CFalse
        if "doNonBondingTerms"       in keywordArguments:
            if keywordArguments.pop ( "doNonBondingTerms"       ): self.doNonBondingTerms                = True
            else:                                                  self.doNonBondingTerms                = False
        if "fixCavityRadius"         in keywordArguments:
            if keywordArguments.pop ( "fixCavityRadius"         ): self.cObject2.fixCavityRadius         = CTrue
            else:                                                  self.cObject2.fixCavityRadius         = CFalse
        if "maximumL"                in keywordArguments: self.cObject2.maximumL                = keywordArguments.pop ( "maximumL"                )
        if "cavityRadius"            in keywordArguments: self.cObject2.cavityRadius            = keywordArguments.pop ( "cavityRadius"            )
        if "cavityRadiusIncrement"   in keywordArguments: self.cObject2.cavityRadiusIncrement   = keywordArguments.pop ( "cavityRadiusIncrement"   )
        if "dielectricInside"        in keywordArguments: self.cObject2.dielectricInside        = keywordArguments.pop ( "dielectricInside"        )
        if "dielectricOutside"       in keywordArguments: self.cObject2.dielectricOutside       = keywordArguments.pop ( "dielectricOutside"       )
        if "empirical1"              in keywordArguments: self.cObject2.empirical1              = keywordArguments.pop ( "empirical1"              )
        if "empirical2"              in keywordArguments: self.cObject2.empirical2              = keywordArguments.pop ( "empirical2"              )
        if "kirkwoodRadiusIncrement" in keywordArguments: self.cObject2.kirkwoodRadiusIncrement = keywordArguments.pop ( "kirkwoodRadiusIncrement" )
        if "pressure"                in keywordArguments: self.cObject2.pressure                = keywordArguments.pop ( "pressure"                )
        if "surfaceTension"          in keywordArguments: self.cObject2.surfaceTension          = keywordArguments.pop ( "surfaceTension"          )
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )

    def SetUp ( self, MMAtomContainer mmAtoms, QCAtomContainer qcAtoms, LJParameterContainer ljParameters, \
                   LJParameterContainer ljParameters14, Selection fixedatoms, SelfPairList interactions14, \
                                SelfPairList exclusions, symmetry, isolates, configuration, log = logFile ):
        """Set up the energy calculation."""
#        cdef Boolean              doKirkwood
        cdef Coordinates3         coordinates3, gradients3
        cdef Integer2DArray       waterAtomIndices
        cdef QCMMInteractionState qcmmState
        cdef Selection            cavitySelection, radiusSelection
        cdef SSBPModelState       state
        cdef CCoordinates3       *ccoordinates3, *cgradients3
        cdef CInteger2DArray     *cWaterAtomIndices
        cdef CQCAtomContainer    *cQCAtoms
        cdef CReal1DArray        *cQCCharges, *cQCMMPotentials
        cdef CSelection          *cCavitySelection, *cRadiusSelection
        cdef CSymmetricMatrix    *cQCQCPotentials
        # . NB state.
        if self.doNonBondingTerms:
            super ( NBModelSSBP, self ).SetUp ( mmAtoms, qcAtoms, ljParameters, ljParameters14, fixedatoms, interactions14, \
                                                                  exclusions, symmetry, isolates, configuration, log = None )
        # . There must be a configuration.
        if ( configuration is not None ):
            # . If there is an existing state it must (or should be) correct.
            if not hasattr ( configuration, "ssbpState" ):
                # . Get selections.
                cavitySelection  = getattr3 ( configuration, "cavitySelection" , None )
                radiusSelection  = getattr3 ( configuration, "radiusSelection" , None )
                waterAtomIndices = getattr3 ( configuration, "waterAtomIndices", None )
                if cavitySelection  is None: cCavitySelection  = NULL
                else:                        cCavitySelection  = Selection_Clone      ( cavitySelection.cObject  )
                if radiusSelection  is None: cRadiusSelection  = NULL
                else:                        cRadiusSelection  = Selection_Clone      ( radiusSelection.cObject  )
                if waterAtomIndices is None: cWaterAtomIndices = NULL
                else:                        cWaterAtomIndices = Integer2DArray_Clone ( waterAtomIndices.cObject, NULL )
                # . An appropriate QC/MM interaction state should have been set up by
                # . NBModelFull but without qcqcPotentials which needs to be allocated.
                QQC = ( qcAtoms is not None ) and ( qcAtoms.size > 0 )
                # . Create a new QC/MM state if there are QC atoms.
                if QQC:
                    # . An old state exists.
                    if self.doNonBondingTerms:
                        qcmmState = configuration.qcmmstate
                        QCMMInteractionState_AllocateQCQCPotentials ( qcmmState.cObject, NULL )
                    # . Create a new state.
                    else:
                        qcmmState = QCMMInteractionState.WithExtent ( qcAtoms.size, includeQCQC = True )
                        setattr ( configuration, "qcmmstate", qcmmState )
                    cQCAtoms        = qcAtoms.cObject
                    cQCCharges      = qcmmState.cObject.qcCharges
                    cQCMMPotentials = qcmmState.cObject.qcmmPotentials
                    cQCQCPotentials = qcmmState.cObject.qcqcPotentials
                else:
                    cQCAtoms        = NULL
                    cQCCharges      = NULL
                    cQCMMPotentials = NULL
                    cQCQCPotentials = NULL
                # . Create a new SSBP state.
# . Pyrex doesn't seem to like this.
#                if self.cObject2.doKirkwood == CTrue: doKirkword = CTrue
#                else:                                 doKirkwood = CFalse
                state = SSBPModelState.Raw ( )
                if self.cObject2.doKirkwood == CTrue:
                    state.cObject = SSBPModelState_Setup ( CTrue , self.cObject2.maximumL, mmAtoms.cObject, cQCAtoms,
                                                               cCavitySelection, cRadiusSelection, cWaterAtomIndices,
                                                                 cQCCharges, cQCMMPotentials, cQCQCPotentials, NULL )
                else:
                    state.cObject = SSBPModelState_Setup ( CFalse, self.cObject2.maximumL, mmAtoms.cObject, cQCAtoms,
                                                               cCavitySelection, cRadiusSelection, cWaterAtomIndices,
                                                                 cQCCharges, cQCMMPotentials, cQCQCPotentials, NULL )
                state.isOwner = True
                setattr ( configuration, "ssbpState", state )
                # . Summary.
                state.Summary ( log = log )
            # . Get the state.
            state = configuration.ssbpState
            # . Initialize the NB state (including gradients).
            coordinates3 = getattr3 ( configuration, "coordinates3", None )
            gradients3   = getattr3 ( configuration, "gradients3"  , None )
            if coordinates3 is None: ccoordinates3 = NULL
            else:                    ccoordinates3 = coordinates3.cObject
            if gradients3   is None: cgradients3   = NULL
            else:                    cgradients3   = gradients3.cObject
            SSBPModelState_Initialize ( state.cObject, ccoordinates3, cgradients3 )

    def Summary ( self, log = logFile ):
        """Summary."""
        super ( NBModelSSBP, self ).Summary ( log = log )
        if LogFileActive ( log ):
            summary = log.GetSummary ( pageWidth = 90 )
            summary.Start ( "SSBP Model Summary" )
            summary.Entry ( "Angular Potential"      , "{!r}".format ( self.cObject2.doAngularPotential      == CTrue ) )
            summary.Entry ( "Cavity Potential"       , "{!r}".format ( self.cObject2.doCavityPotential       == CTrue ) )
            summary.Entry ( "Empirical Correction"   , "{!r}".format ( self.cObject2.doEmpiricalCorrection   == CTrue ) )
            summary.Entry ( "Hard Sphere Restriction", "{!r}".format ( self.cObject2.doHardSphereRestriction == CTrue ) )
            summary.Entry ( "Kirkwood"               , "{!r}".format ( self.cObject2.doKirkwood              == CTrue ) )
            summary.Entry ( "Non-Bonding Terms"      , "{!r}".format ( self.doNonBondingTerms                         ) )
            summary.Entry ( "Fix Cavity Radius"      , "{!r}".format ( self.cObject2.fixCavityRadius         == CTrue ) )
            if self.cObject2.doCavityPotential == CTrue:
                summary.Entry ( "Cavity Radius"            , "{:14.4g}".format ( self.cObject2.cavityRadius            ) )
                summary.Entry ( "Cavity Radius Increment"  , "{:14.4g}".format ( self.cObject2.cavityRadiusIncrement   ) )
            if self.cObject2.doEmpiricalCorrection == CTrue:
                summary.Entry ( "Empirical Parameter 1"    , "{:14.4g}".format ( self.cObject2.empirical1              ) )
                summary.Entry ( "Empirical Parameter 2"    , "{:14.4g}".format ( self.cObject2.empirical2              ) )
            if self.cObject2.doHardSphereRestriction == CTrue:
                summary.Entry ( "Pressure"                 , "{:14.4g}".format ( self.cObject2.pressure                ) )
                summary.Entry ( "Surface Tension"          , "{:14.4g}".format ( self.cObject2.surfaceTension          ) )
            if self.cObject2.doKirkwood == CTrue:
                summary.Entry ( "Dielectric Inside"        , "{:14.4g}".format ( self.cObject2.dielectricInside        ) )
                summary.Entry ( "Dielectric Outside"       , "{:14.4g}".format ( self.cObject2.dielectricOutside       ) )
                summary.Entry ( "Maximum L"                , "{:d}"    .format ( self.cObject2.maximumL                ) )
                summary.Entry ( "Kirkwood Radius Increment", "{:14.4g}".format ( self.cObject2.kirkwoodRadiusIncrement ) )
            summary.Stop ( )
