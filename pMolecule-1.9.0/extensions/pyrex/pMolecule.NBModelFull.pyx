#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelFull.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines a simple full NB model."""

from pCore import logFile, LogFileActive

from QCMMLinkAtomCouplingOptions import QCMMLinkAtomCoupling_ToEnum, QCMMLinkAtomCoupling_ToString

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelFull ( NBModel ):
    """Define a full NB model."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            NBModelFull_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModelFull"

    def __getstate__ ( self ):
        """Return the state."""
        return { "dielectric"           : self.cObject.dielectric           ,
                 "electrostaticScale14" : self.cObject.electrostaticscale14 }

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = NBModelFull_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.label   = "Full"

    def Energy ( self, configuration ):
        """Energy and gradients."""
        cdef CNBModelFullState *cnbState
        cdef NBModelFullState   nbState
        # . Initialization.
        energies = []
        # . Get the state object.
        if hasattr ( configuration, "nbState" ):
            nbState  = configuration.nbState
            cnbState = nbState.cObject
            # . Calculate the energies and gradients.
            NBModelFull_MMMMEnergy   ( self.cObject, cnbState )
            NBModelFull_QCMMEnergyLJ ( self.cObject, cnbState )
            nbState.GetEnergies ( energies )
        return energies

    def QCMMGradients ( self, configuration ):
        """Calculate the QC/MM electrostatic gradients."""
        cdef NBModelFullState nbState
        # . Get the state object.
        if hasattr ( configuration, "nbState" ):
            nbState = configuration.nbState
            NBModelFull_QCMMGradients ( self.cObject, nbState.cObject )

    def QCMMPotentials ( self, configuration ):
        """Calculate the QC/MM electrostatic potentials."""
        cdef NBModelFullState nbState
        # . Get the state object.
        if hasattr ( configuration, "nbState" ):
            nbState = configuration.nbState
            NBModelFull_QCMMPotentials ( self.cObject, nbState.cObject )

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        if "dielectric"           in keywordArguments: self.cObject.dielectric           = keywordArguments.pop ( "dielectric" )
        if "electrostaticScale14" in keywordArguments: self.cObject.electrostaticscale14 = keywordArguments.pop ( "electrostaticScale14" )
        if "qcmmCoupling"         in keywordArguments:
            coupling = keywordArguments.pop ( "qcmmCoupling" )
            option   = QCMMLinkAtomCoupling_ToEnum.get ( coupling, None )
            if option is None: raise TypeError ( "Unrecognized QC/MM coupling option: " + coupling + "." )
            else:              self.cObject.qcmmcoupling = option
        if not keywordArguments.pop ( "skipInvalidOptionCheck", False ):
            if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )
        return keywordArguments

    def SetUp ( self, MMAtomContainer mmAtoms, QCAtomContainer qcAtoms, LJParameterContainer ljParameters, LJParameterContainer ljParameters14, Selection fixedatoms, SelfPairList interactions14, SelfPairList exclusions, symmetry, isolates, configuration, log = logFile ):
        """Set up the energy calculation."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3             gradients3
        cdef NBModelFullState       nbState
        cdef QCMMInteractionState   qcmmstate
        cdef CCoordinates3         *ccoordinates3, *cgradients3
        cdef CLJParameterContainer *cljParameters, *cljParameters14
        cdef CPairList             *cexclusions, *cinteractions14
        cdef CQCAtomContainer      *cqcAtoms
        cdef CReal1DArray          *cQCCharges, *cQCMMPotentials
        cdef CSelection            *cfixedatoms
        # . Symmetry handling.
        if symmetry is not None: raise TypeError ( "NBModelFull class unable to handle symmetry." )
        # . There must be a configuration.
        if ( configuration is not None ):
            # . If there is an existing state it must (or should be) correct.
            # . qcmmstate should also be present in the case of QC/MM interactions.
            # . Create a new state.
            if not hasattr ( configuration, "nbState" ):
                # . Check for QC atoms.
                QQC = ( qcAtoms is not None ) and ( qcAtoms.size > 0 )
                # . Create a new QC/MM state if there are QC atoms.
                if QQC:
                    qcmmstate = QCMMInteractionState.WithExtent ( qcAtoms.size )
                    setattr ( configuration, "qcmmstate", qcmmstate )
                    cQCCharges      = qcmmstate.cObject.qcCharges
                    cQCMMPotentials = qcmmstate.cObject.qcmmPotentials
                else:
                    cQCCharges      = NULL
                    cQCMMPotentials = NULL
                # . Create a new NB state.
                if exclusions is None:     cexclusions     = NULL
                else:                      cexclusions     = exclusions.cObject
                if fixedatoms is None:     cfixedatoms     = NULL
                else:                      cfixedatoms     = fixedatoms.cObject
                if interactions14 is None: cinteractions14 = NULL
                else:                      cinteractions14 = interactions14.cObject
                if ljParameters   is None: cljParameters   = NULL
                else:                      cljParameters   = ljParameters.cObject
                if ljParameters14 is None: cljParameters14 = NULL
                else:                      cljParameters14 = ljParameters14.cObject
                if qcAtoms is None:        cqcAtoms        = NULL
                else:                      cqcAtoms        = qcAtoms.cObject
                nbState         = NBModelFullState.Raw ( )
                nbState.cObject = NBModelFullState_Setup ( mmAtoms.cObject, cqcAtoms, cfixedatoms, cexclusions, cinteractions14, cljParameters, cljParameters14, cQCCharges, cQCMMPotentials )
                nbState.isOwner = True
                setattr ( configuration, "nbState", nbState )
                # . Summary.
                nbState.Summary ( log = log )
            # . Get the state.
            nbState = configuration.nbState
            # . Initialize the NB state (including gradients).
            if hasattr ( configuration, "coordinates3" ): coordinates3 = configuration.coordinates3
            else:                                         coordinates3 = None
            if hasattr ( configuration, "gradients3"   ): gradients3   = configuration.gradients3
            else:                                         gradients3   = None
            if coordinates3 is None: ccoordinates3 = NULL
            else:                    ccoordinates3 = coordinates3.cObject
            if gradients3   is None: cgradients3   = NULL
            else:                    cgradients3   = gradients3.cObject
            NBModelFullState_Initialize ( nbState.cObject, ccoordinates3, cgradients3 )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( valueWidth = 11 )
            summary.Start ( "Full NB Model Summary" )
            summary.Entry ( "Dielectric",      "{:.6f}".format ( self.cObject.dielectric           )  )
            summary.Entry ( "El. 1-4 Scaling", "{:.6f}".format ( self.cObject.electrostaticscale14 )  )
            summary.Entry ( "QC/MM Coupling", QCMMLinkAtomCoupling_ToString[self.cObject.qcmmcoupling] )
            summary.Stop ( )

    # . Properties - readonly.
    property dielectric:
        def __get__ ( self ): return self.cObject.dielectric
    property electrostaticScale14:
        def __get__ ( self ): return self.cObject.electrostaticscale14
