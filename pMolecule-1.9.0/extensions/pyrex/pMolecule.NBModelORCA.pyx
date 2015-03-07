#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelORCA.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines a full NB model compatible with the ORCA program."""

from pCore import logFile, LogFileActive

from QCMMLinkAtomCouplingOptions import QCMMLinkAtomCoupling_ToEnum, QCMMLinkAtomCoupling_ToString

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelORCA ( NBModel ):
    """Defines a full NB model compatible with the ORCA program."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner: NBModelORCA_Deallocate ( &self.cObject )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModelORCA"

    def __getstate__ ( self ):
        """Return the state."""
        return { "electrostaticScale14" : self.cObject.electrostaticscale14 }

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = NBModelORCA_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.label   = "ORCA"

    def Energy ( self, configuration ):
        """Energy and gradients."""
        cdef CNBModelOrcaState *cnbState
        cdef NBModelORCAState   nbState
        # . Initialization.
        energies = []
        # . Get the state object.
        if hasattr ( configuration, "nbState" ):
            nbState  = configuration.nbState
            cnbState = nbState.cObject
            # . Calculate the energies and gradients.
            NBModelORCA_MMMMEnergy   ( self.cObject, cnbState )
            NBModelORCA_QCMMEnergyLJ ( self.cObject, cnbState )
            nbState.GetEnergies ( energies )
        return energies

    def QCMMGradients ( self, configuration ):
        """Calculate the QC/MM electrostatic gradients."""
        cdef CNBModelOrcaState *cnbState
        cdef NBModelORCAState        nbState
        # . Get the state object.
        if hasattr ( configuration, "nbState" ):
            nbState  = configuration.nbState
            cnbState = nbState.cObject
            # . Correct the MM gradients.
            NBModelORCAState_Finalize ( cnbState )

    def QCMMPotentials ( self, configuration ):
        """Calculate the QC/MM electrostatic potentials."""
        pass

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        if "electrostaticScale14" in keywordArguments: self.cObject.electrostaticscale14 = keywordArguments.pop ( "electrostaticScale14" )
        if "qcmmCoupling"         in keywordArguments:
            coupling = keywordArguments.pop ( "qcmmCoupling" )
            option   = QCMMLinkAtomCoupling_ToEnum.get ( coupling, None )
            if ( option is None ) or ( option == QCMMLinkAtomCoupling_MM ): raise TypeError ( "Invalid QC/MM coupling option: " + coupling + "." )
            else:                                                           self.cObject.qcmmcoupling = option
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )

    def SetUp ( self, MMAtomContainer mmAtoms, QCAtomContainer qcAtoms, LJParameterContainer ljParameters, LJParameterContainer ljParameters14, Selection fixedatoms, SelfPairList interactions14, SelfPairList exclusions, symmetry, isolates, configuration, log = logFile ):
        """Set up the energy calculation."""
        cdef Coordinates3           coordinates3, gradients3
        cdef NBModelORCAState       nbState
        cdef CCoordinates3         *ccoordinates3, *cgradients3
        cdef CLJParameterContainer *cljParameters, *cljParameters14
        cdef CPairList             *cexclusions, *cinteractions14
        cdef CQCAtomContainer      *cqcAtoms
        cdef CSelection            *cfixedatoms
        # . Symmetry handling.
        if symmetry is not None: raise TypeError ( "NBModelORCA class unable to handle symmetry." )
        # . There must be a configuration.
        if ( configuration is not None ):
            # . If there is an existing state it must (or should be) correct.
            # . Create a new state.
            if not hasattr ( configuration, "nbState" ):
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
                nbState         = NBModelORCAState.Raw ( )
                nbState.cObject = NBModelORCAState_Setup ( mmAtoms.cObject, cqcAtoms, cfixedatoms, cexclusions, cinteractions14, cljParameters, cljParameters14, self.cObject.qcmmcoupling )
                nbState.isOwner = True
                if nbState.cObject == NULL: raise ValueError ( "Error creating NB state." )
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
            NBModelORCAState_Initialize ( nbState.cObject, ccoordinates3, cgradients3 )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "ORCA NB Model Summary" )
            summary.Entry ( "El. 1-4 Scaling", "{:.6f}".format ( self.cObject.electrostaticscale14 )  )
            summary.Entry ( "QC/MM Coupling",  QCMMLinkAtomCoupling_ToString[self.cObject.qcmmcoupling] )
            summary.Stop ( )

    # . Properties - readonly.
    property electrostaticScale14:
        def __get__ ( self ): return self.cObject.electrostaticscale14
