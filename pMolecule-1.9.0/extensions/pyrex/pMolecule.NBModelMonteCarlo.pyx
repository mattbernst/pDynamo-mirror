#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModelMonteCarlo.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines a simple Monte Carlo NB model."""

from pCore import logFile, LogFileActive

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModelMonteCarlo ( NBModel ):
    """Define a Monte Carlo NB model."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            NBModelMonteCarlo_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.NBModelMonteCarlo"

    def __getstate__ ( self ):
        """Return the state."""
        return { "buffer"                 : self.cObject.buffer      ,
                 "cutoff"                 : self.cObject.cutoff      ,
                 "dielectric"             : self.cObject.dielectric  ,
                 "electrostaticUnderflow" : self.cObject.underflowel ,
                 "lennardJonesUnderflow"  : self.cObject.underflowlj }

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = NBModelMonteCarlo_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.label   = "Monte Carlo"

    def Energy ( self, configuration ):
        """Energy and gradients."""
        cdef CNBModelMonteCarloState *cnbState
        cdef  NBModelMonteCarloState   nbState
        # . Initialization.
        energies = []
        # . Get the state object.
        if hasattr ( configuration, "nbState" ):
            nbState  = configuration.nbState
            cnbState = nbState.cObject
            # . Calculate the energies.
            NBModelMonteCarlo_MMMMEnergyFull ( self.cObject, cnbState )
            nbState.GetEnergies ( energies )
        return energies

    def IsolateInteractionEnergy ( self, Integer isolate, NBModelMonteCarloState nbState ):
        """Calculate the interaction energy between an isolate and the rest."""
        if ( isolate < 0 ) or ( isolate >= nbState.cObject.nisolates ): raise IndexError ( "Isolate index ({:d}) out of range [0,{:d}).".format ( isolate, nbState.cObject.nisolates ) )
        NBModelMonteCarlo_MMMMEnergySingle ( self.cObject, isolate, nbState.cObject )
        return nbState.GetInteractionEnergies ( )

    def QCMMGradients ( self, configuration ):
        """Calculate the QC/MM electrostatic gradients."""
        raise ValueError ( "QC/MM gradients not available." )

    def QCMMPotentials ( self, configuration ):
        """Calculate the QC/MM electrostatic potentials."""
        raise ValueError ( "QC/MM interactions not available." )

    def SetOptions ( self, **keywordArguments ):
        """Set options for the model."""
        if "buffer"                 in keywordArguments: self.cObject.buffer      = keywordArguments.pop ( "buffer"                 )
        if "cutoff"                 in keywordArguments: self.cObject.cutoff      = keywordArguments.pop ( "cutoff"                 )
        if "dielectric"             in keywordArguments: self.cObject.dielectric  = keywordArguments.pop ( "dielectric"             )
        if "electrostaticUnderflow" in keywordArguments: self.cObject.underflowel = keywordArguments.pop ( "electrostaticUnderflow" )
        if "lennardJonesUnderflow"  in keywordArguments: self.cObject.underflowlj = keywordArguments.pop ( "lennardJonesUnderflow"  )
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )

    def SetUp ( self, MMAtomContainer mmAtoms, QCAtomContainer qcAtoms, LJParameterContainer ljParameters, LJParameterContainer ljParameters14, Selection fixedatoms, SelfPairList interactions14, SelfPairList exclusions, symmetry, isolates, configuration, log = logFile ):
        """Set up the energy calculation."""
        cdef NBModelMonteCarloState     nbState
        cdef SelectionContainer         pisolates
        cdef CLJParameterContainer *cljParameters
        cdef CSelectionContainer   *cisolates
        cdef CSelection            *cfixedatoms
        # . QC atoms.
        if qcAtoms is not None: raise TypeError ( "The Monte Carlo NB module cannot handle QC atoms." )
        # . There must be a configuration.
        if ( configuration is not None ):
            # . Check for gradients.
            if hasattr ( configuration, "gradients3" ) and ( configuration.gradients3 is not None ): raise TypeError ( "The Monte Carlo NB module cannot calculate derivatives." )
            # . Check for symmetry.
            if ( symmetry is None ) or not symmetry.crystalClass.IsOrthogonal ( ): raise TypeError ( "The Monte Carlo NB module only works with crystal classes that are orthogonal." )
            # . If there is an existing state it must (or should be) correct.
            nbState = getattr3 ( configuration, "nbState", None )
            # . Create a new state.
            if nbState is None:
                # . Isolates.
                if isinstance ( isolates, SelectionContainer ):
                    pisolates = isolates
                    cisolates = pisolates.cObject
                else:
                    cisolates = NULL
                # . Other data.
                if fixedatoms     is None: cfixedatoms   = NULL
                else:                      cfixedatoms   = fixedatoms.cObject
                if ljParameters   is None: cljParameters = NULL
                else:                      cljParameters = ljParameters.cObject
                # . Create the state.
                nbState         = NBModelMonteCarloState.Raw ( )
                nbState.cObject = NBModelMonteCarloState_Setup ( cisolates, mmAtoms.cObject, cljParameters, cfixedatoms )
                nbState.isOwner = True
                setattr ( configuration, "nbState", nbState )
                # . Summary.
                nbState.Summary ( log = log )
            # . Get the state.
            nbState = configuration.nbState
            # . Initialize the NB state.
            nbState.Initialize ( configuration )
            # . Check that the minimum image convention is satisfied.
            if ( SymmetryParameters_IsMinimumImageConventionSatisfied ( nbState.cObject.symmetryParameters, self.cObject.cutoff ) == CFalse ): raise ValueError ( "The Monte Carlo NB module requires that the minimum image convention be valid." )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( )
            summary.Start ( "Monte Carlo NB Model Summary" )
            summary.Entry ( "Buffer",     "{:.6f}".format ( self.cObject.buffer     , ) )
            summary.Entry ( "Cutoff",     "{:.6f}".format ( self.cObject.cutoff     , ) )
            summary.Entry ( "Dielectric", "{:.6f}".format ( self.cObject.dielectric , ) )
            summary.Stop ( )

    # . Properties - readonly.
    property buffer:
        def __get__ ( self ): return self.cObject.buffer
    property cutoff:
        def __get__ ( self ): return self.cObject.cutoff
    property dielectric:
        def __get__ ( self ): return self.cObject.dielectric
