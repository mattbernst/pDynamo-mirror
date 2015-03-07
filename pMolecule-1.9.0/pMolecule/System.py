#-------------------------------------------------------------------------------
# . File      : System.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""The System class is the essential class for simulating molecular systems.
It combines information about the atomic composition of the system along with
the ways in which its potential energy can be calculated."""

from pCore                      import Clone, Coordinates3, logFile, LogFileActive, RawObjectConstructor, Real1DArray, \
                                       Selection, Transformation3Container, Transformation3Container_Identity, Vector3
from Atom                       import AtomContainer
from ChargeConstraintContainer  import ChargeConstraintContainer
from Configuration              import Configuration
from Connectivity               import Connectivity
from CrystalClass               import CrystalClass
from ElectronicState            import ElectronicState
from EnergyModel                import EnergyModel
from EnergyTerms                import EnergyTerms
from HardConstraintContainer    import HardConstraintContainer
from MMModel                    import MMModel
from NBModel                    import NBModel
from QCAtomContainer            import QCAtomContainer_FromAtomContainer
from QCModel                    import QCModel
from QCParameters               import QCParameters_Define
from Sequence                   import Sequence
from SoftConstraintContainer    import SoftConstraintContainer
from Symmetry                   import Symmetry
from SymmetryParameters         import SymmetryParameters
from SymmetryParameterGradients import SymmetryParameterGradients

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class System ( object ):
    """Define an atomic or molecular system and its energy models."""

    # . Default attributes.
    defaultAttributes =  { "_atoms"           : None ,
                           "_configuration"   : None ,
                           "_connectivity"    : None ,
                           "_electronicState" : None ,
                           "_energyModel"     : None ,
                           "_hardConstraints" : None ,
                           "_label"           : None ,
                           "_sequence"        : None ,
                           "_symmetry"        : None }

    # . Methods.
    def __getstate__ ( self ):
        """Get the state of the object."""
        # . Initialization.
        state = {}
        # . Sequence.
        if self.sequence is not None: state["sequence"] = self.sequence.ToMapping ( )
        # . Connectivity.
        state.update ( self.connectivity.ToMapping ( ) )
        # . Other data.
        for key in self.__class__.defaultAttributes:
            label = key[1:]
            if label not in ( "atoms", "connectivity", "sequence" ):
                value = self.__dict__.get ( key, None )
                if value is not None: state[label] = value
        return state

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state of the object."""
        # . Sequence is present.
        if "sequence" in state:
            sequence = Sequence.FromMapping ( state["sequence"] )
            atoms    = sequence.GatherAtoms ( )
            self.__dict__["_atoms"   ] = AtomContainer.FromIterable ( atoms, attributes = state.get ( "atoms", None ) )
            self.__dict__["_sequence"] = sequence
        # . Sequence is absent.
        else:
            self.__dict__["_atoms"   ] = AtomContainer.FromMapping ( state["atoms"] )
        # . Connectivity.
        self.__dict__["_connectivity"] = Connectivity.FromAtomContainer ( self.atoms, bonds    = state.get ( "bonds"   , None ) , 
                                                                                      ringSets = state.get ( "ringSets", None ) )
        # . Other attributes.
        for key in self.__class__.defaultAttributes:
            label = key[1:]
            if label not in ( "atoms", "connectivity", "sequence" ):
                value = state.get ( label, None )
                if value is not None: self.__dict__[key] = value

    def AtomicCharges ( self, qcChargeModel = None, qcAtomsOnly = False, spinDensities = False ):
        """Calculate the partial charges for the atoms."""
        charges = None
        if self.energyModel is not None:
            em = self.energyModel
            if ( em.mmModel is not None ) and ( not qcAtomsOnly ):
                if spinDensities:
                    charges = Real1DArray.WithExtent ( em.mmAtoms.size )
                    charges.Set ( 0.0e+00 )
                else:
                    charges = em.mmAtoms.AtomicCharges ( )
            if em.qcModel is not None:
                qccharges = em.qcModel.AtomicCharges ( self.configuration, chargeModel = qcChargeModel, spinDensities = spinDensities )
                if charges is None:
                    charges = qccharges
                else:
                    for ( index, q ) in zip ( em.qcAtoms, qccharges ): charges[index] += q
        return charges

    def BondsFromCoordinates3 ( self, coordinates = None, radii = None, safety = 0.45 ):
        """Estimate bonds from coordinates using a distance search."""
        if coordinates is None:
            if hasattr ( self, "coordinates3" ): coordinates = self.coordinates3
        if radii is None: radii = self.atoms.GetItemAttributes ( "covalentRadius" )
        if self.connectivity is None: self.__dict__["_connectivity"] = Connectivity.FromAtomContainer ( self.atoms )
        self.connectivity.BondsFromCoordinates ( coordinates, radii, safety )

    def DefineChargeConstraints ( self, constraintList ):
        """Define charge constraints for the system."""
        # . Reset the hard constraints.
        hc = self.hardConstraints
        if hc is not None: hc.ClearChargeConstraints ( )
        # . Get the constraints.
        constraints = None
        if constraintList is not None:
            constraints = ChargeConstraintContainer.FromList ( constraintList )
        # . Set the constraints.
        if constraints is None:
            if ( hc is not None ) and hc.IsEmpty ( ): hc = None
        else:
            if constraints.HighestIndex ( ) >= len ( self.atoms ): raise ValueError ( "Invalid atom index in charge constraint." )
            if hc is None: hc = HardConstraintContainer ( )
            hc.DefineChargeConstraints ( constraints )
        # . Finish up.
        self.__dict__["_hardConstraints"] = hc

    def DefineFixedAtoms ( self, selection ):
        """Define fixed atoms for a system."""
        # . Reset the energy model.
        em = self.energyModel
        if em is not None: em.ClearFixedAtoms ( self.configuration )
        # . Reset the hard constraints.
        hc = self.hardConstraints
        if hc is not None: hc.ClearFixedAtoms ( )
        # . Treat the argument.
        if selection is None:
            if ( hc is not None ) and hc.IsEmpty ( ): hc = None
        elif isinstance ( selection, Selection ):
            if len ( selection ) > 0:
                if hc is None: hc = HardConstraintContainer ( )
                hc.DefineFixedAtoms ( selection )
                if em is not None: em.DeactivateFixedAtomMMTerms ( selection )
            else:
                if ( hc is not None ) and hc.IsEmpty ( ): hc = None
        else:
            raise TypeError ( "Invalid argument - None or a selection required." )
        # . Finish up.
        self.__dict__["_hardConstraints"] = hc

    def DefineMMModel ( self, mmModel, buildModel = True, log = logFile ):
        """Define a molecular mechanical model for the system."""
        if isinstance ( mmModel, MMModel ):
            # . Get an energy model in the correct state.
            if self.energyModel is None: self.__dict__["_energyModel"] = EnergyModel ( )
            else:                        self.energyModel.ClearMMModel ( self.configuration )
            em = self.energyModel
            # . Define the MM model and initialize mmTerms.
            em.mmModel = mmModel
            em.MakeLabel ( )
            # . Construct the MM data structures.
            if buildModel: mmModel.BuildModel ( self.connectivity, em, self.sequence, log = log )
            # . Deactivate MM terms if fixed atoms are defined.
            if ( self.hardConstraints is not None ): em.DeactivateFixedAtomMMTerms ( self.hardConstraints.fixedAtoms )
            # . Check the total MM charge.
            em.CheckActiveMMAtomTotalCharge ( )

    def DefineNBModel ( self, nbModel ):
        """Define a non-bonding model for the system."""
        if isinstance ( nbModel, NBModel ):
            # . Get an energy model in the correct state.
            if self.energyModel is None: self.__dict__["_energyModel"] = EnergyModel ( )
            else:                        self.energyModel.ClearNBModel ( self.configuration )
            em = self.energyModel
            # . Save the new model.
            em.nbModel = nbModel
            # . Check for an MM model and make sure some options are correctly set.
            # . Need something more general here (method in NBModel class - searches attributes for matches).
            if em.mmModel is not None:
                if hasattr ( em.mmModel, "electrostaticScale14" ):
                    try:    em.nbModel.SetOptions ( electrostaticScale14 = float ( em.mmModel.__dict__["electrostaticScale14"] ) )
                    except: pass
            em.MakeLabel ( )

    def DefineQCModel ( self, qcModel, log = logFile, qcSelection = None ):
        """Define a quantum chemical model for the system."""
        if isinstance ( qcModel, QCModel ):
            # . Get an energy model in the correct state.
            if self.energyModel is None: self.__dict__["_energyModel"] = EnergyModel ( )
            else:                        self.energyModel.ClearQCModel ( self.configuration, fixedAtoms = None )
            em = self.energyModel
            # . Assign a default electronic state if necessary.
            if self.electronicState is None: self.electronicState = ElectronicState ( )
            # . Get the QC atom selection and the boundary atoms.
            if ( qcSelection is None ) or ( qcSelection.size == len ( self.atoms ) ):
                atomindices   = Selection.FromIterable ( range ( len ( self.atoms ) ) )
                boundaryatoms = {}
            else:
                atomindices   = qcSelection
                boundaryatoms = em.IdentifyBoundaryAtoms ( atomindices )
            # . Set up the QC atom data.
            qcAtoms = QCAtomContainer_FromAtomContainer ( self.atoms, atomindices, boundaryatoms, getattr ( qcModel, "linkAtomRatio", True ) )
            # . Get the qcParameters.
            anlist       = qcAtoms.UniqueAtomicNumbers ( )
            qcParameters = qcModel.DefineParameters ( anlist, log = log )
            # . Fill in the remaining qcatom data.
            qcAtoms.AssignParameterData ( anlist, qcParameters )
            # . Save the data.
            em.qcAtoms = qcAtoms
            em.qcModel = qcModel
            em.__dict__["qcParameters"] = qcParameters
            em.MakeLabel ( )
            # . Modify the MM terms if appropriate.
            em.DeactivateQCAtomMMTerms ( )
            # . Deactivate MM terms if fixed atoms are defined.
            if ( self.hardConstraints is not None ): em.DeactivateFixedAtomMMTerms ( self.hardConstraints.fixedAtoms )
            # . Check the total MM charge.
            em.CheckActiveMMAtomTotalCharge ( )

    def DefineSoftConstraints ( self, softConstraints ):
        """Assign soft constraints to a system."""
        # . Get an energy model in the correct state.
        if self.energyModel is None: self.__dict__["_energyModel"] = EnergyModel ( )
        em = self.energyModel
        # . Set or remove the soft constraints as required.
        if isinstance ( softConstraints, SoftConstraintContainer ): em.softConstraints = softConstraints
        else:                                                       del em.softConstraints

    def DefineSymmetry ( self, a = None, alpha = None, b = None, beta = None, c = None, crystalClass = None, gamma = None, transformations = None ):
        """Define symmetry for a system."""
        # . Initialization.
        self.__dict__["_symmetry"] = None
        self.symmetryParameters    = None            
        # . A crystal class is present.
        if crystalClass is not None:
            if not isinstance ( crystalClass, CrystalClass ): raise TypeError ( "Invalid crystal class." )
            if ( transformations is None ): transformations = Transformation3Container_Identity ( )
            elif not isinstance ( transformations, Transformation3Container ): raise TypeError ( "Invalid symmetry transformations." )
            # . Create the symmetry and symmetry parameters.
            symmetry = Symmetry ( )
            symmetry.crystalClass      = crystalClass
            symmetry.transformations   = transformations
            self.__dict__["_symmetry"] = symmetry
            self.symmetryParameters    = crystalClass.CreateSymmetryParameters ( a = a, alpha = alpha, b = b, beta = beta, c = c, gamma = gamma )

    def DipoleMoment ( self, center = None ):
        """Calculate the dipole moment vector of the system."""
        dipole = None
        if self.energyModel is not None:
            em = self.energyModel
            if em.mmModel is not None: dipole = em.mmAtoms.DipoleMoment ( self.configuration, center )
            else:                      dipole = Vector3.Null ( )
            if em.qcModel is not None: dipole.AddScaledVector3 ( 1.0, em.qcModel.DipoleMoment ( self.configuration, center ) )
        return dipole

    # . This should probably be split up to better allow reuse.
    # . E.g. EnergyInitialize, EnergyCalculation, EnergyFinalize.
    def Energy ( self, log = logFile, doGradients = False ):
        """Calculate the energy and, optionally, the gradients for a system."""
        # . Do nothing if coordinates and an energy model are absent.
        if ( self.coordinates3 is None ) or ( self.energyModel is None ): return
        # . Set up configuration.
        self.configuration.ClearTemporaryAttributes ( )
        self.configuration.SetTemporaryAttribute ( "energyTerms", EnergyTerms ( ) )
        # . Make this clearer here (no need to reallocate each time).
        if doGradients:
            g = Coordinates3.WithExtent ( len ( self.atoms ) )
            g.Set ( 0.0 )
            self.configuration.SetTemporaryAttribute ( "gradients3", g )
            if self.symmetry is not None:
                self.configuration.SetTemporaryAttribute ( "symmetryParameterGradients", SymmetryParameterGradients ( ) )
        # . Alias various data structures.
        em  = self.energyModel
        et  = self.configuration.energyTerms
        fa  = None
        crd = self.configuration.coordinates3
        grd = getattr ( self.configuration, "gradients3", None )
        if self.hardConstraints is not None: fa = self.hardConstraints.fixedAtoms
        # . Calculate MM terms.
        if ( em.mmModel is not None ) and ( em.mmTerms is not None ):
            for mmterm in em.mmTerms: et.Append ( mmterm.Energy ( crd, grd ) )
        # . Calculate NB terms.
        if em.nbModel is not None:
            em.nbModel.SetUp ( em.mmAtoms, em.qcAtoms, em.ljParameters, em.ljParameters14, fa, em.interactions14, em.exclusions, self.symmetry, self.connectivity.isolates, self.configuration, log = log )
            et.Extend ( em.nbModel.Energy ( self.configuration ) )
        # . Calculate QC and QC/MM electrostatic terms.
        if em.qcModel is not None:
            if ( em.nbModel is not None ):                 em.nbModel.QCMMPotentials ( self.configuration )
            et.Extend ( em.qcModel.Energy ( em.qcAtoms, em.qcParameters, self.electronicState, self.configuration, log = log ) )
            if ( em.nbModel is not None ) and doGradients: em.nbModel.QCMMGradients  ( self.configuration )
        # . Calculate the soft constraint terms.
        if em.softConstraints is not None:
            ( scEnergy, scState ) = em.softConstraints.Energy ( crd, grd )
            self.configuration.SetTemporaryAttribute ( "scState", scState )
            et.Append ( scEnergy )
        # . Calculate the total energy.
        et.Total ( )
        # . Zero out unnecessary gradient terms.
        if ( fa is not None ) and doGradients: grd.SetRowSelection ( fa, 0.0 )
        # . Print the energy terms if necessary.
        if LogFileActive ( log ):
            if doGradients: et.RMSGradient ( grd.RMSValue ( ) )
            et.Summary ( log = log )
        return et.potentialEnergy

    @classmethod
    def FromAtoms ( selfClass, atoms, bonds = None, withSequence = False ):
        """Constructor given a list of atoms."""
        self = selfClass ( )
        self.__dict__["_atoms"       ] = AtomContainer.FromIterable     ( atoms )
        self.__dict__["_connectivity"] = Connectivity.FromAtomContainer ( self.atoms, bonds = bonds )
        if withSequence: self.__dict__["_sequence"] = Sequence.FromAtomContainer ( self.atoms )
        return self

    @classmethod
    def FromSequence ( selfClass, sequence, bonds = None ):
        """Constructor given a sequence."""
        atoms = sequence.GatherAtoms ( )
        self  = selfClass ( )
        self.__dict__["_atoms"       ] = AtomContainer.FromIterable     ( atoms )
        self.__dict__["_connectivity"] = Connectivity.FromAtomContainer ( self.atoms, bonds = bonds )
        self.__dict__["_sequence"    ] = sequence
        return self

    # . Rationalize Merge and Prune - especially for HardConstraints.
    # . Merge would be better using a constructor with a list of merge items (including possibly None).

    def Merge ( self, others, information = {} ):
        """Merging."""
        # . Initialization.
        toMergeItems  = [ self ] + list ( others )
        numberToMerge = len ( toMergeItems )
        # . Basic system.
        merged = System ( )
        # . Attributes.
        attributes = self.__class__.defaultAttributes.keys ( )
        attributes.sort ( )
        for attribute in attributes:
            items = []
            for toMerge in toMergeItems:
                item = getattr ( toMerge, attribute, None )
                if item is not None: items.append ( item )
            # . All objects being merged must have the attribute.
            if len ( items ) == numberToMerge:
                if hasattr ( items[0], "Merge" ):
                    mergedItem = items[0].Merge ( items[1:], information = information )
                    if mergedItem is not None: merged.__dict__[attribute] = mergedItem
        # . Charge constraints.
        chargeConstraints = []
        increments        = []
        items             = [ getattr ( item, "hardConstraints" ) for item in toMergeItems ]
        for ( hardConstraints, increment ) in zip ( items, information["atomIncrements"] ):
            if ( hardConstraints is not None ) and ( len ( hardConstraints.chargeConstraints ) > 0 ):
                chargeConstraints.append ( hardConstraints.chargeConstraints )
                increments.append ( increment )
        if len ( chargeConstraints ) > 0:
            merged.DefineChargeConstraints ( ChargeConstraintContainer.Merge ( chargeConstraints, information = { "atomIncrements" : increments } ) )
        # . Fixed atoms.
        fixedAtoms = []
        increments = []
        items      = [ getattr ( item, "hardConstraints" ) for item in toMergeItems ]
        for ( hardConstraints, increment ) in zip ( items, information["indexIncrements"] ):
            if ( hardConstraints is not None ) and ( len ( hardConstraints.fixedAtoms ) > 0 ):
                fixedAtoms.append ( hardConstraints.fixedAtoms )
                increments.append ( increment )
        if len ( fixedAtoms ) > 0:
            merged.DefineFixedAtoms ( fixedAtoms[0].Merge ( fixedAtoms[1:], { "indexIncrements" : increments } ) )
        # . Label.
        merged.label = "Merged System"
        # . Finish up.
        return merged

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        # . Basic system.
        pruned = self.__class__ ( )
        # . Attributes.
        attributes  = self.__class__.defaultAttributes.keys ( )
        attributes.sort ( )
        for attribute in attributes:
            item = self.__dict__.get ( attribute, None )
            if ( item is not None ) and hasattr ( item, "Prune" ):
                prunedItem = item.Prune ( selection, information = information )
                if ( prunedItem is not None ): pruned.__dict__[attribute] = prunedItem
        # . Hard constraints.
        if self.hardConstraints is not None:
            if ( self.hardConstraints.chargeConstraints is not None ) and ( len ( self.hardConstraints.chargeConstraints ) > 0 ):
                pruned.DefineChargeConstraints ( self.hardConstraints.chargeConstraints.Prune ( selection, information = information ) )
            if ( self.hardConstraints.fixedAtoms is not None ) and ( len ( self.hardConstraints.fixedAtoms ) > 0 ):
                pruned.DefineFixedAtoms ( self.hardConstraints.fixedAtoms.Prune ( selection ) )
        # . Label.
        pruned.label = "Pruned System"
        # . Finish up.
        return pruned

    def PruneToQCRegion ( self ):
        """Return a system consisting of the QC region, including boundary atoms."""
        # . Initialization.
        qcRegion = None
        # . Do nothing if there is no QC region.
        if ( self.energyModel is not None ) and ( self.energyModel.qcModel is not None ):
            # . Get the atoms to retain.
            toKeep = self.energyModel.qcAtoms.GetFullSelection ( )
            # . Prune the system.
            qcRegion = self.Prune ( toKeep )
            # . Other attributes.
            qcRegion.coordinates3 = self.energyModel.qcAtoms.GetCoordinates3 ( self.coordinates3 )
            if self.label is None: qcRegion.label = "QC Region"
            else:                  qcRegion.label = self.label + " - QC Region"
            # . Electronic state and QC model.
            qcRegion.electronicState = Clone ( self.electronicState )
            qcRegion.DefineQCModel ( self.energyModel.qcModel )
        # . Finish up.
        return qcRegion

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass ( )
        return self

    def Summary ( self, log = logFile ):
        """Summarizing."""
        if LogFileActive ( log ):
            if self.label is None: title = self.__class__.__name__ + " Summary"
            else:                  title = "Summary for " + self.__class__.__name__ + " \"" + self.label + "\""
            log.Heading ( title, includeBlankLine = True )
            attributes = self.__class__.defaultAttributes.keys ( )
            attributes.sort ( )
            for attribute in attributes:
                item = self.__dict__.get ( attribute, None )
                if ( item is not None ) and hasattr ( item, "Summary" ): item.Summary ( log = log )

    # . Best place?
    def VerifyFixedAtomCoordinates ( self, trial, log = logFile, reference = None, tolerance = 1.0e-03 ):
        """Verify that the coordinates of fixed atoms in a trial set are the same as in a reference set."""
        # . Initialization.
        fixedAtoms = None
        isOK       = False
        if self.hardConstraints is not None: fixedAtoms = self.hardConstraints.fixedAtoms
        if reference            is     None: reference  = self.coordinates3
        # . Basic checks.
        if   reference is None: message = "Missing reference coordinate set."
        elif trial     is None: message = "Missing trial coordinate set."
        elif ( trail.rows != len ( system.atoms ) ) or ( len ( reference ) != len ( trial ) ):
            message = "Coordinate sets of incompatible size."
        elif fixedAtoms is None:
            isOK    = True
            message = "No fixed atoms are defined."
        # . Check the coordinates.
        else:
            set1 = Coordinates3 ( len ( fixedAtoms ) ) ; set1.Gather ( reference, fixedAtoms )
            set2 = Coordinates3 ( len ( fixedAtoms ) ) ; set2.Gather ( other    , fixedAtoms )
            set1.AddScaledMatrix ( -1.0, set2 )
            biggest = set1.AbsoluteMaximum ( )
            if biggest >= tolerance:
                message = "The maximum difference in the fixed atom coordinates, {:.6g}, is greater than the tolerance, {:.6g}.".format ( biggest, tolerance )
            else:
                isOK    = True
                message = "The fixed atom coordinates are compatible."
        # . Printing.
        if LogFileActive ( log ): log.Paragraph ( message )
        # . Finish up.
        return isOK

    # . Properties.
    # . System attribute getters.
    @property
    def atoms           ( self ): return self.__dict__ [ "_atoms"           ]
    @property
    def configuration   ( self ): return self.__dict__ [ "_configuration"   ]
    @property
    def connectivity    ( self ): return self.__dict__ [ "_connectivity"    ]
    @property
    def electronicState ( self ): return self.__dict__ [ "_electronicState" ]
    @property
    def energyModel     ( self ): return self.__dict__ [ "_energyModel"     ]
    @property
    def hardConstraints ( self ): return self.__dict__ [ "_hardConstraints" ]
    @property
    def label           ( self ): return self.__dict__ [ "_label"           ]
    @property
    def sequence        ( self ): return self.__dict__ [ "_sequence"        ]
    @property
    def symmetry        ( self ): return self.__dict__ [ "_symmetry"        ]

    # . Configuration attribute getters.
    @property
    def coordinates3 ( self ):
        if ( self.configuration is not None ): return getattr ( self.configuration, "coordinates3"      , None )
        else:                                  None
    @property
    def symmetryParameters ( self ):
        if ( self.configuration is not None ): return getattr ( self.configuration, "symmetryParameters", None )
        else:                                  None

    # . System attribute setters.
    # . Connectivity and sequence here too?
    @electronicState.setter
    def electronicState ( self, value ):
        if ( value is None ) or isinstance ( value, ElectronicState ):
            # . Reset the state.
            oldValue = self.__dict__.get ( "_electronicState", None )
            self.__dict__["_electronicState"] = value
            # . Decide what to do with a qcState if it exists.
            qcState = None
            if self.configuration is not None: qcState = getattr ( self.configuration, "qcState", None )
            if qcState is not None:
                removeState = False
                # . Remove previous electronic state.
                if ( value is None ) or ( oldValue is None ): removeState = True
                # . The electronic state has changed so try and modify qcState if possible.
                elif ( value.charge != oldValue.charge ) or ( value.multiplicity != oldValue.multiplicity ):
                    if hasattr ( qcState, "ModifyElectronicState" ): qcState.ModifyElectronicState ( value )
                    else: removeState = True
                # . Remove the state.
                if removeState: delattr ( self.configuration, "qcState" )
        else:
            raise AttributeError ( "Invalid electronicState attribute (wrong type)." )
    @label.setter
    def label ( self, value ):
        if ( value is None ) or isinstance ( value, basestring ): self._label = value

    # . Configuration attribute setters.
    @coordinates3.setter
    def coordinates3 ( self, value ):
        if ( value is None ) or ( isinstance ( value, Coordinates3 ) and ( value.rows == len ( self.atoms ) ) ):
            if ( self.configuration is None ): self.__dict__["_configuration"] = Configuration ( )
            self.configuration.coordinates3 = value
        else: raise AttributeError ( "Invalid coordinates3 attribute (wrong dimension or type)." )
    @symmetryParameters.setter
    def symmetryParameters ( self, value ):
        # . Check for crystal class compatibility here?
        if ( value is None ) or isinstance ( value, SymmetryParameters ):
            if ( self.configuration is None ): self.__dict__["_configuration"] = Configuration ( )
            self.configuration.symmetryParameters = value
        else: raise AttributeError ( "Invalid symmetryParameters attribute (wrong crystal class or type)." )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
