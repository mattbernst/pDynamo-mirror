#-------------------------------------------------------------------------------
# . File      : pMolecule.MMAtomContainer.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A container for MM atoms."""

import copy

from pCore import Clone, logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MMAtomContainer:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef MMAtomContainer new
        new           = self.__class__.Raw ( )
        new.cObject   = MMAtomContainer_Clone ( self.cObject )
        new.atomTypes = Clone ( self.atomTypes )
        new.isOwner   = True
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            MMAtomContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.MMAtomContainer"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i
        atomFields = [ "isActive", "atomType", "ljType", "charge" ]
        atoms      = []
        for i from 0 <= i < self.cObject.natoms:
            data = [ self.cObject.data[i].QACTIVE == CTrue ,
                     self.cObject.data[i].atomtype        ,
                     self.cObject.data[i].ljtype          ,
                     self.cObject.data[i].charge          ]
            atoms.append ( data )
        return { "atomFields" : atomFields, "atoms" : atoms, "atomTypes" : self.atomTypes }

    def __init__ ( self, numberOfAtoms ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfAtoms )

    def __len__ ( self ):
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer i
        # . Basic construction.
        atoms          = state["atoms"    ]
        self.atomTypes = state["atomTypes"]
        self._Allocate ( len ( atoms ) )
        # . Fill the object.
        for ( i, ( q, aType, ljType, charge ) ) in enumerate ( atoms ):
            if q: self.cObject.data[i].QACTIVE = CTrue
            else: self.cObject.data[i].QACTIVE = CFalse
            self.cObject.data[i].atomtype = aType
            self.cObject.data[i].charge   = charge
            self.cObject.data[i].ljtype   = ljType

    def _Allocate ( self, numberOfAtoms ):
        """Allocation."""
        self.cObject = MMAtomContainer_Allocate ( numberOfAtoms )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.atomTypes = None
        self.cObject   = NULL
        self.isOwner   = False

    def ActivateAtoms ( self ):
        """Activate all atoms."""
        cdef Integer i
        for i from 0 <= i < self.cObject.natoms: self.cObject.data[i].QACTIVE = CTrue

    def AtomicCharges ( self ):
        """Get the atomic charges for the active atoms."""
        cdef Integer     i
        cdef Real1DArray charges
        charges = Real1DArray.WithExtent ( self.size )
        charges.Set ( 0.0 )
        for i from 0 <= i < self.cObject.natoms:
            if self.cObject.data[i].QACTIVE == CTrue: charges.cObject.data[i] = self.cObject.data[i].charge
        return charges

    def AtomTypes ( self ):
        """Get the types of the atoms."""
        cdef Integer i
        types = []
        for i from 0 <= i < self.cObject.natoms:
            types.append ( self.atomTypes[self.cObject.data[i].atomtype] )
        return types

    def DeactivateQCAtoms ( self, Selection qcAtoms, Selection boundaryatoms ):
        """Deactivate QC atoms."""
        # . All QC atoms are deactivated except for boundary atoms.
        cdef Integer i
        if qcAtoms is not None:
            for i in qcAtoms: self.cObject.data[i].QACTIVE = CFalse
        if boundaryatoms is not None:
            for i in boundaryatoms: self.cObject.data[i].QACTIVE = CTrue

    def DipoleMoment ( self, configuration, Vector3 center = None ):
        """Dipole moment."""
        cdef Coordinates3  pcoordinates3
        cdef Vector3       dipole
        cdef CVector3     *ccenter
        pcoordinates3 = configuration.coordinates3
        if center    is None: ccenter   = NULL
        else:                 ccenter   = center.cObject
        dipole         = Vector3.Raw ( )
        dipole.cObject = MMAtomContainer_DipoleMoment ( self.cObject, pcoordinates3.cObject, ccenter )
        dipole.isOwner = True
        return dipole

    def LennardJonesRadii ( self, LJParameterContainer ljParameters ):
        """Get the Lennard-Jones radii for the atoms."""
        cdef Integer     i
        cdef Real1DArray radii
        # . Get the LJ radii.
        state  = ljParameters.__getstate__ ( )
        sigmas = state.get ( "sigmas", None )
        if sigmas is None: raise ValueError ( "Unable to extract radii from Lennard-Jones container." )
        # . Get the atom radii.
        radii = Real1DArray.WithExtent ( self.size )
        for i from 0 <= i < self.cObject.natoms:
            radii.cObject.data[i] = sigmas[self.cObject.data[i].ljtype]
        return radii

    def LennardJonesWellDepths ( self, LJParameterContainer ljParameters ):
        """Get the Lennard-Jones well-depths for the atoms."""
        cdef Integer     i
        cdef Real1DArray wellDepths
        # . Get the LJ well-depths.
        state    = ljParameters.__getstate__ ( )
        epsilons = state.get ( "epsilons", None )
        if epsilons is None: raise ValueError ( "Unable to extract well-depths from Lennard-Jones container." )
        # . Get the atom well-depths.
        wellDepths = Real1DArray.WithExtent ( self.size )
        for i from 0 <= i < self.cObject.natoms:
            wellDepths.cObject.data[i] = epsilons[self.cObject.data[i].ljtype]
        return wellDepths

    # . Note that the construction of a new set of types can disrupt correspondence with other MM objects (e.g. the LJ parameters).
    # . However, this should not affect the energy, etc. If labels are required use LJParameterContainer.parameterKeys.
    def Merge ( self, items, information = {} ):
        """Merging - slow version."""
        cdef Integer i, iOld
        cdef MMAtomContainer new
        # . Initialization.
        new = None
        # . Get the selection increments.
        increments = information.get ( "ljTypeIncrements" , None )
        mapping    = information.get ( "ljTypeMapping"    , None )
        # . Do the merge.
        if ( increments is not None ) and ( mapping is not None ):
            # . Get the old indices.
            atoms     = []
            oldTypes  = []
            ntypes    = 0
            for ( increment, item ) in zip ( increments, [ self ] + items ):
                state  = item.__getstate__ ( )
                atoms0 = state["atoms"]
                for ( q, aType, ljType, charge ) in atoms0:
                    atoms.append ( ( q, aType + ntypes, mapping[ljType + increment], charge ) )
                oldTypes.extend ( state["atomTypes"] )
                ntypes = len ( oldTypes )
            # . Construct new types.
            newTypes = list ( set ( oldTypes ) )
            oldToNew = []
            newTypes.sort ( )
            for oldType in oldTypes:
                oldToNew.append ( newTypes.index ( oldType ) )
            # . Construct the object.
            state = { "atoms" : atoms, "atomTypes" : newTypes }
            new   = self.__class__.Raw ( )
            new.__setstate__ ( state )
            # . Reindex the types.
            for i from 0 <= i < new.cObject.natoms:
                iOld = new.cObject.data[i].atomtype
                new.cObject.data[i].atomtype = oldToNew[iOld]
        return new

    def NumberOfAtoms ( self ):
        """Return the number of atoms."""
        if self.cObject == NULL: return 0
        else:                    return self.cObject.natoms

    def NumberOfAtomTypes ( self ):
        """Return the number of atom types."""
        if self.atomTypes is None: return 0
        else:                      return len ( self.atomTypes )

    def NumberOfInactiveAtoms ( self ):
        """Return the number of inactive atoms."""
        cdef Integer i, n
        n = 0
        for i from 0 <= i < self.cObject.natoms:
            if self.cObject.data[i].QACTIVE == CFalse: n = n + 1
        return n

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef Integer i
        cdef MMAtomContainer new
        # . Prune atom data.
        new         = MMAtomContainer.Raw ( )
        new.cObject = MMAtomContainer_Prune ( self.cObject, selection.cObject )
        # . Prune atom types.
        newTypes = set ( )
        for i from 0 <= i < new.cObject.natoms:
            newTypes.add ( self.atomTypes[new.cObject.data[i].atomtype] )
        newTypes = list ( newTypes )
        newTypes.sort ( )
        # . Get type index.
        oldToNew = {}
        for ( i, newType ) in enumerate ( newTypes ):
            oldToNew[self.atomTypes.index ( newType )] = i
        # . Assign type data.
        new.atomTypes = newTypes
        for i from 0 <= i < new.cObject.natoms:
            iOld = new.cObject.data[i].atomtype
            new.cObject.data[i].atomtype = oldToNew[iOld]
        # . Finish up.
        return new

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetAtomicCharges ( self, charges ):
        """Set the atomic charges for the active atoms."""
        cdef Integer i
        if len ( charges ) == len ( self ):
            for i from 0 <= i < self.cObject.natoms:
                if self.cObject.data[i].QACTIVE == CTrue: self.cObject.data[i].charge = charges[i]
        else:
            raise TypeError ( "Invalid charge array." )

    def SummaryEntry ( self, object summary ):
        """Summary entry."""
        if summary is not None:
            n = self.NumberOfInactiveAtoms ( )
            summary.Entry ( "Number of MM Atoms"     , "{:d}".format ( self.NumberOfAtoms     ( ) ) )
            summary.Entry ( "Number of MM Atom Types", "{:d}".format ( self.NumberOfAtomTypes ( ) ) )
            if n > 0: summary.Entry ( "Number of Inactive MM Atoms", "{:d}".format ( n ) )
            summary.Entry ( "Total MM Charge", "{:.2f}".format ( self.TotalCharge ( ) )  )

    def TotalCharge ( self ):
        """Get the atomic charges for the active atoms."""
        cdef Real q
        cdef Integer    i
        q = 0.0
        for i from 0 <= i < self.cObject.natoms:
            if self.cObject.data[i].QACTIVE == CTrue: q = q + self.cObject.data[i].charge
        return q

    # . Properties.
    property size:
        def __get__ ( self ):
            return MMAtomContainer_Size ( self.cObject )
