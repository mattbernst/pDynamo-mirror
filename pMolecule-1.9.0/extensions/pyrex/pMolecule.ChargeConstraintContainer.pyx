#-------------------------------------------------------------------------------
# . File      : pMolecule.ChargeConstraintContainer.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A container for charge constraints."""

from pCore import CLibraryError, logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ChargeConstraintContainer:
    """A container for charge constraints."""

    def __copy__ ( self ):
        """Copying."""
        constraints = self.__getstate__ ( )
        new         = self.__class__.FromList ( constraints )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            ChargeConstraintContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.ChargeConstraintContainer"

    def __getstate__ ( self ):
        """Return the state as a list."""
        cdef Integer c, t
        cdef CChargeConstraint *constraint
        # . Always return a list for a constraint charge or spin component even if empty.
        state = []
        if self.cObject != NULL:
            for c from 0 <= c < self.cObject.numberOfConstraints:
                constraint = self.cObject.constraints[c]
                charges    = []
                spins      = []
                for t from 0 <= t < constraint.numberOfCharges:
                    charges.append ( [ Integer1DArray_GetItem ( constraint.chargeIndices, t, NULL ), Real1DArray_GetItem ( constraint.chargeWeights, t, NULL ) ] )
                for t from 0 <= t < constraint.numberOfSpins:
                    spins.append ( [ Integer1DArray_GetItem ( constraint.spinIndices, t, NULL ), Real1DArray_GetItem ( constraint.spinWeights, t, NULL ) ] )
                state.append ( ( charges, spins, constraint.target ) )
        return state

    def __init__ ( self, constraints ):
        """Constructor with options."""
        self._Initialize ( )
        self.BuildFromList ( constraints )

    def __len__ ( self ):
        """Length."""
        return ChargeConstraintContainer_Size ( self.cObject )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self.BuildFromList ( state )

    def _Allocate ( self, extent ):
        """Allocation."""
        self.cObject = ChargeConstraintContainer_Allocate ( extent, NULL )
        self.isOwner = True
        if self.cObject == NULL: raise CLibraryError ( "Error allocating C object." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def BuildFromList ( self, constraints ):
        """Build a container from a list."""
        cdef Integer c, i, q, t, s
        cdef Real    w
        cdef CChargeConstraint *constraint
        try:
            n = len ( constraints )
            if n > 0:
                self._Allocate ( n )
                for c from 0 <= c < n:
                    ( charges, spins, target ) = constraints[c]
                    if charges is None: q = 0
                    else:               q = len ( charges )
                    if spins   is None: s = 0
                    else:               s = len ( spins   )
                    constraint = ChargeConstraint_Allocate ( q, s, NULL )
                    if constraint == NULL: raise CLibraryError ( "Out of memory." )
                    for ( t, ( i, w ) ) in enumerate ( charges ):
                        Integer1DArray_SetItem ( constraint.chargeIndices, t, i, NULL )
                        Real1DArray_SetItem    ( constraint.chargeWeights, t, w, NULL )
                    for ( t, ( i, w ) ) in enumerate ( spins   ):
                        Integer1DArray_SetItem ( constraint.spinIndices  , t, i, NULL )
                        Real1DArray_SetItem    ( constraint.spinWeights  , t, w, NULL )
                    constraint.target = target
                    ChargeConstraintContainer_SetItem ( self.cObject, c, &constraint, NULL )
        except:
            raise CLibraryError ( "Error building container from list." )

    def Deviations ( self, Real1DArray charges, Real1DArray spins, Real1DArray deviations ):
        """Calculate the deviations from the constraints given arrays of charges and spins."""
        cdef CReal1DArray *cCharges, *cDeviations, *cSpins
        if charges    is None: cCharges    = NULL
        else:                  cCharges    = charges.cObject
        if deviations is None: cDeviations = NULL
        else:                  cDeviations = deviations.cObject
        if spins      is None: cSpins      = NULL
        else:                  cSpins      = spins.cObject
        ChargeConstraintContainer_Deviations ( self.cObject, cCharges, cSpins, cDeviations, NULL )

    @classmethod
    def FromList ( selfClass, constraints ):
        """Constructor from list."""
        return selfClass ( constraints )

    def HighestIndex ( self ):
        """Return the highest index of the container."""
        return ChargeConstraintContainer_HighestIndex ( self.cObject )

    @classmethod
    def Merge ( selfClass, items, information = {} ):
        """Constructor by merging several existing containers."""
        cdef Integer i, iOld
        cdef ChargeConstraintContainer self
        # . Initialization.
        self = None
        # . Get the selection increments.
        increments = information.get ( "atomIncrements", None )
        # . Do the merge.
        if increments is not None:
            # . Initialization.
            constraints = []
            # . Gather parameters and terms.
            for ( increment, item ) in zip ( increments, items ):
                if item is not None:
                    localState = item.__getstate__ ( )
                    for constraint in localState:
                        for term in constraint[0]: term[0] += increment # . Charges.
                        for term in constraint[1]: term[0] += increment # . Spins.
                        constraints.append ( constraint )
            # . Construct the object.
            self = selfClass.FromList ( constraints )
        return self

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef ChargeConstraintContainer new
        new         = self.__class__.Raw ( )
        new.cObject = ChargeConstraintContainer_Prune ( self.cObject, selection.cObject, NULL )
        new.isOwner = True
        if new.cObject == NULL: return None
        else:                   return new

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    property hasCharges:
        def __get__ ( self ): return ( ( self.cObject != NULL ) and ( self.cObject.hasCharges == CTrue ) )

    property hasSpins:
        def __get__ ( self ): return ( ( self.cObject != NULL ) and ( self.cObject.hasSpins   == CTrue ) )
