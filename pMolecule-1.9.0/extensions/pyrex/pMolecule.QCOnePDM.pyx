#-------------------------------------------------------------------------------
# . File      : pMolecule.QCOnePDM.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle quantum chemical one-particle densities and associated data."""

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Dictionaries for interconverting between C and Python representations.
#===================================================================================================================================
# . Density types.
QCOnePDMDensityType_ToEnum     = { "Alpha" : QCOnePDMDensityType_Alpha , "Beta" : QCOnePDMDensityType_Beta  , "Spin" : QCOnePDMDensityType_Spin  , "Total" : QCOnePDMDensityType_Total }
QCOnePDMDensityType_ToString   = { QCOnePDMDensityType_Alpha : "Alpha" , QCOnePDMDensityType_Beta  : "Beta" , QCOnePDMDensityType_Spin  : "Spin" , QCOnePDMDensityType_Total : "Total" }

# . Occupancy types.
QCOnePDMOccupancyType_ToEnum   = { "Cardinal" : QCOnePDMOccupancyType_Cardinal, "Fixed Fractional" : QCOnePDMOccupancyType_FractionalFixed, "Variable Fractional" : QCOnePDMOccupancyType_FractionalVariable }
QCOnePDMOccupancyType_ToString = { QCOnePDMOccupancyType_Cardinal : "Cardinal", QCOnePDMOccupancyType_FractionalFixed : "Fixed Fractional", QCOnePDMOccupancyType_FractionalVariable : "Variable Fractional" }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCOnePDM:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef QCOnePDM new
        new         = self.__class__.Raw ( )
        new.cObject = QCOnePDM_Clone ( self.cObject, NULL )
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        if self.isOwner:
            QCOnePDM_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.QCOnePDM"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Real1DArray     occupancies
        cdef SymmetricMatrix density
        state = {}
        if self.cObject != NULL:
            state.update ( { "isValid"         : self.cObject.isValid == CTrue, \
                             "numberOccupied"  : self.cObject.numberOccupied  , \
                             "numberOrbitals"  : self.cObject.numberOrbitals  , \
                             "fermiBroadening" : self.cObject.fermiBroadening , \
                             "occupancyFactor" : self.cObject.occupancyFactor , \
                             "totalCharge"     : self.cObject.totalCharge     , \
                             "densityType"     : QCOnePDMDensityType_ToString  [self.cObject.densityType  ], \
                             "occupancyType"   : QCOnePDMOccupancyType_ToString[self.cObject.occupancyType]  } )
            if ( self.cObject.density != NULL ):
                density              = SymmetricMatrix.Raw ( )
                density.cObject      = SymmetricMatrix_Clone ( self.cObject.density )
                density.isOwner      = True
                state["density"]     = density
            if ( self.cObject.occupancies != NULL ):
                occupancies          = Real1DArray.Raw ( )
                occupancies.cObject  = Real1DArray_Clone ( self.cObject.occupancies, NULL )
                occupancies.isOwner  = True
                state["occupancies"] = occupancies
        return state

    def __init__ ( self, size ):
        """Constructor with size."""
        self._Initialize ( )
        self._Allocate ( size )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( 0 )
        self.SetOptions ( **state )

    def _Allocate ( self, size ):
        """Allocation."""
        cdef Integer i
        cdef Status  status
        status       = Status_Continue
        self.cObject = QCOnePDM_Allocate ( size, &status )
        self.isOwner = True
        if status   != Status_Continue: raise ValueError ( "Object allocation error." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = True
        self.owner   = None

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def ResetDensityFromDensity ( self, SymmetricMatrix density, SymmetricMatrix overlap = None ):
        """Reset the density using another density.

        The overlap is required for non-orthogonal basis sets.
        """
        cdef Boolean           isOK
        cdef CSymmetricMatrix *cOverlap
        if overlap is None: cOverlap = NULL
        else:               cOverlap = overlap.cObject
        isOK = QCOnePDM_ResetDensityFromDensity  ( self.cObject, density.cObject, cOverlap, NULL )
        if isOK == False: raise ValueError ( "Trying to reset density with an invalid density." )

    def ResetDensityFromOrbitals ( self, Real2DArray orbitals ):
        """Reset the density using a set of orbitals."""
        cdef Boolean isOK
        isOK = QCOnePDM_ResetDensityFromOrbitals ( self.cObject, orbitals.cObject, NULL )
        if isOK == False: raise ValueError ( "Trying to reset density with an invalid set of orbitals." )

    def SetOptions ( self, **keywordArguments ):
        """Set options for the density."""
        cdef Real1DArray     occupancies
        cdef SymmetricMatrix density
        if self.cObject != NULL:
            if "isValid"         in keywordArguments :
                if keywordArguments.pop ( "isValid" ): self.cObject.isValid = CTrue
                else:                                  self.cObject.isValid = CFalse
            if "numberOccupied"  in keywordArguments : self.cObject.numberOccupied  = keywordArguments.pop ( "numberOccupied"  )
            if "numberOrbitals"  in keywordArguments : self.cObject.numberOrbitals  = keywordArguments.pop ( "numberOrbitals"  )
            if "fermiBroadening" in keywordArguments : self.cObject.fermiBroadening = keywordArguments.pop ( "fermiBroadening" )
            if "occupancyFactor" in keywordArguments : self.cObject.occupancyFactor = keywordArguments.pop ( "occupancyFactor" )
            if "totalCharge"     in keywordArguments : self.cObject.totalCharge     = keywordArguments.pop ( "totalCharge"     )
            if "densityType"     in keywordArguments : self.cObject.densityType     = QCOnePDMDensityType_ToEnum  [keywordArguments.pop ( "densityType"   )]
            if "occupancyType"   in keywordArguments : self.cObject.occupancyType   = QCOnePDMOccupancyType_ToEnum[keywordArguments.pop ( "occupancyType" )]
            if "density" in keywordArguments:
                density = keywordArguments.pop ( "density" )
                SymmetricMatrix_Deallocate ( &(self.cObject.density) )
                self.cObject.density = SymmetricMatrix_Clone ( density.cObject )
            if "occupancies" in keywordArguments:
                occupancies = keywordArguments.pop ( "occupancies" )
                Real1DArray_Deallocate ( &(self.cObject.occupancies) )
                self.cObject.occupancies = Real1DArray_Clone ( occupancies.cObject, NULL )
        if len ( keywordArguments ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( keywordArguments.keys ( ) ) ) + "." )

    # . Properties.
    property density:
        def __get__ ( self ):
            cdef SymmetricMatrix item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.density != NULL ):
                item = SymmetricMatrix.Raw ( )
                item.cObject = self.cObject.density
                item.isOwner = False
            return item
    property energies:
        def __get__ ( self ):
            cdef Real1DArray item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.energies != NULL ):
                item = Real1DArray.Raw ( )
                item.cObject = self.cObject.energies
                item.isOwner = False
                item.owner   = self
            return item
    property fock:
        def __get__ ( self ):
            cdef SymmetricMatrix item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.fock != NULL ):
                item = SymmetricMatrix.Raw ( )
                item.cObject = self.cObject.fock
                item.isOwner = False
            return item
    property occupancies:
        def __get__ ( self ):
            cdef Real1DArray item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.occupancies != NULL ):
                item = Real1DArray.Raw ( )
                item.cObject = self.cObject.occupancies
                item.isOwner = False
                item.owner   = self
            return item
    property orbitals:
        def __get__ ( self ):
            cdef Real2DArray item
            item = None
            if ( self.cObject != NULL ) and ( self.cObject.orbitals != NULL ):
                item = Real2DArray.Raw ( )
                item.cObject = self.cObject.orbitals
                item.isOwner = False
                item.owner   = self
            return item
