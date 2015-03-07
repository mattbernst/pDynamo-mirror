#-------------------------------------------------------------------------------
# . File      : pCore.Transformation3Container.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A container for transformations of dimension 3."""

from Clone         import Clone
from LogFileWriter import logFile, LogFileActive
from Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Transformation3Container:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef Transformation3Container new
        items = Clone ( self.items )
        new   = self.__class__.WithTransformations ( items )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Transformation3Container_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getitem__ ( self, Integer index ):
        """Get an item."""
        if ( index < 0 ) or ( index >= len ( self ) ): raise IndexError
        else:                                          return self.items[index]

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pCore.Transformation3Container"

    def __getstate__ ( self ):
        """Return the state."""
        return { "items" : self.items }

    def __init__ ( self, transformations ):
        """Constructor given a set of transformations."""
        self._Initialize ( )
        self.CreateTransformationList ( transformations )
        self._Allocate ( len ( self.items ) )
        self.FillCObject ( )

    def __len__ ( self ):
        """Return the number of transformations."""
        return len ( self.items )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self.items = state["items"]
        self._Allocate ( len ( self.items ) )
        self.FillCObject ( )

    def _Allocate ( self, size ):
        """Allocation."""
        self.cObject = Transformation3Container_Allocate ( size )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.items   = None

    def CreateTransformationList ( self, transformations ):
        """Create a list of unique transformations."""
        cdef Transformation3 new, old
        self.items = []
        for new in transformations:
            isUnique = True
            for old in self.items:
                if ( Transformation3_IsEqual ( new.cObject, old.cObject ) == CTrue ):
                    isUnique = False
                    break
            if isUnique: self.items.append ( new )

    def FillCObject ( self ):
        """Fill the C object from items."""
        cdef Integer             i
        cdef Transformation3 t
        # . Loop over the transformations.
        for ( i, t ) in enumerate ( self.items ):
            # . Copy the pointer.
            self.cObject.items[i] = t.cObject
        # . Find the identity.
        Transformation3Container_FindIdentity ( self.cObject )
        # . Find inverses.
        Transformation3Container_FindInverses ( self.cObject )

    def IsIdentity ( self ):
        """Does this container correspond to the identity?"""
        cdef Transformation3 item
        isIdentity = False
        if len ( self.items ) == 1:
            item       = self.items[0]
            isIdentity = ( Transformation3_IsIdentity ( item.cObject ) == CTrue )
        return isIdentity

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        cdef Integer i
        if LogFileActive ( log ):
            nt3 = len ( self.items )
            # . Inverse information.
            npair = 0
            nself = 0
            for i from 0 <= i < self.cObject.nitems:
                if   self.cObject.inverses[i] == i: nself = nself + 1
                elif self.cObject.inverses[i] >= 0: npair = npair + 1
            # . Rotation information - not necessarily valid unless have symmetryParameters.
#            nimproper = 0
#            nproper   = 0
#            for t3 in self.items:
#                if   t3.rotation.IsImproperRotation ( ): nimproper = nimproper + 1
#                elif t3.rotation.IsProperRotation   ( ): nproper   = nproper   + 1
            summary = log.GetSummary ( )
            summary.Start ( "Transformation3 Container Summary" )
            summary.Entry ( "Transformations"    , "{:d}".format ( nt3                       ) )
            summary.Entry ( "Identity"           , "{:d}".format ( self.cObject.identity     ) )
            summary.Entry ( "Inverse Pairs"      , "{:d}".format ( npair // 2                ) )
            summary.Entry ( "Self Inverses"      , "{:d}".format ( nself                     ) )
            summary.Entry ( "Absent Inverses"    , "{:d}".format ( nt3 - npair - nself       ) )
#            summary.Entry ( "Improper Rotations" , "{:d}".format ( nimproper                 ) )
#            summary.Entry ( "Proper Rotations"   , "{:d}".format ( nproper                   ) )
#            summary.Entry ( "Invalid Rotations"  , "{:d}".format ( nt3 - nimproper - nproper ) )
            summary.Stop ( )

    @classmethod
    def WithTransformations ( selfClass, transformations ):
        """Constructor given a set of transformations."""
        return selfClass ( transformations )

#===================================================================================================================================
# . Class methods.
#===================================================================================================================================
def Transformation3Container_Identity ( ):
    """Construct a container containing only the identity transformation."""
    cdef Transformation3Container self
    identity = Transformation3 ( )
    identity.rotation.Set    ( 0.0 )
    identity.translation.Set ( 0.0 )
    for i in range ( 3 ):
        identity.rotation[i,i] = 1.0
    self = Transformation3Container.WithTransformations ( [ identity ] )
    return self
