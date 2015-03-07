#-------------------------------------------------------------------------------
# . File      : pMolecule.QCAtomContainer.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle the data necessary to define quantum chemical atoms."""

from pCore   import logFile, LogFileActive, RawObjectConstructor
from Element import PeriodicTable

import math

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Shell labels.
_ShellLabels = [ "s", "p", "d", "f", "g", "h", "i", "j", "k" ]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCAtomContainer:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef QCAtomContainer new
        new         = self.__class__.Raw ( )
        new.cObject = QCAtomContainer_Clone ( self.cObject )
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            QCAtomContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getitem__ ( self, Integer index ):
        """Get an item."""
        if ( index < 0 ) or ( index >= len ( self ) ): raise IndexError
        else:                                          return self.cObject.data[index].index

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.QCAtomContainer"

    def __getstate__ ( self ):
        """Return the state."""
        cdef i, m
        # . Atom data.
        atomFields = [ "isBoundary"           ,
                       "linkFactor"           ,
                       "widthExponent"        ,
                       "widthNormalization"   ,
                       "atomicNumber"         ,
                       "center"               ,
                       "index"                ,
                       "qcPartner"            ,
                       "densityStart"         ,
                       "fittingStart"         ,
                       "orbitalStart"         ,
                       "poissonStart"         ,
                       "densityFunctions"     ,
                       "fittingFunctions"     ,
                       "orbitalFunctions"     ,
                       "poissonFunctions"     ,
                       "densityWorkStart"     , 
                       "fittingWorkStart"     , 
                       "orbitalWorkStart"     , 
                       "poissonWorkStart"     , 
                       "densityWorkFunctions" , 
                       "fittingWorkFunctions" , 
                       "orbitalWorkFunctions" , 
                       "poissonWorkFunctions" ]
        atoms      = []
        mmPartners = {}
        for i from 0 <= i < self.cObject.natoms:
            isBoundary = ( self.cObject.data[i].QBOUNDARY == CTrue )
            atoms.append ( [ isBoundary                        ,
                             self.cObject.data[i].linkfactor   ,
                             self.cObject.data[i].widthe       ,
                             self.cObject.data[i].widthn       ,
                             self.cObject.data[i].atomicNumber ,
                             self.cObject.data[i].center       ,
                             self.cObject.data[i].index        ,
                             self.cObject.data[i].qcpartner    ,
                             self.cObject.data[i].dstart       ,
                             self.cObject.data[i].fstart       ,
                             self.cObject.data[i].ostart       ,
                             self.cObject.data[i].pstart       ,
                             self.cObject.data[i].ndbasis      ,
                             self.cObject.data[i].nfbasis      ,
                             self.cObject.data[i].nobasis      ,
                             self.cObject.data[i].npbasis      ,
                             self.cObject.data[i].dstartw      ,
                             self.cObject.data[i].fstartw      ,
                             self.cObject.data[i].ostartw      ,
                             self.cObject.data[i].pstartw      ,
                             self.cObject.data[i].ndbasisw     ,
                             self.cObject.data[i].nfbasisw     ,
                             self.cObject.data[i].nobasisw     ,
                             self.cObject.data[i].npbasisw   ] )
            if self.cObject.data[i].nmmpartners > 0:
                partners      = []
                for m from 0 <= m < self.cObject.data[i].nmmpartners:
                    partners.append ( self.cObject.data[i].mmpartners[m] )
                mmPartners[i] = partners
        # . Other data.
        if self.cObject.QLINKRATIO   == CTrue: QLINKRATIO   = True
        else:                                  QLINKRATIO   = False
        if self.cObject.QTOSPHERICAL == CTrue: QTOSPHERICAL = True
        else:                                  QTOSPHERICAL = False
        # . Finish up.
        return { "atomFields"            : atomFields                  ,
                 "atoms"                 : atoms                       ,
                 "boundaryAtoms"         : self.cObject.nboundary      ,
                 "densityFunctions"      : self.cObject.ndbasis        ,
                 "densityWorkFunctions"  : self.cObject.ndbasisw       ,
                 "energyBaseLine"        : self.cObject.energybaseline ,
                 "fittingFunctions"      : self.cObject.nfbasis        ,
                 "fittingWorkFunctions"  : self.cObject.nfbasisw       ,
                 "mmPartners"            : mmPartners                  ,
                 "nuclearCharge"         : self.cObject.nuclearCharge  ,
                 "orbitalFunctions"      : self.cObject.nobasis        ,
                 "orbitalWorkFunctions"  : self.cObject.nobasisw       ,
                 "poissonFunctions"      : self.cObject.npbasis        ,
                 "poissonWorkFunctions"  : self.cObject.npbasisw       ,
                 "useLinkRatio"          : QLINKRATIO                  ,
                 "useSphericalFunctions" : QTOSPHERICAL                }

    def __init__ ( self, numberOfAtoms ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfAtoms )

    def __len__ ( self ):
        """Length."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef i, m, n
        # . Allocate the object.
        atoms      = state["atoms"     ]
        mmPartners = state["mmPartners"]
        self._Allocate ( len ( atoms ) )
        # . Fill the object.
        # . Atom data.
        for ( i, data ) in enumerate ( atoms ):
            if data[0]: self.cObject.data[i].QBOUNDARY = CTrue
            else:       self.cObject.data[i].QBOUNDARY = CFalse
            self.cObject.data[i].linkfactor   = data[ 1]
            self.cObject.data[i].widthe       = data[ 2]
            self.cObject.data[i].widthn       = data[ 3]
            self.cObject.data[i].atomicNumber = data[ 4]
            self.cObject.data[i].center       = data[ 5]
            self.cObject.data[i].index        = data[ 6]
            self.cObject.data[i].qcpartner    = data[ 7]
            self.cObject.data[i].dstart       = data[ 8]
            self.cObject.data[i].fstart       = data[ 9]
            self.cObject.data[i].ostart       = data[10]
            self.cObject.data[i].pstart       = data[11]
            self.cObject.data[i].ndbasis      = data[12]
            self.cObject.data[i].nfbasis      = data[13]
            self.cObject.data[i].nobasis      = data[14]
            self.cObject.data[i].npbasis      = data[15]
            self.cObject.data[i].dstartw      = data[16]
            self.cObject.data[i].fstartw      = data[17]
            self.cObject.data[i].ostartw      = data[18]
            self.cObject.data[i].pstartw      = data[19]
            self.cObject.data[i].ndbasisw     = data[20]
            self.cObject.data[i].nfbasisw     = data[21]
            self.cObject.data[i].nobasisw     = data[22]
            self.cObject.data[i].npbasisw     = data[23]
            partners = mmPartners.get ( i, {} )
            n        = len ( partners )
            if n > 0:
                self.cObject.data[i].nmmpartners = n
                self.cObject.data[i].mmpartners  = Memory_Allocate_Array_Integer ( n )
                for m from 0 <= m < n: self.cObject.data[i].mmpartners[m] = partners[m]
        # . Scalars.
        self.cObject.energybaseline = state["energyBaseLine"      ]
        self.cObject.nboundary      = state["boundaryAtoms"       ]
        self.cObject.ndbasis        = state["densityFunctions"    ]
        self.cObject.nfbasis        = state["fittingFunctions"    ]
        self.cObject.nobasis        = state["orbitalFunctions"    ]
        self.cObject.npbasis        = state["poissonFunctions"    ]
        self.cObject.ndbasisw       = state["densityWorkFunctions"]
        self.cObject.nfbasisw       = state["fittingWorkFunctions"]
        self.cObject.nobasisw       = state["orbitalWorkFunctions"]
        self.cObject.npbasisw       = state["poissonWorkFunctions"]
        self.cObject.nuclearCharge  = state["nuclearCharge"       ]
        if state["useLinkRatio"         ]: self.cObject.QLINKRATIO   = CTrue
        else:                              self.cObject.QLINKRATIO   = CFalse
        if state["useSphericalFunctions"]: self.cObject.QTOSPHERICAL = CTrue
        else:                              self.cObject.QTOSPHERICAL = CFalse

    def _Allocate ( self, numberOfAtoms ):
        """Allocation."""
        self.cObject = QCAtomContainer_Allocate ( numberOfAtoms )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def AssignParameterData ( self, object atomicNumberlist, QCParameters qcParameters ):
        """Assign parameter data to the atoms."""
        # . Declarations.
        cdef Integer ic, nd, nf, no, np
        # . Assign parameter data to the atoms.
        for i from 0 <= i < self.cObject.natoms:
            ic = atomicNumberlist.index ( self.cObject.data[i].atomicNumber )
            self.cObject.data[i].center = ic
            # . QC Parameters exist.
            if qcParameters is not None:
                if qcParameters.cObject.centers[ic].densitybasis == NULL:
                    nd  = 0
                    ndw = 0
                else:
                    nd  = qcParameters.cObject.centers[ic].densitybasis.nbasis
                    ndw = qcParameters.cObject.centers[ic].densitybasis.nbasisw
                if qcParameters.cObject.centers[ic].orbitalbasis == NULL:
                    no  = 0
                    now = 0
                else:
                    no  = qcParameters.cObject.centers[ic].orbitalbasis.nbasis
                    now = qcParameters.cObject.centers[ic].orbitalbasis.nbasisw
                if qcParameters.cObject.centers[ic].poissonbasis == NULL:
                    np  = 0
                    npw = 0
                else:
                    np  = qcParameters.cObject.centers[ic].poissonbasis.nbasis
                    npw = qcParameters.cObject.centers[ic].poissonbasis.nbasisw
                nf  = nd  + np
                nfw = ndw + npw
                self.cObject.data[i].dstart   = self.cObject.nfbasis
                self.cObject.data[i].fstart   = self.cObject.nfbasis
                self.cObject.data[i].ostart   = self.cObject.nobasis
                self.cObject.data[i].pstart   = self.cObject.nfbasis + nd
                self.cObject.data[i].ndbasis  = nd
                self.cObject.data[i].nfbasis  = nf
                self.cObject.data[i].nobasis  = no
                self.cObject.data[i].npbasis  = np
                self.cObject.data[i].dstartw  = self.cObject.nfbasisw
                self.cObject.data[i].fstartw  = self.cObject.nfbasisw
                self.cObject.data[i].ostartw  = self.cObject.nobasisw
                self.cObject.data[i].pstartw  = self.cObject.nfbasisw + ndw
                self.cObject.data[i].ndbasisw = ndw
                self.cObject.data[i].nfbasisw = nfw
                self.cObject.data[i].nobasisw = now
                self.cObject.data[i].npbasisw = npw
                self.cObject.ndbasis  = self.cObject.ndbasis  + nd
                self.cObject.nfbasis  = self.cObject.nfbasis  + nf
                self.cObject.nobasis  = self.cObject.nobasis  + no
                self.cObject.npbasis  = self.cObject.npbasis  + np
                self.cObject.ndbasisw = self.cObject.ndbasisw + ndw
                self.cObject.nfbasisw = self.cObject.nfbasisw + nfw
                self.cObject.nobasisw = self.cObject.nobasisw + now
                self.cObject.npbasisw = self.cObject.npbasisw + npw
                if qcParameters.cObject.centers[ic].mndoparameters == NULL:
                    self.cObject.nuclearCharge  = self.cObject.nuclearCharge  + qcParameters.cObject.centers[ic].atomicNumber
                else:
                    self.cObject.nuclearCharge  = self.cObject.nuclearCharge  + int ( round ( qcParameters.cObject.centers[ic].mndoparameters.zcore ) )
                    self.cObject.energybaseline = self.cObject.energybaseline + ( qcParameters.cObject.centers[ic].mndoparameters.eheat - qcParameters.cObject.centers[ic].mndoparameters.eisol )
            # . No QC parameters.
            else:
                self.cObject.nuclearCharge = self.cObject.nuclearCharge + self.cObject.data[i].atomicNumber

    def BoundaryAtomQCPartnerPairs ( self ):
        """Return a list of tuples consisting of the indices of the boundary atoms and their QC partners."""
        data = []
        for i from 0 <= i < self.cObject.natoms:
            if self.cObject.data[i].QBOUNDARY == CTrue:
                data.append ( ( self.cObject.data[i].index, self.cObject.data[i].qcpartner ) )
        return data

    def BoundaryAtomSelection ( self ):
        """Return a selection with the boundary atom indices."""
        data = []
        for i from 0 <= i < self.cObject.natoms:
            if self.cObject.data[i].QBOUNDARY == CTrue: data.append ( self.cObject.data[i].index )
        return Selection.FromIterable ( data )

    def GetAtomicNumbers ( self ):
        """Return the atomic numbers of the QC atoms (including boundary atoms)."""
        cdef Integer i
        atomicNumbers = []
        for i from 0 <= i < self.cObject.natoms:
            atomicNumbers.append ( self.cObject.data[i].atomicNumber )
        return atomicNumbers

    def FillCoordinates3 ( QCAtomContainer self, Coordinates3 coordinates3, Coordinates3 qcCoordinates3, toBohrs = False ):
        """Fill up the coordinates for the quantum atoms in Angstroms or atomic units."""
        cdef Boolean toInternalUnits
        if toBohrs: toInternalUnits = CTrue
        else:       toInternalUnits = CFalse
        QCAtomContainer_GetCoordinates3 ( self.cObject, coordinates3.cObject, toInternalUnits, &(qcCoordinates3.cObject) )

    def GetCoordinates3 ( QCAtomContainer self, Coordinates3 coordinates3, toBohrs = False ):
        """Return the coordinates for the quantum atoms in Angstroms or atomic units."""
        cdef Boolean         toInternalUnits
        cdef Coordinates3 qcCoordinates3
        if toBohrs: toInternalUnits = CTrue
        else:       toInternalUnits = CFalse
        qcCoordinates3         = Coordinates3.Raw ( )
        qcCoordinates3.isOwner = True
        QCAtomContainer_GetCoordinates3 ( self.cObject, coordinates3.cObject, toInternalUnits, &(qcCoordinates3.cObject) )
        return qcCoordinates3

    def GetFullSelection ( self ):
        """Return a selection with the indices of the QC atoms (including boundary atoms)."""
        cdef Selection selection
        selection         = Selection.Raw ( )
        selection.cObject = QCAtomContainer_MakeFullSelection ( self.cObject )
        selection.isOwner = True
        return selection

    def GetOrbitalBasisFunctionLabels ( self, QCParameters qcParameters, workBasis = False ):
        """Get labels for the orbital basis functions."""
        cdef CGaussianBasis *basis
        cdef Integer             i, ic, ishell, l
        labels = None
        if self.cObject != NULL:
            # . Number field length.
            nlength = len ( "{:d}".format ( self.cObject.natoms ) )
            # . Basis function and element labels.
            bfs     = []
            slength =  0
            for i from 0 <= i < self.cObject.natoms:
                ic          = self.cObject.data[i].center
                basis       = qcParameters.cObject.centers[ic].orbitalbasis
                isSpherical = ( basis.QSPHERICAL == CTrue ) and ( not workBasis )
                symbol      = PeriodicTable.Symbol ( self.cObject.data[i].atomicNumber )
                slength     = max ( len ( symbol ), slength )
                for ishell from 0 <= ishell < basis.nshells:
                    for l from basis.shells[ishell].type.angularmomentum_low <= l <= basis.shells[ishell].type.angularmomentum_high:
                        shellsymbol = _ShellLabels[l]
                        if l == 0:
                            bfs.append ( ( symbol, shellsymbol ) )
                        else:
                            if isSpherical:
                                bfs.append ( ( symbol, "{:1s}{:1d}{:1d}".format ( shellsymbol, l, 0 ) ) )
                                for m in range ( 1, l + 1 ):
                                    bfs.append ( ( symbol, "{:1s}{:1d}{:1d}" .format ( shellsymbol, l, m ) ) )
                                    bfs.append ( ( symbol, "{:1s}{:1d}{:1d}'".format ( shellsymbol, l, m ) ) )
                            else:
                                for z in range ( 0, l + 1 ):
                                    for y in range ( 0, l - z + 1 ):
                                        x = l - y - z
                                        bfs.append ( ( symbol, "".join ( [ shellsymbol ] + x * [ "x" ] + y * [ "y" ] + z * [ "z" ] ) ) )
            # . Finish up the labels.
            labels = []
            for ( n, ( s, b ) ) in enumerate ( bfs ):
                labels.append ( "{:d}".format ( n ).rjust ( nlength ) + " " + s.ljust ( slength ) + " " + b )
        return labels

    def NumberOfBoundaryAtoms ( self ):
        """Return the number of boundary atoms."""
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nboundary

    def QCAtomSelection ( self ):
        """Return a selection with the QC atom indices."""
        cdef Integer i
        data = []
        for i from 0 <= i < self.cObject.natoms: data.append ( self.cObject.data[i].index )
        data.sort ( )
        return Selection.FromIterable ( data )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetGradients3 ( QCAtomContainer self, Coordinates3 coordinates3, Coordinates3 qcgradients3, Coordinates3 gradients3, toInternalUnits = False ):
        """Convert QC gradients3 to gradients3 with optional unit conversion (atomic units -> pMolecule units)."""
        cdef Boolean QC
        if toInternalUnits: QC = CTrue
        else:               QC = CFalse
        QCAtomContainer_SetGradients3 ( self.cObject, coordinates3.cObject, qcgradients3.cObject, QC, &(gradients3.cObject) )

    def SummaryEntry ( self, object summary ):
        """Summary entry."""
        if summary is not None:
            summary.Entry ( "Number of QC Atoms", "{:d}"  .format ( self.cObject.natoms         ) )
            summary.Entry ( "Boundary Atoms"    , "{:d}"  .format ( self.cObject.nboundary      ) )
            summary.Entry ( "Nuclear Charge"    , "{:d}"  .format ( self.cObject.nuclearCharge  ) )
            summary.Entry ( "Orbital Functions" , "{:d}"  .format ( self.cObject.nobasis        ) )
            summary.Entry ( "Fitting Functions" , "{:d}"  .format ( self.cObject.nfbasis        ) )
            summary.Entry ( "Energy Base Line"  , "{:.5f}".format ( self.cObject.energybaseline ) )

    def UniqueAtomicNumbers ( self ):
        """Return the unique atomic numbers for the atoms as an ordered list."""
        anlist = []
        for i from 0 <= i < self.cObject.natoms:
            if self.cObject.data[i].atomicNumber not in anlist: anlist.append ( self.cObject.data[i].atomicNumber )
        anlist.sort ( )
        return anlist

    # . Properties.
# . Keep?
    property natoms:
        def __get__ ( QCAtomContainer self ): return self.cObject.natoms
    property nobasis:
        def __get__ ( QCAtomContainer self ): return self.cObject.nobasis
    property nuclearCharge:
        def __get__ ( QCAtomContainer self ): return self.cObject.nuclearCharge
    property size:
        def __get__ ( self ):
            if self.cObject == NULL: return 0
            else:                    return self.cObject.natoms

#===================================================================================================================================
# . Class methods.
#===================================================================================================================================
def QCAtomContainer_FromAtomContainer ( atomcontainer, Selection qcSelection, boundaryatoms, QLINKRATIO ):
    """Constructor from an AtomContainer.

    |atomcontainer| is the atom container
    |qcSelection|   is the QC atom selection.
    |boundaryatoms| is a dictionary of MM boundary atoms with their boundary, MM and QC atom partners.
    """
    cdef i, iatom, m, n
    cdef QCAtomContainer self
    self = None
    # . The presence of boundary atoms makes no sense if there are no QC atoms.
    if ( atomcontainer is not None ) and ( qcSelection is not None ) and ( qcSelection.size > 0 ):
        # . Convert qcAtoms to a list.
        qcAtoms = list ( qcSelection )
        # . Process boundary atoms.
        # . Boundary atoms with multiple QC partners are not permitted.
        for ( b, partners ) in boundaryatoms.items ( ):
            if len ( partners[2] ) > 1: raise ValueError ( "A QC/MM boundary atom - {:d} - has multiple QC partners - {!r}.".format ( b, partners[2] ) )
            qcAtoms.append ( b )
        # . Sort the list.
        qcAtoms.sort ( )
        # . Allocate the container.
        self = QCAtomContainer.Raw ( )
        self._Allocate ( len ( qcAtoms ) )
        # . Set some other options.
        self.cObject.nboundary  = len ( boundaryatoms )
        if QLINKRATIO: self.cObject.QLINKRATIO = CTrue
        else:          self.cObject.QLINKRATIO = CFalse
        # . Assign the QC atom data.
        for ( i, iatom ) in enumerate ( qcAtoms ):
            self.cObject.data[i].index = iatom
            # . Boundary atom.
            if iatom in boundaryatoms:
                # . QC partner.
                partners  = boundaryatoms[iatom]
                qcpartner = partners[2].pop ( )
                self.cObject.data[i].QBOUNDARY    = CTrue
                self.cObject.data[i].atomicNumber = 1 # . Link atoms are hydrogens.
                self.cObject.data[i].qcpartner    = qcpartner
                # . MM partners.
                mmpartners = list ( partners[0].union ( partners[1] ) )
                mmpartners.sort ( )
                n = len ( mmpartners )
                if n > 0:
                    self.cObject.data[i].nmmpartners = n
                    self.cObject.data[i].mmpartners  = Memory_Allocate_Array_Integer ( n )
                    for m from 0 <= m < n: self.cObject.data[i].mmpartners[m] = mmpartners[m]
                # . Set the link factor.
                nqc = atomcontainer[qcpartner].atomicNumber
                eqc = PeriodicTable.Element ( nqc )
                rql = eqc.GetSingleBondDistance ( 1 )
                if rql is None: raise ValueError ( "Unknown link atom distance for elements with atomic numbers {:d} and 1.".format ( nqc ) )
                if QLINKRATIO:
                    nmm = atomcontainer[iatom].atomicNumber
                    rqm = eqc.GetSingleBondDistance ( nmm )
                    if rqm is None: raise ValueError ( "Unknown link atom distance for elements with atomic numbers {:d} and {:d}.".format ( nqc, nmm ) )
                    self.cObject.data[i].linkfactor = rql / rqm
                else:
                    self.cObject.data[i].linkfactor = rql
            # . Pure QC atom.
            else:
                self.cObject.data[i].atomicNumber = int ( atomcontainer[iatom].atomicNumber )
    return self
