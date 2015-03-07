#-------------------------------------------------------------------------------
# . File      : pMolecule.MNDOParameters.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Handle the parameter data necessary for a MNDO NDDO-type calculation."""

from pCore import RawObjectConstructor, UNITS_LENGTH_BOHRS_TO_ANGSTROMS

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Conversion factors.
cdef Real _BohrsToAngstroms
cdef Real _MNDOHartreesToElectronVolts
cdef Real _MNDOHartreesToKCal

_BohrsToAngstroms            = UNITS_LENGTH_BOHRS_TO_ANGSTROMS
_MNDOHartreesToElectronVolts = 27.2113834e+00                           # . New value.
_MNDOHartreesToKCal          = _MNDOHartreesToElectronVolts * 23.060529 # . New value.

# . The YAML tag.
#_YAMLTag = "!MNDOParameters"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOParameters:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef MNDOParameters new
        new         = self.__class__.Raw ( )
        new.cObject = MNDOParameters_Clone ( self.cObject )
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        if self.isOwner:
            MNDOParameters_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.MNDOParameters"

    def __getstate__ ( self ):
        """Return the state."""
        state = {}
        # . Scalars.
        parameters = [ [ "atomicNumber" , self.cObject.atomicNumber , "none"          ] ,
                       [ "orbitals"     , self.cObject.norbitals    , "none"          ] ,
                       [ "iii"          , self.cObject.iii          , "none"          ] ,
                       [ "iiid"         , self.cObject.iiid         , "none"          ] ,
                       [ "ir016"        , self.cObject.ir016        , "none"          ] ,
                       [ "ir066"        , self.cObject.ir066        , "none"          ] ,
                       [ "ir244"        , self.cObject.ir244        , "none"          ] ,
                       [ "ir266"        , self.cObject.ir266        , "none"          ] ,
                       [ "ir466"        , self.cObject.ir466        , "none"          ] ,
                       [ "qnd"          , self.cObject.qnd          , "none"          ] ,
                       [ "qnp"          , self.cObject.qnp          , "none"          ] ,
                       [ "qns"          , self.cObject.qns          , "none"          ] ,
                       [ "ad"           , self.cObject.ad0          , "atomic"        ] ,  
                       [ "alp"          , self.cObject.alp0         , "A^-1"          ] ,
                       [ "am"           , self.cObject.am0          , "atomic"        ] ,
                       [ "aq"           , self.cObject.aq0          , "atomic"        ] ,
                       [ "betad"        , self.cObject.betad0       , "eV"            ] ,
                       [ "betap"        , self.cObject.betap0       , "eV"            ] ,
                       [ "betas"        , self.cObject.betas0       , "eV"            ] ,
                       [ "dd"           , self.cObject.dd0          , "atomic"        ] ,
                       [ "eheat"        , self.cObject.eheat0       , "kcal/mole"     ] ,
                       [ "eisol"        , self.cObject.eisol0       , "eV"            ] ,
                       [ "f0sd"         , self.cObject.f0sd0        , "eV"            ] ,
                       [ "gphot"        , self.cObject.gphot0       , "dimensionless" ] ,
                       [ "gpp"          , self.cObject.gpp0         , "eV"            ] ,
                       [ "gp2"          , self.cObject.gp20         , "eV"            ] ,
                       [ "gsp"          , self.cObject.gsp0         , "eV"            ] ,
                       [ "gss"          , self.cObject.gss0         , "eV"            ] ,
                       [ "g2sd"         , self.cObject.g2sd0        , "eV"            ] ,
                       [ "hsp"          , self.cObject.hsp0         , "eV"            ] ,
                       [ "pcore"        , self.cObject.pcore0       , "atomic"        ] ,
                       [ "qq"           , self.cObject.qq0          , "atomic"        ] ,
                       [ "udd"          , self.cObject.udd0         , "eV"            ] ,
                       [ "upp"          , self.cObject.upp0         , "eV"            ] ,
                       [ "uss"          , self.cObject.uss0         , "eV"            ] ,
                       [ "zcore"        , self.cObject.zcore0       , "atomic"        ] ,
                       [ "zetad"        , self.cObject.zetad0       , "atomic"        ] ,
                       [ "zetap"        , self.cObject.zetap0       , "atomic"        ] ,
                       [ "zetas"        , self.cObject.zetas0       , "atomic"        ] ,
                       [ "zdn"          , self.cObject.zdn0         , "atomic"        ] ,
                       [ "zpn"          , self.cObject.zpn0         , "atomic"        ] ,
                       [ "zsn"          , self.cObject.zsn0         , "atomic"        ] ]
        state["Scalar Parameter Fields"] = [ "Label", "Value", "Units" ]
        state["Scalar Parameters"      ] = parameters
        # . Arrays.
        if self.cObject.nam1pm3g > 0:
            parameters = []
            for i from 0 <= i < self.cObject.nam1pm3g:
                parameters.append ( [ self.cObject.fn10[i] ,
                                      self.cObject.fn20[i] ,
                                      self.cObject.fn30[i] ] )
            state["AM1/PM3 Gaussian Parameter Fields"] = [ "fn1"  , "fn2"  , "fn3" ]
            state["AM1/PM3 Gaussian Parameter Units" ] = [ "A eV" , "A^-2" , "A"   ]
            state["AM1/PM3 Gaussian Parameters"      ] = parameters
        if self.cObject.ndiatomic > 0:
            parameters = []
            for i from 0 <= i < self.cObject.ndiatomic:
                if self.cObject.QDIATOMICFLAGS[i] == CTrue:
                    parameters.append ( [ i, self.cObject.diatomicx0[i], self.cObject.diatomica0[i] ] )
            state["Diatomic Parameter Fields"] = [ "atomicNumber" , "coefficient"   , "exponent" ]
            state["Diatomic Parameter Units" ] = [ "none"         , "dimensionless" , "A^-1"     ]
            state["Diatomic Parameters"      ] = parameters
        if self.cObject.npddg > 0:
            parameters = []
            for i from 0 <= i < self.cObject.npddg:
                parameters.append ( [ self.cObject.pddgc0[i], self.cObject.pddge0[i] ] )
            state["PDDG Gaussian Parameter Fields"] = [ "coefficient", "exponentDistance" ]
            state["PDDG Gaussian Parameter Units" ] = [ "eV"         , "A"                ]
            state["PDDG Gaussian Parameters"      ] = parameters
        return state

    def __init__ ( self, atomicNumber, numberOfOrbitals ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )
        self.cObject.atomicNumber = atomicNumber
        self.cObject.norbitals    = numberOfOrbitals

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        # . Allocate the object.
        self._Allocate ( )
        # . Fill the object.
        # . Scalars.
        parameters = state["Scalar Parameters"]
        for ( label, value, units ) in parameters:
            if   label == "atomicNumber" : self.cObject.atomicNumber = value
            elif label == "orbitals"     : self.cObject.norbitals    = value
            elif label == "iii"          : self.cObject.iii          = value
            elif label == "iiid"         : self.cObject.iiid         = value
            elif label == "ir016"        : self.cObject.ir016        = value
            elif label == "ir066"        : self.cObject.ir066        = value
            elif label == "ir244"        : self.cObject.ir244        = value
            elif label == "ir266"        : self.cObject.ir266        = value
            elif label == "ir466"        : self.cObject.ir466        = value
            elif label == "qnd"          : self.cObject.qnd          = value
            elif label == "qnp"          : self.cObject.qnp          = value
            elif label == "qns"          : self.cObject.qns          = value
            elif label == "ad"           : self.cObject.ad0          = value
            elif label == "alp"          : self.cObject.alp0         = value
            elif label == "am"           : self.cObject.am0          = value
            elif label == "aq"           : self.cObject.aq0          = value
            elif label == "betad"        : self.cObject.betad0       = value
            elif label == "betap"        : self.cObject.betap0       = value
            elif label == "betas"        : self.cObject.betas0       = value
            elif label == "dd"           : self.cObject.dd0          = value
            elif label == "eheat"        : self.cObject.eheat0       = value
            elif label == "eisol"        : self.cObject.eisol0       = value
            elif label == "f0sd"         : self.cObject.f0sd0        = value
            elif label == "gphot"        : self.cObject.gphot0       = value
            elif label == "gpp"          : self.cObject.gpp0         = value
            elif label == "gp2"          : self.cObject.gp20         = value
            elif label == "gsp"          : self.cObject.gsp0         = value
            elif label == "gss"          : self.cObject.gss0         = value
            elif label == "g2sd"         : self.cObject.g2sd0        = value
            elif label == "hsp"          : self.cObject.hsp0         = value
            elif label == "pcore"        : self.cObject.pcore0       = value
            elif label == "qq"           : self.cObject.qq0          = value
            elif label == "udd"          : self.cObject.udd0         = value
            elif label == "upp"          : self.cObject.upp0         = value
            elif label == "uss"          : self.cObject.uss0         = value
            elif label == "zcore"        : self.cObject.zcore0       = value
            elif label == "zetad"        : self.cObject.zetad0       = value
            elif label == "zetap"        : self.cObject.zetap0       = value
            elif label == "zetas"        : self.cObject.zetas0       = value
            elif label == "zdn"          : self.cObject.zdn0         = value
            elif label == "zpn"          : self.cObject.zpn0         = value
            elif label == "zsn"          : self.cObject.zsn0         = value
        # . Arrays.
        self.FillBetaUspd ( )
        # . AM1/PM3.
        parameters = state.get ( "AM1/PM3 Gaussian Parameters", [] )
        size       = len ( parameters )
        if size > 0: self.FillFN123 ( size, parameters )
        # . Diatomics.
        parameters = state.get ( "Diatomic Parameters", [] )
        size       = len ( parameters )
        if size > 0: self.FillDiatomic ( size, parameters )
        # . PDDG.
        parameters = state.get ( "PDDG Gaussian Parameters", [] )
        size       = len ( parameters )
        if size > 0: self.FillPDDG ( size, parameters )
        # . Finish up.
        self.ToAtomicUnits ( )
        # . Remaining data.
        MNDOParameters_CalculateOneCenterTEIs ( self.cObject )

    def _Allocate ( self ):
        """Allocation."""
        if self.cObject != NULL: MNDOParameters_Deallocate ( &self.cObject )
        self.cObject = MNDOParameters_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def DeallocateInternalArrays ( self ):
        """Deallocate the arrays containing data in atomic units."""
        if self.cObject != NULL:
            Memory_Deallocate_Real ( &self.cObject.beta      )
            Memory_Deallocate_Real ( &self.cObject.diatomica )
            Memory_Deallocate_Real ( &self.cObject.diatomicx )
            Memory_Deallocate_Real ( &self.cObject.beta      )
            Memory_Deallocate_Real ( &self.cObject.fn1       )
            Memory_Deallocate_Real ( &self.cObject.fn2       )
            Memory_Deallocate_Real ( &self.cObject.fn3       )
            Memory_Deallocate_Real ( &self.cObject.pddgc     )
            Memory_Deallocate_Real ( &self.cObject.pddge     )
            Memory_Deallocate_Real ( &self.cObject.uspd      )

    def FillBetaUspd ( self ):
        """Fill beta and uspd."""
        if self.cObject.norbitals > 0:
            self.cObject.beta0 = Memory_Allocate_Array_Real ( self.cObject.norbitals )
            self.cObject.uspd0 = Memory_Allocate_Array_Real ( self.cObject.norbitals )
            # . s.
            self.cObject.beta0[0] = self.cObject.betas0
            self.cObject.uspd0[0] = self.cObject.uss0
            # . p.
            if self.cObject.norbitals >= 4:
                for i from 1 <= i < 4:
                    self.cObject.beta0[i] = self.cObject.betap0
                    self.cObject.uspd0[i] = self.cObject.upp0
# #ifdef MNDODORBITALS
            # . d.
            if self.cObject.norbitals >= 9:
                for i from 4 <= i < 9:
                    self.cObject.beta0[i] = self.cObject.betad0
                    self.cObject.uspd0[i] = self.cObject.udd0
# #else /*MNDODORBITALS*/
#            # . d.
#            if self.cObject.norbitals >= 9:
#                raise ValueError ( "MNDO d-orbitals unavailable for element with atomic number " + "{:d}".format ( self.cobject.atomicnumber ) + "." )
# #endif /*MNDODORBITALS*/

    def FillDiatomic ( self, Integer nterms, object data ):
        """Fill QDIATOMICFLAGS, diatomica, diatomicx."""
        cdef Real a, x
        cdef Integer    j
        if nterms > 0:
            # . Find maximum atomic number.
            data.sort ( )
            maximumAtomicNumber = data[-1][0]
            # . Allocate and fill object.
            self.cObject.QDIATOMIC = CTrue
            self.cObject.ndiatomic = maximumAtomicNumber+1
            self.cObject.QDIATOMICFLAGS = Memory_Allocate_Array_Boolean_Initialize ( self.cObject.ndiatomic, CFalse   )
            self.cObject.diatomica0     = Memory_Allocate_Array_Real_Initialize  ( self.cObject.ndiatomic, 0.0e+00 )
            self.cObject.diatomicx0     = Memory_Allocate_Array_Real_Initialize  ( self.cObject.ndiatomic, 0.0e+00 )
            for ( j, x, a ) in data:
                self.cObject.QDIATOMICFLAGS[j] = CTrue
                self.cObject.diatomica0    [j] = a
                self.cObject.diatomicx0    [j] = x

    def FillFN123 ( self, Integer nterms, object data ):
        """Fill fn1, fn2 and fn3."""
        cdef Integer i
        if nterms > 0:
            self.cObject.nam1pm3g = nterms
            self.cObject.fn10 = Memory_Allocate_Array_Real ( self.cObject.nam1pm3g )
            self.cObject.fn20 = Memory_Allocate_Array_Real ( self.cObject.nam1pm3g )
            self.cObject.fn30 = Memory_Allocate_Array_Real ( self.cObject.nam1pm3g )
            for ( i, ( fn1, fn2, fn3 ) ) in enumerate ( data ):
                self.cObject.fn10[i] = fn1
                self.cObject.fn20[i] = fn2
                self.cObject.fn30[i] = fn3

    def FillPDDG ( self, Integer nterms, object data ):
        """Fill pddgc and pddge."""
        cdef Integer i
        if nterms > 0:
            self.cObject.npddg  = nterms
            self.cObject.pddgc0 = Memory_Allocate_Array_Real ( self.cObject.npddg )
            self.cObject.pddge0 = Memory_Allocate_Array_Real ( self.cObject.npddg )
            for ( i, ( c, e ) ) in enumerate ( data ):
                self.cObject.pddgc0[i] = c
                self.cObject.pddge0[i] = e

    @classmethod
    def Raw ( selfClass ):
        """Constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def ToAtomicUnits ( self ):
        """Convert parameters to atomic units from MNDO standard units (A, eV)."""
        cdef Integer i
        self.DeallocateInternalArrays ( )
        # . ad, am, aq, dd, qq, zcore, zetap and zetas are already in atomic units.
        self.cObject.ad    = self.cObject.ad0
        self.cObject.alp   = self.cObject.alp0   * _BohrsToAngstroms # . A^-1.
        self.cObject.am    = self.cObject.am0
        self.cObject.aq    = self.cObject.aq0
        self.cObject.betad = self.cObject.betad0 / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.betap = self.cObject.betap0 / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.betas = self.cObject.betas0 / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.dd    = self.cObject.dd0
        self.cObject.eheat = self.cObject.eheat0 / _MNDOHartreesToKCal   # . kcal mol^-1.
        self.cObject.eisol = self.cObject.eisol0 / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.f0sd  = self.cObject.f0sd0  / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.gphot = self.cObject.gphot0
        self.cObject.gpp   = self.cObject.gpp0   / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.gp2   = self.cObject.gp20   / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.gsp   = self.cObject.gsp0   / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.gss   = self.cObject.gss0   / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.g2sd  = self.cObject.g2sd0  / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.hsp   = self.cObject.hsp0   / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.pcore = self.cObject.pcore0
        self.cObject.qq    = self.cObject.qq0
        self.cObject.udd   = self.cObject.udd0   / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.upp   = self.cObject.upp0   / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.uss   = self.cObject.uss0   / _MNDOHartreesToElectronVolts     # . eV.
        self.cObject.zcore = self.cObject.zcore0
        self.cObject.zetad = self.cObject.zetad0
        self.cObject.zetap = self.cObject.zetap0
        self.cObject.zetas = self.cObject.zetas0
        self.cObject.zdn   = self.cObject.zdn0
        self.cObject.zpn   = self.cObject.zpn0
        self.cObject.zsn   = self.cObject.zsn0
        if self.cObject.norbitals > 0:
            self.cObject.beta = Memory_Allocate_Array_Real ( self.cObject.norbitals )
            self.cObject.uspd = Memory_Allocate_Array_Real ( self.cObject.norbitals )
            for i from 0 <= i < self.cObject.norbitals:
                self.cObject.beta[i] = self.cObject.beta0[i] / _MNDOHartreesToElectronVolts # . eV.
                self.cObject.uspd[i] = self.cObject.uspd0[i] / _MNDOHartreesToElectronVolts # . eV.
        if self.cObject.nam1pm3g > 0:
            self.cObject.fn1 = Memory_Allocate_Array_Real ( self.cObject.nam1pm3g )
            self.cObject.fn2 = Memory_Allocate_Array_Real ( self.cObject.nam1pm3g )
            self.cObject.fn3 = Memory_Allocate_Array_Real ( self.cObject.nam1pm3g )
            for i from 0 <= i < self.cObject.nam1pm3g:
                self.cObject.fn1[i] = self.cObject.fn10[i] / ( _BohrsToAngstroms * _MNDOHartreesToElectronVolts )  # . A eV.
                self.cObject.fn2[i] = self.cObject.fn20[i] * _BohrsToAngstroms * _BohrsToAngstroms # . A^-2.
                self.cObject.fn3[i] = self.cObject.fn30[i] / _BohrsToAngstroms                     # . A.
        if self.cObject.ndiatomic > 0:
            self.cObject.diatomica = Memory_Allocate_Array_Real_Initialize ( self.cObject.ndiatomic, 0.0e+00 )
            self.cObject.diatomicx = Memory_Allocate_Array_Real_Initialize ( self.cObject.ndiatomic, 0.0e+00 )
            for i from 0 <= i < self.cObject.ndiatomic:
                self.cObject.diatomica[i] = self.cObject.diatomica0[i] * _BohrsToAngstroms # . A^-1.
                self.cObject.diatomicx[i] = self.cObject.diatomicx0[i]                     # . Dimensionless.
        if self.cObject.npddg > 0:
            self.cObject.pddgc = Memory_Allocate_Array_Real ( self.cObject.npddg )
            self.cObject.pddge = Memory_Allocate_Array_Real ( self.cObject.npddg )
            for i from 0 <= i < self.cObject.npddg:
                self.cObject.pddgc[i] = self.cObject.pddgc0[i] / _MNDOHartreesToElectronVolts # . eV.
                self.cObject.pddge[i] = self.cObject.pddge0[i] / _BohrsToAngstroms            # . A.

    @classmethod
    def Uninitialized ( selfClass, atomicNumber, numberOfOrbitals ):
        """Constructor."""
        return selfClass ( atomicNumber, numberOfOrbitals )
