#-------------------------------------------------------------------------------
# . File      : pMolecule.LJParameterContainer.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A container for Lennard-Jones parameters."""

from pCore import Clone, logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Keywords for analytic forms.
_AMBERKEYWORD = "AMBER"
_OPLSKEYWORD  = "OPLS"

# . The set of analytic keywords.
_ANALYTICFORMS = [ _AMBERKEYWORD, _OPLSKEYWORD ]

# . Keyword indicating that epsilon/sigma <-> table conversion is not possible.
_TABLEKEYWORD = "TABLE"

# . All keywords.
_FORMKEYWORDS = _ANALYTICFORMS + [ _TABLEKEYWORD ]

# . The default form.
_DEFAULTFORM = _OPLSKEYWORD

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class LJParameterContainer:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef LJParameterContainer new
        new               = self.__class__.Raw ( )
        new.analyticForm  = self.analyticForm
        new.cObject       = LJParameterContainer_Clone ( self.cObject )
        new.isOwner       = True
        new.label         = self.label
        new.parameterKeys = Clone ( self.parameterKeys )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            LJParameterContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.LJParameterContainer"

    def __getstate__ ( self ):
        """Return the state."""
        # . Initialization.
        ntypes = self.cObject.ntypes
        a      = []
        b      = []
        # . Get the data.
        state = { "analyticForm"  : self.analyticForm ,
                  "numberOfTypes" : ntypes            }
        if self.parameterKeys is not None:
            state["parameterKeys"] = self.parameterKeys
        if self.analyticForm == _TABLEKEYWORD:
            for i from 0 <= i < ( ntypes * ( ntypes + 1 ) ) // 2:
                a.append ( self.cObject.tableA[i] )
                b.append ( self.cObject.tableB[i] )
            state["aCoefficients"] = a
            state["bCoefficients"] = b
        else:
            for i from 0 <= i < ntypes:
                a.append ( self.cObject.epsilon[i] )
                b.append ( self.cObject.sigma[i]   )
            state["epsilons"] = a
            state["sigmas"  ] = b
        if self.label is not None: state["label"] = self.label
        return state

    def __init__ ( self, analyticForm = _DEFAULTFORM, label = "LJ Parameters" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate   ( )
        self.label   = label
        if analyticForm.upper ( ) in _FORMKEYWORDS: self.analyticForm = analyticForm.upper ( )
        else: raise ValueError ( "Unrecognized analytic form: {:s}.".format ( analyticForm ) )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        # . Allocate the object.
        ntypes             = state["numberOfTypes"]
        self.analyticForm  = state["analyticForm" ]
        self.cObject       = LJParameterContainer_Allocate ( ntypes )
        self.isOwner       = True
        self.parameterKeys = state.get ( "parameterKeys", None )
        # . Fill the object.
        if self.analyticForm == _TABLEKEYWORD:
            a = state["aCoefficients"]
            b = state["bCoefficients"]
            for i from 0 <= i < ( ntypes * ( ntypes + 1 ) ) // 2:
                self.cObject.tableA[i] = a[i]
                self.cObject.tableB[i] = b[i]
            self.MakeEpsilonSigma ( )
        else:
            a = state["epsilons"]
            b = state["sigmas"  ]
            for i from 0 <= i < ntypes:
                self.cObject.epsilon[i] = a[i]
                self.cObject.sigma[i]   = b[i]
            self.MakeTable ( )
        if "label" in state: self.label = state["label"]

    def _Allocate ( self, numberOfTypes ):
        """Allocation."""
        self.cObject = LJParameterContainer_Allocate ( numberOfTypes )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.analyticForm  = _DEFAULTFORM
        self.cObject       = NULL
        self.isOwner       = False
        self.label         = "LJ Parameters"
        self.parameterKeys = None

    def MakeEpsilonSigma ( self ):
        """Create the epsilon and sigma representation from the table representation."""
        if   self.analyticForm == _AMBERKEYWORD: LJParameterContainer_MakeEpsilonSigmaAMBER ( self.cObject )
        elif self.analyticForm == _OPLSKEYWORD:  LJParameterContainer_MakeEpsilonSigmaOPLS  ( self.cObject )

    def MakeTable ( self ):
        """Create the table representation from the epsilon and sigma representation."""
        if   self.analyticForm == _AMBERKEYWORD: LJParameterContainer_MakeTableAMBER ( self.cObject )
        elif self.analyticForm == _OPLSKEYWORD:  LJParameterContainer_MakeTableOPLS  ( self.cObject )

    def Merge ( self, items, information = {} ):
        """Merging - slow version.

	Merging can only be done (for the moment) if all items have the same non-table analytic form.
	"""
        # . Initialization.
        increments = None
        new        = None
        # . Check the labels and analytic forms.
        QOK = ( self.analyticForm != _TABLEKEYWORD )
        for item in items:
            if ( item.analyticForm != self.analyticForm ) or ( item.label != self.label ):
                QOK = False
                break
        # . Do the merge.
        if QOK:
            # . Get the old indices.
            doKeys        = True
            epsilons      = []
            increments    = []
            n             = 0
            parameterKeys = []
            sigmas        = []
            for item in [ self ] + items:
                state = item.__getstate__ ( )
                n0    = state["numberOfTypes"]
                epsilons.extend   ( state["epsilons"] )
                sigmas.extend     ( state["sigmas"  ] )
                increments.append ( n )
                n = n + n0
                parameterKeys0 = state.get ( "parameterKeys", None )
                if parameterKeys0 is None: doKeys = False
                else: parameterKeys.extend ( parameterKeys0 )
            # . Construct the new keys and parameters.
#            print "Old>", parameterKeys, epsilons, sigmas, len ( epsilons )
            if doKeys:
                ( doKeys, oldToNew, parameterKeys, epsilons, sigmas ) = self.MergeKeys ( parameterKeys, epsilons, sigmas )
            else:
                parameterKeys = None
#            print "New>", parameterKeys, epsilons, sigmas, len ( epsilons )
            # . Allocate the object.
            new = self.__class__.Raw ( )
            new.__setstate__ ( { "analyticForm"  : self.analyticForm ,
                                 "epsilons"      : epsilons          , 
                                 "label"         : self.label        ,
                                 "numberOfTypes" : len ( epsilons )  , 
                                 "parameterKeys" : parameterKeys     ,
                                 "sigmas"        : sigmas          } )
            # . Save data.
            if oldToNew is None: oldToNew = [ i for i in range ( n ) ]
            information["ljTypeIncrements"] = increments
            information["ljTypeMapping"   ] = oldToNew
#            print "Information>", increments, oldToNew
        return new

    def MergeKeys ( self, parameterKeys, epsilons, sigmas ):
        """Merge parameter keys."""
        newKeys = set ( parameterKeys )
        if len ( newKeys ) < len ( parameterKeys ):
            # . Keys.
            newKeys = list ( newKeys )
            newKeys.sort ( )
            oldToNew = []
            for oldKey in parameterKeys:
                oldToNew.append ( newKeys.index ( oldKey ) )
            # . Parameters.
            newEpsilons = []
            newSigmas   = []
            for newKey in newKeys:
                index = parameterKeys.index ( newKey )
                newEpsilons.append ( epsilons[index] )
                newSigmas.append   ( sigmas  [index] )
            # . Finish up.
            return ( True , oldToNew, newKeys      , newEpsilons, newSigmas )
        else:
            return ( False, None    , parameterKeys, epsilons   , sigmas    )

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        return self.__copy__ ( )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Scale ( self, Real scale ):
        """Scaling."""
        LJParameterContainer_Scale ( self.cObject, scale )

    def SummaryEntry ( self, summary ):
        """Summary entry."""
        if summary is not None:
            summary.Entry ( self.label + " Form"  ,  self.analyticForm          )
            summary.Entry ( self.label + " Types" , "{:d}".format ( self.size ) )

    # . Properties.
    property size:
        def __get__ ( self ):
            if self.cObject == NULL: return 0
            else:                    return self.cObject.ntypes

#===================================================================================================================================
# . Class methods.
#===================================================================================================================================
def LJParameterContainer_FromEpsilonSigma ( ntypes, epsilons, sigmas, analyticForm = _DEFAULTFORM, parameterKeys = None ):
    """Constructor from epsilon and sigma.

    |ntypes| is the number of types.
    |epsilons| and |sigmas| are the energy and distance parameters, respectively.
    """
    cdef Integer i
    cdef LJParameterContainer self
    self  = None
    if ( ntypes > 0 ) and ( len ( epsilons ) == ntypes ) and ( len ( sigmas ) == ntypes ):
        self = LJParameterContainer.Raw ( )
        self.analyticForm = analyticForm
        self._Allocate ( ntypes )
        for i from 0 <= i < ntypes:
            self.cObject.epsilon[i] = epsilons[i]
            self.cObject.sigma[i]   = sigmas[i]
        self.MakeTable ( )
        self.parameterKeys = parameterKeys
    return self

def LJParameterContainer_FromTableCoefficients ( ntypes, acoefficients, bcoefficients, analyticForm = _DEFAULTFORM, parameterKeys = None ):
    """Constructor from table coefficients.

    |ntypes| is the number of types.
    |acoefficients| and |bcoefficients| are the coefficients for the r^12 and r^6 terms, respectively.
    """
    cdef Integer i
    cdef LJParameterContainer self
    self  = None
    nsize = ( ntypes * ( ntypes + 1 ) ) // 2
    if ( ntypes > 0 ) and ( len ( acoefficients ) == nsize ) and ( len ( bcoefficients ) == nsize ):
        self = LJParameterContainer.Raw ( )
        self.analyticForm = analyticForm
        self._Allocate ( ntypes )
        for i from 0 <= i < nsize:
            self.cObject.tableA[i] = acoefficients[i]
            self.cObject.tableB[i] = bcoefficients[i]
        self.MakeEpsilonSigma ( )
        self.parameterKeys = parameterKeys
    return self
