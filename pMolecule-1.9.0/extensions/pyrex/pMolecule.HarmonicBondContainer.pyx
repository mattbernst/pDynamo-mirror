#-------------------------------------------------------------------------------
# . File      : pMolecule.HarmonicBondContainer.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A container for harmonic bonds."""

from pCore import Clone, logFile, LogFileActive

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class HarmonicBondContainer ( MMTerm ):

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef HarmonicBondContainer new
        new               = self.__class__.Raw ( )
        new.label         = self.label
        new.cObject       = HarmonicBondContainer_Clone ( self.cObject )
        new.isOwner       = True
        new.parameterKeys = Clone ( self.parameterKeys )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            HarmonicBondContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.HarmonicBondContainer"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i
        # . Parameters.
        parameterFields = [ "equilibriumValue", "forceConstant" ]
        parameters      = []
        for i from 0 <= i < self.cObject.nparameters:
            parameters.append ( [ self.cObject.parameters[i].eq ,
                                  self.cObject.parameters[i].fc ] )
        # . Terms.
        termFields = [ "atom1", "atom2", "parameter", "isActive" ]
        terms      = []
        for i from 0 <= i < self.cObject.nterms:
            terms.append ( [ self.cObject.terms[i].atom1 ,
                             self.cObject.terms[i].atom2 ,
                             self.cObject.terms[i].type  ,
                             self.cObject.terms[i].QACTIVE == CTrue ] )
        # . State.
        state = { "is12Interaction" : self.is12Interaction ,
                  "label"           : self.label           ,
                  "parameterFields" : parameterFields      ,
                  "parameters"      : parameters           ,
                  "termFields"      : termFields           ,
                  "terms"           : terms                }
        if self.parameterKeys is not None: state["parameterKeys"] = self.parameterKeys
        return state

    def __init__ ( self, numberOfParameters, numberOfTerms, is12Interaction = True, label = "Harmonic Bond" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfParameters, numberOfTerms )
        self.is12Interaction = is12Interaction
        self.label           = label

    def _Initialize ( self ):
        """Initialization."""
        self.cObject         = NULL
        self.isOwner         = False
        self.is12Interaction = True
        self.label           = "Harmonic Bond"
        self.parameterKeys   = None

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer i
        # . Allocate the object.
        parameters = state["parameters"]
        terms      = state["terms"     ]
        self._Allocate ( len ( parameters ), len ( terms ) )
        self.is12Interaction = state.get ( "is12Interaction" , True            )
        self.label           = state.get ( "label"           , "Harmonic Bond" )
        self.parameterKeys   = state.get ( "parameterKeys"   , None            )
        # . Fill the object.
        # . Parameters.
        for ( i, ( eq, fc ) ) in enumerate ( parameters ):
            self.cObject.parameters[i].eq = eq
            self.cObject.parameters[i].fc = fc
        # . Terms.
        for ( i, ( a1, a2, p, q ) ) in enumerate ( terms ):
            self.cObject.terms[i].atom1 = a1
            self.cObject.terms[i].atom2 = a2
            self.cObject.terms[i].type  = p
            if q: self.cObject.terms[i].QACTIVE = CTrue
            else: self.cObject.terms[i].QACTIVE = CFalse

    def _Allocate ( self, numberOfParameters, numberOfTerms ):
        """Allocation."""
        self.cObject = HarmonicBondContainer_Allocate ( numberOfTerms, numberOfParameters )
        self.isOwner = True

    def ActivateTerms ( self ):
        """Activate terms."""
        HarmonicBondContainer_ActivateTerms ( self.cObject )

    def DeactivateFixedAtomTerms ( self, Selection fixedatoms ):
        """Deactivate terms involving fixed atoms."""
        if fixedatoms is not None: HarmonicBondContainer_DeactivateFixedAtomTerms ( self.cObject, fixedatoms.cObject )

    def DeactivateQCAtomTerms ( self, Selection qcAtoms, Selection boundaryatoms ):
        """Deactivate terms involving QC atoms."""
        if qcAtoms is not None:
            if boundaryatoms is None: HarmonicBondContainer_DeactivateQCAtomTerms ( self.cObject, qcAtoms.cObject, NULL                  )
            else:                     HarmonicBondContainer_DeactivateQCAtomTerms ( self.cObject, qcAtoms.cObject, boundaryatoms.cObject )

    def Energy ( self, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Energy and gradients."""
        if gradients3 is None: energy = HarmonicBondContainer_Energy ( self.cObject, coordinates3.cObject, NULL               )
        else:                  energy = HarmonicBondContainer_Energy ( self.cObject, coordinates3.cObject, gradients3.cObject )
        return ( self.label, energy )

    def Get12Indices ( self ):
        """Return the 1-2 indices of the bonds."""
        cdef Integer i
        indices = []
        for i from 0 <= i < self.cObject.nterms:
            indices.append ( self.cObject.terms[i].atom1 )
            indices.append ( self.cObject.terms[i].atom2 )
        return indices

    def GetUniqueTermIndices ( self ):
        """Return the term indices."""
        cdef Integer i
        indices = []
        for i from 0 <= i < self.cObject.nterms:
            indices.append ( self.cObject.terms[i].atom1 )
            indices.append ( self.cObject.terms[i].atom2 )
        return ( self.cObject.nterms, indices )

#    def GetUsedParameterKeys ( self, parameterKeys ):
#        """Get used parameter keys."""
#        cdef Integer i
#        usedKeys = set ( )
#        for i from 0 <= i < self.cObject.nterms:
#            usedKeys.add ( parameterKeys[self.cObject.terms[i].type] )
#        return usedKeys

    def IdentifyBoundaryAtoms ( self, Selection qcAtoms, results ):
        """Return a dictionary of MM boundary atoms and their MM and QC partners."""
        cdef atom1, atom2, i, j, n, t
        # . Initialization.
        if self.is12Interaction and ( len ( qcAtoms ) > 0 ):
            # . Make the flags representation of the selection.
            n = HarmonicBondContainer_UpperBound ( self.cObject )
            Selection_MakeFlags ( qcAtoms.cObject, n )
            # . Loop to find the boundary atoms.
            boundaryatoms = set ( )
            for t from 0 <= t < self.cObject.nterms:
                i = self.cObject.terms[t].atom1
                j = self.cObject.terms[t].atom2
                if   ( ( qcAtoms.cObject.flags[i] == CTrue ) and ( qcAtoms.cObject.flags[j] == CFalse ) ): boundaryatoms.add ( j )
                elif ( ( qcAtoms.cObject.flags[j] == CTrue ) and ( qcAtoms.cObject.flags[i] == CFalse ) ): boundaryatoms.add ( i )
            # . Boundary atoms exist.
            if len ( boundaryatoms ) > 0:
                # . Find the partners of the boundary atoms - slow version.
                for t from 0 <= t < self.cObject.nterms:
                    atom1 = self.cObject.terms[t].atom1
                    atom2 = self.cObject.terms[t].atom2
                    # . Each term is independent.
                    for ( i, j ) in ( ( atom1, atom2 ), ( atom2, atom1 ) ):
                        if ( i in boundaryatoms ):
                            data = results.get ( i, [ set ( ), set ( ), set ( ) ] )
                            if   ( j in boundaryatoms ):                data[0].add ( j )
                            elif ( qcAtoms.cObject.flags[j] == CFalse ): data[1].add ( j )
                            else:                                       data[2].add ( j )
                            results[i] = data

    def Merge ( self, items, information = {} ):
        """Merging - slow version."""
        cdef Integer i, iOld
        cdef HarmonicBondContainer new
        # . Initialization.
        new = None
        # . Get the selection increments.
        increments = information.get ( "atomIncrements", None )
        # . Do the merge.
        if increments is not None:
            # . Initialization.
            doKeys          = True
            is12Interaction = True
            np              = 0
            parameterKeys   = []
            parameters      = []
            terms           = []
            # . Gather parameters and terms.
            for ( increment, item ) in zip ( increments, [ self ] + items ):
                state = item.__getstate__ ( )
                if np == 0:
                    parameterFields = state["parameterFields"]
                    termFields      = state["termFields"     ]
                is12Interaction = is12Interaction and state.get ( "is12Interaction", True )
                parameterKeys0  = state.get ( "parameterKeys", None )
                parameters0     = state["parameters"]
                terms0          = state["terms"     ]
                np0             = len ( parameters0 )
                parameters.extend ( parameters0 )
                for ( a1, a2, p, q ) in terms0:
                    terms.append ( ( a1 + increment ,
                                     a2 + increment ,
                                     p  + np , q ) )
                if parameterKeys0 is None: doKeys = False
                else: parameterKeys.extend ( parameterKeys0 )
                np = np + np0
            # . Construct the new keys and parameters.
            if doKeys:
                ( doKeys, oldToNew, parameterKeys, parameters ) = self.MergeKeys ( parameterKeys, parameters )
            else:
                parameterKeys = None
            # . Construct the object.
            state = { "is12Interaction" : is12Interaction ,
                      "label"           : self.label      ,
                      "parameterFields" : parameterFields ,
                      "parameterKeys"   : parameterKeys   ,
                      "parameters"      : parameters      ,
                      "termFields"      : termFields      ,
                      "terms"           : terms           }
            new = self.__class__.Raw ( )
            new.__setstate__ ( state )
            # . Reindex.
            if doKeys: new.ResetTermTypes ( oldToNew )
        return new

    def NumberOfParameters ( self ):
        """Return the number of parameters."""
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nparameters

    def NumberOfTerms ( self ):
        """Return the number of terms."""
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nterms

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef HarmonicBondContainer new
        new                 = HarmonicBondContainer.Raw ( )
        new.cObject         = HarmonicBondContainer_Prune ( self.cObject, selection.cObject )
        new.is12Interaction = self.is12Interaction
        new.isOwner         = True
        new.label           = self.label
        new.parameterKeys   = self.parameterKeys
        if new.cObject == NULL: return None
        else:                   return new

# . Need to prune parameters in cObject if also prune keys!
#    def PruneKeys ( self, new ):
#        """Prune parameter keys."""
#        if self.parameterKeys is not None:
#            # . Get new keys.
#            print "Old & New Keys Before Prune>", self.parameterKeys, new.parameterKeys, self.label
#            newKeys = new.GetUsedParameterKeys ( self.parameterKeys )
#            if len ( newKeys ) < len ( self.parameterKeys ):
#                newKeys = list ( newKeys )
#                newKeys.sort ( )
#                print "New Keys> ", newKeys
#                # . Get type index.
#                oldToNew = {}
#                for ( i, newKey ) in enumerate ( newKeys ):
#                    oldToNew[self.parameterKeys.index ( newKey )] = i
#                print "Old->New Mapping>", oldToNew
#                # . Assign type data.
#                new.parameterKeys = newKeys
#                new.ResetTermTypes ( oldToNew )

    def ResetTermTypes ( self, mapping ):
        """Reset term types."""
        cdef Integer i
        for i from 0 <= i < self.cObject.nterms:
            iOld = self.cObject.terms[i].type
            self.cObject.terms[i].type = mapping[iOld]

    def Sort ( self ):
        """Sorting."""
        HarmonicBondContainer_Sort ( self.cObject )

    def SummaryEntry ( self, summary ):
        """Summary entry."""
        cdef Integer n
        if summary is not None:
            n = HarmonicBondContainer_NumberOfInactiveTerms ( self.cObject )
            summary.Entry ( self.label + " Terms",      "{:d}".format ( self.cObject.nterms      ) )
            summary.Entry ( self.label + " Parameters", "{:d}".format ( self.cObject.nparameters ) )
            if n > 0: summary.Entry ( self.label + " Inactive", "{:d}".format ( n ) )

    def UpperBound ( self ):
        return HarmonicBondContainer_UpperBound ( self.cObject )
