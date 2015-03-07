#-------------------------------------------------------------------------------
# . File      : pMolecule.CMAPDihedralContainer.pyx
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A container for CMAP dihedral energy terms."""

from pCore import Clone, logFile, LogFileActive

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CMAPDihedralContainer ( MMTerm ):

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef CMAPDihedralContainer new
        new               = self.__class__.Raw ( )
        new.label         = self.label
        new.cObject       = CMAPDihedralContainer_Clone ( self.cObject )
        new.isOwner       = True
        new.parameterKeys = Clone ( self.parameterKeys )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            CMAPDihedralContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.CMAPDihedralContainer"

    def __getstate__ ( self ):
        """Return the state."""
        cdef Integer i, ix, iy, lengthx, lengthy
        # . Parameters.
        parameterFields = [ "x", "y", "map" ]
        parameters      = []
        for i from 0 <= i < self.cObject.nparameters:
            x = []
            y = []
            f = []
            lengthx = self.cObject.parameters[i].lengthx
            lengthy = self.cObject.parameters[i].lengthy
            for ix from 0 <= ix < lengthx:
                x.append ( Real1DArray_GetItem ( self.cObject.parameters[i].x, ix, NULL ) )
            for iy from 0 <= iy < lengthy:
                y.append ( Real1DArray_GetItem ( self.cObject.parameters[i].y, iy, NULL ) )
            for ix from 0 <= ix < lengthx:
                for iy from 0 <= iy < lengthy:
                   f.append ( Real2DArray_GetItem ( self.cObject.parameters[i].f, ix, iy, NULL ) )
            parameters.append ( [ x, y, f ] )
        # . Terms.
        termFields = [ "atom1", "atom2", "atom3", "atom4", "atom5", "atom6", "atom7", "atom8", "parameter", "isActive" ]
        terms      = []
        for i from 0 <= i < self.cObject.nterms:
            terms.append ( [ self.cObject.terms[i].atom1 ,
                             self.cObject.terms[i].atom2 ,
                             self.cObject.terms[i].atom3 ,
                             self.cObject.terms[i].atom4 ,
                             self.cObject.terms[i].atom5 ,
                             self.cObject.terms[i].atom6 ,
                             self.cObject.terms[i].atom7 ,
                             self.cObject.terms[i].atom8 ,
                             self.cObject.terms[i].type  ,
                             self.cObject.terms[i].isActive == CTrue ] )
        # . State.
        state = { "label"           : self.label      ,
                  "parameterFields" : parameterFields ,
                  "parameters"      : parameters      ,
                  "termFields"      : termFields      ,
                  "terms"           : terms           }
        if self.parameterKeys is not None: state["parameterKeys"] = self.parameterKeys
        return state

    def __init__ ( self, numberOfParameters, numberOfTerms, label = "CMAP Dihedral" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfParameters, numberOfTerms )
        self.label = label

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef Integer       i, ix, iy, lengthx, lengthy
        cdef CReal1DArray *x, *y
        cdef CReal2DArray *f
        # . Allocate the object.
        parameters = state["parameters"]
        terms      = state["terms"     ]
        self._Allocate ( len ( parameters ), len ( terms ) )
        if "label" in state: self.label = state["label"]
        self.parameterKeys = state.get ( "parameterKeys", None )
        # . Fill the object.
        # . Parameters.
        for ( i, ( xP, yP, fP ) ) in enumerate ( parameters ):
            lengthx = len ( xP )
            lengthy = len ( yP )
            x = Real1DArray_Allocate ( lengthx, NULL )
            y = Real1DArray_Allocate ( lengthy, NULL )
            f = Real2DArray_Allocate ( lengthx, lengthy, NULL )
            for ix from 0 <= ix < lengthx:
                Real1DArray_SetItem ( x, ix, xP[ix], NULL )
            for iy from 0 <= iy < lengthy:
                Real1DArray_SetItem ( y, iy, yP[iy], NULL )
            n = 0
            for ix from 0 <= ix < lengthx:
                for iy from 0 <= iy < lengthy:
                    Real2DArray_SetItem ( f, ix, iy, fP[n], NULL )
                    n += 1
            self.cObject.parameters[i] = BicubicSpline_MakeFromReal2DArray ( &x, &y, &f, BicubicSplineType_Periodic, NULL )
        # . Terms.
        for ( i, ( a1, a2, a3, a4, a5, a6, a7, a8, p, q ) ) in enumerate ( terms ):
            self.cObject.terms[i].atom1 = a1
            self.cObject.terms[i].atom2 = a2
            self.cObject.terms[i].atom3 = a3
            self.cObject.terms[i].atom4 = a4
            self.cObject.terms[i].atom5 = a5
            self.cObject.terms[i].atom6 = a6
            self.cObject.terms[i].atom7 = a7
            self.cObject.terms[i].atom8 = a8
            self.cObject.terms[i].type  = p
            if q: self.cObject.terms[i].isActive = CTrue
            else: self.cObject.terms[i].isActive = CFalse

    def _Allocate ( self, numberOfParameters, numberOfTerms ):
        """Allocation."""
        self.cObject = CMAPDihedralContainer_Allocate ( numberOfTerms, numberOfParameters )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject       = NULL
        self.isOwner       = False
        self.label         = "CMAP Dihedral"
        self.parameterKeys = None

    def ActivateTerms ( self ):
        """Activate terms."""
        CMAPDihedralContainer_ActivateTerms ( self.cObject )

    def DeactivateFixedAtomTerms ( self, Selection fixedatoms ):
        """Deactivate terms involving fixed atoms."""
        if fixedatoms is not None: CMAPDihedralContainer_DeactivateFixedAtomTerms ( self.cObject, fixedatoms.cObject )

    def DeactivateQCAtomTerms ( self, Selection qcAtoms, Selection boundaryatoms ):
        """Deactivate terms involving QC atoms."""
        if qcAtoms is not None: CMAPDihedralContainer_DeactivateQCAtomTerms ( self.cObject, qcAtoms.cObject, NULL )

    def Energy ( self, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Energy and gradients."""
        if gradients3 is None: energy = CMAPDihedralContainer_Energy ( self.cObject, coordinates3.cObject, NULL               )
        else:                  energy = CMAPDihedralContainer_Energy ( self.cObject, coordinates3.cObject, gradients3.cObject )
        return ( self.label, energy )

    def GetUniqueTermIndices ( self ):
        """Return the term indices."""
        cdef Integer i
        indices = []
        for i from 0 <= i < self.cObject.nterms:
            indices.append ( self.cObject.terms[i].atom1 )
            indices.append ( self.cObject.terms[i].atom2 )
            indices.append ( self.cObject.terms[i].atom3 )
            indices.append ( self.cObject.terms[i].atom4 )
            indices.append ( self.cObject.terms[i].atom5 )
            indices.append ( self.cObject.terms[i].atom6 )
            indices.append ( self.cObject.terms[i].atom7 )
            indices.append ( self.cObject.terms[i].atom8 )
        return ( self.cObject.nterms, indices )

    def Merge ( self, items, information = {} ):
        """Merging - slow version."""
        # . Initialization.
        new = None
        # . Get the selection increments.
        increments = information.get ( "atomIncrements", None )
        # . Do the merge.
        if increments is not None:
            # . Gather terms.
            np         = 0
            parameters = []
            terms      = []
            for ( increment, item ) in zip ( increments, [ self ] + items ):
                state       = item.__getstate__ ( )
                if np == 0:
                    parameterFields = state["parameterFields"]
                    termFields      = state["termFields"     ]
                parameters0 = state["parameters"]
                terms0      = state["terms"     ]
                np0         = len ( parameters0 )
                parameters.extend ( parameters0 )
                for ( a1, a2, a3, a4, a5, a6, a7, a8, p, q ) in terms0:
                    terms.append ( ( a1 + increment ,
                                     a2 + increment ,
                                     a3 + increment ,
                                     a4 + increment ,
                                     a5 + increment ,
                                     a6 + increment ,
                                     a7 + increment ,
                                     a8 + increment ,
                                     p  + np , q ) )
                np = np + np0
            # . Construct the object.
            state = { "label"           : self.label      ,
                      "parameterFields" : parameterFields ,
                      "parameters"      : parameters      ,
                      "termFields"      : termFields      ,
                      "terms"           : terms           }
            new = self.__class__.Raw ( )
            new.__setstate__ ( state )
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
        cdef CMAPDihedralContainer new
        new         = CMAPDihedralContainer.Raw ( )
        new.cObject = CMAPDihedralContainer_Prune ( self.cObject, selection.cObject )
        new.isOwner = True
        new.label   = self.label
        if new.cObject == NULL: return None
        else:                   return new

    def Sort ( self ):
        """Sorting."""
        CMAPDihedralContainer_Sort ( self.cObject )

    def SummaryEntry ( self, summary ):
        """Summary entry."""
        cdef Integer n
        if summary is not None:
            n = CMAPDihedralContainer_NumberOfInactiveTerms ( self.cObject )
            summary.Entry ( self.label + " Terms",      "{:d}".format ( self.cObject.nterms      ) )
            summary.Entry ( self.label + " Parameters", "{:d}".format ( self.cObject.nparameters ) )
            if n > 0: summary.Entry ( self.label + " Inactive", "{:d}".format ( n ) )

    def UpperBound ( self ):
        return CMAPDihedralContainer_UpperBound ( self.cObject )
