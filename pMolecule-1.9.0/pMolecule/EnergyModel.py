#-------------------------------------------------------------------------------
# . File      : EnergyModel.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the energy model for a system.

The model contains all the data and definitions necessary for an energy
calculation apart from those that are stored in the system's configuration.

The order of definition of models is MM before QC before NB. Softconstraints
are independent and can be defined at any time.
"""

from pCore                   import Clone, logFile, LogFileActive, pObjectProperty, RawObjectConstructor, SelfPairList, SelfPairList_FromIntegerPairs
from LJParameterContainer    import LJParameterContainer
from MMAtomContainer         import MMAtomContainer
from MMModel                 import MMModel
from MMModelError            import MMModelError
from MMTermContainer         import MMTermContainer
from NBModel                 import NBModel
from QCAtomContainer         import QCAtomContainer
from QCModel                 import QCModel
from QCParameters            import QCParameters
from SoftConstraintContainer import SoftConstraintContainer

import math

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The tolerance for ensuring that the total MM active charge is integral.
_DefaultChargeTolerance = 1.0e-4

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class EnergyModel ( object ):
    """Define an energy model."""

    # . Different attribute sets.
    mmAttributes   = [ "mmModel"        ,
                       "mmAtoms"        ,
                       "mmTerms"        ,
                       "exclusions"     ,
                       "interactions14" ,
                       "ljParameters"   ,
                       "ljParameters14" ]
    modelAttributes = [ "mmModel"       ,
                        "qcModel"       ,
                        "nbModel"       ]
    qcAttributes    = [ "qcModel"       ,
                        "qcAtoms"       ,
                        "qcParameters"  ]

    # . All attributes.
    defaultAttributes = [ "label" ] + mmAttributes + qcAttributes + [ "nbModel", "softConstraints" ]

    def __getstate__ ( self ):
        state = {}
        for key in self.__class__.defaultAttributes:
            value = self.__dict__.get ( key, None )
            if value is not None: state[key] = value
        return state

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        for key in self.__class__.defaultAttributes:
            value = state.get ( key, None )
            if value is not None: self.__dict__[key] = value

    def ActivateMMTerms ( self ):
        """Activate all MM atoms and terms."""
        if self.mmAtoms is not None: self.mmAtoms.ActivateAtoms ( )
        if ( self.mmTerms is not None ) and ( len ( self.mmTerms ) > 0 ):
            for mmterm in self.mmTerms:
                if hasattr ( mmterm, "ActivateTerms" ): mmterm.ActivateTerms ( )

    def CheckActiveMMAtomTotalCharge ( self, tolerance = _DefaultChargeTolerance ):
        """Throw an error if the total active MM charge is not integral."""
        if hasattr ( self, "mmAtoms" ) and ( self.mmAtoms is not None ):
            totalcharge    = sum ( self.mmAtoms.AtomicCharges ( ) )
            nearestinteger = round ( totalcharge )
            if math.fabs ( totalcharge - nearestinteger ) > tolerance:
                raise MMModelError ( "Total active MM charge is not integral (or zero): {:.3f}.".format ( totalcharge ), totalcharge )

    def ClearFixedAtoms ( self, configuration ):
        """Clear data pertaining to fixed atoms."""
        self.ActivateMMTerms         ( )
        self.DeactivateQCAtomMMTerms ( )
        if self.nbModel is not None: self.nbModel.Clear ( configuration )

    def ClearMMModel ( self, configuration ):
        """Clear an MM model."""
        self.ClearNBModel ( configuration )
        self.ClearQCModel ( configuration )
        if self.mmModel is not None:
            for attribute in self.__class__.mmAttributes: delattr ( self, attribute )
        self.MakeLabel ( )

    def ClearNBModel ( self, configuration ):
        """Clear an NB model."""
        if self.nbModel is not None:
            self.nbModel.Clear ( configuration )
            del self.nbModel
        self.MakeLabel ( )

    def ClearQCModel ( self, configuration, fixedAtoms = None ):
        """Clear a QC model."""
        self.ClearNBModel ( configuration )
        if self.qcModel is not None:
            self.qcModel.Clear ( configuration )
            for attribute in self.__class__.qcAttributes: delattr ( self, attribute )
            if self.mmModel is not None:
                self.ActivateMMTerms ( )
                self.DeactivateFixedAtomMMTerms ( fixedAtoms )
        self.MakeLabel ( )

    def DeactivateFixedAtomMMTerms ( self, fixedAtoms ):
        """Deactivate the MM terms between fixed atoms."""
        if ( self.mmTerms is not None ) and ( len ( self.mmTerms ) > 0 ) and ( fixedAtoms is not None ):
            for mmterm in self.mmTerms:
                if hasattr ( mmterm, "DeactivateFixedAtomTerms" ): mmterm.DeactivateFixedAtomTerms ( fixedAtoms )

    def DeactivateQCAtomMMTerms ( self ):
        """Deactivate MM atoms that are purely QC and the necessary MM terms between QC atoms."""
        if ( self.mmAtoms is not None ) and ( self.qcAtoms is not None ):
            # . Get selections for the qc (all) and boundary atoms.
            qcSelection = self.qcAtoms.QCAtomSelection ( )
            if self.qcAtoms.NumberOfBoundaryAtoms ( ) > 0: baSelection = self.qcAtoms.BoundaryAtomSelection ( )
            else:                                          baSelection = None
            # . Deactivate terms.
            self.mmAtoms.DeactivateQCAtoms ( qcSelection, baSelection )
            if ( self.mmTerms is not None ) and ( len ( self.mmTerms ) > 0 ):
                for mmterm in self.mmTerms:
                    if hasattr ( mmterm, "DeactivateQCAtomTerms" ): mmterm.DeactivateQCAtomTerms ( qcSelection, baSelection )

    def Get12Exclusions ( self ):
        """Get a 1-2 exclusion list from the MM terms."""
        exclusions = None
        if ( self.mmTerms is not None ) and ( len ( self.mmTerms ) > 0 ):
            indices = []
            for mmterm in self.mmTerms:
                if hasattr ( mmterm, "Get12Indices" ): indices.extend ( mmterm.Get12Indices ( ) )
            exclusions = SelfPairList_FromIntegerPairs ( indices )
        return exclusions

    def IdentifyBoundaryAtoms ( self, selection ):
        """Identify atoms in |selection| that are covalently bound to atoms outside of |selection|."""
        results = {}
        if ( self.mmTerms is not None ) and ( len ( self.mmTerms ) > 0 ) and ( selection is not None ):
            for mmterm in self.mmTerms:
                if hasattr ( mmterm, "IdentifyBoundaryAtoms" ):
                    mmterm.IdentifyBoundaryAtoms ( selection, results )
        return results

    def MakeLabel ( self ):
        """Make a label for the model."""
        labels = []
        for attribute in self.__class__.modelAttributes:
            item = getattr ( self, attribute )
            if ( item is not None ) and ( item.label is not None ): labels.append ( item.label )
        if len ( labels ) > 0: self.label = "/".join ( labels )
        else:                  del self.label

    def Merge ( self, others, information = {} ):
        """Merging."""
        # . Initialization.
        merged        = None
        toMergeItems  = [ self ] + list ( others )
        numberToMerge = len ( toMergeItems )
        # . Need increments.
        increments = information.get ( "atomIncrements", None )
        if ( increments is not None ) and ( len ( increments ) == numberToMerge ):
            # . MM and NB models - all must be of the same class.
            models = {}
            for attribute in ( "mmModel", "nbModel" ):
                items = [ getattr ( item, attribute ) for item in ( toMergeItems ) ]
                model = items.pop ( 0 )
                for item in items:
                    if not isinstance ( model, item.__class__ ):
                        model = None
                        break
                models[attribute] = model
            # . QC models.
            hasQCModels = False
            for toMerge in toMergeItems:
                item = getattr ( item, "qcModel", None )
                if item is not None:
                    hasQCModels = True
                    break
            # . A uniform MM model was found.
            if models["mmModel"] is not None:
                # . Basic object.
                merged = self.__class__ ( )
                merged.label   = Clone ( self.label )
                merged.mmModel = Clone ( models["mmModel"] )
                if ( models["nbModel"] is not None ) and ( not hasQCModels ): merged.nbModel = Clone ( models["nbModel"] )
                # . Attributes which must be defined for all merging objects.
                # . "ljParameters" must be before "mmAtoms".
                for attribute in ( "ljParameters", "ljParameters14", "mmAtoms" ):
                    items = []
                    for toMerge in toMergeItems:
                        item = getattr ( toMerge, attribute, None )
                        if item is not None: items.append ( item )
                    if len ( items ) == numberToMerge: setattr ( merged, attribute, items[0].Merge ( items[1:], information = information ) )
                # . Attributes which can be optional.
                for attribute in ( "exclusions", "interactions14", "mmTerms", "softConstraints" ):
                    items         = []
                    newIncrements = []
                    for ( toMerge, increment ) in zip ( toMergeItems, increments ):
                        item = getattr ( toMerge, attribute, None )
                        if item is not None:
                            items.append         ( item      )
                            newIncrements.append ( increment )
                    if len ( items ) > 0:
                        mergedItem = items[0].Merge ( items[1:], information = { "atomIncrements" : newIncrements, "indexIncrements" : newIncrements } )
                        if mergedItem is not None: setattr ( merged, attribute, mergedItem )
                # . Finish up.
                merged.ActivateMMTerms ( )
        return merged

    def MMObjects ( self ):
        """Return a list of MM objects defined by the model."""
        items = []
        for attribute in self.__class__.mmAttributes:
            if attribute == "mmTerms":
                if self.mmTerms is not None:
                    for item in self.mmTerms: items.append ( item )
            else:
                item = getattr ( self, attribute )
                if item is not None: items.append ( item )
        return items

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        pruned = None
        if self.mmModel is not None:
            # . Basic data.
            pruned = self.__class__ ( )
            if self.label   is not None: pruned.label   = Clone ( self.label   )
            if self.mmModel is not None: pruned.mmModel = Clone ( self.mmModel )
            if ( self.nbModel is not None ) and ( self.qcModel is None ): pruned.nbModel = Clone ( self.nbModel )
            # . Other attributes.
            for attribute in ( "exclusions", "interactions14", "ljParameters", "ljParameters14", "mmAtoms", "mmTerms", "softConstraints" ):
                item = getattr ( self, attribute )
                if hasattr ( item, "Prune" ):
                    prunedItem = item.Prune ( selection, information = information )
                    if prunedItem is not None: setattr ( pruned, attribute, prunedItem )
            # . Finish up.
            pruned.ActivateMMTerms ( )
        return pruned

    def QCObjects ( self ):
        """Return a list of QC objects defined by the model."""
        items = []
        for attribute in self.__class__.qcAttributes:
            item = getattr ( self, attribute )
            if item is not None: items.append ( item )
        return items

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ) and ( ( self.mmModel is not None ) or ( self.qcModel is not None ) ):
            if self.label is None: title = "Energy Model Summary"
            else:                  title = "Summary for Energy Model \"" + self.label + "\""
            log.Heading ( title, includeBlankLine = True )
            if self.mmModel is not None:
                summary = log.GetSummary ( pageWidth = 100 )
                if self.mmModel.label is None: summary.Start ( "MM Model Summary" )
                else:                          summary.Start ( "Summary for MM Model \"" + self.mmModel.label + "\"" )
                for item in self.MMObjects ( ):
                    if hasattr ( item, "SummaryEntry" ): item.SummaryEntry ( summary )
                summary.Stop ( )
            if self.qcModel is not None:
                summary = log.GetSummary ( )
                if self.qcModel.label is None: summary.Start ( "QC Model Summary" )
                else:                          summary.Start ( "Summary for QC Model \"" + self.qcModel.label + "\"" )
                for item in self.QCObjects ( ):
                    if hasattr ( item, "SummaryEntry" ): item.SummaryEntry ( summary )
                summary.Stop ( )
            if self.nbModel         is not None: self.nbModel.Summary         ( log )
            if self.softConstraints is not None: self.softConstraints.Summary ( log )

    # . Properties.
    label           = pObjectProperty ( "label", basestring, default = None, doc = "Label." )

    mmModel         = pObjectProperty ( "mmModel",         MMModel,                 doc = "MM Model."                   )
    mmAtoms         = pObjectProperty ( "mmAtoms",         MMAtomContainer,         doc = "MM Atom Container."          )
    mmTerms         = pObjectProperty ( "mmTerms",         MMTermContainer,         doc = "MM Term Container."          )
    exclusions      = pObjectProperty ( "exclusions",      SelfPairList,            doc = "Exclusions."                 )
    interactions14  = pObjectProperty ( "interactions14",  SelfPairList,            doc = "1-4 Interactions."           )
    ljParameters    = pObjectProperty ( "ljParameters",    LJParameterContainer,    doc = "LJ Parameter Container."     )
    ljParameters14  = pObjectProperty ( "ljParameters14",  LJParameterContainer,    doc = "1-4 LJ Parameter Container." )

    qcModel         = pObjectProperty ( "qcModel",         QCModel,                 doc = "QC Model."                   )
    qcAtoms         = pObjectProperty ( "qcAtoms",         QCAtomContainer,         doc = "QC Atom Container."          )
    qcParameters    = pObjectProperty ( "qcParameters",    QCParameters,            doc = "QC Parameters."              )

    nbModel         = pObjectProperty ( "nbModel",         NBModel,                 doc = "NB Model."                   )

    softConstraints = pObjectProperty ( "softConstraints", SoftConstraintContainer, doc = "Soft Constraint Container."  )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
