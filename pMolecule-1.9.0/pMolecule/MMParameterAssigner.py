#-------------------------------------------------------------------------------
# . File      : MMParameterAssigner.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MM parameter assigner."""

import glob, os, os.path

from pCore                    import Clone, logFile, LogFileActive, YAMLMappingFile_ToObject, YAMLPickleFileExtension

from FourierDihedralContainer import FourierDihedralContainer
from HarmonicAngleContainer   import HarmonicAngleContainer
from HarmonicBondContainer    import HarmonicBondContainer
from LJParameterContainer     import LJParameterContainer
from MMAtomContainer          import MMAtomContainer
from MMModelError             import MMModelError
from MMTermContainer          import MMTermContainer

from CMAPDihedralParameterContainer       import CMAPDihedralParameterContainer
from FourierDihedralParameterContainer    import FourierDihedralParameterContainer
from FourierOutOfPlaneParameterContainer  import FourierOutOfPlaneParameterContainer
from HarmonicAngleParameterContainer      import HarmonicAngleParameterContainer
from HarmonicBondParameterContainer       import HarmonicBondParameterContainer
from HarmonicOutOfPlaneParameterContainer import HarmonicOutOfPlaneParameterContainer
from LennardJonesParameterContainer       import LennardJonesParameterContainer
from UreyBradleyParameterContainer        import UreyBradleyParameterContainer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
# . Mapping from files names to classes.
_PathClassMapping = { "cmapDihedralParameters"       : CMAPDihedralParameterContainer       ,
                      "fourierDihedralParameters"    : FourierDihedralParameterContainer    ,
                      "fourierOutOfPlaneParameters"  : FourierOutOfPlaneParameterContainer  ,
                      "harmonicAngleParameters"      : HarmonicAngleParameterContainer      ,
                      "harmonicBondParameters"       : HarmonicBondParameterContainer       ,
                      "harmonicOutOfPlaneParameters" : HarmonicOutOfPlaneParameterContainer ,
                      "lennardJones14Parameters"     : LennardJonesParameterContainer       ,
                      "lennardJonesParameters"       : LennardJonesParameterContainer       ,
                      "ureyBradleyParameters"        : UreyBradleyParameterContainer        }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMParameterAssigner ( object ):
    """Assign MM parameters."""

    defaultAttributes = { "dataPath"                 : None,
                          "lennardJonesParameters"   : None,
                          "lennardJones14Parameters" : None,
                          "lennardJonesScale14"      :  1.0,
                          "parameterContainers"      : None }

    def __init__ ( self, dataPath, **keywordArguments ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        for ( key, value ) in                 keywordArguments.iteritems ( ): setattr ( self, key, value )
        self.dataPath = dataPath
        self.ReadData ( )

    def AssignParameters ( self, connectivity, atomTypes, atomCharges, energyModel, log ):
        """Assign parameters to the model."""
        # . Initialization.
        missingParameters = set ( )
        mmTerms           = []
        uniqueAtomTypes   = self.FindUniqueAtomTypes ( atomTypes )
        # . Atom data.
        energyModel.mmAtoms = self.MakeMMAtomContainer ( atomTypes, atomCharges, uniqueAtomTypes )
        # . LJ parameters.
        for ( localAttribute, modelAttribute ) in ( ( "lennardJonesParameters"  , "ljParameters"   ),
                                                    ( "lennardJones14Parameters", "ljParameters14" ) ):
            container = getattr ( self, localAttribute, None )
            if container is not None:
                ( parameters, localMissingParameters ) = container.MakeParameterContainer ( uniqueAtomTypes )
                missingParameters.update ( localMissingParameters )
                if parameters is not None: setattr ( energyModel, modelAttribute, parameters )
        # . Exclusions.
        ( energyModel.exclusions, energyModel.interactions14 ) = connectivity.Make1234And14PairLists ( )
        energyModel.exclusions.label     = "Exclusions"
        energyModel.interactions14.label = "1-4 Interactions"
        # . MM terms.
        for container in self.parameterContainers:
            ( localMMTerms, localMissingParameters ) = container.MakeMMTermsFromConnectivity ( atomTypes, connectivity )
            missingParameters.update ( localMissingParameters )
            mmTerms.extend           ( localMMTerms           )
        if len ( mmTerms ) > 0: energyModel.mmTerms = MMTermContainer ( mmTerms )
        energyModel.ActivateMMTerms ( )
        self.CheckMissingParameters ( missingParameters, log )

    def CheckMissingParameters ( self, missingParameters, log ):
        """Print any missing parameters and raise an error."""
        numberMissing = len ( missingParameters )
        if numberMissing > 0:
            if LogFileActive ( log ):
                # . Sort.
                missingParameters = list ( missingParameters )
                missingParameters.sort ( )
                # . Find label data.
                labelCount  = 0
                labelLength = 0
                tagLength   = 0
                for ( tag, labels ) in missingParameters:
                    labelCount = max ( len ( labels ), labelCount )
                    tagLength  = max ( len ( tag    ), tagLength  )
                    for label in labels:
                        labelLength = max ( len ( label ), labelLength )
                # . Output.
                table = log.GetTable ( columns = [ tagLength + 2 ] + labelCount * [ max ( 10, labelLength + 2 ) ] )
                table.Start  ( )
                table.Title  ( "Missing Force Field Parameters" )
                for ( tag, labels ) in missingParameters:
                    table.Entry ( tag, alignment = "l" )
                    for label in labels: table.Entry ( label )
                    if len ( labels ) < labelCount: table.EndRow ( )
                table.Stop ( )
            raise MMModelError ( "There are {:d} missing force field parameters.".format ( numberMissing ) )

    def FindUniqueAtomTypes ( self, atomTypes ):
        """Find a sorted list of unique atom type labels."""
        temporary       = set  ( atomTypes )
        uniqueAtomTypes = list ( temporary )
        uniqueAtomTypes.sort ( )
        return uniqueAtomTypes

    def MakeLennardJones14Parameters ( self ):
        """Make an appropriate container of LJ Parameters for the 1-4 interactions."""
        lj    = self.lennardJonesParameters
        scale = self.lennardJonesScale14
        # . Remove 1-4 parameters if the scaling is zero.
        if scale == 0.0:
            lj14 = None
        # . Processing only if LJs are present too.
        elif lj is not None:
            # . Clone and scale.
            lj14           = Clone ( lj )
            lj14.termLabel = "1-4 Lennard-Jones"
            lj14.ScaleEnergies ( scale )
            # . Update by any predefined 1-4 values.
            oldLJ14 = self.lennardJones14Parameters
            if oldLJ14 is not None:
                lj14.UpdateParameters ( oldLJ14 )
        # . Reset the 14 parameters.
        self.lennardJones14Parameters = lj14

    def MakeMMAtomContainer ( self, atomTypes, atomCharges, uniqueAtomTypes ):
        """Make an MM atom container."""
        # . Type label index.
        labelIndex = {}
        for ( i, label ) in enumerate ( uniqueAtomTypes ):
            labelIndex[label] = i
        # . MM atom properties.
        atoms = []
        for ( charge, label ) in zip ( atomCharges, atomTypes ):
            index = labelIndex[label]
            atoms.append ( ( True, index, index, charge ) )
        # . Create the object.
        state = { "atoms"     : atoms           ,
                  "atomTypes" : uniqueAtomTypes }
        mm = MMAtomContainer.Raw ( )
        mm.__setstate__ ( state )
        return mm

    def ReadData ( self ):
        """Read the parameter data."""
        if self.dataPath is not None:
            self.parameterContainers = []
            paths                    = glob.glob ( os.path.join ( self.dataPath, "*Parameters" + YAMLPickleFileExtension ) )
            for path in paths:
                ( head, tail ) = os.path.split ( path )
                container      = YAMLMappingFile_ToObject ( path, _PathClassMapping[os.path.splitext ( tail )[0]] )
                if isinstance ( container, LennardJonesParameterContainer ):
                    if path.find ( "14Parameters" ) >= 0: self.lennardJones14Parameters = container
                    else:                                 self.lennardJonesParameters   = container
                else: self.parameterContainers.append ( container )
            self.MakeLennardJones14Parameters ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
