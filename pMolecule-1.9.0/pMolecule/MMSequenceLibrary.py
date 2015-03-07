#-------------------------------------------------------------------------------
# . File      : MMSequenceLibrary.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for handling an MM sequence library."""

import os, os.path

from pCore      import logFile, LogFileActive, YAMLMappingFile_ToObject, YAMLPickleFileExtension
from MMSequence import MMSequenceComponent, MMSequenceLink, MMSequenceVariant
from Sequence   import Sequence

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Paths.
_ComponentPath = "components"
_LinkPath      = "links"
_VariantPath   = "variants"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMSequenceLibrary ( object ):
    """Library of MM sequence components, links and variants."""

    defaultAttributes = { "componentPaths" :  None ,
                          "components"     :  None ,
                          "isSetup"        : False ,
                          "libraryPath"    :  None ,
                          "linkPaths"      :  None ,
                          "links"          :  None ,
                          "variantPaths"   :  None ,
                          "variants"       :  None }

    def __init__ ( self, path ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )
        self.DefineLibraryPath ( path )
        self.SetupPaths ( )
        self.ResetCache ( )

    def CheckMissingItems ( self, missingItems, log  ):
        """Check missing items."""
        # . Check for missing items.
        if len ( missingItems ) > 0:
            if LogFileActive ( log ):
                missingItems = list ( missingItems )
                missingItems.sort ( )
                table = log.GetTable ( columns = [ 12, 24 ] )
                table.Start ( )
                table.Title ( "Missing MM Sequence Components" )
                table.Heading ( "Type"  )
                table.Heading ( "Label" )
                for ( label, key ) in missingItems:
                    table.Entry ( label, alignment = "l" )
                    table.Entry ( key,   alignment = "l" )
                table.Stop ( )
            raise ValueError ( "There are {:d} missing items.".format ( len ( missingItems ) ) )

    def DefineLibraryPath ( self, path ):
        """Define the library path."""
        self.libraryPath = path

    def GetComponent ( self, componentLabel, missingItems = None ):
        """Get a component."""
        key = MMSequenceComponent.MakeKey ( componentLabel )
        if key not in self.components:
            for componentPath in self.componentPaths:
                fileName = os.path.join ( componentPath, key + YAMLPickleFileExtension )
                if os.path.exists ( fileName ):
                    self.components[key] = YAMLMappingFile_ToObject ( fileName, MMSequenceComponent )
                    break
        component = self.components.get ( key, None )
        if ( component is None ) and ( missingItems is not None ): missingItems.add ( ( "Component", key ) )
        return component

    def GetLink ( self, linkLabel, leftComponentLabel, rightComponentLabel, missingItems = None ):
        """Get a link."""
        keys = MMSequenceLink.MakeKeys ( linkLabel, leftComponentLabel, rightComponentLabel )
        key  = keys[0]
        if key not in self.links:
            found = False
            for tag in keys:
                if found: break
                for linkPath in self.linkPaths:
                    fileName = os.path.join ( linkPath, tag + YAMLPickleFileExtension )
                    if os.path.exists ( fileName ):
                        self.links[key] = YAMLMappingFile_ToObject ( fileName, MMSequenceLink )
                        found = True
                        break
        link = self.links.get ( key, None )
        if ( link is None ) and ( missingItems is not None ): missingItems.add ( ( "Link", key ) )
        return link

    def GetVariant ( self, variantLabel, componentLabel, missingItems = None ):
        """Get a variant."""
        keys = MMSequenceVariant.MakeKeys ( variantLabel, componentLabel )
        key  = keys[0]
        if key not in self.variants:
            found = False
            for tag in keys:
                if found: break
                for variantPath in self.variantPaths:
                    fileName = os.path.join ( variantPath, tag + YAMLPickleFileExtension )
                    if os.path.exists ( fileName ):
                        self.variants[key] = YAMLMappingFile_ToObject ( fileName, MMSequenceVariant )
                        found = True
                        break
        variant = self.variants.get ( key, None )
        if ( variant is None ) and ( missingItems is not None ): missingItems.add ( ( "Variant", key ) )
        return variant

    def ResetCache ( self ):
        """Clear stored components, links and variants."""
        self.components = {}
        self.links      = {}
        self.variants   = {}

    def SetupPaths ( self ):
        """Set up all paths."""
        # . Define paths.
        for ( attribute, tailPath ) in ( ( "componentPaths", _ComponentPath ) ,
                                         ( "linkPaths"     , _LinkPath      ) ,
                                         ( "variantPaths"  , _VariantPath   ) ):
            localPaths = []
            for rootPath in ( self.libraryPath, ):
                path = os.path.join ( rootPath, tailPath )
                localPaths.append ( path )
            setattr ( self, attribute, localPaths )

    @classmethod
    def TestForLibrary ( selfClass, path ):
        """Test for a library and return if present."""
        self   = selfClass ( path )
        exists = False
        for path in self.componentPaths:
            if os.path.exists ( path ):
                exists = True
                break
        if not exists: self = None
        return self

    def TypeSequence ( self, sequence, atomTypes, atomCharges, untypedAtoms, log = logFile ):
        """Type the atoms of a sequence."""
        # . Initialization.
        componentsToType = {}
        missingItems     = set ( )
        variantsToType   = {}
        self.ResetCache ( )
        # . Get all MM sequence objects.
        # . Components.
        for entity in sequence.children:
            for component in entity.children:
                mmComponent = self.GetComponent ( component.genericLabel )
                if mmComponent is not None: componentsToType[component] = mmComponent
        # . Variants.
        for variant in sequence.variants:
            component = variant.component
            if component in componentsToType:
                mmVariant = self.GetVariant ( variant.label, component.genericLabel, missingItems = missingItems )
                if mmVariant is not None:
                    mmVariants = variantsToType.get ( component, None )
                    if mmVariants is None:
                        mmVariants = []
                        variantsToType[component] = mmVariants
                    mmVariants.append ( mmVariant )
        # . Links.
        for link in sequence.links:
            leftComponent  = link.leftComponent
            rightComponent = link.rightComponent
            if ( leftComponent in componentsToType ) and ( rightComponent in componentsToType ):
                mmLink = self.GetLink ( link.label, leftComponent.genericLabel, rightComponent.genericLabel, missingItems = missingItems )
                if mmLink is not None:
                    for ( component, variant ) in ( ( leftComponent, mmLink.leftVariant ), ( rightComponent, mmLink.rightVariant ) ):
                        mmVariants = variantsToType.get ( component, None )
                        if mmVariants is None:
                            mmVariants = []
                            variantsToType[component] = mmVariants
                        mmVariants.append ( variant )
        # . Type the atoms.
        for ( component, mmComponent ) in componentsToType.iteritems ( ):
            ( componentTypes, componentCharges ) = mmComponent.TypeSequenceComponentAtoms ( component )
            mmVariants = variantsToType.pop ( component, [] )
            for mmVariant in mmVariants:
                mmVariant.TypeSequenceComponentAtoms ( component, atomCharges = componentCharges, atomTypes = componentTypes )
            if len ( componentTypes ) == len ( component.children ):
                indices = componentTypes.keys ( )
                for index in indices:
                    atomCharges[index] = componentCharges.get ( index, 0.0 )
                    atomTypes  [index] = componentTypes[index]
                untypedAtoms.Exclude ( indices )
        # . Check for untyped components due to links.
        for component in variantsToType:
            missingItems.append ( ( "Component", component.genericLabel ) )
        # . Check for missing items.
        self.CheckMissingItems ( missingItems, log )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
