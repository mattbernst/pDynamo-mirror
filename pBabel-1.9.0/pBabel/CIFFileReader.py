#-------------------------------------------------------------------------------
# . File      : CIFFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Read data from a CIF file."""

import os.path

from pCore          import Coordinates3, logFile, LogFileActive, Matrix33, Selection, TextFileReader , \
                           Transformation3_FromSymmetryOperationString,  Transformation3Container    , \
                           Vector3
from ExportImport   import _Importer
from pMolecule      import CrystalClassCubic, CrystalClassHexagonal, CrystalClassMonoclinic, CrystalClassOrthorhombic, \
                           CrystalClassRhombohedral, CrystalClassTetragonal, CrystalClassTriclinic, System
from STARFileReader import STARFileReader, STARFileTable

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Crystal classes.
_CrystalClasses = { "cubic"        : CrystalClassCubic        ( ) ,
                    "hexagonal"    : CrystalClassHexagonal    ( ) ,
                    "monoclinic"   : CrystalClassMonoclinic   ( ) ,
                    "orthorhombic" : CrystalClassOrthorhombic ( ) ,
                    "rhombohedral" : CrystalClassRhombohedral ( ) ,
                    "tetragonal"   : CrystalClassTetragonal   ( ) ,
                    "triclinic"    : CrystalClassTriclinic    ( ) ,
                    "trigonal"     : CrystalClassHexagonal    ( ) }

# . Default disorder group.
_DefaultDisorderGroup = "1"

# . Undefined character.
_UndefinedCharacter = "."

#===================================================================================================================================
# . File reader class.
#===================================================================================================================================
class CIFFileReader ( STARFileReader ):
    """A class for reading CIF files."""

    defaultattributes = { }
    defaultattributes.update ( STARFileReader.defaultattributes )

    def AssignCoordinates ( self, dataBlock, system, xyzC, xyzF ):
        """Assign coordinates to a system."""
        # . Simplest treatment done.
        # . Assignment done in the following order:
        #   1. Cartesians if present.
        #   2. Fractionals transformed by pDynamo transformation.
        #   3. Fractionals transformed by data block transformation.
        if system is not None:
            if xyzC is not None:
                system.coordinates3 = xyzC
            elif xyzF is not None:
                sp = system.symmetryParameters
                if sp is not None:
                    rotation    = sp.M
                    translation = None
                else:
                    ( rotation, translation ) = self.ExtractTransformation ( dataBlock, "Cartn" )
                if rotation is not None:
                    xyzF.Rotate ( rotation )
                    if translation is not None: xyzF.Translate ( translation )
                    system.coordinates3 = xyzF

    def FileType ( self ): return "CIF"

    def ExtractAtomData ( self, dataBlock, blockCode ):
        """Extract atom data from a data block."""
        system = xyzC = xyzF = None
        table  = dataBlock.get ( "atom_site", None )
        if table is not None:
            # . Atom labels and symbols.
            labels  = table.ColumnValues ( "atom_site_label"       )
            symbols = table.ColumnValues ( "atom_site_type_symbol" )
            # . Make the symbols from the labels if necessary.
            if symbols is None:
                if labels is not None:
                    symbols = []
                    for label in labels:
                        characters = []
                        for c in label:
                            if c.isalpha ( ): characters.append ( c )
                            else:             break
                        symbols.append ( "".join ( characters ) )
            # . Make the basic system.
            if symbols is not None:
                system = System.FromAtoms ( symbols )
                system.label = blockCode
            # . Use the labels in the data block.
            if labels is not None:
                for ( atom, label ) in zip ( system.atoms, labels ):
                    atom.label = label
            # . Coordinates.
            xyzC = self.ExtractCoordinates ( table, "atom_site_Cartn" )
            xyzF = self.ExtractCoordinates ( table, "atom_site_fract" )
        return ( system, xyzC, xyzF )

    def ExtractCellData ( self, dataBlock ):
        """Extract crystal cell data from a data block."""
        a     = self.ExtractReal ( dataBlock, "cell_length_a"    )
        b     = self.ExtractReal ( dataBlock, "cell_length_b"    )
        c     = self.ExtractReal ( dataBlock, "cell_length_c"    )
        alpha = self.ExtractReal ( dataBlock, "cell_angle_alpha" )
        beta  = self.ExtractReal ( dataBlock, "cell_angle_beta"  )
        gamma = self.ExtractReal ( dataBlock, "cell_angle_gamma" )
        return ( a, b, c, alpha, beta, gamma )

    def ExtractCoordinates ( self, table, tag ):
        """Extract coordinates from a table."""
        x = table.RealColumnValues ( tag + "_x" )
        y = table.RealColumnValues ( tag + "_y" )
        z = table.RealColumnValues ( tag + "_z" )
        if ( x is not None ) and ( y is not None ) and ( z is not None ):
            xyz = Coordinates3.WithExtent ( len ( x ) )
            for ( i, ( u, v, w ) ) in enumerate ( zip ( x, y, z ) ):
                if ( u is None ) or ( v is None ) or ( w is None ):
                    xyz.FlagCoordinateAsUndefined ( i )
                else:
                    xyz[i,0] = u
                    xyz[i,1] = v
                    xyz[i,2] = w
        else: xyz = None
        return xyz

    def ExtractCrystalClass ( self, dataBlock ):
        """Extract the crystal class from a data block."""
        crystalSystem = None
        for tag in ( "space_group_crystal_system", "symmetry_cell_setting" ):
            result = self.ExtractString ( dataBlock, tag )
            if result is not None:
                crystalSystem = result.lower ( )
                break
        return _CrystalClasses.get ( crystalSystem, None )

    def ExtractCrystalData ( self, dataBlock, system ):
        """Extract crystal data from a data block."""
        if system is not None:
            # . Get the basic data.
            ( a, b, c, alpha, beta, gamma ) = self.ExtractCellData               ( dataBlock )
            crystalClass                    = self.ExtractCrystalClass           ( dataBlock )
            transformations                 = self.ExtractCrystalTransformations ( dataBlock )
            # . Check that all items exist.
            isOK = True
            for item in ( a, b, c, alpha, beta, gamma, crystalClass, transformations ):
                isOK = isOK and ( item is not None )
            # . Assign symmetry to the system.
            if isOK:
                system.DefineSymmetry ( a               = a               ,
                                        b               = b               ,
                                        c               = c               ,
                                        alpha           = alpha           ,
                                        beta            = beta            ,
                                        gamma           = gamma           ,
                                        crystalClass    = crystalClass    ,
                                        transformations = transformations )

    def ExtractCrystalTransformations ( self, dataBlock ):
        """Extract the crystal transformations from a data block."""
        # . Get the symmetry operations.
        operations = self.ExtractSymmetryOperations ( dataBlock, "space_group", "space_group_symop_operation_xyz" )
        if operations is None:
            operations = self.ExtractSymmetryOperations ( dataBlock, "symmetry_equiv_pos", "symmetry_equiv_pos_as_xyz" )
        if operations is None:
            return None
        else:
            transformations = []
            for operation in operations:
                transformations.append ( Transformation3_FromSymmetryOperationString ( operation ) )
            return Transformation3Container.WithTransformations ( transformations )

    def ExtractDisorderData ( self, dataBlock, system, disorderGroup = _DefaultDisorderGroup ):
        """Extract atom data from a data block."""
        result = system
        table  = dataBlock.get ( "atom_site", None )
        if table is not None:
            # . Disorder and multiplicity information.
            assemblies     = table.ColumnValues ( "atom_site_disorder_assembly"     )
            groups         = table.ColumnValues ( "atom_site_disorder_group"        )
            multiplicities = table.ColumnValues ( "atom_site_symmetry_multiplicity" )
            # . Check for the disorder groups to include.
            if groups is not None:
                uniqueGroups = set ( groups )
                uniqueGroups.discard (       disorderGroup )
                uniqueGroups.discard ( "-" + disorderGroup )
                uniqueGroups.discard ( _UndefinedCharacter )
                if len ( uniqueGroups ) > 0:
                    toInclude = set ( [ _UndefinedCharacter, disorderGroup, "-" + disorderGroup ] )
                    indices   = []
                    for ( i, group ) in enumerate ( groups ):
                        if group in toInclude: indices.append ( i )
                    if len ( indices ) < len ( groups ):
                        result       = system.Prune ( Selection ( indices ) )
                        result.label = system.label
        return result

    def ExtractSymmetryOperations ( self, dataBlock, shortTag, longTag ):
        """Extract symmetry operations from a data block."""
        operations = None
        # . Composite table.
        table = dataBlock.get ( shortTag, None )
        if isinstance ( table, STARFileTable ):
            operations = table.ColumnValues ( longTag )
        # . Single element table or a data value.
        if operations is None:
            item = dataBlock.get ( longTag, None )
            if item is not None:
                if isinstance ( item, STARFileTable ):
                    operations = item.ColumnValues ( longTag )
                else:
                    operations = [ item ]
        return operations
 
    def ExtractTransformation ( self, dataBlock, tag ):
        """Extract a transformation from a data block."""
        # . tag is "fract" for Cartesians -> fractionals and "Cartn" for Cartesians -> fractionals.
        isOK        = True
        rotation    = Matrix33 ( ) ; rotation.Set    ( 0.0 )
        translation = Vector3  ( ) ; translation.Set ( 0.0 )
        for i in range ( 3 ):
            for j in range ( 3 ):
                c = self.ExtractReal ( dataBlock, "atom_sites_{:s}_tran_matrix_{:d}{:d}".format ( tag, i, j ) )
                if c is None: isOK = False
                else:         rotation[i,j] = c
            c = self.ExtractReal ( dataBlock, "atom_sites_{:s}_tran_vector_{:d}".format ( tag, i ) )
            if c is None: isOK = False
            else:         translation[i] = c
        if isOK: return ( rotation, translation )
        else:    return ( None    , None        )

    def ToSystem ( self, blockCode = None, disorderGroup = _DefaultDisorderGroup ):
        """Extract a system from the file."""
        # . Initialization.
        system = None
        # . Get the data block.
        ( blockCode, dataBlock ) = self.GetDataBlock ( blockCode )
        if dataBlock is not None:
            # . Basic system.
            ( system, xyzC, xyzF ) = self.ExtractAtomData ( dataBlock, blockCode )
            # . Crystal data.
            self.ExtractCrystalData ( dataBlock, system )
            # . Assign coordinates.
            self.AssignCoordinates ( dataBlock, system, xyzC, xyzF )
            # . Prune system if necessary due to disorder.
            system = self.ExtractDisorderData ( dataBlock, system, disorderGroup = disorderGroup )
        return system

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def CIFFile_ToSystem ( path, blockCode = None, disorderGroup = _DefaultDisorderGroup, log = logFile ):
    """Helper function that reads a system from a CIF file."""
    inFile = CIFFileReader ( path )
    inFile.Parse   ( )
    inFile.Summary ( log = log )
    system = inFile.ToSystem ( blockCode = blockCode, disorderGroup = disorderGroup )
    return system

def CIFFile_ToSystems ( path, disorderGroup = _DefaultDisorderGroup, log = logFile ):
    """Helper function that reads all the systems from a CIF file."""
    inFile = CIFFileReader ( path )
    inFile.Parse   ( )
    inFile.Summary ( log = log )
    systems = []
    for key in sorted ( inFile.dataBlocks.keys ( ) ):
        systems.append ( inFile.ToSystem ( blockCode = key, disorderGroup = disorderGroup ) )
    return systems

# . Importer definitions.
_Importer.AddHandler ( { System : CIFFile_ToSystem } , [ "cif", "CIF" ], "Crystallographic Information File", defaultFunction = CIFFile_ToSystem )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
