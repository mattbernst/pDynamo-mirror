#-------------------------------------------------------------------------------
# . File      : mmCIFModel.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Ancillary data structures needed when building mmCIF systems and reading mmCIF files."""

from pCore     import Coordinates3, logFile, LogFileActive, Selection
from pMolecule import Atom, PeriodicTable, Sequence, SequenceLinearPolymer, System

# . This module should be merged into PDBModel as ultimately they are doing the same thing.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Default entity type.
_DEFAULTENTITYTYPE = "non-polymer"

# . Default label.
_DefaultLabel = ""

# . Default model number.
_DEFAULTMODELNUMBER = "1"

# . Default occupancy.
_DEFAULTOCCUPANCY = 1.0

# . Entity types.
_ENTITYTYPES = ( "non-polymer", "polymer", "water" )

# . Separators for keys and labels.
_MAJORSEPARATOR = "/"
_MINORSEPARATOR = "|"

# . Undefined character.
_UNDEFINEDCHARACTER = "."

# . Undefined float.
_UNDEFINEDFLOAT = 9999.0

# . Undefined integer.
_UNDEFINEDINTEGER = 9999

# . Unknown character.
_UNKNOWNCHARACTER = "?"

# . The wildcard character.
_WILDCARD = "*"

#===================================================================================================================================
# . Base class.
#===================================================================================================================================
class mmCIFModelItem ( object ):
    """Base class for all mmCIF model items."""

    defaultattributes = { }

    def __init__ ( self, *args, **kwargs ):
        """Constructor."""
        self.__dict__.update ( self.__class__.defaultattributes )
        for ( key, value ) in kwargs.iteritems ( ): setattr ( self, key, value )
        self.InitializeContainerAttributes ( )

    @staticmethod
    def FetchFromKeywordArguments ( kwargs, key, default = _UNDEFINEDCHARACTER ):
        """Get a value from kwargs."""
        value = kwargs.get ( key, default )
        if value == _UNKNOWNCHARACTER: value = default
        return value

    @staticmethod
    def GetIndex ( item ):
        """Get the value of the item index."""
        return getattr ( item, "index", -1 )

    def InitializeContainerAttributes ( self ): pass

    @staticmethod
    def MakeKey ( *args ): return _MINORSEPARATOR.join ( args )

    @staticmethod
    def ParseKey ( key ): return tuple ( key.split ( _MINORSEPARATOR ) )

    @staticmethod
    def ToFloat ( value, default = _UNDEFINEDFLOAT ):
        """Convert a value to a float."""
        if value == _UNDEFINEDCHARACTER: return default
        else:                            return float ( value )

    @staticmethod
    def ToInteger ( value, default = _UNDEFINEDINTEGER ):
        """Convert a value to an integer."""
        if value == _UNDEFINEDCHARACTER: return default
        else:                            return int ( value )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFModel ( mmCIFModelItem ):
    """Hold data pertaining to a mmCIF model."""

    # . Attributes that should be defined.
    defaultattributes = { "label"           : None,    \
                          "log"             : logFile, \
                          "maximumwarnings" : 100,     \
                          "nfatal"          : 0,       \
                          "nwarnings"       : 0,       \
                          "QFINALIZED"      : False,   \
                          "warningtable"    : None     }
    defaultattributes.update ( mmCIFModelItem.defaultattributes )

    def AddAsymmetricUnitDefinition ( self, **kwargs ):
        """Add an asymmetric unit definition."""
        asymmetricunitdefinition = None
        index                    = len ( self.asymmetricunitdefinitions )
        if self.QFINALIZED:
            raise ValueError ( "Cannot add an asymmetric unit definition to a finalized model." )
        else:
            # . Get the entity.
            entitylabel      = mmCIFModel.FetchFromKeywordArguments ( kwargs, "entity_id" )
            entitydefinition = self.entitydefinitions.get ( entitylabel, None )
            if entitydefinition is None: entitydefinition = self.AddEntityDefinition ( id = entitylabel, type = _DEFAULTENTITYTYPE )
            # . Get the asymmetric unit definition.
            if entitydefinition is not None:
                label = mmCIFModel.FetchFromKeywordArguments ( kwargs, "id" )
                if label in self.asymmetricunitdefinitions.keys ( ):
                    self.Warning ( "Asymmetric Unit Definition {:d}".format ( index ), "Adding duplicate asymmetric unit definition (" + label + ") to model.", False )
                else:
                    asymmetricunitdefinition = mmCIFModelAsymmetricUnitDefinition ( entitydefinition = entitydefinition, index = index, label = label )
                    self.asymmetricunitdefinitions[label] = asymmetricunitdefinition
        return asymmetricunitdefinition

    def AddAtomSite ( self, **kwargs ):
        """Add an atom site."""
        atom = None
        if self.QFINALIZED:
            raise ValueError ( "Cannot add an atom site to a finalized model." )
        else:
            # . Initialization.
            warningtag = "Atom Site {:d}".format ( len ( self.atoms ) )
            # . Get the entity.
            auth_asym_id    = mmCIFModel.FetchFromKeywordArguments ( kwargs, "auth_asym_id"    )
            label_asym_id   = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_asym_id"   )
            label_entity_id = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_entity_id" )
            entitykey       = mmCIFModelEntity.MakeKey ( label_entity_id, label_asym_id )
            entity          = self.entities.get ( entitykey, None )
            if entity is None: entity = self.AddEntity ( entitykey )
            # . Check the alternate name.
            if entity is not None:
                if   entity.alternateasymlabel is None: entity.alternateasymlabel = auth_asym_id
                elif entity.alternateasymlabel != auth_asym_id: self.Warning ( warningtag, "Mismatch in alternate names for entity asymmetric unit - " + entity.alternateasymlabel + " and " + auth_asym_id + ".", False )
                # . Get the component.
                auth_comp_id   = mmCIFModel.FetchFromKeywordArguments ( kwargs, "auth_comp_id"      )
                auth_seq_id    = mmCIFModel.FetchFromKeywordArguments ( kwargs, "auth_seq_id"       )
                insertioncode  = mmCIFModel.FetchFromKeywordArguments ( kwargs, "pdbx_pdb_ins_code" )
                name           = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_comp_id"     )
                sequencenumber = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_seq_id"      )
                if name == _UNDEFINEDCHARACTER:
                    name         = auth_comp_id
                    auth_comp_id = None
                if sequencenumber == _UNDEFINEDCHARACTER:
                    sequencenumber = auth_seq_id
                    auth_seq_id    = None
                componentkey      = mmCIFModelComponent.MakeKey ( name, sequencenumber, insertioncode )
                component         = entity.components.get ( componentkey, None )
                if component is None:
                    component = entity.AddComponent ( componentkey )
                    self.components.append ( component )
                # . Check alternate names.
                if component is not None:
                    if auth_comp_id is not None:
                        if   component.alternatename     is None         : component.alternatename     = auth_comp_id
                        elif component.alternatename     != auth_comp_id : self.Warning ( warningtag, "Mismatch in alternate names for component - "     + component.alternatename     + " and " + auth_comp_id + ".", False )
                    if auth_seq_id  is not None:
                        if   component.alternatesequence is None         : component.alternatesequence = auth_seq_id
                        elif component.alternatesequence != auth_seq_id  : self.Warning ( warningtag, "Mismatch in alternate sequences for component - " + component.alternatesequence + " and " + auth_seq_id  + ".", False )
                    # . Get the atom data.
                    atomicNumber = PeriodicTable.AtomicNumber ( mmCIFModel.FetchFromKeywordArguments ( kwargs, "type_symbol" ) )
                    auth_atom_id = mmCIFModel.FetchFromKeywordArguments ( kwargs, "auth_atom_id"  )
                    name         = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_atom_id" )
                    if name == _UNDEFINEDCHARACTER:
                        name         = auth_atom_id
                        auth_atom_id = None
                    atom    = component.atoms.get ( name, None )
                    atomKey = name
                    if atom is None:
                        atom = component.AddAtom ( name, atomicNumber = atomicNumber, entity = entity )
                        self.atoms.append ( atom )
                    # . Check alternate labels.
                    if atom is not None:
                        if auth_comp_id is not None:
                            if   atom.alternatelabel is None         : atom.alternatelabel = auth_atom_id
                            elif atom.alternatelabel != auth_atom_id : self.Warning ( warningtag, "Mismatch in alternate labels for atom - " + atom.alternatelabel + " and " + auth_atom_id + ".", False )
                        # . Get the atom data.
                        label_alt_id       = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_alt_id" )
                        pdbx_pdb_model_num = mmCIFModel.FetchFromKeywordArguments ( kwargs, "pdbx_pdb_model_num", default = _DEFAULTMODELNUMBER )
                        self.altLocs.add      ( label_alt_id       )
                        self.modelNumbers.add ( pdbx_pdb_model_num )
                        datakey            = mmCIFModelAtomData.MakeKey ( label_alt_id, pdbx_pdb_model_num )
                        data               = atom.data.get ( datakey, None )
                        if data is None:
                            try:
                                x         = mmCIFModel.ToFloat ( mmCIFModel.FetchFromKeywordArguments ( kwargs, "cartn_x" ) )
                                y         = mmCIFModel.ToFloat ( mmCIFModel.FetchFromKeywordArguments ( kwargs, "cartn_y" ) )
                                z         = mmCIFModel.ToFloat ( mmCIFModel.FetchFromKeywordArguments ( kwargs, "cartn_z" ) )
                                occupancy = mmCIFModel.ToFloat ( mmCIFModel.FetchFromKeywordArguments ( kwargs, "occupancy" ), default = _DEFAULTOCCUPANCY )
                            except:
                                self.Warning ( warningtag, "Unable to convert atom site floating point data.", False )
                            atom.data[datakey] = mmCIFModelAtomData ( altLoc = label_alt_id, atom = atom, model = pdbx_pdb_model_num, x = x, y = y, z = z, occupancy = occupancy )
                        else:
                            self.Warning ( warningtag, "Adding duplicate atomsite (" + _MAJORSEPARATOR.join ( [ datakey, atomKey, componentkey, entitylabel ] ) + ") to model.", False )
        return atom

    def AddConnection ( self, **kwargs ):
        """Add a connection."""
        connection = None
        index      = len ( self.connections )
        if self.QFINALIZED:
            raise ValueError ( "Cannot add a connection to a finalized model." )
        else:
            # . Identify the atoms in the connection.
            definitions = {}
            for ( key, value ) in kwargs.iteritems ( ):
                if "ptnr1_" in key: definitions[key.replace ( "ptnr1_", "" )] = value
            ( entity1, component1, atom1, atomdata1 ) = self.GetAtomObjects ( **definitions )
            definitions = {}
            for ( key, value ) in kwargs.iteritems ( ):
                if "ptnr2_" in key: definitions[key.replace ( "ptnr2_", "" )] = value
            ( entity2, component2, atom2, atomdata2 ) = self.GetAtomObjects ( **definitions )
            # . Get the id and type.
            id   = mmCIFModel.FetchFromKeywordArguments ( kwargs, "id"      )
            type = mmCIFModel.FetchFromKeywordArguments ( kwargs, "type_id" )
            # . Add the connection.
            if ( atom1 is not None ) and ( atom2 is not None ) and ( atomdata1 is not None ) and ( atomdata2 is not None ) and ( entity1 is not None ) and ( entity2 is not None ):
                symmetry1  = mmCIFModel.FetchFromKeywordArguments ( kwargs, "ptnr1_symmetry" )
                symmetry2  = mmCIFModel.FetchFromKeywordArguments ( kwargs, "ptnr2_symmetry" )
                connection = mmCIFModelConnection ( atom1 = atom1, atom2 = atom2, entity1 = entity1, entity2 = entity2, symmetry1 = symmetry1, symmetry2 = symmetry2, type = type )
                self.connections.append  ( connection )
                self.connectiontypes.add ( type       )
            else:
                if id == _UNDEFINEDCHARACTER: self.Warning ( "Connection {:d}".format ( index ), "Unable to identify atoms in connection.",                      False )
                else:                         self.Warning ( "Connection {:d}".format ( index ), "Unable to identify atoms in connection with id - " + id + ".", False )
        return connection

    def AddEntity ( self, key ):
        """Add an entity."""
        ( label, asymlabel ) = mmCIFModelEntity.ParseKey ( key )
        # . Find the type.
        if label in self.entitydefinitions: type = self.entitydefinitions[label].type
        else:                               type = _DEFAULTENTITYTYPE
        # . Create the entity.
        index = len ( self.entities )
        if   type == "non-polymer" : entity = mmCIFModelNonPolymerEntity ( asymlabel = asymlabel, index = index, label = label )
        elif type == "polymer"     : entity = mmCIFModelPolymerEntity    ( asymlabel = asymlabel, index = index, label = label )
        elif type == "water"       : entity = mmCIFModelWaterEntity      ( asymlabel = asymlabel, index = index, label = label )
        self.entities[key] = entity
        return entity

    def AddEntityDefinition ( self, **kwargs ):
        """Add an entity definition."""
        entitydefinition = None
        if self.QFINALIZED:
            raise ValueError ( "Cannot add an entity definition to a finalized model." )
        else:
            index = len ( self.entitydefinitions )
            label = mmCIFModel.FetchFromKeywordArguments ( kwargs, "id" )
            if label in self.entitydefinitions.keys ( ):
                self.Warning ( "Entity Definition {:d}".format ( index ), "Adding duplicate entity definition (" + label + ") to model.", False )
            else:
                type = mmCIFModel.FetchFromKeywordArguments ( kwargs, "type" ).lower ( )
                if type not in _ENTITYTYPES:
                    if type != _UNDEFINEDCHARACTER: self.Warning ( "Entity Definition {:d}".format ( index ), "Unknown entity definition type - " + type + ".", False )
                    type = _DEFAULTENTITYTYPE
                entitydefinition = mmCIFModelEntityDefinition ( label = label, type = type )
                self.entitydefinitions[label] = entitydefinition
        return entitydefinition

    def Finalize ( self ):
        """Finalize construction of the model."""
        self.WarningStop ( )
        for entity in self.entities.values ( ): entity.Finalize ( )
        self.QFINALIZED = True

    def GetAtomObjects ( self, **kwargs ):
        """Get the objects pertaining to an atom."""
        # . Initialization.
        atom      = None
        atomdata  = None
        component = None
        entity    = None
        # . Get the entity.
        # . Be optimistic and try specific match first.
        label_asym_id   = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_asym_id"   )
        label_entity_id = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_entity_id" )
        entitykey       = mmCIFModelEntity.MakeKey ( label_entity_id, label_asym_id )
        entity          = self.entities.get ( entitykey, None )
        # . Failed so search for label_asym_id.
        if entity is None:
            for ( key, value ) in self.entities.iteritems ( ):
                if key.endswith ( _MINORSEPARATOR + label_asym_id ):
                    entity = value
                    break
        # . Get the component.
        if entity is not None:
            insertioncode  = mmCIFModel.FetchFromKeywordArguments ( kwargs, "pdbx_pdb_ins_code" )
            name           = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_comp_id"     )
            sequencenumber = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_seq_id"      )
            if name           == _UNDEFINEDCHARACTER: name           = mmCIFModel.FetchFromKeywordArguments ( kwargs, "auth_comp_id" )
            if sequencenumber == _UNDEFINEDCHARACTER: sequencenumber = mmCIFModel.FetchFromKeywordArguments ( kwargs, "auth_seq_id"  )
            componentkey      = mmCIFModelComponent.MakeKey ( name, sequencenumber, insertioncode )
            component         = entity.components.get ( componentkey, None )
            # . Get the atom.
            if component is not None:
                name = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_atom_id" )
                if name == _UNDEFINEDCHARACTER: name = mmCIFModel.FetchFromKeywordArguments ( kwargs, "auth_atom_id" )
                atom = component.atoms.get ( name, None )
                # . Get the atom data.
                if atom is not None:
                    label_alt_id       = mmCIFModel.FetchFromKeywordArguments ( kwargs, "label_alt_id" )
                    pdbx_pdb_model_num = mmCIFModel.FetchFromKeywordArguments ( kwargs, "pdbx_pdb_model_num", default = _DEFAULTMODELNUMBER )
                    datakey            = mmCIFModelAtomData.MakeKey ( label_alt_id, pdbx_pdb_model_num )
                    atomdata           = atom.data.get ( datakey, None )
        return ( entity, component, atom, atomdata )

    def InitializeContainerAttributes ( self ):
        """Initialize container attributes."""
        self.altLocs                   = set ( )
        self.asymmetricunitdefinitions = {}
        self.atoms                     = []
        self.components                = []
        self.connections               = []
        self.connectiontypes           = set ( )
        self.entities                  = {}
        self.entitydefinitions         = {}
        self.modelNumbers              = set ( )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            summary = log.GetSummary ( pageWidth = 100 )
            if self.label is None: summary.Start ( "mmCIF Model Summary" )
            else:                  summary.Start ( "Summary for mmCIF Model " + self.label )
            summary.Entry ( "Atoms",                      "{:d}".format ( len ( self.atoms                     ) ) )
            summary.Entry ( "Components",                 "{:d}".format ( len ( self.components                ) ) )
            summary.Entry ( "Entities",                   "{:d}".format ( len ( self.entities                  ) ) )
            summary.Entry ( "AltLocs",                    "{:d}".format ( len ( self.altLocs                   ) ) )
            summary.Entry ( "Asymmetric Unit Definitons", "{:d}".format ( len ( self.asymmetricunitdefinitions ) ) )
            summary.Entry ( "Entity Definitions",         "{:d}".format ( len ( self.entitydefinitions         ) ) )
            summary.Entry ( "Connections",                "{:d}".format ( len ( self.connections               ) ) )
            summary.Entry ( "Connection Types",           "{:d}".format ( len ( self.connectiontypes           ) ) )
            summary.Entry ( "Models",                     "{:d}".format ( len ( self.modelNumbers              ) ) )
            summary.Entry ( "Warnings",                   "{:d}".format (       self.nwarnings                   ) )
            summary.Stop ( )

    def ToSystem ( self, altLoc = _UNDEFINEDCHARACTER, modelNumber = _DEFAULTMODELNUMBER, embeddedHydrogens = False ):
        """Return a system from the model."""
        system = None
        if self.QFINALIZED:
            # . Initialization.
            index          =  0
            atomPaths      = []
            majorSeparator = Sequence.defaultAttributes["labelSeparator"]
            minorSeparator = Sequence.defaultAttributes["fieldSeparator"]
            systemAtoms    = []
            undefined      = []
            xyz            = []
            # . Loop over entities.
            entities = list ( self.entities.values ( ) )
            entities.sort ( key = mmCIFModelEntity.GetIndex )
            entityLabels = {}
            for entity in entities:
                # . Entity label.
                if entity.asymlabel == _UNDEFINEDCHARACTER: entityLabel = _DefaultLabel
                else:                                       entityLabel = entity.asymlabel
                if ( entity.asymlabel != entity.label ) and ( entity.label != _UNDEFINEDCHARACTER ): entityLabel += ( minorSeparator + entity.label )
                entityLabels[entity] = entityLabel
                # . Loop over components.
                components = list ( entity.components.values ( ) )
                components.sort ( key = mmCIFModelComponent.GetIndex )
                for component in components:
                    # . Component path.
                    fields = []
                    for item in ( component.name, component.sequencenumber, component.insertioncode ):
                        if item == _UNDEFINEDCHARACTER: fields.append ( _DefaultLabel )
                        else:                           fields.append ( item          )
                    n = 3
                    for item in reversed ( fields ):
                        if item == _DefaultLabel: n -= 1
                        else: break
                    componentLabel = minorSeparator.join ( fields[0:max(n,1)] )
                    # . Loop over atoms - putting hydrogens at the end if necessary.
                    atoms = list ( component.atoms.values ( ) )
                    if embeddedHydrogens: atoms.sort ( key = mmCIFModelAtom.GetIndex                    )
                    else:                 atoms.sort ( key = mmCIFModelAtom.GetNonEmbeddedHydrogenIndex )
                    # . Generate atoms and associated data.
                    for atom in atoms:
                        systemAtoms.append ( Atom ( atomicNumber = atom.atomicNumber, formalCharge = atom.formalCharge ) )
                        if atom.label == _UNDEFINEDCHARACTER: atomLabel = _DefaultLabel
                        else:                                 atomLabel = atom.label
                        atomPaths.append ( entityLabel + majorSeparator + componentLabel + majorSeparator + atomLabel )
                        ( QUNDEFINED, x, y, z ) = atom.GetCoordinateData ( altLoc, modelNumber )
                        if QUNDEFINED: undefined.append ( index )
                        xyz.append ( ( x, y, z ) )
                        index += 1
            # . Make the sequence.
            sequence = Sequence.FromAtomPaths ( atomPaths, atoms = systemAtoms )
            # . Make the system.
            system = System.FromSequence ( sequence )
            if self.label is not None: system.label = self.label
            # . Coordinates.
            coordinates3 = Coordinates3.WithExtent ( len ( systemAtoms ) )
            for ( i, ( x, y, z ) ) in enumerate ( xyz ):
                coordinates3[i,0] = x
                coordinates3[i,1] = y
                coordinates3[i,2] = z
            system.coordinates3 = coordinates3
            # . Undefined coordinates.
            if len ( undefined ) > 0:
                for i in undefined:
                    system.coordinates3.FlagCoordinateAsUndefined ( i )
            # . Define polymer data.
            for oldEntity in entities:
                if isinstance ( oldEntity, mmCIFModelPolymerEntity ):
                    newEntity = system.sequence.childIndex.get ( entityLabels[oldEntity] )
                    system.sequence.linearPolymers.append ( SequenceLinearPolymer ( isCyclic               = oldEntity.QCYCLIC      , \
                                                                                    leftTerminalComponent  = newEntity.children[ 0] , \
                                                                                    rightTerminalComponent = newEntity.children[-1] ) )
            # . Finish up.
            return system

    def Warning ( self, tag, message, QFATAL ):
        """Print a warning."""
        if ( self.log is not None ) and ( self.maximumwarnings > 0 ):
            if self.nwarnings == 0: self.WarningStart ( )
            self.nwarnings += 1
            if self.nwarnings <= self.maximumwarnings:
                self.warningtable.Entry ( tag, alignment = "l" )
                self.warningtable.Entry ( message )
                if QFATAL:
                    self.warningtable.Entry ( "Fatal" )
                    self.nfatal += 1
                else:
                    self.warningtable.Entry ( "Warning" )

    def WarningStart ( self ):
        """Start warning printing."""
        if ( self.log is not None ):
            self.warningtable = self.log.GetTable ( columns = [ 20, 80, 10 ] )
            self.warningtable.Start ( )
            self.warningtable.Title ( "mmCIF Model Building Warnings" )
            self.warningtable.Heading ( "Tag"      )
            self.warningtable.Heading ( "Message"  )
            self.warningtable.Heading ( "Severity" )

    def WarningStop ( self ):
        """Terminate warning printing."""
        if ( self.warningtable is not None ):
            self.warningtable.Stop ( )
            self.warningtable = None
            if self.nfatal > 0:
                self.log.Paragraph ( "There have been fatal errors!" )
                raise ValueError ( "Fatal errors generating model." )
            else:
                self.log.Paragraph ( "There have been warnings. Proceed with caution!" )
            self.log.LineBreak ( )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFModelAsymmetricUnitDefinition ( mmCIFModelItem ):
    """Define a mmCIFModel asymmetric unit definition."""

    # . Attributes that should be defined.
    defaultattributes = { "entitydefinition" : None, \
                          "label"            : _UNDEFINEDCHARACTER }
    defaultattributes.update ( mmCIFModelItem.defaultattributes )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFModelAtom ( mmCIFModelItem ):
    """Define a mmCIFModel atom."""

    # . Attributes that should be defined.
    defaultattributes = { "alternatelabel" : None, \
                          "atomicNumber"   :   -1, \
                          "component"      : None, \
                          "entity"         : None, \
                          "formalCharge"   :    0, \
                          "index"          :   -1, \
                          "label"          : _UNDEFINEDCHARACTER }
    defaultattributes.update ( mmCIFModelItem.defaultattributes )

    def GetCoordinateData ( self, altLoc, modelNumber ):
        """Get coordinate data given an altLoc and modelNumber."""
        # . Initialization.
        data = None
        x    = _UNDEFINEDFLOAT
        y    = _UNDEFINEDFLOAT
        z    = _UNDEFINEDFLOAT
        # . Try the specific altLoc if this exists.
        if altLoc != _UNDEFINEDCHARACTER:
            datakey = mmCIFModelAtomData.MakeKey ( altLoc, modelNumber )
            data    = self.data.get ( datakey, None )
        # . Try for undefined altLoc.
        if data is None:
            datakey = mmCIFModelAtomData.MakeKey ( _UNDEFINEDCHARACTER, modelNumber )
            data    = self.data.get ( datakey, None )
        # . Get the highest-occupancy data for the model.
        if data is None:
            for ( key, value ) in self.data.iteritems ( ):
                if key.endswith ( _MINORSEPARATOR + modelNumber ):
                    if ( data is None ) or ( data.occupancy < value.occupancy ): data = value
        # . Data found.
        if data is not None:
            x = data.x
            y = data.y
            z = data.z
        # . Finish up.
        QUNDEFINED = ( x == _UNDEFINEDFLOAT ) or ( y == _UNDEFINEDFLOAT ) or ( z == _UNDEFINEDFLOAT )
        return ( QUNDEFINED, x, y, z )

    @staticmethod
    def GetNonEmbeddedHydrogenIndex ( item ):
        """Generate an index which ensures hydrogens have lower priority than other atoms."""
        if item.atomicNumber == 1: weight = 1
        else:                      weight = 0
        return ( weight, item.index )

    def InitializeContainerAttributes ( self ):
        """Initialize container attributes."""
        self.data = {}

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFModelAtomData ( mmCIFModelItem ):
    """Define coordinate data for a mmCIFModel atom."""

    # . Attributes that should be defined.
    defaultattributes = { "altLoc"    : _UNDEFINEDCHARACTER, \
                          "atom"      : None, \
                          "model"     : _UNDEFINEDCHARACTER, \
                          "occupancy" : _UNDEFINEDFLOAT,     \
                          "x"         : _UNDEFINEDFLOAT,     \
                          "y"         : _UNDEFINEDFLOAT,     \
                          "z"         : _UNDEFINEDFLOAT      }
    defaultattributes.update ( mmCIFModelItem.defaultattributes )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFModelConnection ( mmCIFModelItem ):
    """Define a mmCIFModel connection."""

    # . Attributes that should be defined.
    defaultattributes = { "atom1"     : None, \
                          "atom2"     : None, \
                          "entity1"   : None, \
                          "entity2"   : None, \
                          "symmetry1" : _UNDEFINEDCHARACTER, \
                          "symmetry2" : _UNDEFINEDCHARACTER, \
                          "type"      : _UNDEFINEDCHARACTER  }
    defaultattributes.update ( mmCIFModelItem.defaultattributes )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFModelComponent ( mmCIFModelItem ):
    """Hold data for a mmCIFModel component."""

    # . Attributes that should be defined.
    defaultattributes = { "alternatename"     : None, \
                          "alternatesequence" : None, \
                          "entity"            : None, \
                          "index"             :   -1, \
                          "insertioncode"     : _UNDEFINEDCHARACTER, \
                          "name"              : _UNDEFINEDCHARACTER, \
                          "sequencenumber"    : _UNDEFINEDCHARACTER  }
    defaultattributes.update ( mmCIFModelItem.defaultattributes )

    def AddAtom ( self, key, **kwargs ):
        """Add an atom."""
        atom = mmCIFModelAtom ( atomicNumber = kwargs.get ( "atomicNumber", -1 ), component = self, entity = kwargs.get ( "entity", None ), index = len ( self.atoms ), label = key )
        self.atoms[key] = atom
        return atom

    def InitializeContainerAttributes ( self ):
        """Initialize container attributes."""
        self.atoms = {}

#===================================================================================================================================
# . Classes for entities.
#===================================================================================================================================
# . Basic entity.
class mmCIFModelEntity ( mmCIFModelItem ):
    """Base class for mmCIFModel entities."""

    # . Attributes that should be defined.
    defaultattributes = { "alternateasymlabel" : None,                \
                          "asymlabel"          : _UNDEFINEDCHARACTER, \
                          "components"         : None,                \
                          "index"              : -1,                  \
                          "label"              : _UNDEFINEDCHARACTER  }
    defaultattributes.update ( mmCIFModelItem.defaultattributes )

    def AddComponent ( self, key ):
        """Add a component."""
        ( name, sequencenumber, insertioncode ) = mmCIFModelComponent.ParseKey ( key )
        component = mmCIFModelComponent ( entity = self, index = len ( self.components ), insertioncode = insertioncode, name = name, sequencenumber = sequencenumber )
        self.components[key] = component
        return component

    def InitializeContainerAttributes ( self ):
        """Initialize container attributes."""
        self.components = {}

    def Finalize ( self ):
        """Finalize a model."""
        pass

    def Type ( self ): return "Base Class"

# . Non-polymer entity.
class mmCIFModelNonPolymerEntity ( mmCIFModelEntity ):
    """Container for mmCIFModel non-polymer entities."""

    def Type ( self ): return "Non-Polymer"

# . Polymer entity.
class mmCIFModelPolymerEntity ( mmCIFModelEntity ):
    """Container for mmCIFModel polymer entities."""

    # . Attributes that should be defined.
    defaultattributes = { "leftterminalcomponent"  :  None, \
                          "QCYCLIC"                : False, \
                          "rightterminalcomponent" :  None  }
    defaultattributes.update ( mmCIFModelEntity.defaultattributes )

    def Finalize ( self ):
        """Finalize a model."""
        # . Determine the left and right terminii of the polymer.
        if len ( self.components ) > 0:
            components = list ( self.components.values ( ) )
            components.sort ( key = mmCIFModelComponent.GetIndex )
            self.leftterminalcomponent  = components[ 0]
            self.rightterminalcomponent = components[-1]

    def Type ( self ): return "Polymer"

# . Water entity.
class mmCIFModelWaterEntity ( mmCIFModelEntity ):
    """Container for mmCIFModel water entities."""

    def Type ( self ): return "Water"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFModelEntityDefinition ( mmCIFModelItem ):
    """Define a mmCIFModel entity definition."""

    # . Attributes that should be defined.
    defaultattributes = { "label" : _UNDEFINEDCHARACTER, \
                          "type"  : _DEFAULTENTITYTYPE   }
    defaultattributes.update ( mmCIFModelItem.defaultattributes )

#===================================================================================================================================
# . Test the module.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
