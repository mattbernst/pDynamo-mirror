"""Utilities for compilation and installation."""

import glob, os, os.path, shutil, sys

# add this to detect the Windows platform in line 305
import platform

from distutils.core import setup, Extension
from distutils.util import get_platform

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Directory names.
CINCLUDEDIRECTORY   = "cinclude"
CLIBRARYDIRECTORY   = "clibrary"
EXTENSIONSDIRECTORY = "extensions"
PSOURCEDIRECTORY    = "psource"
PYREXDIRECTORY      = "pyrex"

# . Environment variables.
PDYNAMOROOT = "PDYNAMO_ROOT"

# . Extension names.
CEXTENSION            = ".c"
LIBRARYEXTENSION      = ".a"
PYREXEXTENSION        = ".pyx"
SHAREDOBJECTEXTENSION = ".so"

# . File names.
DEPENDENCYFILE = "dependencies"
LIBRARYFILE    = "libraries"

# . Libraries.
# . May need to add "irc" to libraries for Intel C compiler.
SYSTEMLIBRARIES    = [ "m" ]
SYSTEMLIBRARYPATHS = [ "/usr/local/lib", "/usr/lib" ]

#===================================================================================================================================
# . Class to allow compilation of C-libraries with no optimization. Only necessary for dlamch in clapack!
# . Ideally a cleaner and more general solution could be found.
#===================================================================================================================================
# . Specific imports.
from distutils                    import log
from distutils.command.build_clib import build_clib
from distutils.errors             import DistutilsSetupError
from types                        import ListType, TupleType

# . Class definition.
class BuildClibNoOpt ( build_clib ):
    """Subclass of distutils.command.build_clib which enables compiler options to be played with."""

    def build_libraries ( self, libraries ):
        """Build libraries."""
        # . Loop over libraries.
        for ( lib_name, build_info ) in libraries:
            # . Get sources.
            sources = build_info.get ( 'sources' )
            if ( sources is None ) or ( type ( sources ) not in ( ListType, TupleType ) ):
                raise DistutilsSetupError ( "in 'libraries' option (library '%s'), 'sources' must be present and must be a list of source filenames") % lib_name
            sources = list ( sources )
            # . Remove optimization strings from compiler options.
            new = []
            for token in self.compiler.compiler_so:
                if not token.startswith ( "-O" ): new.append ( token )
            self.compiler.compiler_so = new
            # . Do the work.
            log.info ( "building '%s' library", lib_name )
            # . Compile the source code to object files in the library directory.
            extra_preargs = build_info.get ( "extra_preargs" )
            include_dirs  = build_info.get ( "include_dirs"  )
            macros        = build_info.get ( "macros"        )
            objects       = self.compiler.compile ( sources, output_dir = self.build_temp, macros = macros, include_dirs = include_dirs, debug = self.debug, extra_preargs = extra_preargs )
            #  Create a static library.
            self.compiler.create_static_lib ( objects, lib_name, output_dir = self.build_clib, debug = self.debug )

#===================================================================================================================================
# . Directory install class.
#===================================================================================================================================
class DirectoryToInstall ( object ):
    """Class for a directory that is to be installed."""

    def __init__ ( self, **kwargs ):
        """Constructor."""
        for ( key, value ) in kwargs.iteritems ( ): setattr ( self, key, value )

    @staticmethod
    def Compare ( self, other ):
        """Comparison function."""
        # . Cross-dependence.
        QOTHERINSELF = ( other.name in self.dependencynames )
        QSELFINOTHER = ( self.name in other.dependencynames )
        if QOTHERINSELF and QSELFINOTHER: raise ValueError ( "The directories " + self.name + " and " + other.name + " are mutually dependent." )
        elif QOTHERINSELF: return  1
        elif QSELFINOTHER: return -1
        else:
            # . Number of dependencies.
            nself  = len (  self.dependencynames )
            nother = len ( other.dependencynames )
            if   nself > nother: return  1
            elif nself < nother: return -1
            else:
                # . Names.
                if   self.name > other.name: return  1
                elif self.name < other.name: return -1
                else:                        return  0

    def CheckForCLibraries ( self ):
        """Check for C libraries."""
        # . Read the library file.
        items  = []
        target = os.path.join ( self.extensionspath, LIBRARYFILE )
        if os.path.exists ( target ):
            tfile = file ( target, "r" )
            for line in tfile:
                tokens = line.split ( )
#                print tokens
                if ( len ( tokens ) > 0 ) and ( not tokens[0].startswith ( "#" ) ):
                    directoryname   = tokens[0]
                    libraryname     = directoryname
                    QNOOPTIMIZATION = False
                    QTHIRDPARTY     = False
                    if len ( tokens ) > 1:
                        libraryname = tokens[1]
                        if len ( tokens ) > 2:
                            QNOOPTIMIZATION = ( "NoOptimization" in tokens[2:] )
                            QTHIRDPARTY     = ( "ThirdParty"     in tokens[2:] )
                    items.append ( ( directoryname, libraryname, QNOOPTIMIZATION, QTHIRDPARTY ) )
            tfile.close ( )
        setattr ( self, "clibraries", items )

    def CheckForDependencies ( self ):
        """Check for directory dependencies."""
        items  = []
        target = os.path.join ( self.path, DEPENDENCYFILE )
        if os.path.exists ( target ):
            tfile = file ( target, "r" )
            for line in tfile:
                line = line.strip ( )
                if len ( line ) > 0: items.append ( line )
            tfile.close ( )
        setattr ( self, "dependencynames", items )

    def CheckForExtensions ( self ):
        """Check for directory extensions."""
        # . Extension directory.
        target = os.path.join ( self.path, EXTENSIONSDIRECTORY )
        if os.path.exists ( target ):
            setattr ( self, "extensionspath", target )
            # . Check for a C include directory.
            target = os.path.join ( self.extensionspath, CINCLUDEDIRECTORY )
            if os.path.exists ( target ): setattr ( self, "cincludepath", target )
            # . Check for libraries.
            self.CheckForCLibraries ( )
            # . Check for a C library directory and make it if necessary.
            if len ( self.clibraries ) > 0:
                target = os.path.join ( self.extensionspath, CLIBRARYDIRECTORY )
                if not os.path.exists ( target ): os.mkdir ( target )
                setattr ( self, "clibrarypath", target )
                setattr ( self, "QCLIBRARIES",  True   )
            # . Check for a Pyrex directory with files.
            target = os.path.join ( self.extensionspath, PYREXDIRECTORY )
            if os.path.exists ( target ):
                pyrexfiles = glob.glob ( os.path.join ( target, self.name + ".*" + PYREXEXTENSION ) )
                if len ( pyrexfiles ) > 0:
                    setattr ( self, "pyrexpath", target )
                    # . Process the file names.
                    pyrexroots = []
                    for pyrexfile in pyrexfiles:
                        ( head, tailext ) = os.path.split    ( pyrexfile )
                        ( tail, ext     ) = os.path.splitext ( tailext   )
                        pyrexroots.append ( tail )
                    pyrexroots.sort ( )
                    setattr ( self, "pyrexfileroots", pyrexroots )
                    setattr ( self, "QPYREXFILES",    True       )
                    # . Check for a psource directory.
                    target = os.path.join ( self.extensionspath, PSOURCEDIRECTORY )
                    if not os.path.exists ( target ): os.mkdir ( target )
                    setattr ( self, "psourcepath", target )

    def CheckPyrexFiles ( self, QVERBOSE = True ):
        """Check that the Pyrex files have been converted."""
        if getattr ( self, "QPYREXFILES", False ):
            errors = []
            for name in self.pyrexfileroots:
                ppath = os.path.join ( self.pyrexpath,   name + PYREXEXTENSION )
                cpath = os.path.join ( self.psourcepath, name + CEXTENSION     )
                if os.path.exists ( ppath ):
# . Could check for modification time here but no good if just unpacked distribution as these are not kept.
#                    QOK = os.path.exists ( cpath ) and ( os.path.getmtime ( cpath ) > os.path.getmtime ( ppath ) )
                    if not os.path.exists ( cpath ): errors.append ( name + PYREXEXTENSION )
            if len ( errors ) > 0:
                if QVERBOSE:
                    print "\nInvalid Pyrex-derived C files:\n"
                    for name in errors: print name
                    print
                raise ValueError ( "There are missing or out-of-date Pyrex-derived C files." )

    def CompileCLibraries ( self, QTHIRDPARTY = True, QVERBOSE = True ):
        """Compile the C-libraries."""
        if getattr ( self, "QCLIBRARIES", False ):
            # . Get the build directory name.
            buildpath = os.path.join ( "build", "temp" + ".%s-%s" % ( get_platform ( ), sys.version[0:3] ) )
            # . Get the include directories.
            includedirectories = [ self.cincludepath ]
            for item in reversed ( getattr ( self, "dependencyobjects", [] ) ):
                path = getattr ( item, "cincludepath", None )
                if path is not None: includedirectories.append ( path )
            # . Loop over libraries.
            for ( directoryname, libraryname, QNOOPTIMIZATION, QISTHIRDPARTY ) in self.clibraries:
                # . Check for compilation of a thirdparty library - avoids having to recompile these everytime something else is changed.
                if ( not QISTHIRDPARTY ) or ( QISTHIRDPARTY and QTHIRDPARTY ):
                    # . Get the source file list.
                    sourcefiles = glob.glob ( os.path.join ( self.extensionspath, directoryname, "*" + CEXTENSION ) )
                    if len ( sourcefiles ) > 0:
                        # . Make build_info.
                        build_info = { }
                        build_info["sources"]      = sourcefiles
                        build_info["include_dirs"] = includedirectories
                        build_info["macros"]       = None
                        # . Compile the library.
                        if QNOOPTIMIZATION:
                            build_info["extra_preargs"] = [ "-O0" ]
                            setup ( name = libraryname, libraries = [ ( libraryname, build_info ) ], script_args = [ "BuildClibNoOpt" ], cmdclass = { "BuildClibNoOpt" : BuildClibNoOpt } )
                        else:
                            setup ( name = libraryname, libraries = [ ( libraryname, build_info ) ], script_args = [ "build_clib"     ] )
                        # . Move the library to the appropriate place.
                        os.rename ( os.path.join ( buildpath, "lib" + libraryname + LIBRARYEXTENSION ), os.path.join ( self.clibrarypath, "lib" + libraryname + LIBRARYEXTENSION ) )

    def ConvertPyrexToC ( self, QCYTHON = False, QVERBOSE = True ):
        """Convert Pyrex to C files."""
        # . Import the correct modules.
        if QCYTHON: tag = "Cython"
        else:       tag = "Pyrex"
#        try:
        if QCYTHON: import Cython.Compiler.Errors as Errors, Cython.Compiler.Main as Main, Cython.Compiler.Version as Version
        else:       import Pyrex.Compiler.Errors  as Errors, Pyrex.Compiler.Main  as Main, Pyrex.Compiler.Version  as Version
        QPYREX = True
#        except:
#            QPYREX = False
        # . Compile.
        if QPYREX:
            if getattr ( self, "QPYREXFILES", False ):
                # . Get the pxddirectories (self and then dependencies in reverse order).
                pxddirectories = [ self.pyrexpath ]
                for item in reversed ( getattr ( self, "dependencyobjects", [] ) ):
                    path = getattr ( item, "pyrexpath", None )
                    if path is not None: pxddirectories.append ( path )
                # . Convert the files.
                if QVERBOSE: print "\nConverting files in " + self.pyrexpath + " with " + tag + " version " + Version.version + ":\n"
                for name in self.pyrexfileroots:
                    if QVERBOSE: print " -> " + name + PYREXEXTENSION
                    PyrexCompile ( Errors, Main, os.path.join ( self.pyrexpath, name + PYREXEXTENSION ), pxddirectories = pxddirectories )
                    os.rename ( os.path.join ( self.pyrexpath, name + CEXTENSION ), os.path.join ( self.psourcepath, name + CEXTENSION ) )
                if QVERBOSE: print
        else:
            raise ValueError ( "Error locating the " + tag + " compiler." )

    @classmethod
    def FromPathName ( selfclass, path ):
        """Constructor from path name."""
        ( head, tail ) = os.path.split ( path )
        name           = tail.rsplit ( "-", 1 )[0]
        self = selfclass ( name = name, path = path )
        return self

    def GetCLibraryNames ( self ):
        """Return a list of C library names in reverse order."""
        names = []
        for data in reversed ( getattr ( self, "clibraries", [] ) ): names.append ( data[1] )
        return names

    def HasExtensions ( self ):
        """Does this directory have extensions?"""
        return ( getattr ( self, "QCLIBRARIES", False ) or getattr ( self, "QPYREXFILES", False ) )

    def MakeExtensions ( self, QVERBOSE = True ):
        """Make extensions."""
        if getattr ( self, "QPYREXFILES", False ):
            # . Get the build and destination path names.
            buildpath       = os.path.join ( "build", "lib" + ".%s-%s" % ( get_platform ( ), sys.version[0:3] ), self.name )
            destinationpath = os.path.join ( self.path, self.name )
            # . Get the include directories.
            includedirectories = [ self.cincludepath ]
            for item in reversed ( getattr ( self, "dependencyobjects", [] ) ):
                path = getattr ( item, "cincludepath", None )
                if path is not None: includedirectories.append ( path )
            # . Get the library directories.
            librarydirectories = [ self.clibrarypath ]
            for item in reversed ( getattr ( self, "dependencyobjects", [] ) ):
                path = getattr ( item, "clibrarypath", None )
                if path is not None: librarydirectories.append ( path )
            librarydirectories.extend ( SYSTEMLIBRARYPATHS )
            # . Get the library names.
            clibraries = self.GetCLibraryNames ( )
            for item in reversed ( getattr ( self, "dependencyobjects", [] ) ):
                clibraries.extend ( item.GetCLibraryNames ( ) )
            clibraries.extend ( SYSTEMLIBRARIES )
            # . Make the list of extension modules.
            extensions = []
            for name in self.pyrexfileroots:
                extensions.append ( Extension ( name, [ os.path.join ( self.psourcepath, name + CEXTENSION ) ], include_dirs = includedirectories, libraries = clibraries, library_dirs = librarydirectories, runtime_library_dirs = librarydirectories ) )
            # . Compile the extension modules.
            setup ( name = self.name, ext_modules = extensions, script_args = [ "build_ext" ] )
            # . Move the files to the appropriate place
            # wew, provide support for Windows with MinGW and Cygwin
            osname = platform.system().lower()
            if osname.find("windows") > -1:
                SHAREDOBJECTEXTENSION = ".pyd"
            elif osname.find("cygwin") > -1:
                SHAREDOBJECTEXTENSION = ".dll"
                
            for name in self.pyrexfileroots:
                tail = name.split ( "." )[-1]
                os.rename ( os.path.join ( buildpath, tail + SHAREDOBJECTEXTENSION ), os.path.join ( destinationpath, tail + SHAREDOBJECTEXTENSION ) )

    def ResolveDependencies ( self, itemdictionary ):
        """Resolve the directory dependencies."""
        names = getattr ( self, "dependencynames", [] )
        items = []
        for name in names:
            try:    items.append ( itemdictionary[name] )
            except: raise ValueError ( "Unable to resolve dependency \"" + name + "\" for directory " + self.name + "." )
        setattr ( self, "dependencyobjects", items )

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def FindRootDirectory ( ):
    """Find or guess a root directory."""
    # . Is the variable defined?
    ROOTDIRECTORY = os.getenv ( PDYNAMOROOT )

    # . Guess a value from the current directory.
    if ROOTDIRECTORY is None:
        ( head, tail ) = os.path.split ( os.getcwd ( ) )
        while len ( tail ) > 0:
            if "dynamo" in tail.lower ( ):
                os.environ[PDYNAMOROOT] = os.path.join ( head, tail )
                break
            ( head, tail ) = os.path.split ( head )

def FindSubdirectory ( root, pattern ):
    """Find a subdirectory of root matching a particular pattern."""
    paths = glob.glob ( os.path.join ( root, pattern ) )
    if len ( paths ) == 1:
        ( head, tail ) = os.path.split ( paths[0] )
        return tail
    else:
        raise ValueError ( "Multiple subdirectories found matching " + pattern + "." )

def GetRootDirectories ( names ):
    """Get the root directories in increasing order of dependency.

    All directories are returned by default otherwise only those specified by |names|.
    """

    # . Get the root directory name.
    ROOTDIRECTORY = os.getenv ( PDYNAMOROOT )
    if ROOTDIRECTORY is None: raise ValueError ( "Unable to find pDynamo root directory." )

    # . Get all directories in root.
    toprocess = {}
    for candidate in glob.glob ( os.path.join ( ROOTDIRECTORY, "*" ) ):
        if os.path.isdir ( candidate ):
            directory = DirectoryToInstall.FromPathName ( path = candidate )
            directory.CheckForDependencies ( )
            directory.CheckForExtensions   ( )
            toprocess[directory.name] = directory

    # . Satisfy the dependencies of the directories.
    for directory in toprocess.itervalues ( ):
        directory.ResolveDependencies ( toprocess )

    # . Create the sorted directory list.
    directories = list ( toprocess.values ( ) )
    directories.sort ( cmp = DirectoryToInstall.Compare )

    # . Prune the list if necessary.
    if ( names is not None ) and ( len ( names ) > 0 ):
        items = directories
        directories = []
        for item in items:
            if item.name in names: directories.append ( item )

    # . Finish up.
    return directories

def PyrexCompile ( Errors, Main, source, pxddirectories = None ):
    """Convert a Pyrex file to C."""
    options = Main.default_options
    # . Change for Pyrex version 0.9.6: options is now a dictionary not a CompilationOptions object so convert to a CompilationOptions object.
    if isinstance ( options, dict ): options = Main.CompilationOptions ( options )
    if pxddirectories is not None: options.include_path.extend ( pxddirectories )
    context = Main.Context ( options.include_path )
    try:
        result   = context.compile ( source, options )
        QFAILURE = ( result.num_errors > 0 )
    except Errors.PyrexError, e:
        print >>sys.stderr, e
        QFAILURE = True
    if QFAILURE: raise ValueError ( "There was a Pyrex compiler error." )

def RemoveBuildDirectory ( ):
    """Remove the build directory if it exists."""
    if os.path.exists ( "build" ): shutil.rmtree ( "build" )

def WriteShellFile ( inpath, outpath, variables ):
    """Write a shell file."""
    infile   = file ( inpath, "r" )
    instring = infile.read ( ) % variables
    infile.close  ( )
    outfile  = file ( outpath, "w" )
    outfile.write ( instring )
    outfile.close ( )

def WriteShellFiles ( scratch = None ):
    """Write shell files containing the environment variables needed by pDynamo."""
    # . Get all required paths.
    # . Root.
    FindRootDirectory ( )
    pdynamoroot = os.getenv ( PDYNAMOROOT )
    # . Package and module directories.
    pbabel           = FindSubdirectory ( pdynamoroot, "pBabel-*"           )
    pcore            = FindSubdirectory ( pdynamoroot, "pCore-*"            )
    pmolecule        = FindSubdirectory ( pdynamoroot, "pMolecule-*"        )
    pmoleculescripts = FindSubdirectory ( pdynamoroot, "pMoleculeScripts-*" )
    # . Scratch.
    if scratch is None: scratch = "$PDYNAMO_ROOT/scratch"
    # . Get the new Python path.
    pythonpath = os.getenv ( "PYTHONPATH" )
    if pythonpath is None: pythonpath = ""
    oldtokens = pythonpath.split ( ":" )
    newtokens = []
    for token in oldtokens:
        if ( token.find ( pbabel ) < 0 ) and ( token.find ( pcore ) < 0 ) and ( token.find ( pmolecule ) < 0 ) and ( token.find ( pmoleculescripts ) < 0 ): newtokens.append ( token )
    newtokens.extend ( [ "$PDYNAMO_ROOT/" + pbabel, "$PDYNAMO_ROOT/" + pcore, "$PDYNAMO_ROOT/" + pmolecule, "$PDYNAMO_ROOT/" + pmoleculescripts ] )
    pythonpath = ":".join ( newtokens )
    # . Write the shell files.
    variables = { "pbabel" : pbabel, "pcore" : pcore, "pdynamoroot" : pdynamoroot, "pmolecule" : pmolecule, "pmoleculescripts" : pmoleculescripts, "pythonpath" : pythonpath, "scratch" : scratch }
    WriteShellFile ( os.path.join ( "templates", "environment_bash.tmpl"   ), "environment_bash.com",   variables )
    WriteShellFile ( os.path.join ( "templates", "environment_cshell.tmpl" ), "environment_cshell.com", variables )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
