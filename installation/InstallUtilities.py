"""Utilities for compilation and installation."""

import functools, glob, os, os.path, shutil, sys

from distutils.core import setup, Extension
from distutils.util import get_platform

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Directory names.
_CIncludeDirectory   = "cinclude"
_CLibraryDirectory   = "clibrary"
_ExtensionsDirectory = "extensions"
_PSourceDirectory    = "psource"
_PyrexDirectory      = "pyrex"

# . Environment variables.
_PDynamoRoot = "PDYNAMO_ROOT"

# . Extension names.
_CExtension            = ".c"
_LibraryExtension      = ".a"
_PyrexExtension        = ".pyx"
_SharedObjectExtension = ".so"

# . File names.
_DependencyFile = "dependencies"
_LibraryFile    = "libraries"

# . Libraries.
# . May need to add "irc" to libraries for Intel C compiler.
# . Default library data.
_SystemLibraries    = [ "m" ]
_SystemLibraryPaths = [ "/usr/local/lib" , "/usr/lib" ]

# . Atlas library data.
# . This will depend on compiler and operating system so it needs to be generalized.
_AtlasCLibrariesToOmit  = [ "dcblas", "df2c", "df2cblas", "df2cdlamch", "df2clapack" ]
_AtlasLibraryPaths      = [ "/usr/local/lib/gcc/x86_64-linux-gnu/4.6" , "/usr/local/lib" , "/usr/lib" ]
_AtlasSerialLibraries   = [ "lapack"   , "f77blas"   , "cblas"   , "atlas" , "m" , "gfortran" ]
_AtlasThreadedLibraries = [ "ptlapack" , "ptf77blas" , "ptcblas" , "atlas" , "m" , "gfortran" ]

# . Compiler flags and options.
_FullOptimizationFlag   = "-O3"
_NoOptimizationFlag     = "-O0"
_NumberOfThreadsFlag    = "-DMAXIMUMNUMBEROFTHREADS={:d}"
_OpenMPFlag             = "-fopenmp"
_OpenMPPreprocessorFlag = "-DUSEOPENMP"
_OptimizationFlagPrefix = "-O"

#===================================================================================================================================
# . Class to allow compilation of C-libraries with additional compiler options.
# . Originally designed for removing all optimization (necessary for dlamch in clapack!).
#===================================================================================================================================
# . Specific imports.
from distutils                    import log
from distutils.command.build_clib import build_clib
from distutils.errors             import DistutilsSetupError
from types                        import ListType, TupleType

# . Class definition.
class BuildClibWithCompileOptions ( build_clib ):
    """Subclass of distutils.command.build_clib which enables compiler options to be played with."""

    def build_libraries ( self, libraries ):
        """Build libraries."""
        # . Loop over libraries.
        for ( lib_name, build_info ) in libraries:
            # . Get sources.
            sources = build_info.get ( 'sources' )
            if ( sources is None ) or ( type ( sources ) not in ( ListType, TupleType ) ):
                raise DistutilsSetupError ( "in 'libraries' option (library '{:s}') 'sources' must be present and must be a list of source filenames".format ( lib_name ) )
            sources = list ( sources )
            # . Remove optimization strings from compiler options if requested.
            if build_info.get ( "removeOptimizationFlags", False ):
                new = []
                for token in self.compiler.compiler_so:
                    if not token.startswith ( _OptimizationFlagPrefix ): new.append ( token )
                self.compiler.compiler_so = new
            # . Do the work.
            log.info ( "building '{:s}' library".format ( lib_name ) )
            # . Compile the source code to object files in the library directory.
            extra_preargs = build_info.get ( "extra_preargs" )
            include_dirs  = build_info.get ( "include_dirs"  )
            macros        = build_info.get ( "macros"        )
            objects       = self.compiler.compile ( sources, output_dir = self.build_temp, macros = macros, include_dirs = include_dirs, debug = self.debug, extra_preargs = extra_preargs )
            #  Create a static library.
            self.compiler.create_static_lib ( objects, lib_name, output_dir = self.build_clib, debug = self.debug )

#===================================================================================================================================
# . Installation options class.
#===================================================================================================================================
class InstallationOptions ( object ):
    """Class for determining installation options."""

    def __init__ ( self ):
        """Constructor."""
        pass

    @classmethod
    def FromCommandLine ( selfClass ):
        """Constructor by getting arguments from a command line."""
        # . Initialization.
        self = selfClass ( )
        # . Parse the command line.
        from argparse import ArgumentParser
        parser = ArgumentParser ( epilog = "The package names are optional as all packages are installed by default." )
        parser.add_argument (        "--atlas"        , action = "store_true"  , dest = "useSerialAtlas"   , default = False , help = "link with the serial ATLAS libraries"           )
        parser.add_argument ( "-c" , "--cLibraries"   , action = "store_true"  , dest = "doCLibraries"     , default = False , help = "install C libraries                  [default]" )
        parser.add_argument ( "-e" , "--extensions"   , action = "store_true"  , dest = "doExtensions"     , default = False , help = "install Cython/Pyrex extensions      [default]" )
        parser.add_argument ( "-f" , "--full"         , action = "store_true"  , dest = "doFull"           , default = False , help = "do a full installation"                         )
        parser.add_argument (        "--noClearUp"    , action = "store_false" , dest = "doClearUp"        , default = True  , help = "do not clear up after installation"             )
        parser.add_argument (        "--noThirdParty" , action = "store_false" , dest = "doThirdParty"     , default = True  , help = "do not install thirdparty libraries"            )
        parser.add_argument (        "--openMP"       , action = "store_true"  , dest = "useOpenMP"        , default = False , help = "compile with OpenMP"                            )
        parser.add_argument (        "packageName"    , nargs  = "*"           ,                                               help = "the name of the package to install"             )
        parser.add_argument (        "--ptAtlas"      , action = "store_true"  , dest = "useThreadedAtlas" , default = False , help = "link with the threaded ATLAS libraries"         )
        parser.add_argument ( "-p" , "--pyrex"        , action = "store_true"  , dest = "doPyrex"          , default = False , help = "convert Cython/Pyrex to C"                      )
        parser.add_argument ( "-q" , "--quiet"        , action = "store_false" , dest = "verbose"          ,                   help = "quiet mode"                                     )
        parser.add_argument ( "-s" , "--shell"        , action = "store_true"  , dest = "doShellFiles"     , default = False , help = "write shell files                    [default]" )
        parser.add_argument (        "--threads"      , type   = int           , dest = "threads"          , default = -1    , help = "the number of threads"                          )
        parser.add_argument ( "-v" , "--verbose"      , action = "store_true"  , dest = "verbose"          , default = True  , help = "verbose mode                         [default]" )
        arguments = parser.parse_args ( )
        # . Check what is to be done.
        # . A full installation overrides everything else.
        if arguments.doFull:
            doCLibraries = doExtensions = doPyrex = doShellFiles = True
        # . Get flags.
        else:
            doCLibraries = arguments.doCLibraries
            doExtensions = arguments.doExtensions
            doPyrex      = arguments.doPyrex
            doShellFiles = arguments.doShellFiles
            # . If all are False do a standard installation.
            if ( not doCLibraries ) and ( not doExtensions ) and ( not doPyrex ) and ( not doShellFiles ): doCLibraries = doExtensions = doShellFiles = True
        # . Get the maximum number of threads for the current processor.
        maximumThreads = selfClass.MaximumNumberOfThreads ( )
        # . Remaining arguments.
        # . Threaded Atlas overrides serial Atlas on threaded machines.
        doClearUp        = arguments.doClearUp
        doThirdParty     = arguments.doThirdParty
        numberOfThreads  = 0
        packageNames     = arguments.packageName
        useOpenMP        = arguments.useOpenMP        and ( maximumThreads > 1 )
        useThreadedAtlas = arguments.useThreadedAtlas and ( maximumThreads > 1 )
        useSerialAtlas   = arguments.useSerialAtlas and ( not useThreadedAtlas )
        verbose          = arguments.verbose
        # . Process the number of threads when OpenMP is to be used.
        # . Is this needed?
        if useOpenMP:
            if maximumThreads > 2: numberOfThreads = max ( min ( arguments.threads, maximumThreads ), 2 )
            else:                  numberOfThreads = 2
        # . Set all attributes.
        self.doClearUp        = doClearUp
        self.doCLibraries     = doCLibraries
        self.doExtensions     = doExtensions
        self.doPyrex          = doPyrex
        self.doShellFiles     = doShellFiles
        self.doThirdParty     = doThirdParty
        self.numberOfThreads  = numberOfThreads
        self.packageNames     = packageNames
        self.useOpenMP        = useOpenMP
        self.useSerialAtlas   = useSerialAtlas
        self.useThreadedAtlas = useThreadedAtlas
        self.verbose          = verbose
        # . Finish up.
        return self

    def IsCLibraryToBeCompiled ( self, libraryName, isThirdParty ):
        """Check whether a C library is to be compiled."""
        doCompilation = ( not isThirdParty ) or ( isThirdParty and self.doThirdParty )
        if doCompilation and ( self.useSerialAtlas or self.useThreadedAtlas ):
            doCompilation = ( libraryName not in _AtlasCLibrariesToOmit )
        return doCompilation

    def IsCLibraryToBeLinked ( self, libraryName ):
        """Check whether a C library is to be linked."""
        if ( self.useSerialAtlas or self.useThreadedAtlas ):
            doLink = ( libraryName not in _AtlasCLibrariesToOmit )
        else:
            doLink = True
        return doLink

    @staticmethod
    def MaximumNumberOfThreads ( ):
        """Find the maximum number of threads."""
        try:
            import multiprocessing
            return multiprocessing.cpu_count ( )
        except:
            return 1

    # . Properties.
    @property
    def compilerFlags ( self ):
        item = self.__dict__.get ( "_compilerFlags", None )
        if item is None:
            if self.useOpenMP: item = [ _NumberOfThreadsFlag.format ( self.numberOfThreads ), _OpenMPPreprocessorFlag, _OpenMPFlag ]
            else: item = []
            setattr ( self, "_compilerFlags", item )
        return item
    @property
    def linkerFlags ( self ):
        item = self.__dict__.get ( "_linkerFlags", None )
        if item is None:
            if self.useOpenMP: item = [ _OpenMPFlag ]
            else: item = []
            setattr ( self, "_linkerFlags", item )
        return item
    @property
    def systemLibraries ( self ):
        item = self.__dict__.get ( "_systemLibraries", None )
        if item is None:
            if   self.useSerialAtlas   : item = _AtlasSerialLibraries
            elif self.useThreadedAtlas : item = _AtlasThreadedLibraries
            else                       : item = _SystemLibraries 
            setattr ( self, "_systemLibraries", item )
        return item
    @property
    def systemLibraryPaths ( self ):
        item = self.__dict__.get ( "_systemLibraryPaths", None )
        if item is None:
            if   self.useSerialAtlas   : item = _AtlasLibraryPaths
            elif self.useThreadedAtlas : item = _AtlasLibraryPaths
            else                       : item = _SystemLibraryPaths 
            setattr ( self, "_systemLibraryPaths", item )
        return item

#===================================================================================================================================
# . Directory install class.
#===================================================================================================================================
class PackageToInstall ( object ):
    """Class for a package that is to be installed."""

    def __init__ ( self, **kwargs ):
        """Constructor."""
        for ( key, value ) in kwargs.iteritems ( ): setattr ( self, key, value )

    @staticmethod
    def Compare ( self, other ):
        """Comparison function."""
        # . Cross-dependence.
        otherIsInSelf = ( other.name in self.dependencyNames )
        selfIsInOther = ( self.name in other.dependencyNames )
        if otherIsInSelf and selfIsInOther: raise ValueError ( "The directories " + self.name + " and " + other.name + " are mutually dependent." )
        elif otherIsInSelf: return  1
        elif selfIsInOther: return -1
        else:
            # . Number of dependencies.
            nself  = len (  self.dependencyNames )
            nother = len ( other.dependencyNames )
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
        target = os.path.join ( self.extensionsPath, _LibraryFile )
        if os.path.exists ( target ):
            tFile = open ( target, "r" )
            for line in tFile:
                tokens = line.split ( )
                if ( len ( tokens ) > 0 ) and ( not tokens[0].startswith ( "#" ) ):
                    directoryName    = tokens[0]
                    libraryName      = directoryName
                    isThirdParty     = False
                    optimizationFlag = None
                    if len ( tokens ) > 1:
                        libraryName      = tokens[1]
                        if len ( tokens ) > 2:
                            isThirdParty = ( "ThirdParty" in tokens[2:] )
                            if   "FullOptimization" in tokens[2:]: optimizationFlag = _FullOptimizationFlag
                            elif "NoOptimization"   in tokens[2:]: optimizationFlag = _NoOptimizationFlag
                    items.append ( ( directoryName, libraryName, isThirdParty, optimizationFlag ) )
            tFile.close ( )
        setattr ( self, "cLibraries", items )

    def CheckForDependencies ( self ):
        """Check for package dependencies."""
        items  = []
        target = os.path.join ( self.path, _DependencyFile )
        if os.path.exists ( target ):
            tFile = open ( target, "r" )
            for line in tFile:
                line = line.strip ( )
                if len ( line ) > 0: items.append ( line )
            tFile.close ( )
        setattr ( self, "dependencyNames", items )

    def CheckForExtensions ( self ):
        """Check for package extensions."""
        # . Extension directory.
        target = os.path.join ( self.path, _ExtensionsDirectory )
        if os.path.exists ( target ):
            setattr ( self, "extensionsPath", target )
            # . Check for a C include directory.
            target = os.path.join ( self.extensionsPath, _CIncludeDirectory )
            if os.path.exists ( target ): setattr ( self, "cIncludePath", target )
            # . Check for libraries.
            self.CheckForCLibraries ( )
            # . Check for a C library directory and make it if necessary.
            if len ( self.cLibraries ) > 0:
                target = os.path.join ( self.extensionsPath, _CLibraryDirectory )
                if not os.path.exists ( target ): os.mkdir ( target )
                setattr ( self, "cLibraryPath" , target )
                setattr ( self, "hasCLibraries", True   )
            # . Check for a Pyrex directory with files.
            target = os.path.join ( self.extensionsPath, _PyrexDirectory )
            if os.path.exists ( target ):
                pyrexFiles = glob.glob ( os.path.join ( target, self.name + ".*" + _PyrexExtension ) )
                if len ( pyrexFiles ) > 0:
                    setattr ( self, "pyrexPath", target )
                    # . Process the file names.
                    pyrexRoots = []
                    for pyrexFile in pyrexFiles:
                        ( head, tailext ) = os.path.split    ( pyrexFile )
                        ( tail, ext     ) = os.path.splitext ( tailext   )
                        pyrexRoots.append ( tail )
                    pyrexRoots.sort ( )
                    setattr ( self, "pyrexFileRoots", pyrexRoots )
                    setattr ( self, "hasPyrexFiles",  True       )
                    # . Check for a psource directory.
                    target = os.path.join ( self.extensionsPath, _PSourceDirectory )
                    if not os.path.exists ( target ): os.mkdir ( target )
                    setattr ( self, "pSourcePath", target )

    def CheckPyrexFiles ( self ):
        """Check that the Pyrex files have been converted."""
        if getattr ( self, "hasPyrexFiles", False ):
            errors = []
            for name in self.pyrexFileRoots:
                pPath = os.path.join ( self.pyrexPath,   name + _PyrexExtension )
                cPath = os.path.join ( self.pSourcePath, name + _CExtension     )
                if os.path.exists ( pPath ):
# . Could check for modification time here but no good if just unpacked distribution as these are not kept.
#                    isOK = os.path.exists ( cPath ) and ( os.path.getmtime ( cPath ) > os.path.getmtime ( pPath ) )
                    if not os.path.exists ( cPath ): errors.append ( name + _PyrexExtension )
            if len ( errors ) > 0:
                if self.options.verbose:
                    print ( "\nInvalid Pyrex-derived C files:\n" )
                    for name in errors: print ( name )
                    print ( "" )
                raise ValueError ( "There are missing or out-of-date Pyrex-derived C files." )

    def CompileCLibraries ( self ):
        """Compile the C-libraries."""
        if getattr ( self, "hasCLibraries", False ):
            # . Get the build directory name.
            buildPath = os.path.join ( "build", "temp" + ".{:s}-{:s}".format ( get_platform ( ), sys.version[0:3] ) )
            # . Get the include directories.
            includeDirectories = [ self.cIncludePath ]
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                path = getattr ( item, "cIncludePath", None )
                if path is not None: includeDirectories.append ( path )
            # . Loop over libraries.
            for ( directoryName, libraryName, isThirdParty, optimizationFlag ) in self.cLibraries:
                # . Check to see whether this library is to be compiled.
                if self.options.IsCLibraryToBeCompiled ( libraryName, isThirdParty ):
                    # . Get the source file list.
                    sourceFiles = glob.glob ( os.path.join ( self.extensionsPath, directoryName, "*" + _CExtension ) )
                    if len ( sourceFiles ) > 0:
                        # . Make build_info.
                        build_info = { }
                        build_info["sources"]       = sourceFiles
                        build_info["extra_preargs"] = list ( self.options.compilerFlags )
                        build_info["include_dirs"]  = includeDirectories
                        build_info["macros"]        = None
                        if optimizationFlag is not None:
                            build_info["removeOptimizationFlags"] = True
                            build_info["extra_preargs"          ].append ( optimizationFlag )
                        # . Compile the library.
                        setup ( name        = libraryName ,
                                libraries   = [ ( libraryName, build_info )   ] ,
                                script_args = [ "BuildClibWithCompileOptions" ] ,
                                cmdclass    = { "BuildClibWithCompileOptions" : BuildClibWithCompileOptions } )
                        self.report["C Files"    ] += len ( sourceFiles )
                        self.report["C Libraries"] += 1
                        # . Move the library to the appropriate place.
                        os.rename ( os.path.join ( buildPath, "lib" + libraryName + _LibraryExtension ), os.path.join ( self.cLibraryPath, "lib" + libraryName + _LibraryExtension ) )

    def ConvertPyrexToC ( self ):
        """Convert Pyrex to C files."""
         # . Check for Pyrex files.
        if getattr ( self, "hasPyrexFiles", False ):
            # . Import the correct modules.
            tag = "Cython"
            try:
                import Cython.Compiler.Errors as Errors, Cython.Compiler.Main as Main, Cython.Compiler.Version as Version
            except:
                raise ValueError ( "Error locating the " + tag + " compiler." )
            # . Compile.
            # . Get the pxdDirectories (self and then dependencies in reverse order).
            pxdDirectories = [ self.pyrexPath ]
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                path = getattr ( item, "pyrexPath", None )
                if path is not None: pxdDirectories.append ( path )
            # . Convert the files.
            if self.options.verbose: print ( "\nConverting files in " + self.pyrexPath + " with " + tag + " version " + Version.version + ":\n" )
            for name in self.pyrexFileRoots:
                if self.options.verbose: print ( " -> " + name + _PyrexExtension )
                PyrexCompile ( Errors, Main, os.path.join ( self.pyrexPath, name + _PyrexExtension ), pxdDirectories = pxdDirectories )
                os.rename ( os.path.join ( self.pyrexPath, name + _CExtension ), os.path.join ( self.pSourcePath, name + _CExtension ) )
                self.report["Pyrex Files"] += 1
            if self.options.verbose: print ( "" )

    @classmethod
    def FromPathName ( selfclass, path ):
        """Constructor from path name."""
        ( head, tail ) = os.path.split ( path )
        name           = tail.rsplit ( "-", 1 )[0]
        self = selfclass ( name = name, path = path )
        return self

    def GetCLibraryNames ( self, options ):
        """Return a list of C library names in reverse order."""
        names = []
        for data in reversed ( getattr ( self, "cLibraries", [] ) ):
            name = data[1]
            if options.IsCLibraryToBeLinked ( name ): names.append ( name )
        return names

    def HasExtensions ( self ):
        """Does this package have extensions?"""
        return ( getattr ( self, "hasCLibraries", False ) or getattr ( self, "hasPyrexFiles", False ) )

    def MakeExtensions ( self ):
        """Make extensions."""
        if getattr ( self, "hasPyrexFiles", False ):
            # . Get the build and destination path names.
            buildPath       = os.path.join ( "build", "lib" + ".{:s}-{:s}".format ( get_platform ( ), sys.version[0:3] ), self.name )
            destinationPath = os.path.join ( self.path, self.name )
            # . Get the include directories.
            includeDirectories = [ self.cIncludePath ]
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                path = getattr ( item, "cIncludePath", None )
                if path is not None: includeDirectories.append ( path )
            # . Get the library directories.
            libraryDirectories = [ self.cLibraryPath ]
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                path = getattr ( item, "cLibraryPath", None )
                if path is not None: libraryDirectories.append ( path )
            libraryDirectories.extend ( self.options.systemLibraryPaths )
            # . Get the library names.
            cLibraries = self.GetCLibraryNames ( self.options )
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                cLibraries.extend ( item.GetCLibraryNames ( self.options ) )
            cLibraries.extend ( self.options.systemLibraries )
            # . Compile and link arguments.
            compileArguments = list ( self.options.compilerFlags )
            linkArguments    = list ( self.options.linkerFlags   )
            # . Make the list of extension modules.
            extensions = []
            for name in self.pyrexFileRoots:
                extensions.append ( Extension ( name, [ os.path.join ( self.pSourcePath, name + _CExtension ) ] ,
                                                                     extra_compile_args   = compileArguments   ,
                                                                     extra_link_args      = linkArguments      ,
                                                                     include_dirs         = includeDirectories ,
                                                                     libraries            = cLibraries         ,
                                                                     library_dirs         = libraryDirectories ,
                                                                     runtime_library_dirs = libraryDirectories ) )
            # . Compile the extension modules.
            setup ( name = self.name, ext_modules = extensions, script_args = [ "build_ext" ] )
            self.report["Extension Files"] += len ( extensions )
            # . Move the files to the appropriate place.
            for name in self.pyrexFileRoots:
                tail = name.split ( "." )[-1]
                os.rename ( os.path.join ( buildPath, tail + _SharedObjectExtension ), os.path.join ( destinationPath, tail + _SharedObjectExtension ) )

    def ProcessWithOptions ( self, options, report ):
        """Process the package given a set of installation options."""
        # . Reporting.
        report["Packages"] += 1
        # . Assign options and report.
        self.options = options
        self.report  = report
        # . C libraries.
        if options.doCLibraries: self.CompileCLibraries ( )
        # . Handle pyrex.
        if options.doPyrex: self.ConvertPyrexToC ( )
        # . Create the extension modules.
        if options.doExtensions:
            if not options.doPyrex: self.CheckPyrexFiles ( )
            self.MakeExtensions ( )

    def ResolveDependencies ( self, itemDictionary ):
        """Resolve the package dependencies."""
        names = getattr ( self, "dependencyNames", [] )
        items = []
        for name in names:
            try:    items.append ( itemDictionary[name] )
            except: raise ValueError ( "Unable to resolve dependency \"" + name + "\" for package " + self.name + "." )
        setattr ( self, "dependencyObjects", items )

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def FindRootDirectory ( ):
    """Find or guess a root directory."""
    # . Initialization.
    rootDirectory = None
    # . Guess a value from the current directory.
    ( head, tail ) = os.path.split ( os.getcwd ( ) )
    while len ( tail ) > 0:
        if "dynamo" in tail.lower ( ):
            rootDirectory = os.path.join ( head, tail )
            os.environ[_PDynamoRoot] = rootDirectory
            break
        ( head, tail ) = os.path.split ( head )
    # . Use the environment variable in case of failure.
    if rootDirectory is None:
        rootDirectory = os.getenv ( _PDynamoRoot )
    # . Finish up.
    if rootDirectory is None: raise ValueError ( "Unable to find pDynamo root directory." )
    return rootDirectory
 
def FindSubdirectory ( root, pattern, isFatal ):
    """Find a subdirectory of root matching a particular pattern."""
    paths = glob.glob ( os.path.join ( root, pattern ) )
    if   len ( paths ) == 0:
        if isFatal: raise ValueError ( "No subdirectory found matching " + pattern + "." )
        return None
    elif len ( paths ) == 1:
        ( head, tail ) = os.path.split ( paths[0] )
        return tail
    else:
        raise ValueError ( "Multiple subdirectories found matching " + pattern + "." )

def GetRootDirectories ( names ):
    """Get the root directories in increasing order of dependency.

    All directories are returned by default otherwise only those specified by |names|.
    """
    # . Get the root directory name.
    rootDirectory = os.getenv ( _PDynamoRoot )
    # . Get all directories in root.
    toProcess = {}
    for candidate in glob.glob ( os.path.join ( rootDirectory, "*" ) ):
        if os.path.isdir ( candidate ):
            directory = PackageToInstall.FromPathName ( path = candidate )
            directory.CheckForDependencies ( )
            directory.CheckForExtensions   ( )
            toProcess[directory.name] = directory
    # . Satisfy the dependencies of the directories.
    for directory in toProcess.values ( ):
        directory.ResolveDependencies ( toProcess )
    # . Create the sorted directory list.
    directories = list ( toProcess.values ( ) )
    directories.sort ( key = functools.cmp_to_key ( PackageToInstall.Compare ) )
    # . Prune the list if necessary.
    if ( names is not None ) and ( len ( names ) > 0 ):
        items = directories
        directories = []
        for item in items:
            if item.name in names: directories.append ( item )
    # . Finish up.
    return directories

def InstallationSummary ( report, totalTime, message ):
    """Summarize the results of the installation."""
    print ( "\n---------------------\nInstallation Summary:\n---------------------" )
    if len ( report ) > 0:
        for key in sorted ( report.keys ( ) ):
            print ( "{:<17s} {:5d}".format ( key, report[key] ) )
    print ( "Processing Time   {:s}".format ( TimeToString ( totalTime ) ) )
    print ( "\n{:s}".format ( message ) )

def PyrexCompile ( Errors, Main, source, pxdDirectories = None ):
    """Convert a Pyrex file to C."""
    options = Main.default_options
    # . Change for Pyrex version 0.9.6: options is now a dictionary not a CompilationOptions object so convert to a CompilationOptions object.
    if isinstance ( options, dict ): options = Main.CompilationOptions ( options )
    if pxdDirectories is not None: options.include_path.extend ( pxdDirectories )
    context = Main.Context ( options.include_path, {} )
    try:
        result = Main.compile ( source, options )
        failed = ( result.num_errors > 0 )
    except Errors.PyrexError as e:
        print ( e )
        failed = True
    if failed: raise ValueError ( "There was a Pyrex compiler error." )

def RemoveBuildDirectory ( ):
    """Remove the build directory if it exists."""
    if os.path.exists ( "build" ): shutil.rmtree ( "build" )

def TimeToString ( time ):
    """Convert a floating point time to a string."""
    value  = time
    fields = []
    for ( size, tag ) in ( ( 60*60*24, "d" ), ( 60*60, "h" ), ( 60, "m" ) ):
        ( newf, value ) = divmod ( value, size )
        newi = int ( round ( newf ) )
        if not ( ( newi == 0 ) and ( len ( fields ) == 0 ) ):
            fields.append ( "{:2d}{:1s}".format ( newi, tag ) )
    tag = "s"
    fields.append ( "{:5.3f}{:1s}".format ( value, tag ) )
    return " ".join ( fields )

def WriteShellFile ( inPath, outPath, variables, report ):
    """Write a shell file."""
    inFile   = open ( inPath, "r" )
    inString = inFile.read ( ).format ( **variables )
    inFile.close  ( )
    outFile  = open ( outPath, "w" )
    outFile.write ( inString )
    outFile.close ( )
    report["Shell Files"] += 1

def WriteShellFiles ( report, parameters = None, scratch = None ):
    """Write shell files containing the environment variables needed by pDynamo."""
    # . Get all required paths.
    # . Root.
    FindRootDirectory ( )
    pDynamoRoot = os.getenv ( _PDynamoRoot )
    # . Package directories.
    pBabel             = FindSubdirectory ( pDynamoRoot, "pBabel-*"          , True  )
    pCore              = FindSubdirectory ( pDynamoRoot, "pCore-*"           , True  )
    pGraph             = FindSubdirectory ( pDynamoRoot, "pGraph-*"          , False )
    pMolecule          = FindSubdirectory ( pDynamoRoot, "pMolecule-*"       , True  )
    pMoleculeScripts   = FindSubdirectory ( pDynamoRoot, "pMoleculeScripts-*", True  )
    pSandBox           = FindSubdirectory ( pDynamoRoot, "pSandBox"          , False )
    packageDirectories = [ pBabel, pCore ]
    if pGraph   is not None: packageDirectories.append ( pGraph   )
    packageDirectories.extend ( [ pMolecule, pMoleculeScripts ] )
    if pSandBox is not None: packageDirectories.append ( pSandBox )
    # . Parameters.
    if parameters is None: parameters = "$PDYNAMO_ROOT/parameters"
    # . Scratch.
    if scratch is None: scratch = "$PDYNAMO_ROOT/scratch"
    # . Get the new Python path.
    pythonPath = os.getenv ( "PYTHONPATH" )
    if pythonPath is None: pythonPath = ""
    oldTokens = pythonPath.split ( ":" )
    newTokens = []
    for token in oldTokens:
        found = False
        for label in packageDirectories:
            if token.find ( label ) >= 0:
                found = True
                break
        if not found: newTokens.append ( token )
    newTokens.extend ( [ "$PDYNAMO_ROOT/" + pBabel   , "$PDYNAMO_ROOT/" + pCore            ] )
    if pGraph   is not None: newTokens.append ( "$PDYNAMO_ROOT/" + pGraph   )
    newTokens.extend ( [ "$PDYNAMO_ROOT/" + pMolecule, "$PDYNAMO_ROOT/" + pMoleculeScripts ] )
    if pSandBox is not None: newTokens.append ( "$PDYNAMO_ROOT/" + pSandBox )
    pythonPath = ":".join ( newTokens )
    # . Write the shell files.
    variables = { "parameters"       : parameters       ,
                  "pBabel"           : pBabel           ,
                  "pCore"            : pCore            ,
                  "pDynamoRoot"      : pDynamoRoot      ,
                  "pMolecule"        : pMolecule        ,
                  "pMoleculeScripts" : pMoleculeScripts ,
                  "pythonPath"       : pythonPath       ,
                  "scratch"          : scratch          }
    templatePath = os.path.join ( pDynamoRoot, "installation", "templates" )
    WriteShellFile ( os.path.join ( templatePath, "environment_bash.tmpl"   ), "environment_bash.com",   variables, report )
    WriteShellFile ( os.path.join ( templatePath, "environment_cshell.tmpl" ), "environment_cshell.com", variables, report )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
