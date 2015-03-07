"""pDynamo installation script."""

import time

from collections      import defaultdict
from InstallUtilities import FindRootDirectory    , \
                             GetRootDirectories   , \
                             InstallationOptions  , \
                             InstallationSummary  , \
                             RemoveBuildDirectory , \
                             WriteShellFiles

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def InstallpDynamo ( ):
    """Main function."""

    # . Initialization.
    report    = defaultdict ( int )
    startTime = time.time ( )

    # . Get the installation options.
    options = InstallationOptions.FromCommandLine ( )

    # . Processing.
    try:

        # . Find the root directory.
        FindRootDirectory ( )

        # . Check for package processing.
        if options.doCLibraries or options.doExtensions or options.doPyrex:

            # . Get the packages to install in the correct order.
            packages = GetRootDirectories ( options.packageNames )

            # . Process all packages with extensions.
            for package in packages:
                if package.HasExtensions ( ):
                    package.ProcessWithOptions ( options, report )

            # . Clear up.
            if options.doClearUp: RemoveBuildDirectory ( )

        # . Shell.
        if options.doShellFiles: WriteShellFiles ( report )

        # . Terminating message.
        message = "Installation terminated successfully."

    except ValueError as e:
        if hasattr ( e, "args" ) and ( len ( e.args ) > 0 ): message = "Exiting: " + e.args[0]
        else:                                                message = "Exiting with error."

    # . Finish up.
    totalTime = time.time ( ) - startTime
    InstallationSummary ( report, totalTime, message )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    # . Run the script.
    InstallpDynamo ( )
