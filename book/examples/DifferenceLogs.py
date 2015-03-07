"""Difference generated logs with those supplied with the distribution (text only)."""

import glob, os.path, subprocess

from Definitions import logExtension, logPath, rootDirectory

# . Get the reference log directory.
logReference = os.path.join ( rootDirectory, "book", "logs" )

# . Get the file sets (tails only).
files1 = glob.glob ( os.path.join ( logPath     , "*" + logExtension ) )
files2 = glob.glob ( os.path.join ( logReference, "*" + logExtension ) )

set1 = set ( )
set2 = set ( )
for file in files1: set1.add ( os.path.split ( file )[-1] )
for file in files2: set2.add ( os.path.split ( file )[-1] )

# . Get common and uncommon files.
commonFiles = list ( set1 & set2 )
singleFiles = list ( set1 ^ set2 )
commonFiles.sort ( )
singleFiles.sort ( )

# . Print the results.
subprocess.call ( "echo " + 80 * "=", shell = True )
if len ( singleFiles ) > 0:
    subprocess.call ( "echo " + 80 * "=", shell = True )
    print ( "Single files:" )
    for file in singleFiles: print ( file )

if len ( commonFiles ) > 0:
    for file in commonFiles:

        file1 = os.path.join ( logPath     , file )
        file2 = os.path.join ( logReference, file )

        subprocess.call ( "echo " + 80 * "=",                 shell = True )
        subprocess.call ( "echo Differencing file = " + file, shell = True )
        subprocess.call ( "diff " + file1 + " " + file2,      shell = True )

subprocess.call ( "echo " + 80 * "=", shell = True )
print ( "Common files = {:2d}".format ( len ( commonFiles ) ) )
print ( "Single files = {:2d}".format ( len ( singleFiles ) ) )
