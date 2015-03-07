-------------------------
Installation Instructions
-------------------------

------------
Requirements
------------

Installation has been tested on Linux and MAC OS-X machines. Details about
installation on Windows machines may be found on the web site.

It is necessary to have Python 2.7 and a C compiler.

------------
Installation
------------

If there are older verions of pDynamo installed on your machine, the environment
variables and Python path definitions corresponding to these versions should be
removed. This can be done by undefining the appropriate variables or by
overwriting them with the definitions appropriate to the new version (see the
next section on how to do this).

To install, unpack the distribution, go to the installation subdirectory and
then type:

python Install.py

Installation generates a lot of output and could take a few minutes to complete
depending on the machine.

Installation is performed in the pDynamo distribution directories. The Python
distribution itself is NOT modified. This is because the current pDynamo
distribution should be regarded as beta (not fully stable) and could change
substantially in the next few months.

By default, installation does NOT process the Cython files that link the Python
and C parts of pDynamo together. Instead, it uses the preprocessed files that
are included in the distribution. This may cause problems on some machines. If
this is the case, install the current version of Cython and repeat the 
installation using the command:

python Install.py -f

A version of Cython is included in the distribution in the thirdParty subdirectory.
Alternatively, the Cython website is:

"http://www.cython.org".

--------------------
Third Party Software
--------------------

The pDynamo distribution includes some third party software in the thirdParty
subdirectory, including the Python packages Cython and PyYAML. PyYAML is essential
for the proper working of pDynamo and should be installed, whereas Cython is useful
for the installation of pDynamo on certain machines and with certain operating
systems. Both packages can be installed by going to the relevant directory
and typing "python setup.py install". This command may require superuser
privileges.

---------------
Running pDynamo
---------------

To run pDynamo, it is necessary to define a few environment variables and to add
directories to the Python path.

To aid users, the installation script writes files containing these definitions
for certain Unix shells. These may be found in the installation directory, after
installation, with the names "environment_shell.com". Users should choose the
file that is appropriate for the shell that they employ, modify it (if required)
and then add it to their shell configuration file (e.g. ".bash_profile" or
".cshrc").

The environment variables that need to be defined are:

PDYNAMO_ROOT       should point to the distribution directory.

PDYNAMO_PARAMETERS should be given the value "$PDYNAMO_ROOT/parameters".

PDYNAMO_SCRATCH    should point to a scratch directory that can be used by the
                   examples provided with the distribution.

PDYNAMO_STYLE      should point to the location of the CSS style file that is
                   used when generating pDynamo output in XHTML format. A
                   basic stylefile is given in
                   "$PDYNAMO_ROOT/parameters/cssStyleSheets/defaultStyle.css".

The directories that need to be added to the PYTHONPATH are:

$PDYNAMO_ROOT/pBabel-x.y
$PDYNAMO_ROOT/pCore-x.y
$PDYNAMO_ROOT/pMolecule-x.y
$PDYNAMO_ROOT/pMoleculeScripts-x.y

-------
Testing
-------

To test the distribution using the examples from the book, go to the examples
directory in the book directory and type:

python RunExamples.py

This will run a selection of the shorter examples. To run the full suite of
examples type:

python RunExamples.py --all

All examples take several hours to run and requires approximately one gigabyte
of disk space. Please note that the order in which some of the examples are run
is important. Thus, for example, Example26Q should be run before Example26LJ as
the former generates a file that is required by the latter.

Differences between the generated outputs and the reference outputs provided in
the distribution can be compared by typing:

python DifferenceOutputs.py

There are likely to be many numerical differences but these should be quite
small.

The sets of package tests can be run from the installation directory (once the
environment variables are set up) by typing:

python RunPackageTests.py
