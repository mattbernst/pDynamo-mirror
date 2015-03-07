# pDynamo-mirror
*This is a mirror of the pDynamo computational chemistry software, which does not currently use a public version control repository. It is provided to improve visibility of this package and its components. I, Matt Ernst, have no relationship to the original developers.*

This software is governed by the [CeCILL license](http://www.cecill.info/index.en.html)

https://sites.google.com/site/pdynamomodeling/

pDynamo is an open source program library that has been designed for the simulation of molecular systems using quantum chemical (QC), molecular mechanical (MM) and hybrid QC/MM potential energy functions. pDynamo is written in a mixture of Python and C and complements the existing version of Dynamo, now denoted fDynamo, that is written in Fortrans 90 and 95.

The current release version of pDynamo has the following capabilities:

* A density functional theory QC method employing Gaussian basis sets.
* Semi-empirical QC methods of the MNDO type, including AM1, MNDO, PDDG, PM3 and RM1.
* Support for some standard MM force fields, including AMBER, CHARMM and OPLS-AA.
* Hybrid QC/MM methods using any combination of the QC and MM potentials implemented in the library.
* Coupling to third-party QC programs.
* Energy calculations.
* Geometry optimizations.
* Transition state searches.
* Reaction path calculations.
* Normal mode analyses.
* Property calculations, such as charges and dipoles.
* Molecular dynamics simulations.
* Monte Carlo simulations.
* Various "soft" geometrical constraints.
* The ability to handle various common molecular file formats.
* Miscellaneous analysis tools.
