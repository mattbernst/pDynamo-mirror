------------------------------
Timings of pDynamo Benchmarks:
------------------------------

DHFR is a standard test from the joint AMBER-CHARMM set of benchmarks. The
pDynamo test is not exactly comparable as a different NB model is used (ABFS
instead of PME) and no SHAKE constraints are applied.

CS and PKA are pDynamo-specific tests. CS was provided by Alexey Aleksandrov
and PKA is from the pDynamo PDB file tutorial.

Currently DHFR is used as an MM benchmark and CS and PKA as QC/MM benchmarks.
CS and PKA give very similar results so only CS is enabled by default.

For the Langevin/Velocity Verlet and System/SystemWithTimings options, there is
little difference in timings so all benchmarks are done with Langevin and
SystemWithTimings.

--------------------------------------------------
Installation         Citrate Synthase         DHFR
--------------------------------------------------
Serial                        56m 30s      11m
Serial Atlas                  36m          11m
Threaded Atlas 4              41m 20s      11m
OMP4                          46m 30s       4m 30s
OMP8                          45m 50s       3m 20s
OMP4 + Serial Atlas           25m 40s       4m 30s
OMP8 + Serial Atlas           24m 25s       3m 20s
--------------------------------------------------

The conclusions are fairly clear. MM simulations benefit from OpenMP whereas
QC and QC/MM simulations benefit from both the use of Atlas and OpenMP. The
parallel efficiency with OpenMP is reasonable for the parts of the code that
use it.

-------------------------------------------------
Timings of pDynamo Examples, Tests and Tutorials:
-------------------------------------------------

Jobs run with a standard Linux installation in serial (unless otherwise noted):

Book Examples (all) ~  3h
Package Tests       ~  30m
Tutorials:
    aSmallProtein   ~  30m
    bAlaPhiPsiPMF   ~   1h
    dftQCModels     <  20s
    orca            <  20s
    pdbFiles        ~   6h 30m - serial
                        2h 30m - OMP 4 + serial Atlas
                        2h     - OMP 8 + serial Atlas

*** Full outputs can be found in the log files. ***
