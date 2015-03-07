#-------------------------------------------------------------------------------
# . File      : ParallelizationOptions.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Basic parallelization options."""

ParallelizationOptions = set ( [ "Serial" ] )

# . Determine the various parallelization options.
try:
    import concurrent.futures
    ParallelizationOptions.add ( "Futures" )
    HasFutures = True
except:
    HasFutures = False

try:
    from mpi4py import MPI
    ParallelizationOptions.add ( "MPI2" )
    HasMPI2 = True
except:
    HasMPI2 = False
    
try:
    import multiprocessing
    ParallelizationOptions.add ( "Multiprocessing" )
    HasMultiprocessing = True
except:
    HasMultiprocessing = False

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    print ( "{!r}".format ( ParallelizationOptions ) )
