/*------------------------------------------------------------------------------
! . File      : ExecutionEnvironment.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _EXECUTIONENVIRONMENT
# define _EXECUTIONENVIRONMENT

# include <stdlib.h>

/*
! . OpenMP:
!
!   The number of threads is 1 if the pragma if clause is not satisfied.
!   Otherwise, it is determined in order of priority via:
!   pragma num_threads clause
!   omp_set_num_threads library function
!   OMP_NUM_THREADS environment variable
!   implementation default - usually the number of CPUs on a node
!
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Options.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Thread options. */
/*# define USEOPENMP*/
/*# define USEPTHREAD*/

/*
# if defined ( USEOPENMP ) || defined ( USEPTHREAD )
# define MAXIMUMNUMBEROFTHREADS 4
# endif
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Includes.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifdef USEOPENMP
# include "omp.h"
# endif

# ifdef USEPTHREAD
# include <pthread.h>
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Thread types. */
/*
typedef enum { ThreadType_OpenMP = 0, ThreadType_PThread = 1, ThreadType_Serial = 2 } ThreadType ;
*/

/* . The execution environment type. */
/*
# if defined ( USEOPENMP )
typedef struct {
          Integer    numberOfThreads ;
    const ThreadType threadType      ;
} ExecutionEnvironment = { MAXIMUMNUMBEROFTHREADS, ThreadType_OpenMP } ;
# elif defined ( USEPTHREAD )
typedef struct {
          Boolean    includeMainThreadInThreadCount ;
          Integer    numberOfThreads ;
    const ThreadType threadType      ;
} ExecutionEnvironment = { MAXIMUMNUMBEROFTHREADS, ThreadType_PThread } ;
# else
typedef struct {
          Integer    numberOfThreads ;
    const ThreadType threadType      ;
} ExecutionEnvironment = { 0, ThreadType_Serial } ;
# endif
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/

# endif
