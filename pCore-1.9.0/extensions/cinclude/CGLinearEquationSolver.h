/*------------------------------------------------------------------------------
! . File      : CGLinearEquationSolver.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _CGLINEAREQUATIONSOLVER
# define _CGLINEAREQUATIONSOLVER

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The report type. */
typedef struct {
    Boolean isConverged     ;
    Integer iterations      ;
    Real    finalResidual   ;
    Real    initialResidual ;
    Real    rhsNorm2        ;
} CGLinearEquationSolverReport ;

/* . The target type. */
typedef struct {
    Real1DArray *rhs      ;
    Real1DArray *solution ;
    void        *object   ;
    void ( * ApplyMatrix         ) ( void *target, Real1DArray *x, Real1DArray *y ) ;
    void ( * ApplyPreconditioner ) ( void *target, Real1DArray *x, Real1DArray *y ) ;
} CGLinearEquationSolverTarget ;

/* . The state type. */
typedef struct {
    Real1DArray *b ;
    Real1DArray *h ;
    Real1DArray *r ;
    Real1DArray *x ;
    CGLinearEquationSolverTarget *target ;
} CGLinearEquationSolverState ;

/* . The solver type. */
typedef struct {
    Integer convergenceMode   ;
    Integer maximumIterations ;
    Real    errorTolerance    ;
} CGLinearEquationSolver ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Report. */
extern void CGLinearEquationSolverReport_Initialize ( CGLinearEquationSolverReport *self ) ;

/* . Solver. */
extern CGLinearEquationSolver      *CGLinearEquationSolver_Allocate   ( void ) ;
extern CGLinearEquationSolver      *CGLinearEquationSolver_Clone      ( const CGLinearEquationSolver  *self ) ;
extern void                         CGLinearEquationSolver_Deallocate (       CGLinearEquationSolver **self ) ;
extern void                         CGLinearEquationSolver_PCGSolver  ( const CGLinearEquationSolver  *self, CGLinearEquationSolverState *state, CGLinearEquationSolverReport *report ) ;

/* . State. */
extern CGLinearEquationSolverState *CGLinearEquationSolverState_Allocate        ( void ) ;
extern void                         CGLinearEquationSolverState_Deallocate      (       CGLinearEquationSolverState **self ) ;
extern CGLinearEquationSolverState *CGLinearEquationSolverState_SetupFromTarget ( CGLinearEquationSolverTarget *target, Status *status ) ;

/* . Target. */
extern void CGLinearEquationSolverTarget_Initialize ( CGLinearEquationSolverTarget *self ) ;

# endif
