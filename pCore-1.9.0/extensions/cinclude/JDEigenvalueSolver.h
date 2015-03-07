/*------------------------------------------------------------------------------
! . File      : JDEigenvalueSolver.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _JDEIGENVALUESOLVER
# define _JDEIGENVALUESOLVER

# include "cprimme.h"
# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The report type. */
typedef struct {
    Boolean isConverged        ;
    Boolean solutionChecked    ;
    Integer convergedPairs     ;
    Integer numberMatrixVectorMultiplications  ;
    Integer returnCode         ;
    Real    eigenvalueError    ;
    Real    eigenvectorError   ;
    Real    normalizationError ;
} JDEigenvalueSolverReport ;

/* . The target type. */
typedef struct {
    Real1DArray *eigenvalues  ;
    Real2DArray *eigenvectors ;
    void        *object       ;
    void ( * ApplyMatrix         ) ( void *x, void *y, int *blockSize, struct primme_params *primme ) ;
    void ( * ApplyPreconditioner ) ( void *x, void *y, int *blockSize, struct primme_params *primme ) ;
} JDEigenvalueSolverTarget ;

/* . The state type. */
typedef struct {
    Integer1DArray           *iWork         ;
    Real1DArray              *eigenvalues   ;
    Real1DArray              *residualNorms ;
    Real1DArray              *rWork         ;
    Real2DArray              *eigenvectors  ;
    JDEigenvalueSolverTarget *target        ;
    struct primme_params      primme        ;
} JDEigenvalueSolverState ;

/* . The solver type. */
typedef struct {
    Boolean usePreconditioning ;
    Integer maximumMatrixVectorMultiplications ;
    Integer printLevel     ;
    Real    errorTolerance ;
} JDEigenvalueSolver ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Report. */
void JDEigenvalueSolverReport_Initialize ( JDEigenvalueSolverReport *self ) ;

/* . Solver. */
extern JDEigenvalueSolver      *JDEigenvalueSolver_Allocate      ( void ) ;
extern void                     JDEigenvalueSolver_CheckSolution ( const JDEigenvalueSolver  *self, JDEigenvalueSolverState *state, Real1DArray *referenceEigenvalues, JDEigenvalueSolverReport *report ) ;
extern JDEigenvalueSolver      *JDEigenvalueSolver_Clone         ( const JDEigenvalueSolver  *self ) ;
extern void                     JDEigenvalueSolver_CopyTo        ( const JDEigenvalueSolver  *self, JDEigenvalueSolver *other ) ;
extern void                     JDEigenvalueSolver_Deallocate    (       JDEigenvalueSolver **self ) ;
extern void                     JDEigenvalueSolver_Initialize    (       JDEigenvalueSolver  *self ) ;
extern void                     JDEigenvalueSolver_Solve         ( const JDEigenvalueSolver  *self, JDEigenvalueSolverState *state, JDEigenvalueSolverReport *report ) ;

/* . State. */
extern JDEigenvalueSolverState *JDEigenvalueSolverState_Allocate        ( void ) ;
extern void                     JDEigenvalueSolverState_Deallocate      ( JDEigenvalueSolverState **self ) ;
extern JDEigenvalueSolverState *JDEigenvalueSolverState_SetupFromTarget ( JDEigenvalueSolverTarget *target, Status *status ) ;

/* . Target. */
void JDEigenvalueSolverTarget_Initialize ( JDEigenvalueSolverTarget *self ) ;

# endif
