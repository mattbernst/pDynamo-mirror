/*------------------------------------------------------------------------------
! . File      : CGLinearEquationSolver.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module implements a preconditioned conjugate-gradient linear equation solver.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "Memory.h"
# include "CGLinearEquationSolver.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Default_ConvergenceMode     1
# define Default_MaximumIterations 500
# define Default_ErrorTolerance    1.0e-10

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Solver. */
static Boolean CGLinearEquationSolver_IsConverged ( const Integer convergenceMode, const Real errorTolerance, const Real rNorm2, const Real hNorm2, const Real bNorm2 ) ;

/*==================================================================================================================================
! . Report.
!=================================================================================================================================*/
void CGLinearEquationSolverReport_Initialize ( CGLinearEquationSolverReport *self )
{
    if ( self != NULL )
    {
        self->finalResidual   = 0.0e+00 ;
        self->initialResidual = 0.0e+00 ;
        self->isConverged     = False   ;
        self->iterations      = 0       ;
        self->rhsNorm2        = 0.0e+00 ;
    }
}

/*==================================================================================================================================
! . Solver.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
CGLinearEquationSolver *CGLinearEquationSolver_Allocate ( void )
{
    CGLinearEquationSolver *self = NULL ;
    self = ( CGLinearEquationSolver * ) Memory_Allocate ( sizeof ( CGLinearEquationSolver ) ) ;
    if ( self != NULL )
    {
        self->convergenceMode   = Default_ConvergenceMode   ;
        self->maximumIterations = Default_MaximumIterations ;
        self->errorTolerance    = Default_ErrorTolerance    ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
CGLinearEquationSolver *CGLinearEquationSolver_Clone ( const CGLinearEquationSolver *self )
{
    CGLinearEquationSolver *new = NULL ;
    if ( self != NULL )
    {
        new = CGLinearEquationSolver_Allocate ( ) ;
        if ( new != NULL )
        {
            new->convergenceMode   = self->convergenceMode   ;
            new->maximumIterations = self->maximumIterations ;
            new->errorTolerance    = self->errorTolerance    ;
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CGLinearEquationSolver_Deallocate ( CGLinearEquationSolver **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for convergence.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean CGLinearEquationSolver_IsConverged ( const Integer convergenceMode, const Real errorTolerance, const Real rNorm2, const Real hNorm2, const Real bNorm2 )
{
    Real stepTest = 0.0e+00 ;
    switch ( convergenceMode )
    {
        case 1: stepTest = rNorm2          ; break ;
        case 2: stepTest = rNorm2 / bNorm2 ; break ;
        case 3: stepTest = hNorm2          ; break ;
        case 4: stepTest = hNorm2 / bNorm2 ; break ;
    }
    return ( stepTest <= errorTolerance ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Solver for SPD matrices with or without preconditioning.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . H can be R if there is no preconditioning. */

void CGLinearEquationSolver_PCGSolver ( const CGLinearEquationSolver *self, CGLinearEquationSolverState *state, CGLinearEquationSolverReport *report )
{
    CGLinearEquationSolverReport_Initialize ( report ) ;
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        auto Boolean      doPreconditioning, isConverged ;
        auto Integer      iterations  ;
        auto Real         alpha, beta, bNorm2, denominator, h0Norm2, hNorm2, oldRH, r0Norm2, rDotH, rNorm2 ;
        auto Real1DArray *b, *h, *r, *x ;
        auto CGLinearEquationSolverTarget *target ;
        auto void                         *object ;

        /* . Aliases. */
        target = state->target ;
        b      = state->b ;
        h      = state->h ;
        r      = state->r ;
        x      = state->x ;
        object = target->object ;

        /* . Initialization. */
        doPreconditioning = ( target->ApplyPreconditioner != NULL ) ;
        bNorm2            = Real1DArray_Norm2 ( b ) ;
        iterations        = 0 ;

        /* . Denominator for convergence checks. */
        denominator = 1.0e+00 ;
             if ( self->convergenceMode == 2 ) denominator = bNorm2 ;
        else if ( self->convergenceMode == 4 )
        {
            if ( doPreconditioning )
            {
                (*target->ApplyPreconditioner) ( object, b, r ) ;
                denominator = Real1DArray_Norm2 ( r ) ;
            }
            else denominator = bNorm2 ;
        }

        /* . Compute initial residual r = b - A*x. */
        (*target->ApplyMatrix)     ( object, x, r ) ;
        Real1DArray_AddScaledArray ( r, -1.0e+00, b, NULL ) ;
        Real1DArray_Scale          ( r, -1.0e+00  ) ;
        r0Norm2 = Real1DArray_Norm2 ( r ) ;
        rNorm2  = r0Norm2 ;

        /* . Preconditioning h = p * r (in b). */
        if ( doPreconditioning )
        {
            (*target->ApplyPreconditioner) ( object, r, b ) ;
            h0Norm2 = Real1DArray_Norm2 ( b ) ;
        }
        else
        {
            Real1DArray_CopyTo  ( r, b, NULL ) ;
            h0Norm2 = r0Norm2 ;
        }

        /* . Initial convergence check. */
        isConverged = CGLinearEquationSolver_IsConverged ( self->convergenceMode, self->errorTolerance, r0Norm2, h0Norm2, denominator );
        if ( ! isConverged )
        {
            rDotH = Real1DArray_Dot ( r, b, NULL ) ;
            for ( iterations = 1 ; iterations <= self->maximumIterations ; iterations++ )
            {
                /* . New x. */
                (*target->ApplyMatrix) ( object, b, h ) ;
                alpha = rDotH / Real1DArray_Dot ( h, b, NULL ) ;
                Real1DArray_AddScaledArray ( x, alpha, b, NULL ) ;

                /* . New r. */
                Real1DArray_AddScaledArray ( r, -alpha, h, NULL ) ;
                rNorm2 = Real1DArray_Norm2 ( r ) ;

                /* . New h. */
                if ( doPreconditioning )
                {
                    (*target->ApplyPreconditioner) ( object, r, h ) ;
                    hNorm2 = Real1DArray_Norm2 ( h ) ;
                }
                else
                {
                    Real1DArray_CopyTo  ( r, h, NULL ) ;
                    hNorm2 = rNorm2 ;
                }

                /* . Check for termination. */
                isConverged = CGLinearEquationSolver_IsConverged ( self->convergenceMode, self->errorTolerance, rNorm2, hNorm2, denominator ) ;
                if ( ( isConverged ) || ( iterations >= self->maximumIterations ) ) break ;

                /* . New p. */
                oldRH = rDotH ;
                rDotH = Real1DArray_Dot ( r, h, NULL ) ;
                beta  = rDotH / oldRH ;
                Real1DArray_Scale ( b, beta ) ;
                Real1DArray_AddScaledArray ( b, 1.0e+00, h, NULL ) ;
            }
        }

        /* . Finish up. */
        if ( report != NULL )
        {
            report->finalResidual   = rNorm2      ;
            report->initialResidual = r0Norm2     ;
            report->isConverged     = isConverged ;
            report->iterations      = iterations  ;
            report->rhsNorm2        = bNorm2      ;
        }
    }
}

/*==================================================================================================================================
! . State.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
CGLinearEquationSolverState *CGLinearEquationSolverState_Allocate ( void )
{
    CGLinearEquationSolverState *self = NULL ;
    self = ( CGLinearEquationSolverState * ) Memory_Allocate ( sizeof ( CGLinearEquationSolverState ) ) ;
    if ( self != NULL )
    {
        self->b      = NULL ;
        self->h      = NULL ;
        self->r      = NULL ;
        self->x      = NULL ;
        self->target = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CGLinearEquationSolverState_Deallocate ( CGLinearEquationSolverState **self )
{
    if ( (*self) != NULL )
    {
        Real1DArray_Deallocate ( &((*self)->h) ) ;
        Real1DArray_Deallocate ( &((*self)->r) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setup a state given a target.
!---------------------------------------------------------------------------------------------------------------------------------*/
CGLinearEquationSolverState *CGLinearEquationSolverState_SetupFromTarget ( CGLinearEquationSolverTarget *target, Status *status )
{
    CGLinearEquationSolverState *self = NULL ;
    if ( target != NULL )
    {
        auto Integer n = 0 ;
        if ( target->rhs != NULL ) n = Real1DArray_Length ( target->rhs ) ;
        if ( n > 0 )
        {
            self = CGLinearEquationSolverState_Allocate ( ) ;
            if ( self != NULL )
            {
                /* . Assignment. */
                self->b      = target->rhs      ;
                self->x      = target->solution ;
                self->target = target           ;
                /* . Allocate space. */
                self->h = Real1DArray_Allocate ( n, status ) ;
                self->r = Real1DArray_Allocate ( n, status ) ;
                if ( ( self->h == NULL ) || ( self->r == NULL ) ) CGLinearEquationSolverState_Deallocate ( &self ) ;
            }
            if ( self == NULL ) Status_Set ( status, Status_MemoryAllocationFailure ) ;
        }
    }
    return self ;
}

/*==================================================================================================================================
! . Target.
!=================================================================================================================================*/
void CGLinearEquationSolverTarget_Initialize ( CGLinearEquationSolverTarget *self )
{
    if ( self != NULL )
    {
        self->rhs                 = NULL ;
        self->solution            = NULL ;
        self->object              = NULL ;
        self->ApplyMatrix         = NULL ;
        self->ApplyPreconditioner = NULL ;
    }
}
