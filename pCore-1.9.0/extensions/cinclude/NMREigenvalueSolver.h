/*------------------------------------------------------------------------------
! . File      : NMREigenvalueSolver.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _NMREIGENVALUESOLVER
# define _NMREIGENVALUESOLVER

# include "Boolean.h"
# include "Integer.h"
# include "Macros.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*# define USENMREIGENVALUESOLVER*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void NMREigenvalueSolver_CheckSolution (       SymmetricMatrix *self                 ,
                                                const Integer          numberOfEigenvalues  ,
                                                      Real1DArray     *eigenvalues          ,
                                                      Real1DArray     *referenceEigenvalues ,
                                                      Real2DArray     *eigenvectors         ,
                                                      Real            *eigenvalueError      ,
                                                      Real            *eigenvectorError     ,
                                                      Real            *normalizationError   ) ;
extern void NMREigenvalueSolver_FullCheck ( SymmetricMatrix *self, Real1DArray *eigenvalues, Real2DArray *eigenvectors, Status *status ) ;
extern void NMREigenvalueSolver_Solve     ( SymmetricMatrix *self, Real1DArray *eigenvalues, Real2DArray *eigenvectors, Status *status ) ;

# endif
