/*------------------------------------------------------------------------------
! . File      : LAPACKEigenvalueSolver.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _LAPACKEIGENVALUESOLVER
# define _LAPACKEIGENVALUESOLVER

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Real1DArray.h"
# include "Real2DArray.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void LAPACKEigenvalueSolver_Diagonalize ( SymmetricMatrix *self, Real1DArray *eigenValues, Real2DArray *eigenVectors, Status *status ) ;

# endif
