/*------------------------------------------------------------------------------
! . File      : LAPACKEigenvalueSolver.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . LAPACK eigenvalue solver.
!=================================================================================================================================*/

# include <stdio.h>

# include "LAPACKEigenvalueSolver.h"
# include "Memory.h"

# include "f2clapack.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Quick diagonalization interface with no extraneous objects.
!---------------------------------------------------------------------------------------------------------------------------------*/
void LAPACKEigenvalueSolver_Diagonalize ( SymmetricMatrix *self, Real1DArray *eigenValues, Real2DArray *eigenVectors, Status *status )
{
    if ( ( self != NULL ) && ( eigenValues != NULL ) )
    {
        /* . Various checks - self and eigenValues need to be compact as does dimension 1 of eigenVectors. */
             if ( ! ( SymmetricMatrix_IsCompact ( self ) && Real1DArray_IsCompact ( eigenValues ) ) ) Status_Set ( status, Status_NonCompactArray     ) ;
        else if ( ( eigenVectors != NULL ) && ( ! Real2DArray_IsCompact ( eigenVectors, 1, NULL ) ) ) Status_Set ( status, Status_NonCompactDimension ) ;
        else
        {
            /* . LAPACK specific. */
            auto integer          iFail = 0, *iWork = NULL, lIWork, lRWork, m, n ;
            auto Real            *rWork = NULL ;
            auto SymmetricMatrix *work  = NULL ;
            /* . Work space. */
            n      = ( integer ) self->dimension ;
            lIWork = 3 + 5 * n ;
            lRWork  = 1 + 6 * n + n * n ;
            iWork  = ( integer * ) calloc ( lIWork, sizeof ( integer ) ) ;
            rWork  = Memory_Allocate_Array_Real ( ( Integer ) lRWork ) ;
            work   = SymmetricMatrix_Clone ( self ) ;
            /* . Diagonalization. */
            if ( ( iWork != NULL ) && ( rWork != NULL ) && ( work != NULL ) )
            {
                if ( eigenVectors == NULL )
                {
                    dspevd_ ( "N", "U", &n, work->data, Real1DArray_Data ( eigenValues ), NULL, NULL, rWork, &lRWork, iWork, &lIWork, &iFail ) ;
                }
                else
                {
                    m = ( integer ) eigenVectors->stride0 ;
                    dspevd_ ( "V", "U", &n, work->data, Real1DArray_Data ( eigenValues ), Real2DArray_Data ( eigenVectors ), &m, rWork, &lRWork, iWork, &lIWork, &iFail ) ;
                    Real2DArray_Transpose ( eigenVectors, status ) ;
                }
                if ( iFail != 0 ) { Status_Set ( status, Status_DiagonalizationFailure ) ; printf ( "\nLAPACK Diagonalization Error = %d\n", iFail ) ; }
            }
            else Status_Set ( status, Status_MemoryAllocationFailure ) ;
            /* . Finish up. */
            Memory_Deallocate          ( iWork ) ;
            Memory_Deallocate          ( rWork ) ;
            SymmetricMatrix_Deallocate ( &work ) ;
        }
    }
}
