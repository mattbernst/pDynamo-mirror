/*------------------------------------------------------------------------------
! . File      : DFTFunctionalModel.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _DFTFUNCTIONALMODEL
# define _DFTFUNCTIONALMODEL

# include "Boolean.h"
# include "DFTIntegratorDataBlock.h"
# include "Integer.h"
# include "Integer1DArray.h"
# include "Status.h"
# include "xc.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The DFT functional model type. */
typedef struct {
    Boolean       hasLaplacian        ;
    Boolean       hasSigma            ;
    Boolean       hasTau              ;
    Boolean       isSpinRestricted    ;
    Integer       numberOfFunctionals ;
    Integer       order               ;
    xc_func_type *functionals         ;
} DFTFunctionalModel ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern DFTFunctionalModel *DFTFunctionalModel_Allocate    ( const Integer numberOfFunctionals, Status *status ) ;
extern DFTFunctionalModel *DFTFunctionalModel_Clone       ( const DFTFunctionalModel  *self, Status *status ) ;
extern void                DFTFunctionalModel_Deallocate  (       DFTFunctionalModel **self ) ;
extern void                DFTFunctionalModel_Evaluate    ( const DFTFunctionalModel  *self, DFTIntegratorDataBlock *data ) ;
extern DFTFunctionalModel *DFTFunctionalModel_MakeFromIDs ( const Integer1DArray *ids, const Boolean isSpinRestricted, Status *status ) ;

# endif
