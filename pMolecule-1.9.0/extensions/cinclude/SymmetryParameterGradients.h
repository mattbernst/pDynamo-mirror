/*------------------------------------------------------------------------------
! . File      : SymmetryParameterGradients.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _SYMMETRYPARAMETERGRADIENTS
# define _SYMMETRYPARAMETERGRADIENTS

# include "Coordinates3.h"
# include "Definitions.h"
# include "Matrix33.h"
# include "SymmetryParameters.h"
# include "Transformation3.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The symmetryparametergradients type. */
typedef struct {
    Real dEda      ;
    Real dEdb      ;
    Real dEdc      ;
    Real dEdalpha  ;
    Real dEdbeta   ;
    Real dEdgamma  ;
    Matrix33 *dEdM ;
} SymmetryParameterGradients ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern SymmetryParameterGradients *SymmetryParameterGradients_Allocate              ( void ) ;
extern void                        SymmetryParameterGradients_Deallocate            ( SymmetryParameterGradients **self ) ;
extern void                        SymmetryParameterGradients_CrystalDerivatives    ( SymmetryParameterGradients  *self, const SymmetryParameters *symmetryParameters ) ;
extern void                        SymmetryParameterGradients_FractionalDerivatives ( SymmetryParameterGradients  *self, const SymmetryParameters *symmetryParameters,
                                                                                                          const Coordinates3 *coordinates3, Coordinates3 *gradients3  ) ;
extern void                        SymmetryParameterGradients_ImageDerivatives      ( SymmetryParameterGradients  *self, const SymmetryParameters *symmetryParameters,
                                                                                                                               const Transformation3 *transformation3,
                                                                                                    const Coordinates3 *coordinates3, const Coordinates3 *gradients3  ) ;

# endif
