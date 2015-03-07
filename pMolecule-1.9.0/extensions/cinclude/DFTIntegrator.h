/*------------------------------------------------------------------------------
! . File      : DFTIntegrator.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _DFTINTEGRATOR
# define _DFTINTEGRATOR

# include "Coordinates3.h"
# include "Definitions.h"
# include "DFTFunctionalModel.h"
# include "DFTGrid.h"
# include "QCAtomContainer.h"
# include "QCOnePDM.h"
# include "QCParameters.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void DFTIntegrator_Integrate ( const DFTFunctionalModel *functionalModel    ,
                                            DFTGrid            *grid               ,
                                      const QCAtomContainer    *qcAtoms            ,
                                      const QCParameter        *qcParameters       ,
                                            Coordinates3       *qcCoordinates      ,
                                      const QCOnePDM           *densityP           ,
                                      const QCOnePDM           *densityQ           ,
                                      const Boolean             inCore             ,
                                      const Boolean             isSpinUnrestricted ,
                                            Real               *eQuad              ,
                                            Real               *rhoQuad            ,
# ifdef USEOPENMP
                                            SymmetricMatrix   **fockA              ,
                                            SymmetricMatrix   **fockB              ,
                                            Coordinates3      **gradients3         ,
# else
                                            SymmetricMatrix    *fockA              ,
                                            SymmetricMatrix    *fockB              ,
                                            Coordinates3       *gradients3         ,
# endif
                                            Status             *status             ) ;
# endif
