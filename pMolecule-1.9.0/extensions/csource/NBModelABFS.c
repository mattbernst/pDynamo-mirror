/*------------------------------------------------------------------------------
! . File      : NBModelABFS.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This module implements atom-based force-switching non-bonding interaction
! . procedures. QC coordinates and gradients are in atomic units.
!=================================================================================================================================*/

/* # define DEBUGPRINTING */
# ifdef DEBUGPRINTING
# include <stdio.h>
# endif

# include "ExecutionEnvironment.h"
# include "IndexedSelection.h"
# include "Memory.h"
# include "NBModelABFS.h"
# include "PairwiseInteraction.h"
# include "Selection.h"
# include "SymmetricMatrix.h"
# include "Transformation3.h"
# include "Units.h"

/* . Defaults for some variables. */
# define DEFAULT_CHECKFORINVERSES     True
# define DEFAULT_DAMPINGCUTOFF         0.5
# define DEFAULT_DIELECTRIC            1.0
# define DEFAULT_ELECTROSTATICSCALE14  1.0
# define DEFAULT_IMAGEEXPANDFACTOR       0
# define DEFAULT_INNERCUTOFF           8.0
# define DEFAULT_LISTCUTOFF           13.5
# define DEFAULT_OUTERCUTOFF          12.0
# define DEFAULT_QCMMCOUPLING         QCMMLinkAtomCoupling_RC
# define DEFAULT_USECENTERING         False

/*==================================================================================================================================
! . Local procedures.
!=================================================================================================================================*/
/* . Updating. */
static Boolean CheckForImageUpdate ( const SymmetryParameters        *symmetryParameters          ,
                                           SymmetryParameters        *referenceParameters         ,
                                           ImageList                 *imageList                   ,
                                     const Real                       listCutoff                  ,
                                     const Real                       outerCutoff                 ,
                                     const Real                       maximumDisplacement         ) ;

static Boolean CheckForUpdate      ( const Coordinates3              *coordinates3                ,
                                           Coordinates3              *reference3                  ,
                                           Selection                 *fixedAtoms                  ,
                                     const Real                       listCutoff                  ,
                                     const Real                       outerCutoff                 ,
                                           Real                      *maximumDisplacement         ) ;

static Status GenerateImageLists   ( const PairListGenerator          *generator                  ,
                                           Selection                  *mmSelection                ,
                                           Selection                  *qcbSelection               ,
                                           Selection                  *qcpSelection               ,
                                           Selection                  *freeSelection              ,
                                     const Coordinates3               *coordinates3               ,
                                     const Coordinates3               *qcCoordinates3             ,
                                     const Coordinates3               *mmCoordinates3             ,
                                     const Transformation3Container   *transformations            ,
                                     const SymmetryParameters         *symmetryParameters         ,
                                           RegularGrid                *mmGrid                     ,
                                           RegularGridOccupancy       *mmOccupancy                ,
                                           RegularGrid                *qcGrid                     ,
                                           RegularGridOccupancy       *qcOccupancy                ,
                                     const QCMMLinkAtomCoupling        qcmmCoupling               ,
                                     const Boolean                     checkForInverses           ,
                                     const Integer                     expandFactor               ,
                                           ImageList                 **inbmmmm                    ,
                                           ImageList                 **inbqcmmel                  ,
                                           ImageList                 **inbqcmmlj                  ,
                                           ImageList                 **inbqcqcel                  ,
                                           ImageList                 **inbqcqclj                  ) ;

static Status GenerateLists        ( const PairListGenerator          *generator                  ,
                                           Selection                  *mmSelection                ,
                                           Selection                  *qcbSelection               ,
                                           Selection                  *qcpSelection               ,
                                           Selection                  *freeSelection              ,
                                           PairList                   *exclusions                 ,
                                           PairList                   *qcmmExclusions             ,
                                     const Coordinates3               *coordinates3               ,
                                     const Coordinates3               *qcCoordinates3             ,
                                     const Coordinates3               *mmCoordinates3             ,
                                           RegularGrid                *mmGrid                     ,
                                           RegularGridOccupancy       *mmOccupancy                ,
                                           RegularGrid                *qcGrid                     ,
                                           RegularGridOccupancy       *qcOccupancy                ,
                                     const QCMMLinkAtomCoupling        qcmmCoupling               ,
                                           PairList                  **nbmmmm                     ,
                                           PairList                  **nbqcmmlj                   ,
                                           PairList                  **nbqcmmel                   ) ;

static Status GenerateLists14      (       Selection                  *mmSelection                ,
                                           Selection                  *qcbSelection               ,
                                           Selection                  *qcpSelection               ,
                                           Selection                  *freeSelection              ,
                                           PairList                   *interactions14             ,
                                     const QCMMLinkAtomCoupling        qcmmCoupling               ,
                                           PairList                  **nbmmmm                     ,
                                           PairList                  **nbqcmmlj                   ,
                                           PairList                  **nbqcmmel                   ) ;

/* . Image energies and gradients. */
static void MMMMImageEnergy        ( const PairwiseInteractionABFS    *mmmmPairwiseInteraction    ,
                                     const MMAtomContainer            *mmAtoms                    ,
                                     const LJParameterContainer       *ljParameters               ,
                                           SymmetryParameters         *symmetryParameters         ,
                                           ImageList                  *imageList                  ,
                                     const Real                        electrostaticScale         ,
                                     const Coordinates3               *coordinates3               ,
                                           Real                       *eElectrostatic             ,
                                           Real                       *eLennardJones              ,
# ifdef USEOPENMP
                                     const Integer                     numberOfThreads            ,
                                           Coordinates3              **gradientsArray             ,
# else
                                           Coordinates3               *gradients3                 ,
# endif
                                           SymmetryParameterGradients *symmetryParameterGradients ) ;

static void QCMMImageGradients     ( const PairwiseInteractionABFS    *qcmmPairwiseInteraction    ,
                                     const Real1DArray                *mmCharges                  ,
                                     const QCAtomContainer            *qcAtoms                    ,
                                           SymmetryParameters         *symmetryParameters         ,
                                           ImageList                  *imageList                  ,
                                     const Real                        electrostaticScale         ,
                                     const Coordinates3               *coordinates3               ,
                                     const Coordinates3               *mmCoordinates3             ,
                                     const Real1DArray                *qcCharges                  ,
                                           Coordinates3               *gradients3                 ,
                                           Coordinates3               *mmGradients3               ,
                                           SymmetryParameterGradients *symmetryParameterGradients ) ;

static void QCMMImagePotentials    ( const PairwiseInteractionABFS    *qcmmPairwiseInteraction    ,
                                     const Real1DArray                *mmCharges                  ,
                                     const QCAtomContainer            *qcAtoms                    ,
                                           SymmetryParameters         *symmetryParameters         ,
                                           ImageList                  *imageList                  ,
                                     const Real                        electrostaticScale         ,
                                     const Coordinates3               *coordinates3               ,
                                     const Coordinates3               *mmCoordinates3             ,
                                           Real1DArray                *potentials                 ) ;

static void QCQCImageGradients     ( const PairwiseInteractionABFS    *qcqcPairwiseInteraction    ,
                                     const QCAtomContainer            *qcAtoms                    ,
                                           SymmetryParameters         *symmetryParameters         ,
                                           ImageList                  *imageList                  ,
                                     const Real                        electrostaticScale         ,
                                     const Coordinates3               *coordinates3               ,
                                     const Coordinates3               *qcCoordinates3             ,
                                     const Real1DArray                *qcCharges                  ,
                                           Coordinates3               *gradients3                 ,
                                           SymmetryParameterGradients *symmetryParameterGradients ) ;

static void QCQCImagePotentials    ( const PairwiseInteractionABFS    *qcqcPairwiseInteraction    ,  
                                           SymmetryParameters         *symmetryParameters         ,  
                                           ImageList                  *imageList                  ,  
                                     const Real                        electrostaticScale         ,  
                                     const Coordinates3               *qcCoordinates3             ,  
                                           SymmetricMatrix            *potentials                 ) ;

/*==================================================================================================================================
! . Basic procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
NBModelABFS *NBModelABFS_Allocate ( void )
{
    NBModelABFS *self = NULL ;
    self = ( NBModelABFS * ) Memory_Allocate ( sizeof ( NBModelABFS ) ) ;
    if ( self != NULL )
    {
        self->checkForInverses               = DEFAULT_CHECKFORINVERSES     ;
        self->dampingCutoff                  = DEFAULT_DAMPINGCUTOFF        ;
        self->dielectric                     = DEFAULT_DIELECTRIC           ;
        self->electrostaticScale14           = DEFAULT_ELECTROSTATICSCALE14 ;
        self->imageExpandFactor              = DEFAULT_IMAGEEXPANDFACTOR    ;
        self->innerCutoff                    = DEFAULT_INNERCUTOFF          ;
        self->listCutoff                     = DEFAULT_LISTCUTOFF           ;
        self->outerCutoff                    = DEFAULT_OUTERCUTOFF          ;
        self->qcmmCoupling                   = DEFAULT_QCMMCOUPLING         ;
        self->useCentering                   = DEFAULT_USECENTERING         ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
NBModelABFS *NBModelABFS_Clone ( const NBModelABFS *self )
{
    NBModelABFS *new = NULL ;
    if ( self != NULL )
    {
        new = NBModelABFS_Allocate ( ) ;
        new->checkForInverses     = self->checkForInverses     ;
        new->dampingCutoff        = self->dampingCutoff        ;
        new->dielectric           = self->dielectric           ;
        new->electrostaticScale14 = self->electrostaticScale14 ;
        new->imageExpandFactor    = self->imageExpandFactor    ;
        new->innerCutoff          = self->innerCutoff          ;
        new->listCutoff           = self->listCutoff           ;
        new->outerCutoff          = self->outerCutoff          ;
        new->qcmmCoupling         = self->qcmmCoupling         ;
        new->useCentering         = self->useCentering         ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFS_Deallocate ( NBModelABFS **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electrostatic and LJ MM/MM energy and gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFS_MMMMEnergy ( const NBModelABFS *self, const PairwiseInteractionABFS *mmmmPairwiseInteraction, NBModelABFSState *nbState )
{
    if ( ( self != NULL ) && ( mmmmPairwiseInteraction != NULL ) && ( nbState != NULL ) )
    {
# ifdef USEOPENMP
        auto Coordinates3 *gradients[1] = { nbState->gradients3 }, **gradientsArray = NULL ;
        if ( nbState->gradients3 != NULL ) gradientsArray = gradients ; 
# endif
        auto Real eScale ;
        /* . Initialization. */
        eScale = 1.0e+00 / self->dielectric ;
        /* . MM/MM Image. */
        MMMMImageEnergy ( mmmmPairwiseInteraction             ,
                          nbState->mmAtoms                    ,
                          nbState->ljParameters               ,
                          nbState->symmetryParameters         ,
                          nbState->inbmmmm                    ,
                          eScale                              ,
                          nbState->coordinates3               ,
                          &(nbState->eimmmel)                 ,
                          &(nbState->eimmmlj)                 ,
# ifdef USEOPENMP
                          nbState->numberOfThreads            ,
                          nbState->gradientsArray             ,
# else
                          nbState->gradients3                 ,
# endif
                          nbState->symmetryParameterGradients ) ;
        /* . MM/MM. */
        PairwiseInteractionABFS_MMMMEnergy ( mmmmPairwiseInteraction  ,
                                             nbState->mmAtoms         ,
                                             nbState->ljParameters    ,
                                             nbState->nbmmmm          ,
                                             eScale                   ,
                                             1.0e+00                  ,
                                             nbState->coordinates3    ,
                                             nbState->coordinates3    ,
                                             &(nbState->emmel)        ,
                                             &(nbState->emmlj)        ,
# ifdef USEOPENMP
                                             nbState->numberOfThreads ,
                                             nbState->gradientsArray  ,
                                             nbState->gradientsArray  ) ;
# else
                                             nbState->gradients3      ,
                                             nbState->gradients3      ) ;
# endif
        /* . MM/MM 1-4. */
        eScale *= self->electrostaticScale14 ;
        PairwiseInteractionABFS_MMMMEnergy ( mmmmPairwiseInteraction ,
                                             nbState->mmAtoms        ,
                                             nbState->ljParameters14 ,
                                             nbState->nbmmmm14       ,
                                             eScale                  ,
                                             1.0e+00                 ,
                                             nbState->coordinates3   ,
                                             nbState->coordinates3   ,
                                             &(nbState->emmel14)     ,
                                             &(nbState->emmlj14)     ,
# ifdef USEOPENMP
                                             1                       ,
                                             gradientsArray          ,
                                             gradientsArray          ) ;
# else
                                             nbState->gradients3     ,
                                             nbState->gradients3     ) ;
# endif

# ifdef USEOPENMP
        /* . Accumulate gradients. */
        NBModelABFSState_AccumulateGradients ( nbState ) ;
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM LJ energy and gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFS_QCMMEnergyLJ ( const NBModelABFS *self, const PairwiseInteractionABFS *mmmmPairwiseInteraction, NBModelABFSState *nbState )
{
    if ( ( self != NULL ) && ( mmmmPairwiseInteraction != NULL ) && ( nbState != NULL ) && ( nbState->qcAtoms != NULL ) )
    {
# ifdef USEOPENMP
        auto Coordinates3 *gradients[1] = { nbState->gradients3 }, **gradientsArray = NULL ;
        if ( nbState->gradients3 != NULL ) gradientsArray = gradients ; 
# endif
        /* . QC/MM. */
        PairwiseInteractionABFS_MMMMEnergy ( mmmmPairwiseInteraction ,
                                             nbState->mmAtoms        ,
                                             nbState->ljParameters   ,
                                             nbState->nbqcmmlj       ,
                                             0.0e+00                 ,
                                             1.0e+00                 , 
                                             nbState->coordinates3   ,
                                             nbState->coordinates3   ,
                                             NULL                    ,
                                             &(nbState->eqcmmlj)     ,
# ifdef USEOPENMP
                                             1                       ,
                                             gradientsArray          ,
                                             gradientsArray          ) ;
# else
                                             nbState->gradients3     ,
                                             nbState->gradients3     ) ;
# endif
        /* . QC/MM 1-4. */
        PairwiseInteractionABFS_MMMMEnergy ( mmmmPairwiseInteraction ,
                                             nbState->mmAtoms        ,
                                             nbState->ljParameters14 ,
                                             nbState->nbqcmmlj14     ,
                                             0.0e+00                 ,
                                             1.0e+00                 ,
                                             nbState->coordinates3   ,
                                             nbState->coordinates3   ,
                                             NULL                    ,
                                             &(nbState->eqcmmlj14)   ,
# ifdef USEOPENMP
                                             1                       ,
                                             gradientsArray          ,
                                             gradientsArray          ) ;
# else
                                             nbState->gradients3     ,
                                             nbState->gradients3     ) ;
# endif
        /* . QC/MM Image. */
        MMMMImageEnergy ( mmmmPairwiseInteraction             ,
                          nbState->mmAtoms                    ,
                          nbState->ljParameters               ,
                          nbState->symmetryParameters         ,
                          nbState->inbqcmmlj                  ,
                          0.0e+00                             ,
                          nbState->coordinates3               ,
                          NULL                                ,
                          &(nbState->eimqcmmlj)               ,
# ifdef USEOPENMP
                          1                                   ,
                          gradientsArray                      ,
# else
                          nbState->gradients3                 ,
# endif
                          nbState->symmetryParameterGradients ) ;
        /* . QC/QC Image. */
        MMMMImageEnergy ( mmmmPairwiseInteraction             ,
                          nbState->mmAtoms                    ,
                          nbState->ljParameters               ,
                          nbState->symmetryParameters         ,
                          nbState->inbqcqclj                  ,
                          0.0e+00                             ,
                          nbState->coordinates3               ,
                          NULL                                ,
                          &(nbState->eimqcqclj)               ,
# ifdef USEOPENMP
                          1                                   ,
                          gradientsArray                      ,
# else
                          nbState->gradients3                 ,
# endif
                          nbState->symmetryParameterGradients ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM electrostatic gradients in regular units.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFS_QCMMGradients ( const NBModelABFS *self, const PairwiseInteractionABFS *qcmmPairwiseInteraction, const PairwiseInteractionABFS *qcqcPairwiseInteraction, NBModelABFSState *nbState )
{
    if ( ( self != NULL ) && ( qcmmPairwiseInteraction != NULL ) && ( qcqcPairwiseInteraction != NULL ) && ( nbState != NULL ) && ( nbState->qcAtoms != NULL ) && ( nbState->gradients3 != NULL ) )
    {
        auto Real eScale ;
        /* . Initialization. */
        eScale = 1.0e+00 / self->dielectric ;
        /* . QC/MM Image. */
        QCMMImageGradients ( qcmmPairwiseInteraction             ,
                             nbState->mmCharges                  ,
                             nbState->qcAtoms                    ,
                             nbState->symmetryParameters         ,
                             nbState->inbqcmmel                  ,
                             eScale                              ,
                             nbState->coordinates3               ,
                             nbState->mmCoordinates3             ,
                             nbState->qcCharges                  ,
                             nbState->gradients3                 ,
                             nbState->mmGradients3               ,
                             nbState->symmetryParameterGradients ) ;
        /* . QC/QC Image. */
        QCQCImageGradients ( qcqcPairwiseInteraction             ,
                             nbState->qcAtoms                    ,
                             nbState->symmetryParameters         ,
                             nbState->inbqcqcel                  ,
                             eScale                              ,
                             nbState->coordinates3               ,
                             nbState->qcCoordinates3             ,
                             nbState->qcCharges                  ,
                             nbState->gradients3                 ,
                             nbState->symmetryParameterGradients ) ;
        /* . QC/MM. */
        PairwiseInteractionABFS_QCMMGradients ( qcmmPairwiseInteraction ,
                                                nbState->mmCharges      ,
                                                nbState->qcAtoms        ,
                                                nbState->nbqcmmel       ,
                                                eScale                  ,
                                                nbState->coordinates3   ,
                                                nbState->mmCoordinates3 ,
                                                nbState->qcCharges      ,
                                                nbState->gradients3     ,
                                                nbState->mmGradients3   ) ;
        /* . QC/MM 1-4. */
        eScale *= self->electrostaticScale14 ;
        PairwiseInteractionABFS_QCMMGradients ( qcmmPairwiseInteraction ,
                                                nbState->mmCharges      ,
                                                nbState->qcAtoms        ,
                                                nbState->nbqcmmel14     ,
                                                eScale                  ,
                                                nbState->coordinates3   ,
                                                nbState->mmCoordinates3 ,
                                                nbState->qcCharges      ,
                                                nbState->gradients3     ,
                                                nbState->mmGradients3   ) ;
        /* . Put MM gradients3 into gradients3. */
        if ( self->qcmmCoupling != QCMMLinkAtomCoupling_MM ) QCAtomContainer_SetMMGradients3 ( nbState->qcAtoms, nbState->mmGradients3, nbState->gradients3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM electrostatic potentials in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
void NBModelABFS_QCMMPotentials  ( const NBModelABFS *self, const PairwiseInteractionABFS *qcmmPairwiseInteraction, const PairwiseInteractionABFS *qcqcPairwiseInteraction, NBModelABFSState *nbState )
{
    if ( ( self != NULL ) && ( qcmmPairwiseInteraction != NULL ) && ( qcqcPairwiseInteraction != NULL ) && ( nbState != NULL ) && ( nbState->qcAtoms != NULL ) )
    {
        auto Real eScale ;
        /* . Initialization. */
        eScale = 1.0e+00 / self->dielectric ;
        /* . QC/MM Image. */
        QCMMImagePotentials ( qcmmPairwiseInteraction     ,
                              nbState->mmCharges          ,
                              nbState->qcAtoms            ,
                              nbState->symmetryParameters ,
                              nbState->inbqcmmel          ,
                              eScale                      ,
                              nbState->coordinates3       ,
                              nbState->mmCoordinates3     ,
                              nbState->qcmmPotentials     ) ;
        /* . QC/QC Image. */
        QCQCImagePotentials ( qcqcPairwiseInteraction     ,
                              nbState->symmetryParameters ,
                              nbState->inbqcqcel          ,
                              eScale                      ,
                              nbState->qcCoordinates3     ,
                              nbState->qcqcPotentials     ) ;
        /* . QC/MM. */
        PairwiseInteractionABFS_QCMMPotentials ( qcmmPairwiseInteraction ,
                                                 nbState->mmCharges      ,
                                                 nbState->qcAtoms        ,
                                                 nbState->nbqcmmel       ,
                                                 eScale                  ,
                                                 nbState->coordinates3   ,
                                                 nbState->mmCoordinates3 ,
                                                 nbState->qcmmPotentials ) ;
        /* . QC/MM 1-4. */
        eScale *= self->electrostaticScale14 ;
        PairwiseInteractionABFS_QCMMPotentials ( qcmmPairwiseInteraction ,
                                                 nbState->mmCharges      ,
                                                 nbState->qcAtoms        ,
                                                 nbState->nbqcmmel14     ,
                                                 eScale                  ,
                                                 nbState->coordinates3   ,
                                                 nbState->mmCoordinates3 ,
                                                 nbState->qcmmPotentials ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Update the non-bonding lists if necessary.
! . Updates are done if:
!   (i)   the cutoffs (or other options) change.
!   (ii)  the state is new (i.e. is empty).
!   (iii) the coordinates or symmetryParameters have changed significantly.
! . 1-4 lists only need to be created when the state is new.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean NBModelABFS_Update  ( const NBModelABFS *self, const PairListGenerator *generator, NBModelABFSState *nbState, Status *status )
{
    Boolean doUpdate = False ;
    if ( ( self != NULL ) && ( generator != NULL ) && ( nbState != NULL ) )
    {
        auto Real   maximumDisplacement = 0.0e+00 ;
        auto Status localStatus ;
        /* . Check for a new state. */
        doUpdate = nbState->isNew ;
        /* . Update the 1-4 lists. Otherwise they remain unchanged. */
        if ( doUpdate )
        {
            localStatus = GenerateLists14 ( nbState->mmSelection    ,
                                            nbState->qcbSelection   ,
                                            nbState->qcpSelection   ,
                                            nbState->freeSelection  ,
                                            nbState->interactions14 ,
                                            self->qcmmCoupling      ,
                                            &(nbState->nbmmmm14)    ,
                                            &(nbState->nbqcmmlj14)  ,
                                            &(nbState->nbqcmmel14)  ) ;
            if ( localStatus != Status_Continue ) Status_SafeSet ( status, localStatus ) ;
        }
        /* . Check for new options. */
        doUpdate = doUpdate || ( self->listCutoff   != nbState->listCutoff  ) ||
                               ( self->outerCutoff  != nbState->outerCutoff ) ;
        if ( doUpdate )
        {
            nbState->listCutoff  = self->listCutoff  ;
            nbState->outerCutoff = self->outerCutoff ;
        }
        /* . Check for significant coordinate changes. */
        doUpdate = doUpdate || CheckForUpdate ( nbState->inputCoordinates3     ,
                                                nbState->referenceCoordinates3 ,
                                                nbState->fixedAtoms            ,
                                                self->listCutoff               ,
                                                nbState->outerCutoff           ,
                                                &maximumDisplacement           ) ;
        /* . Generate coordinate-dependent quantities. */
        NBModelABFSState_InitializeCoordinates3 ( nbState, doUpdate ) ;
        /* . Update the lists and make sure the reference values are up-to-date. */
        if ( doUpdate )
        {
            NBModelABFSState_GridInitialize ( nbState, generator, status ) ;
            localStatus = GenerateLists ( generator               ,
                                          nbState->mmSelection    ,
                                          nbState->qcbSelection   ,
                                          nbState->qcpSelection   ,
                                          nbState->freeSelection  ,
                                          nbState->exclusions     ,
                                          nbState->qcmmExclusions ,
                                          nbState->coordinates3   ,
                                          nbState->qcCoordinates3 ,
                                          nbState->mmCoordinates3 ,
                                          nbState->mmGrid         ,
                                          nbState->mmOccupancy    ,
                                          nbState->qcGrid         ,
                                          nbState->qcOccupancy    ,
                                          self->qcmmCoupling      ,
                                          &(nbState->nbmmmm)      ,
                                          &(nbState->nbqcmmlj)    ,
                                          &(nbState->nbqcmmel)    ) ;
            if ( localStatus != Status_Continue ) Status_SafeSet ( status, localStatus ) ;
            Coordinates3_CopyTo ( nbState->inputCoordinates3, nbState->referenceCoordinates3, NULL ) ;
        }
        /* . Symmetry is present. */
        if ( nbState->transformations != NULL )
        {
            /* . Check for significant symmetry changes. */
            doUpdate = doUpdate || CheckForImageUpdate ( nbState->symmetryParameters          ,
                                                         nbState->referenceSymmetryParameters ,
                                                         nbState->inbmmmm                     ,
                                                         self->listCutoff                     ,
                                                         nbState->outerCutoff                 ,
                                                         maximumDisplacement                  ) ;
            /* . Update the lists and make sure the reference values are up-to-date. */
            if ( doUpdate )
            {
                NBModelABFSState_GridInitialize ( nbState, generator, status ) ;
                localStatus = GenerateImageLists ( generator                   ,
                                                   nbState->mmSelection        ,
                                                   nbState->qcbSelection       ,
                                                   nbState->qcpSelection       ,
                                                   nbState->freeSelection      ,
                                                   nbState->coordinates3       ,
                                                   nbState->qcCoordinates3     ,
                                                   nbState->mmCoordinates3     ,
                                                   nbState->transformations    ,
                                                   nbState->symmetryParameters ,
                                                   nbState->mmGrid             ,
                                                   nbState->mmOccupancy        ,
                                                   nbState->qcGrid             ,
                                                   nbState->qcOccupancy        ,
                                                   self->qcmmCoupling          ,
                                                   self->checkForInverses      ,
                                                   self->imageExpandFactor     ,
                                                   &(nbState->inbmmmm)         ,
                                                   &(nbState->inbqcmmel)       ,
                                                   &(nbState->inbqcmmlj)       ,
                                                   &(nbState->inbqcqcel)       ,
                                                   &(nbState->inbqcqclj)       ) ;
                if ( localStatus != Status_Continue ) Status_SafeSet ( status, localStatus ) ;
                SymmetryParameters_CopyTo ( nbState->symmetryParameters, nbState->referenceSymmetryParameters ) ;
            }
        }
        /* . Accumulate statistics. */
        if ( doUpdate )
        {
            NBModelABFSState_GridFinalize         ( nbState ) ;
            NBModelABFSState_StatisticsAccumulate ( nbState ) ;
        }
        /* . The state is no longer new. */
        nbState->isNew = False ;
    }
    return doUpdate ;
}

/*==================================================================================================================================
! . Private procedures.
! . These depend upon neither NBModelABFS nor NBModelABFSState.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Check to see whether image list updates are needed.
! . An update occurs whenever the difference in a lattice translation added
! . to the maximum displacement of a particle exceeds half the buffer distance
! . (defined as the difference of the list and outer cutoffs).
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean CheckForImageUpdate ( const SymmetryParameters *symmetryParameters  ,
                                           SymmetryParameters *referenceParameters ,
                                           ImageList          *imageList           ,
                                     const Real                listCutoff          ,
                                     const Real                outerCutoff         ,
                                     const Real                maximumDisplacement )
{
    Boolean doUpdate = False ;
    if ( ( symmetryParameters != NULL ) && ( referenceParameters != NULL ) && ( imageList != NULL ) )
    {
        auto Real      buffac, di ;
        auto Image    *image        = NULL ;
        auto Matrix33 *dM           = NULL ;
        auto Vector3  *displacement = NULL ;

        /* . Initialization. */
        buffac = listCutoff - outerCutoff - maximumDisplacement ;

        /* . Allocate space. */
        dM           = Matrix33_Allocate ( ) ;
        displacement = Vector3_Allocate  ( ) ;

        /* . Get the difference in the lattice vector matrices. */
        Matrix33_CopyTo            ( symmetryParameters->M, dM, NULL ) ;
        Matrix33_AddScaledMatrix33 ( dM, - 1.0e+00, referenceParameters->M, NULL ) ;

        /* . Loop over the images in the list. */
        List_Iterate_Initialize ( imageList->images ) ;
        while ( ( image = ImageList_Iterate ( imageList ) ) != NULL )
        {
            /* . Get the transformation. */
            Vector3_CopyTo ( image->transformation3->translation, displacement, NULL ) ;
            displacement->data[0] += ( Real ) image->a ;
            displacement->data[1] += ( Real ) image->b ;
            displacement->data[2] += ( Real ) image->c ;

            /* . Get the displacement and its size squared. */
            Matrix33_ApplyToVector3 ( dM, displacement ) ;
            di = Vector3_Norm2 ( displacement ) ;

            /* . Check for displacements that are too large. */
            if ( di > buffac ) { doUpdate = True ; break ; }
        }

        /* . Finish up. */
        Matrix33_Deallocate ( &dM           ) ;
        Vector3_Deallocate  ( &displacement ) ;
    }
    return doUpdate ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for a regular update.
! . Updates occur whenever a particle has moved by more than half the buffer
! . distance (defined as the difference of the list and outer cutoffs).
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean CheckForUpdate  ( const Coordinates3 *coordinates3        ,
                                       Coordinates3 *reference3          ,
                                       Selection    *fixedAtoms          ,
                                 const Real          listCutoff          ,
                                 const Real          outerCutoff         ,
                                       Real         *maximumDisplacement )
{
    Boolean doUpdate = False ;
    if ( ( coordinates3 != NULL ) && ( reference3 != NULL ) )
    {
        auto Integer i ;
        auto Real    buffac, buffacsq, dx, dy, dz, maxr2 = 0.0e+00, r2, xi, xr, yi, yr, zi, zr ;

        /* . Get the buffer variables. */
        buffac   = 0.5e+00 * ( listCutoff - outerCutoff ) ;
        buffacsq = buffac * buffac ;

        /* . Check for displacements that are greater than |buffac|. */
        if ( fixedAtoms == NULL )
        {
            for ( i = 0 ; i < coordinates3->length0 ; i++ )
            {
                Coordinates3_GetRow ( coordinates3, i, xi, yi, zi ) ;
                Coordinates3_GetRow ( reference3,   i, xr, yr, zr ) ;
                dx    = xi - xr ;
                dy    = yi - yr ;
                dz    = zi - zr ;
                r2    = dx * dx + dy * dy + dz * dz ;
                maxr2 = Maximum ( maxr2, r2 ) ;
                if ( r2 > buffacsq ) { doUpdate = True ; break ; }
            }
        }
        else
        {
            Selection_MakeFlags ( fixedAtoms, coordinates3->length0 ) ;
            for ( i = 0 ; i < coordinates3->length0 ; i++ )
            {
                if ( ! fixedAtoms->flags[i] )
                {
                    Coordinates3_GetRow ( coordinates3, i, xi, yi, zi ) ;
                    Coordinates3_GetRow ( reference3,   i, xr, yr, zr ) ;
                    dx    = xi - xr ;
                    dy    = yi - yr ;
                    dz    = zi - zr ;
                    r2    = dx * dx + dy * dy + dz * dz ;
                    maxr2 = Maximum ( maxr2, r2 ) ;
                    if ( r2 > buffacsq ) { doUpdate = True ; break ; }
                }
            }
        }

        /* . Set the maximum displacement. */
        (*maximumDisplacement) = sqrt ( maxr2 ) ;
    }
    return doUpdate ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate image lists.
! . The lists are MM/MM el. and LJ, QC/MM and QC/QC LJ, QC/MM el. and QC/QC el.
! . It is possible to have all MM, all QC or both MM and QC atoms.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Status GenerateImageLists ( const PairListGenerator        *generator          ,
                                         Selection                *mmSelection        ,
                                         Selection                *qcbSelection       ,
                                         Selection                *qcpSelection       ,
                                         Selection                *freeSelection      ,
                                   const Coordinates3             *coordinates3       ,
                                   const Coordinates3             *qcCoordinates3     ,
                                   const Coordinates3             *mmCoordinates3     ,
                                   const Transformation3Container *transformations    ,
                                   const SymmetryParameters       *symmetryParameters ,
                                         RegularGrid              *mmGrid             ,
                                         RegularGridOccupancy     *mmOccupancy        ,
                                         RegularGrid              *qcGrid             ,
                                         RegularGridOccupancy     *qcOccupancy        ,
                                   const QCMMLinkAtomCoupling      qcmmCoupling       ,
                                   const Boolean                   checkForInverses   ,
                                   const Integer                   expandFactor       ,
                                         ImageList               **inbmmmm            ,
                                         ImageList               **inbqcmmel          ,
                                         ImageList               **inbqcmmlj          ,
                                         ImageList               **inbqcqcel          ,
                                         ImageList               **inbqcqclj          )
{
    Status status = Status_Continue ;
    /* . Deallocate existing lists. */
    ImageList_Deallocate ( inbmmmm   ) ;
    ImageList_Deallocate ( inbqcmmel ) ;
    ImageList_Deallocate ( inbqcmmlj ) ;
    ImageList_Deallocate ( inbqcqcel ) ;
    ImageList_Deallocate ( inbqcqclj ) ;
    /* . Generate the lists. */
    if ( ( generator != NULL ) && ( coordinates3 != NULL ) && ( transformations != NULL ) && ( symmetryParameters != NULL ) && ( generator->cutoff > 0.0 ) )
    {
        auto Boolean          QMM, QOK, QQC, QQCMM, QSKIP, QTSKIP ;
        auto Integer          a, ahigh, ainverse, alow, b, bhigh, binverse, blow, c, chigh, cinverse, clow, t, tinverse ;
	auto Real             defaultscale, scale ;
        auto Coordinates3    *icoordinates3 = NULL, *immCoordinates3 = NULL, *iqcCoordinates3 = NULL ;
        auto PairList        *pairlist = NULL ;
        auto Transformation3 *itransformation3 = NULL ;
        auto Vector3         *displacement = NULL, *ilower = NULL, *iupper = NULL, *lower = NULL, *upper = NULL ;

        /* . Set some flags. */
        QMM   = ( mmSelection == NULL ) || ( ( mmSelection != NULL ) && ( mmSelection->nindices > 0 ) ) ;
        QQC   = ( qcpSelection != NULL ) && ( qcpSelection->nindices > 0 ) ;
        QQCMM = ( QQC && QMM ) ;

        /* . The system must be pure MM or pure QC if checkForInverses is False. */
        /* . QC/MM systems automatically make use of the symmetry to avoid having to create two lists - i.e. QC/MM AND MM/QC for each image. */
        QOK = checkForInverses || ( ! QQCMM ) ;
# ifdef DEBUGPRINTING
printf ( "Into GenerateImageLists %u %u %u\n", QMM, QQC, QQCMM ) ;
# endif
        /* . Allocate space. */
        /* . Image data. */
        icoordinates3    = Coordinates3_Allocate ( coordinates3->length0 ) ;
        itransformation3 = Transformation3_Allocate ( ) ;
        displacement     = Vector3_Allocate         ( ) ;
        ilower           = Vector3_Allocate         ( ) ;
        iupper           = Vector3_Allocate         ( ) ;
        lower            = Vector3_Allocate         ( ) ;
        upper            = Vector3_Allocate         ( ) ;
        QOK = QOK && ( icoordinates3 != NULL ) && ( itransformation3 != NULL ) && ( displacement != NULL ) && ( ilower != NULL ) && ( iupper != NULL ) && ( lower != NULL ) && ( upper != NULL ) ;
        if ( qcmmCoupling != QCMMLinkAtomCoupling_MM )
        {
            immCoordinates3 = Coordinates3_Allocate ( mmCoordinates3->length0 ) ;
            QOK = QOK && ( immCoordinates3 != NULL ) ;
        }

        /* . Check QC-related quantities. */
        if ( QQC )
        {
            if ( qcCoordinates3 != NULL ) iqcCoordinates3 = Coordinates3_Allocate ( qcCoordinates3->length0 ) ;
            QOK = QOK && ( iqcCoordinates3 != NULL ) && ( qcCoordinates3 != NULL ) ;
        }

        /* . Allocate the image lists. */
        if ( QMM   ) { (*inbmmmm)   = ImageList_Allocate ( ) ; QOK = QOK && ( (*inbmmmm  ) != NULL ) ; }
        if ( QQC   ) { (*inbqcqcel) = ImageList_Allocate ( ) ; QOK = QOK && ( (*inbqcqcel) != NULL ) ; }
        if ( QQC   ) { (*inbqcqclj) = ImageList_Allocate ( ) ; QOK = QOK && ( (*inbqcqclj) != NULL ) ; }
        if ( QQCMM ) { (*inbqcmmel) = ImageList_Allocate ( ) ; QOK = QOK && ( (*inbqcmmel) != NULL ) ; }
        if ( QQCMM ) { (*inbqcmmlj) = ImageList_Allocate ( ) ; QOK = QOK && ( (*inbqcmmlj) != NULL ) ; }

        /* . Check memory. */
        if ( QOK )
        {
            /* . Find the orthorhombic box within which interactions are to be sought. */
            /* . Searching is only done on the normal coordinates and not those of the QC atoms - this could be changed but it probably isn't necessary. */
            Coordinates3_EnclosingOrthorhombicBox ( coordinates3, NULL, NULL, lower, upper ) ;
            Vector3_Add ( upper, lower, NULL ) ; /* . To get absolute upper coordinates rather than extents. */
            Vector3_AddScalar ( lower, - generator->cutoff ) ;
            Vector3_AddScalar ( upper,   generator->cutoff ) ;

            /* . Loop over the transformations. */
            for ( t = 0 ; t < transformations->nitems ; t++ )
            {

                /* . Deal with inverses - this can be switched off. */
                tinverse = transformations->inverses[t] ;
                if ( ! checkForInverses ) tinverse = -1 ;

                /* . Set default values. */
                defaultscale = 0.5e+00 ;
                QTSKIP       = False   ;

                /* . Inverses are to be treated. */
                /* . The reasonable assumption is made that the search ranges for the transformation and its inverse are compatible. */
                if ( tinverse >= 0 )
                {
                    /* . Skip this transformation if its inverse will occur later - except for QC/MM cases. */
                    if ( t < tinverse )
                    {
                        QTSKIP = True ;
                        if ( ! QQCMM ) continue ;
                    }

                    /* . Reset the default scale - will work except for some self-inverses (see below). */
                    defaultscale = 1.0e+00 ;
                }

                /* . Generate the image transformation in real space. */
                Transformation3_Copy          ( itransformation3, transformations->items[t] ) ;
                Transformation3_Orthogonalize ( itransformation3, symmetryParameters->M, symmetryParameters->inverseM ) ;

                /* . Generate the coordinates. */
                Coordinates3_CopyTo    ( coordinates3 , icoordinates3   , NULL ) ;
                Coordinates3_Transform ( icoordinates3, itransformation3, NULL ) ;
                if ( qcmmCoupling != QCMMLinkAtomCoupling_MM )
                {
                    Coordinates3_CopyTo    ( mmCoordinates3 , immCoordinates3 , NULL ) ;
                    Coordinates3_Transform ( immCoordinates3, itransformation3, NULL ) ;
                }
                if ( QQC )
                {
                    Coordinates3_CopyTo    ( qcCoordinates3 , iqcCoordinates3 , NULL ) ;
                    Coordinates3_Transform ( iqcCoordinates3, itransformation3, NULL ) ;
                }

                /* . Find the enclosing box. */
                Coordinates3_EnclosingOrthorhombicBox ( icoordinates3, NULL, NULL, ilower, iupper ) ;
                Vector3_Add ( iupper, ilower, NULL ) ; /* . To get absolute upper coordinates rather than extents. */

                /* . Find the limits for the search. */
                SymmetryParameters_FindBoxSearchLimits ( symmetryParameters, lower, upper, ilower, iupper, &alow, &ahigh, &blow, &bhigh, &clow, &chigh ) ;

                /* . Expand the search range (for testing only). */
                if ( expandFactor > 0 )
                {
                    alow  -= expandFactor ;
                    ahigh += expandFactor ;
                    blow  -= expandFactor ;
                    bhigh += expandFactor ;
                    clow  -= expandFactor ;
                    chigh += expandFactor ;
                }

                /* . Loop over the search directions. */
	        for ( a = alow ; a <= ahigh ; a++ )
	        {
	            for ( b = blow ; b <= bhigh ; b++ )
	            {
        	        for ( c = clow ; c <= chigh ; c++ )
		        {
		            /* . Exclude self-interaction. */
		            if ( ( a == 0 ) && ( b == 0 ) && ( c == 0 ) && ( transformations->identity == t ) ) continue ;

                            /* . Set default values. */
                            scale = defaultscale ;
                            QSKIP = QTSKIP       ;

                            /* . Treat self-inverses. */
                            if ( ( tinverse >= 0 ) && ( tinverse == t ) )
                            {
                                /* . Find the inverse translation. */
                                Transformation3Container_FindInverseIntegerTranslation ( transformations, t, a, b, c, displacement, &ainverse, &binverse, &cinverse ) ;

                                /* . The inverse is in the search range. */
                                /* . Note that this is only valid for the self-inverse (otherwise the other transformations search range must be checked). */
                                if ( ( ainverse >= alow ) && ( ainverse <= ahigh ) && ( binverse >= blow ) && ( binverse <= bhigh ) && ( cinverse >= clow ) && ( cinverse <= chigh ) )
                                {
                                    /* . Pure self-inverses only occur once. */
                                    if ( ( a == ainverse ) && ( b == binverse ) && ( c == cinverse ) )
                                    {
                                        scale = 0.5e+00 ;
                                    }
                                    /* . Skip this image if its inverse will occur later - except for QC/MM cases. */
                                    else
                                    {
                                        scale = 1.0e+00 ;
                                        if ( ( a < ainverse ) || ( ( a == ainverse ) && ( b < binverse ) ) || ( ( a == ainverse ) && ( b == binverse ) && ( c < cinverse ) ) )
                                        {
                                            QSKIP = True ;
                                            if ( ! QQCMM ) continue ;
                                        }
                                    }
                                }
                                /* . The inverse is out of the search range but treat it implicitly anyway (this should never occur). */
                                else scale = 1.0e+00 ;
                            }

                            /* . Displace the coordinates. */
                            SymmetryParameters_Displacement ( symmetryParameters, a, b, c, displacement ) ;
                            Coordinates3_Translate ( icoordinates3,   displacement, NULL ) ;
                            Coordinates3_Translate ( immCoordinates3, displacement, NULL ) ;
                            Coordinates3_Translate ( iqcCoordinates3, displacement, NULL ) ;
                            Vector3_Add ( ilower, displacement, NULL ) ;
                            Vector3_Add ( iupper, displacement, NULL ) ;

                            /* . Do a final box check to skip boxes with no overlap. */
                            if ( ( Vector3_Item( ilower, 0 ) <= Vector3_Item( upper, 0 ) ) &&
                                 ( Vector3_Item( ilower, 1 ) <= Vector3_Item( upper, 1 ) ) &&
                                 ( Vector3_Item( ilower, 2 ) <= Vector3_Item( upper, 2 ) ) &&
                                 ( Vector3_Item( iupper, 0 ) >= Vector3_Item( lower, 0 ) ) &&
                                 ( Vector3_Item( iupper, 1 ) >= Vector3_Item( lower, 1 ) ) &&
                                 ( Vector3_Item( iupper, 2 ) >= Vector3_Item( lower, 2 ) ) )
			    {
                                /* . Create the pairlists. */
                                /* . MM/MM and QC/QC - one of an inverse pair only. */
                                if ( ! QSKIP )
                                {
                                    /* . MM/MM. */
                                    if ( QMM )
                                    {
                                        pairlist = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( generator, coordinates3, icoordinates3, NULL, NULL, mmSelection, mmSelection, freeSelection, freeSelection, NULL, mmGrid, mmOccupancy, NULL ) ;
                                        QOK      = ImageList_CreateImage ( (*inbmmmm), a, b, c, scale, transformations->items[t], &pairlist ) ;
                                        if ( ! QOK ) break ;
                                    }
                                    /* . QC/QC el. and LJ. */
                                    if ( QQC )
                                    {
                                        pairlist = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( generator, qcCoordinates3, iqcCoordinates3, NULL, NULL, NULL, NULL, NULL, NULL, NULL, qcGrid, qcOccupancy, NULL ) ;
                                        QOK      = ImageList_CreateImage ( (*inbqcqcel), a, b, c, scale, transformations->items[t], &pairlist ) ;
                                        if ( ! QOK ) break ;
                                        pairlist = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( generator, coordinates3, icoordinates3, NULL, NULL, qcpSelection, qcpSelection, freeSelection, freeSelection, NULL, mmGrid, mmOccupancy, NULL ) ;
                                        QOK      = ImageList_CreateImage ( (*inbqcqclj), a, b, c, scale, transformations->items[t], &pairlist ) ;
                                        if ( ! QOK ) break ;
                                    }
                                }
                                /* . QC/MM el. and LJ - all transformations. */
                                if ( QQCMM )
                                {
                                    if ( qcmmCoupling == QCMMLinkAtomCoupling_MM ) pairlist = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( generator, qcCoordinates3, icoordinates3,   NULL, NULL, NULL, mmSelection, NULL, NULL, NULL, qcGrid, qcOccupancy, NULL ) ;
                                    else                                           pairlist = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( generator, qcCoordinates3, immCoordinates3, NULL, NULL, NULL,        NULL, NULL, NULL, NULL, qcGrid, qcOccupancy, NULL ) ;
                                    QOK      = ImageList_CreateImage ( (*inbqcmmel), a, b, c, scale, transformations->items[t], &pairlist ) ;
                                    if ( ! QOK ) break ;
                                    pairlist = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( generator, coordinates3, icoordinates3, NULL, NULL, qcpSelection, mmSelection, freeSelection, freeSelection, NULL, mmGrid, mmOccupancy, NULL ) ;
                                    QOK      = ImageList_CreateImage ( (*inbqcmmlj), a, b, c, scale, transformations->items[t], &pairlist ) ;
                                    if ( ! QOK ) break ;
                                }
			    }

                            /* . Move back the coordinates. */
                            Vector3_Scale ( displacement, -1.0e+00 ) ;
                            Coordinates3_Translate ( icoordinates3,   displacement, NULL ) ;
                            Coordinates3_Translate ( immCoordinates3, displacement, NULL ) ;
                            Coordinates3_Translate ( iqcCoordinates3, displacement, NULL ) ;
                            Vector3_Add ( ilower, displacement, NULL ) ;
                            Vector3_Add ( iupper, displacement, NULL ) ;
        	        }
                        if ( ! QOK ) break ;
                    }
                    if ( ! QOK ) break ;
	        }
                if ( ! QOK ) break ;
	    }
        }

        /* . Deallocate space. */
        Coordinates3_Deallocate    ( &icoordinates3    ) ;
        Coordinates3_Deallocate    ( &immCoordinates3  ) ;
        Coordinates3_Deallocate    ( &iqcCoordinates3  ) ;
        Transformation3_Deallocate ( &itransformation3 ) ;
        Vector3_Deallocate         ( &displacement     ) ;
        Vector3_Deallocate         ( &ilower           ) ;
        Vector3_Deallocate         ( &iupper           ) ;
        Vector3_Deallocate         ( &lower            ) ;
        Vector3_Deallocate         ( &upper            ) ;

        /* . Remove lists if there is a problem or if any of them are empty. */
        if ( ! QOK || ( ( (*inbmmmm  ) != NULL ) && ( (*inbmmmm  )->nimages <= 0 ) ) ) ImageList_Deallocate ( inbmmmm   ) ;
        if ( ! QOK || ( ( (*inbqcmmel) != NULL ) && ( (*inbqcmmel)->nimages <= 0 ) ) ) ImageList_Deallocate ( inbqcmmel ) ;
        if ( ! QOK || ( ( (*inbqcmmlj) != NULL ) && ( (*inbqcmmlj)->nimages <= 0 ) ) ) ImageList_Deallocate ( inbqcmmlj ) ;
        if ( ! QOK || ( ( (*inbqcqcel) != NULL ) && ( (*inbqcqcel)->nimages <= 0 ) ) ) ImageList_Deallocate ( inbqcqcel ) ;
        if ( ! QOK || ( ( (*inbqcqclj) != NULL ) && ( (*inbqcqclj)->nimages <= 0 ) ) ) ImageList_Deallocate ( inbqcqclj ) ;

        /* . Set the status flag. */
        if ( QOK ) status = Status_Continue    ;
        else       status = Status_OutOfMemory ;
# ifdef DEBUGPRINTING
printf ( "Out of GenerateImageLists %u\n", QOK ) ;
# endif
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate the lists.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Status GenerateLists ( const PairListGenerator     *generator      ,
                                    Selection             *mmSelection    ,
                                    Selection             *qcbSelection   ,
                                    Selection             *qcpSelection   ,
                                    Selection             *freeSelection  ,
                                    PairList              *exclusions     ,
                                    PairList              *qcmmExclusions ,
                              const Coordinates3          *coordinates3   ,
                              const Coordinates3          *qcCoordinates3 ,
                              const Coordinates3          *mmCoordinates3 ,
                                    RegularGrid           *mmGrid         ,
                                    RegularGridOccupancy  *mmOccupancy    ,
                                    RegularGrid           *qcGrid         ,
                                    RegularGridOccupancy  *qcOccupancy    ,
                              const QCMMLinkAtomCoupling   qcmmCoupling   ,
                                    PairList             **nbmmmm         ,
                                    PairList             **nbqcmmlj       ,
                                    PairList             **nbqcmmel       )
{
    Status status = Status_Continue ;
    /* . Deallocate existing lists. */
    PairList_Deallocate ( nbmmmm   ) ;
    PairList_Deallocate ( nbqcmmel ) ;
    PairList_Deallocate ( nbqcmmlj ) ;
    /* . Create new ones. */
    if ( ( generator != NULL ) && ( coordinates3 != NULL ) && ( generator->cutoff > 0.0e+00 ) )
    {
        /* . MM/MM. */
        if ( nbmmmm != NULL )
        {
            (*nbmmmm) = PairListGenerator_SelfPairListFromCoordinates3 ( generator, coordinates3, NULL, mmSelection, freeSelection, exclusions, mmGrid, mmOccupancy, NULL ) ;
            if ( (*nbmmmm) == NULL ) status = Status_OutOfMemory ;
            else if ( (*nbmmmm)->npairs <= 0 ) PairList_Deallocate ( nbmmmm ) ;
        }
        /* . QC/MM electrostatic. */
        if ( qcbSelection != NULL )
        {
            /* . Electrostatic lists. */
            if ( nbqcmmel != NULL )
            {
                if ( qcmmCoupling == QCMMLinkAtomCoupling_MM ) (*nbqcmmel) = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( generator, qcCoordinates3,   coordinates3, NULL, NULL, NULL, mmSelection, NULL, NULL, qcmmExclusions, qcGrid, qcOccupancy, NULL ) ;
                else                                           (*nbqcmmel) = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( generator, qcCoordinates3, mmCoordinates3, NULL, NULL, NULL,        NULL, NULL, NULL,           NULL, qcGrid, qcOccupancy, NULL ) ;
                if ( (*nbqcmmel) == NULL ) status = Status_OutOfMemory ;
                else if ( (*nbqcmmel)->npairs <= 0 ) PairList_Deallocate ( nbqcmmel ) ;
            }
        }
        /* . QC/MM LJ. */
        if ( qcpSelection != NULL )
        {
            if ( nbqcmmlj != NULL )
            {
                (*nbqcmmlj) = PairListGenerator_CrossPairListFromSingleCoordinates3  ( generator, coordinates3, NULL, qcpSelection, mmSelection, freeSelection, exclusions, False, mmGrid, mmOccupancy, NULL ) ;
                if ( (*nbqcmmlj) == NULL ) status = Status_OutOfMemory ;
                else if ( (*nbqcmmlj)->npairs <= 0 ) PairList_Deallocate ( nbqcmmlj ) ;
            }
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate the 1-4 lists.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Status GenerateLists14 ( Selection *mmSelection, Selection *qcbSelection, Selection *qcpSelection, Selection *freeSelection,
                                                                 PairList *interactions14, const QCMMLinkAtomCoupling qcmmCoupling,
                                                                      PairList **nbmmmm, PairList **nbqcmmlj, PairList **nbqcmmel )
{
    Status status = Status_Continue ;
    /* . Deallocate existing lists. */
    PairList_Deallocate ( nbmmmm   ) ;
    PairList_Deallocate ( nbqcmmel ) ;
    PairList_Deallocate ( nbqcmmlj ) ;
    /* . Create new ones. */
    if ( ( interactions14 != NULL ) && ( interactions14->npairs > 0 ) )
    {
        /* . MM/MM. */
        if ( nbmmmm != NULL )
        {
            (*nbmmmm) = SelfPairList_FromSelfPairList ( interactions14, mmSelection, freeSelection, False ) ;
            if ( (*nbmmmm) == NULL ) status = Status_OutOfMemory ;
            else if ( (*nbmmmm)->npairs <= 0 ) PairList_Deallocate ( nbmmmm ) ;
        }
        /* . QC/MM electrostatic. */
        if ( ( qcbSelection != NULL ) && ( qcmmCoupling == QCMMLinkAtomCoupling_MM ) )
        {
            /* . Electrostatic lists. */
            if ( nbqcmmel != NULL )
            {
                (*nbqcmmel) = SelfPairList_ToCrossPairList  ( interactions14, qcbSelection, mmSelection, NULL, True, False, True, False ) ;
                if ( (*nbqcmmel) == NULL ) status = Status_OutOfMemory ;
                else if ( (*nbqcmmel)->npairs <= 0 ) PairList_Deallocate ( nbqcmmel ) ;
            }
        }
        /* . QC/MM LJ. */
        if ( qcpSelection != NULL )
        {
            if ( nbqcmmlj != NULL )
            {
                (*nbqcmmlj) = SelfPairList_ToCrossPairList  ( interactions14, qcpSelection, mmSelection, freeSelection, False, False, False, False ) ;
                if ( (*nbqcmmlj) == NULL ) status = Status_OutOfMemory ;
                else if ( (*nbqcmmlj)->npairs <= 0 ) PairList_Deallocate ( nbqcmmlj ) ;
            }
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . MM/MM image energy.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Tests by Nikolay indicate that parallelization over images rather than interactions is more efficient. */
static void MMMMImageEnergy ( const PairwiseInteractionABFS    *mmmmPairwiseInteraction    ,
                              const MMAtomContainer            *mmAtoms                    ,
                              const LJParameterContainer       *ljParameters               ,
                                    SymmetryParameters         *symmetryParameters         ,
                                    ImageList                  *imageList                  ,
                              const Real                        electrostaticScale         ,
                              const Coordinates3               *coordinates3               ,
                                    Real                       *eElectrostatic             ,
                                    Real                       *eLennardJones              ,
# ifdef USEOPENMP
                              const Integer                     numberOfThreads            ,
                                    Coordinates3              **gradientsArray             ,
# else
                                    Coordinates3               *gradients3                 ,
# endif
                                    SymmetryParameterGradients *symmetryParameterGradients )
{
    Boolean doElectrostatics, doLennardJones ;
    doElectrostatics = ( eElectrostatic != NULL ) ;
    doLennardJones   = ( eLennardJones  != NULL ) ;
    if ( doElectrostatics ) (*eElectrostatic) = 0.0e+00 ;
    if ( doLennardJones   ) (*eLennardJones ) = 0.0e+00 ;
    if ( ( mmmmPairwiseInteraction != NULL ) &&
         ( mmAtoms                 != NULL ) &&
         ( symmetryParameters      != NULL ) &&
         ( imageList               != NULL ) &&
         ( coordinates3            != NULL ) &&
         ( doElectrostatics || doLennardJones ) && ( ( ljParameters != NULL ) || ( ! doLennardJones && ( ljParameters == NULL ) ) ) )
    {
        auto Boolean doGradients ;
        auto Real    eLJ, eQQ ;

        /* . Initialization. */
# ifdef USEOPENMP
        doGradients = ( gradientsArray != NULL ) && ( symmetryParameterGradients != NULL ) ;
# else
        doGradients = ( gradients3     != NULL ) && ( symmetryParameterGradients != NULL ) ;
# endif
        eLJ = 0.0e+00 ;
        eQQ = 0.0e+00 ;
        ImageList_MakeRecords ( imageList ) ;

# ifdef USEOPENMP
	#pragma omp parallel num_threads(numberOfThreads) reduction(+:eLJ,eQQ)
# endif
	{
            auto Integer          r ;
            auto Real             eLocalLJ, eLocalQQ, *eTempLJ = NULL, *eTempQQ = NULL, scale ;
            auto Coordinates3    *icoordinates3    = NULL, *igradients3 = NULL, *localGradients3 = NULL ;
# ifdef USEOPENMP
            auto Coordinates3    *gradients1[1], *gradients2[1], **gradientsArray1 = NULL, **gradientsArray2 = NULL ;
# endif
            auto Image           *image            = NULL ;
            auto Transformation3 *itransformation3 = NULL, *xtransformation3 = NULL ;

            /* . Allocation. */
            icoordinates3    = Coordinates3_Allocate ( coordinates3->length0 ) ;
            itransformation3 = Transformation3_Allocate ( ) ;
            xtransformation3 = Transformation3_Allocate ( ) ;
            if ( doElectrostatics ) eTempQQ = &eLocalQQ ;
            if ( doLennardJones   ) eTempLJ = &eLocalLJ ;
            if ( doGradients )
            {
# ifdef USEOPENMP
                auto Integer thisThread = omp_get_thread_num ( ) ;
# endif
                igradients3     = Coordinates3_Allocate ( coordinates3->length0 ) ;
# ifdef USEOPENMP
                localGradients3 = gradientsArray[thisThread] ;
                gradients1[0]   = localGradients3 ; gradientsArray1 = gradients1 ;
                gradients2[0]   = igradients3     ; gradientsArray2 = gradients2 ;
# else
                localGradients3 = gradients3 ;
# endif
            }

            /* . Loop over the images in the list. */
# ifdef USEOPENMP
	    #pragma omp for schedule(dynamic,1)
# endif
            for ( r = 0 ; r < imageList->numberOfRecords ; r++ )
            {
                /* . Data for the image. */
                image = imageList->records[r] ;
                scale = image->scale ;

                /* . Generate the image transformation in real space. */
                Transformation3_Copy          ( xtransformation3, image->transformation3 ) ;
                xtransformation3->translation->data[0] += ( Real ) image->a ;
                xtransformation3->translation->data[1] += ( Real ) image->b ;
                xtransformation3->translation->data[2] += ( Real ) image->c ;
                Transformation3_Copy          ( itransformation3, xtransformation3 ) ;
                Transformation3_Orthogonalize ( itransformation3, symmetryParameters->M, symmetryParameters->inverseM ) ;

                /* . Generate image coordinates. */
                Coordinates3_CopyTo    ( coordinates3 , icoordinates3   , NULL ) ;
                Coordinates3_Transform ( icoordinates3, itransformation3, NULL ) ;

                /* . Initialize image gradients. */
                if ( doGradients ) Coordinates3_Set ( igradients3, 0.0e+00 ) ;

                /* . Calculate the energy. */
                PairwiseInteractionABFS_MMMMEnergy ( mmmmPairwiseInteraction     ,
                                                     mmAtoms                     ,
                                                     ljParameters                ,
                                                     image->pairlist             ,
                                                     electrostaticScale * scale  ,
                                                     scale                       ,
                                                     coordinates3                ,
                                                     icoordinates3               ,
                                                     eTempQQ                     ,
                                                     eTempLJ                     ,
# ifdef USEOPENMP
                                                     1                           ,
                                                     gradientsArray1             ,
                                                     gradientsArray2             ) ;
# else
                                                     gradients3                  ,
                                                     igradients3                 ) ;
# endif

                /* . Save the energies. */
                if ( doElectrostatics ) eQQ += (*eTempQQ) ;
                if ( doLennardJones   ) eLJ += (*eTempLJ) ;

                /* . Treat the gradients. */
                if ( doGradients )
                {
                    /* . Symmetry parameter gradients - done here before igradients3 is changed. */
# ifdef USEOPENMP
		    #pragma omp critical
# endif
                    SymmetryParameterGradients_ImageDerivatives ( symmetryParameterGradients, symmetryParameters, xtransformation3, coordinates3, igradients3 ) ;
                    /* . Transpose the rotation in place as it is no longer needed. */
                    Matrix33_Transpose ( itransformation3->rotation, NULL ) ;
                    /* . Gradients. */
                    Coordinates3_Rotate         ( igradients3, itransformation3->rotation, NULL ) ;
                    Coordinates3_AddScaledArray ( localGradients3 , 1.0e+00, igradients3 , NULL ) ;
                }
            }

            /* . Finish up. */
            Coordinates3_Deallocate    ( &icoordinates3    ) ;
            Transformation3_Deallocate ( &itransformation3 ) ;
            Transformation3_Deallocate ( &xtransformation3 ) ;
            if ( doGradients ) Coordinates3_Deallocate ( &igradients3 ) ;
        }

        /* . Finish up. */
        if ( doElectrostatics ) (*eElectrostatic) += eQQ ;
        if ( doLennardJones   ) (*eLennardJones ) += eLJ ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM image electrostatic gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void QCMMImageGradients ( const PairwiseInteractionABFS    *qcmmPairwiseInteraction    ,
                                 const Real1DArray                *mmCharges                  ,
                                 const QCAtomContainer            *qcAtoms                    ,
                                       SymmetryParameters         *symmetryParameters         ,
                                       ImageList                  *imageList                  ,
                                 const Real                        electrostaticScale         ,
                                 const Coordinates3               *coordinates3               ,
                                 const Coordinates3               *mmCoordinates3             ,
                                 const Real1DArray                *qcCharges                  ,
                                       Coordinates3               *gradients3                 ,
                                       Coordinates3               *mmGradients3               ,
                                       SymmetryParameterGradients *symmetryParameterGradients )
{
    if ( ( qcmmPairwiseInteraction != NULL ) && ( mmCharges != NULL ) && ( qcAtoms != NULL ) && ( symmetryParameters != NULL ) && ( imageList != NULL ) && ( electrostaticScale != 0.0e+00 ) &&
             ( coordinates3 != NULL )  && ( mmCoordinates3 != NULL ) && ( qcCharges != NULL ) && ( gradients3 != NULL ) && ( mmGradients3 != NULL ) && ( symmetryParameterGradients != NULL ) )
    {
        auto Real scale ;
        auto Coordinates3    *icoordinates3    = NULL, *igradients3 = NULL ;
        auto Image           *image            = NULL ;
        auto Transformation3 *itransformation3 = NULL, *xtransformation3 = NULL ;

        /* . Allocate some empty data structures. */
        icoordinates3    = Coordinates3_Allocate ( mmCoordinates3->length0 ) ;
        igradients3      = Coordinates3_Allocate ( mmCoordinates3->length0 ) ;
        itransformation3 = Transformation3_Allocate ( ) ;
        xtransformation3 = Transformation3_Allocate ( ) ;

        /* . Loop over the images in the list. */
        List_Iterate_Initialize ( imageList->images ) ;
        while ( ( image = ImageList_Iterate ( imageList ) ) != NULL )
        {
            /* . Data for the image. */;
            scale = image->scale ;

            /* . Generate the image transformation in real space. */
            Transformation3_Copy          ( xtransformation3, image->transformation3 ) ;
            xtransformation3->translation->data[0] += ( Real ) image->a ;
            xtransformation3->translation->data[1] += ( Real ) image->b ;
            xtransformation3->translation->data[2] += ( Real ) image->c ;
            Transformation3_Copy          ( itransformation3, xtransformation3 ) ;
            Transformation3_Orthogonalize ( itransformation3, symmetryParameters->M, symmetryParameters->inverseM ) ;

            /* . Generate image coordinates. */
            Coordinates3_CopyTo    ( mmCoordinates3, icoordinates3   , NULL ) ;
            Coordinates3_Transform ( icoordinates3 , itransformation3, NULL ) ;

            /* . Initialize image gradients. */
            Coordinates3_Set ( igradients3, 0.0e+00 ) ;

            /* . Calculate the gradients. */
            PairwiseInteractionABFS_QCMMGradients ( qcmmPairwiseInteraction    ,
                                                    mmCharges                  ,
                                                    qcAtoms                    ,
                                                    image->pairlist            ,
                                                    electrostaticScale * scale ,
                                                    coordinates3               ,
                                                    icoordinates3              ,
                                                    qcCharges                  ,
                                                    gradients3                 ,
                                                    igradients3                ) ;

            /* . Symmetry parameter gradients - done here before igradients3 is changed. */
            SymmetryParameterGradients_ImageDerivatives ( symmetryParameterGradients, symmetryParameters, xtransformation3, mmCoordinates3, igradients3 ) ;

            /* . Transpose the rotation in place as it is no longer needed. */
            Matrix33_Transpose ( itransformation3->rotation, NULL ) ;

            /* . Gradients. */
            Coordinates3_Rotate         ( igradients3 ,  itransformation3->rotation, NULL ) ;
            Coordinates3_AddScaledArray ( mmGradients3, 1.0e+00, igradients3       , NULL ) ;
        }

        /* . Finish up. */
        Coordinates3_Deallocate    ( &icoordinates3    ) ;
        Coordinates3_Deallocate      ( &igradients3      ) ;
        Transformation3_Deallocate ( &itransformation3 ) ;
        Transformation3_Deallocate ( &xtransformation3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM image potentials.
! . Potentials is incremented not reset.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void QCMMImagePotentials ( const PairwiseInteractionABFS *qcmmPairwiseInteraction ,
                                  const Real1DArray             *mmCharges               ,
                                  const QCAtomContainer         *qcAtoms                 ,
                                        SymmetryParameters      *symmetryParameters      ,
                                        ImageList               *imageList               ,
                                  const Real                     electrostaticScale      ,
                                  const Coordinates3            *coordinates3            ,
                                  const Coordinates3            *mmCoordinates3          ,
                                        Real1DArray             *potentials              )
{
    if ( ( qcmmPairwiseInteraction != NULL ) && ( mmCharges != NULL ) && ( qcAtoms != NULL ) && ( symmetryParameters != NULL ) && ( imageList != NULL ) &&
                                   ( electrostaticScale != 0.0e+00 ) && ( coordinates3 != NULL ) && ( mmCoordinates3 != NULL ) && ( potentials != NULL ) )
    {
        auto Real scale ;
        auto Coordinates3    *icoordinates3    = NULL ;
        auto Image           *image            = NULL ;
        auto Transformation3 *itransformation3 = NULL, *xtransformation3 = NULL ;

        /* . Allocate some empty data structures. */
        icoordinates3    = Coordinates3_Allocate ( mmCoordinates3->length0 ) ;
        itransformation3 = Transformation3_Allocate ( ) ;
        xtransformation3 = Transformation3_Allocate ( ) ;

        /*----------------------------------------------------------------------
        ! . Loop over the images in the list.
        !---------------------------------------------------------------------*/
        List_Iterate_Initialize ( imageList->images ) ;
        while ( ( image = ImageList_Iterate ( imageList ) ) != NULL )
        {
            /* . Data for the image. */;
            scale = image->scale ;

            /* . Generate the image transformation in real space. */
            Transformation3_Copy          ( xtransformation3, image->transformation3 ) ;
            xtransformation3->translation->data[0] += ( Real ) image->a ;
            xtransformation3->translation->data[1] += ( Real ) image->b ;
            xtransformation3->translation->data[2] += ( Real ) image->c ;
            Transformation3_Copy          ( itransformation3, xtransformation3 ) ;
            Transformation3_Orthogonalize ( itransformation3, symmetryParameters->M, symmetryParameters->inverseM ) ;

            /* . Generate image coordinates. */
            Coordinates3_CopyTo    ( mmCoordinates3, icoordinates3   , NULL ) ;
            Coordinates3_Transform ( icoordinates3 , itransformation3, NULL ) ;

            /* . Calculate the potentials. */
            PairwiseInteractionABFS_QCMMPotentials ( qcmmPairwiseInteraction    ,
                                                     mmCharges                  ,
                                                     qcAtoms                    ,
                                                     image->pairlist            ,
                                                     electrostaticScale * scale ,
                                                     coordinates3               ,
                                                     icoordinates3              ,
                                                     potentials                 ) ;
        }

        /* . Finish up. */
        Coordinates3_Deallocate    ( &icoordinates3    ) ;
        Transformation3_Deallocate ( &itransformation3 ) ;
        Transformation3_Deallocate ( &xtransformation3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/QC image electrostatic gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void QCQCImageGradients ( const PairwiseInteractionABFS    *qcqcPairwiseInteraction    ,
                                 const QCAtomContainer            *qcAtoms                    ,
                                       SymmetryParameters         *symmetryParameters         ,
                                       ImageList                  *imageList                  ,
                                 const Real                        electrostaticScale         ,
                                 const Coordinates3               *coordinates3               ,
                                 const Coordinates3               *qcCoordinates3             ,
                                 const Real1DArray                *qcCharges                  ,
                                       Coordinates3               *gradients3                 ,
                                       SymmetryParameterGradients *symmetryParameterGradients )
{
    if ( ( qcqcPairwiseInteraction != NULL ) && ( qcAtoms != NULL ) && ( symmetryParameters != NULL ) && ( imageList != NULL ) && ( electrostaticScale != 0.0e+00 ) &&
                 ( coordinates3 != NULL ) && ( qcCoordinates3 != NULL ) && ( qcCharges != NULL ) && ( gradients3 != NULL ) && ( symmetryParameterGradients != NULL ) )
    {
        auto Real scale ;
        auto Coordinates3    *iqcCoordinates3  = NULL, *iqcgradients3 = NULL, *qcgradients3 = NULL ;
        auto Image           *image            = NULL ;
        auto Transformation3 *itransformation3 = NULL, *xtransformation3 = NULL ;

        /* . Allocate some empty data structures. */
        iqcCoordinates3  = Coordinates3_Allocate ( qcCoordinates3->length0 ) ;
        iqcgradients3    = Coordinates3_Allocate ( qcCoordinates3->length0 ) ;
        qcgradients3     = Coordinates3_Allocate ( qcCoordinates3->length0 ) ;
        itransformation3 = Transformation3_Allocate ( ) ;
        xtransformation3 = Transformation3_Allocate ( ) ;

        /* . Initialization. */
        Coordinates3_Set ( qcgradients3, 0.0e+00 ) ;

        /*----------------------------------------------------------------------
        ! . Loop over the images in the list.
        !---------------------------------------------------------------------*/
        List_Iterate_Initialize ( imageList->images ) ;
        while ( ( image = ImageList_Iterate ( imageList ) ) != NULL )
        {
            /* . Data for the image. */;
            scale = image->scale ;

            /* . Generate the image transformation in real space. */
            Transformation3_Copy          ( xtransformation3, image->transformation3 ) ;
            xtransformation3->translation->data[0] += ( Real ) image->a ;
            xtransformation3->translation->data[1] += ( Real ) image->b ;
            xtransformation3->translation->data[2] += ( Real ) image->c ;
            Transformation3_Copy          ( itransformation3, xtransformation3 ) ;
            Transformation3_Orthogonalize ( itransformation3, symmetryParameters->M, symmetryParameters->inverseM ) ;

            /* . Generate image coordinates. */
            Coordinates3_CopyTo    ( qcCoordinates3 , iqcCoordinates3 , NULL ) ;
            Coordinates3_Transform ( iqcCoordinates3, itransformation3, NULL ) ;

            /* . Initialize image gradients. */
            Coordinates3_Set ( iqcgradients3, 0.0e+00 ) ;

            /* . Calculate the gradients. */
            PairwiseInteractionABFS_QCQCGradients ( qcqcPairwiseInteraction    ,
                                                    image->pairlist            ,
                                                    electrostaticScale * scale ,
                                                    qcCoordinates3             ,
                                                    iqcCoordinates3            ,
                                                    qcCharges                  ,
                                                    qcgradients3               ,
                                                    iqcgradients3              ) ;

            /* . Symmetry parameter gradients - done here before igradients3 is changed. */
            SymmetryParameterGradients_ImageDerivatives ( symmetryParameterGradients, symmetryParameters, xtransformation3, qcCoordinates3, iqcgradients3 ) ;

            /* . Transpose the rotation in place as it is no longer needed. */
            Matrix33_Transpose ( itransformation3->rotation, NULL ) ;

            /* . Gradients. */
            Coordinates3_Rotate         ( iqcgradients3, itransformation3->rotation, NULL ) ;
            Coordinates3_AddScaledArray ( qcgradients3 , 1.0e+00, iqcgradients3    , NULL ) ;
        }

        /* . Put the QC gradients into the full gradients. */
        QCAtomContainer_SetGradients3 ( qcAtoms, coordinates3, qcgradients3, False, &gradients3 ) ;

        /* . Finish up. */
        Coordinates3_Deallocate    ( &iqcCoordinates3  ) ;
        Coordinates3_Deallocate    ( &iqcgradients3    ) ;
        Coordinates3_Deallocate    ( &qcgradients3     ) ;
        Transformation3_Deallocate ( &itransformation3 ) ;
        Transformation3_Deallocate ( &xtransformation3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/QC image potentials.
! . Potentials is incremented not reset.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void QCQCImagePotentials ( const PairwiseInteractionABFS *qcqcPairwiseInteraction ,
                                        SymmetryParameters      *symmetryParameters      ,
                                        ImageList               *imageList               ,
                                  const Real                     electrostaticScale      ,
                                  const Coordinates3            *qcCoordinates3          ,
                                        SymmetricMatrix         *potentials              )
{
    if ( ( qcqcPairwiseInteraction != NULL ) && ( symmetryParameters != NULL ) && ( imageList != NULL ) && ( electrostaticScale != 0.0e+00 ) && ( qcCoordinates3 != NULL )  && ( potentials != NULL ) )
    {
        auto Real                scale ;
        auto Coordinates3    *iqcCoordinates3  = NULL ;
        auto Image           *image            = NULL ;
        auto Transformation3 *itransformation3 = NULL, *xtransformation3 = NULL ;

        /* . Allocate some empty data structures. */
        iqcCoordinates3  = Coordinates3_Allocate ( qcCoordinates3->length0 ) ;
        itransformation3 = Transformation3_Allocate ( ) ;
        xtransformation3 = Transformation3_Allocate ( ) ;

        /*----------------------------------------------------------------------
        ! . Loop over the images in the list.
        !---------------------------------------------------------------------*/
        List_Iterate_Initialize ( imageList->images ) ;
        while ( ( image = ImageList_Iterate ( imageList ) ) != NULL )
        {
            /* . Data for the image. */;
            scale = image->scale ;

            /* . Generate the image transformation in real space. */
            Transformation3_Copy          ( xtransformation3, image->transformation3 ) ;
            xtransformation3->translation->data[0] += ( Real ) image->a ;
            xtransformation3->translation->data[1] += ( Real ) image->b ;
            xtransformation3->translation->data[2] += ( Real ) image->c ;
            Transformation3_Copy          ( itransformation3, xtransformation3 ) ;
            Transformation3_Orthogonalize ( itransformation3, symmetryParameters->M, symmetryParameters->inverseM ) ;

            /* . Generate image coordinates. */
            Coordinates3_CopyTo    ( qcCoordinates3 , iqcCoordinates3 , NULL ) ;
            Coordinates3_Transform ( iqcCoordinates3, itransformation3, NULL ) ;

            /* . Calculate the potentials. */
            PairwiseInteractionABFS_QCQCPotentials ( qcqcPairwiseInteraction    ,
                                                     image->pairlist            ,
                                                     electrostaticScale * scale ,
                                                     qcCoordinates3             ,
                                                     iqcCoordinates3            ,
                                                     potentials                 ) ;
        }

        /* . Finish up. */
        Coordinates3_Deallocate    ( &iqcCoordinates3  ) ;
        Transformation3_Deallocate ( &itransformation3 ) ;
        Transformation3_Deallocate ( &xtransformation3 ) ;
    }
}
