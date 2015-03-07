/*------------------------------------------------------------------------------
! . File      : ImageList.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
!=================================================================================================================================*/

# include "Memory.h"
# include "ImageList.h"

/*==================================================================================================================================
! . ImageList procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ImageList *ImageList_Allocate ( void )
{
    ImageList *self = NULL ;
    self = ( ImageList * ) Memory_Allocate ( sizeof ( ImageList ) ) ;
    self->nimages         =  0 ;
    self->npairs          =  0 ;
    self->numberOfRecords = -1 ;
    self->images  = List_Allocate ( ) ;
    self->images->Element_Deallocate = Image_Deallocate ;
    self->records         = NULL ;
   return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create an image.
! . This is only done if the pairlist is non-NULL and non-empty and the
! . transformation is non-NULL.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ImageList_CreateImage ( ImageList *self, const Integer a, const Integer b, const Integer c, const Real scale, Transformation3 *transformation3, PairList **pairlist )
{
    Boolean QSUCCESS = False ;
    if ( ( self != NULL ) && ( (*pairlist) != NULL ) && ( transformation3 != NULL ) )
    {
        /* . So far, so good. */
        QSUCCESS = True ;

        /* . Save the image if there are interactions. */
        if ( (*pairlist)->npairs > 0 )
        {
            auto Image *image = NULL ;

            /* . Allocate an image. */
            image = Image_Allocate ( ) ;
            if ( image == NULL ) QSUCCESS = False;
            else
            {
                /* . Save the data. */
                image->a               = a ;
                image->b               = b ;
                image->c               = c ;
                image->pairlist        = (*pairlist)     ;
                image->scale           = scale           ;
                image->transformation3 = transformation3 ;
                List_Element_Append ( self->images, ( void * ) image ) ;
                self->nimages ++ ;
                self->npairs  += PairList_Length ( (*pairlist) ) ;

                /* . Reset pairlist. */
                (*pairlist) = NULL ;
            }
        }
        else PairList_Deallocate ( pairlist ) ;

    }
    return QSUCCESS ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImageList_Deallocate ( ImageList **self )
{
    if ( (*self) != NULL )
    {
        /* . Representations. */
        free ( (*self)->records ) ;
        (*self)->numberOfRecords =   -1 ;
        (*self)->records         = NULL ;
        /* . Standard data. */
        List_Deallocate ( &((*self)->images) ) ;
        free ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Iterate over the images in a list.
!---------------------------------------------------------------------------------------------------------------------------------*/
Image *ImageList_Iterate ( ImageList *self )
{
    if ( self == NULL ) return NULL ;
    else                return ( Image * ) List_Iterate ( self->images ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the records representation of the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImageList_MakeRecords ( ImageList *self )
{
    if ( ( self != NULL ) && ( self->records == NULL ) )
    {
        auto Integer n ;
        n = ImageList_NumberOfRecords ( self ) ;
        MEMORY_ALLOCATEARRAY ( self->records, n, Image * ) ;
        if ( self->records != NULL )
        {
            auto Image *record ;
            n = 0 ;
            List_Iterate_Initialize ( self->images ) ;
            while ( ( record = ImageList_Iterate ( self ) ) != NULL ) { self->records[n] = record ; n += 1 ; }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return data about the image list.
!---------------------------------------------------------------------------------------------------------------------------------*/
int ImageList_NumberOfImages ( const ImageList *self )
{
    if ( self == NULL ) return 0 ;
    else                return self->nimages ;
}

int ImageList_NumberOfPairs ( const ImageList *self )
{
    if ( self == NULL ) return 0 ;
    else                return self->npairs ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the number of records in the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ImageList_NumberOfRecords ( ImageList *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        if ( self->numberOfRecords >= 0 ) n = self->numberOfRecords ;
        else
        {
            auto Image *record ;
            List_Iterate_Initialize ( self->images ) ;
            while ( ( record = ImageList_Iterate ( self ) ) != NULL ) n += 1 ;
            self->numberOfRecords = n ;
        }
    }
    return n ;
}

/*==================================================================================================================================
! . Image procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate an image.
!---------------------------------------------------------------------------------------------------------------------------------*/
Image *Image_Allocate ( void )
{
    Image *self ;
    self = ( Image * ) malloc ( sizeof ( Image ) ) ;
    self->a               = 0 ;
    self->b               = 0 ;
    self->c               = 0 ;
    self->QOWNER          = False   ;
    self->scale           = 0.0e+00 ;
    self->pairlist        = NULL    ;
    self->transformation3 = NULL    ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocate an image.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Image_Deallocate ( void *vimage )
{
    Image *self ;
    self = ( Image * ) vimage ;
    if ( self != NULL )
    {
        PairList_Deallocate ( &(self->pairlist) ) ;
        if ( self->QOWNER ) Transformation3_Deallocate ( &(self->transformation3) ) ;
        free ( self ) ;
    }
}
