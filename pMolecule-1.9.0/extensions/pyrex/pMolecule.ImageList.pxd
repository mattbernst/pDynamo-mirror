#-------------------------------------------------------------------------------
# . File      : pMolecule.ImageList.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions cimport Boolean, CFalse, CTrue, Integer, Real

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "ImageList.h":

    ctypedef struct CImageList "ImageList":
        pass

    cdef CImageList *ImageList_Allocate       ( )
    cdef void        ImageList_Deallocate     ( CImageList **self )
    cdef Integer     ImageList_NumberOfImages ( CImageList  *self )
    cdef Integer     ImageList_NumberOfPairs  ( CImageList  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ImageList:

    cdef CImageList    *cObject
    cdef public object  isOwner
