#-------------------------------------------------------------------------------
# . File      : pCore.List.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "List.h":

    ctypedef struct CList "List":
        pass

    cdef void List_Iterate_Initialize ( CList *list )
