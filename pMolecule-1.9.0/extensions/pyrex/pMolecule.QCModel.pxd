#-------------------------------------------------------------------------------
# . File      : pMolecule.QCModel.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.Coordinates3         cimport Coordinates3
from pCore.Vector3              cimport Vector3
from pMolecule.DIISSCFConverger cimport DIISSCFConverger
from pMolecule.QCAtomContainer  cimport QCAtomContainer
from pMolecule.QCParameters     cimport QCParameters

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class QCModel:

    cdef public object converger
    cdef public object label
