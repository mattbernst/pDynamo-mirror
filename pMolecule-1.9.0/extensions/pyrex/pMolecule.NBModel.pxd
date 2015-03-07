#-------------------------------------------------------------------------------
# . File      : pMolecule.NBModel.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.Coordinates3             cimport Coordinates3
from pCore.PairList                 cimport CrossPairList, SelfPairList
from pCore.Selection                cimport Selection
from pMolecule.LJParameterContainer cimport LJParameterContainer
from pMolecule.MMAtomContainer      cimport MMAtomContainer
from pMolecule.QCAtomContainer      cimport QCAtomContainer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NBModel:

    pass
