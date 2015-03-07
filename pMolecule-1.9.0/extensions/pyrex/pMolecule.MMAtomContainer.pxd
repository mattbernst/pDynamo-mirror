#-------------------------------------------------------------------------------
# . File      : pMolecule.MMAtomContainer.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions             cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3             cimport CCoordinates3, Coordinates3
from pCore.Real1DArray              cimport Real1DArray
from pCore.Selection                cimport Selection, CSelection
from pCore.Vector3                  cimport Vector3, CVector3
from pMolecule.LJParameterContainer cimport LJParameterContainer

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "MMAtomContainer.h":

    ctypedef struct CMMAtom "MMAtom":
        Boolean QACTIVE
        Integer atomtype
        Integer ljtype
        Real    charge

    ctypedef struct CMMAtomContainer "MMAtomContainer":
       Integer      natoms
       CMMAtom *data

    cdef CMMAtomContainer *MMAtomContainer_Allocate            ( Integer natoms )
    cdef CMMAtomContainer *MMAtomContainer_Clone               ( CMMAtomContainer  *self )
    cdef void              MMAtomContainer_Deallocate          ( CMMAtomContainer **self )
    cdef CVector3         *MMAtomContainer_DipoleMoment        ( CMMAtomContainer  *self, CCoordinates3 *coordinates3, CVector3 *center )
    cdef CMMAtomContainer *MMAtomContainer_Merge               ( CMMAtomContainer  *self, CMMAtomContainer *other, Integer atomtypeincrement,
                                                                						    Integer ljtypeincrement )
    cdef Integer           MMAtomContainer_NumberOfActiveAtoms ( CMMAtomContainer  *self )
    cdef CMMAtomContainer *MMAtomContainer_Prune               ( CMMAtomContainer  *self, CSelection *selection )
    cdef Integer           MMAtomContainer_Size                ( CMMAtomContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MMAtomContainer:

    cdef CMMAtomContainer *cObject
    cdef public object     isOwner
    cdef public object     atomTypes
