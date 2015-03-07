#-------------------------------------------------------------------------------
# . File      : pMolecule.ChargeConstraintContainer.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions   cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Integer1DArray cimport CInteger1DArray, Integer1DArray_GetItem, Integer1DArray_SetItem
from pCore.Real1DArray    cimport CReal1DArray, Real1DArray, Real1DArray_GetItem, Real1DArray_SetItem
from pCore.Selection      cimport CSelection, Selection
from pCore.Status         cimport Status, Status_Continue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "ChargeConstraintContainer.h":

    ctypedef struct CChargeConstraint "ChargeConstraint":
        Integer          numberOfCharges
        Integer          numberOfSpins  
        Real             target         
        CInteger1DArray *chargeIndices  
        CInteger1DArray *spinIndices    
        CReal1DArray    *chargeWeights  
        CReal1DArray    *spinWeights    

    ctypedef struct CChargeConstraintContainer "ChargeConstraintContainer":
        Boolean             hasCharges         
        Boolean             hasSpins           
        Integer             highestIndex       
        Integer             numberOfConstraints
        CChargeConstraint **constraints        

    cdef CChargeConstraintContainer *ChargeConstraintContainer_Allocate     ( Integer numberOfConstraints, Status *status )
    cdef CChargeConstraintContainer *ChargeConstraintContainer_Clone        ( CChargeConstraintContainer  *self, Status *status )
    cdef void                        ChargeConstraintContainer_Deallocate   ( CChargeConstraintContainer **self )
    cdef void                        ChargeConstraintContainer_Deviations   ( CChargeConstraintContainer  *self, CReal1DArray *charges, CReal1DArray *spins, CReal1DArray *deviations, Status *status )
    cdef Integer                     ChargeConstraintContainer_HighestIndex ( CChargeConstraintContainer  *self )
    cdef CChargeConstraintContainer *ChargeConstraintContainer_Prune        ( CChargeConstraintContainer  *self, CSelection *selection, Status *status )
    cdef void                        ChargeConstraintContainer_SetItem      ( CChargeConstraintContainer  *self, Integer index, CChargeConstraint **constraint, Status *status )
    cdef Integer                     ChargeConstraintContainer_Size         ( CChargeConstraintContainer  *self )

    cdef CChargeConstraint *ChargeConstraint_Allocate   ( Integer numberOfCharges, Integer numberOfSpins, Status *status )
    cdef CChargeConstraint *ChargeConstraint_Clone      ( CChargeConstraint  *self, Status *status )
    cdef void               ChargeConstraint_Deallocate ( CChargeConstraint **self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ChargeConstraintContainer:

    cdef CChargeConstraintContainer *cObject
    cdef public object               isOwner

