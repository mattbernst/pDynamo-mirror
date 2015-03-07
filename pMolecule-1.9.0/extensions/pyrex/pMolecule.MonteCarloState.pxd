#-------------------------------------------------------------------------------
# . File      : pMolecule.MonteCarloState.pxd
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions               cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Coordinates3               cimport CCoordinates3, Coordinates3
from pCore.Matrix33                   cimport CMatrix33
from pCore.SelectionContainer         cimport CSelectionContainer
from pCore.Status                     cimport Status
from pCore.Vector3                    cimport CVector3
from pMolecule.NBModelMonteCarlo      cimport CNBModelMonteCarlo
from pMolecule.NBModelMonteCarloState cimport CNBModelMonteCarloState
from pMolecule.SymmetryParameters     cimport CSymmetryParameters

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "MonteCarloState.h":

    ctypedef struct CMonteCarloState "MonteCarloState":
        Integer                  blocks
        Integer                  moves
        Integer                  nreject
        Integer                  nrejectm
        Integer                  nrejectt
        Integer                  nrejectv
        Integer                  ntrym
        Integer                  ntryv
        Real                     beta
        Real                     ecurrent
        Real                     pressure
        Real                     tfact
        Real                     volume
        Real                     acceptanceratio
        Real                     rmax
        Real                     tmax
        Real                     vmax
        Real                     eav
        Real                     eav2
        Real                     etot
        Real                     etot2
        Real                     etotb
        Real                     etotb2
        Real                     hav
        Real                     hav2
        Real                     htot
        Real                     htot2
        Real                     htotb
        Real                     htotb2
        Real                     vav
        Real                     vav2
        Real                     vtot
        Real                     vtot2
        Real                     vtotb
        Real                     vtotb2
        Real                    *random
        CCoordinates3           *oldcoordinates3
        CMatrix33               *rotation
        CSymmetryParameters     *oldsymmetryParameters
        CVector3                *translation
        CCoordinates3           *coordinates3
        CNBModelMonteCarlo      *nbModel
        CNBModelMonteCarloState *nbState
        CSelectionContainer     *isolates
        CSymmetryParameters     *symmetryParameters

    cdef CMonteCarloState *MonteCarloState_Allocate                  ( Integer nparticles, Integer nrandom )
    cdef void              MonteCarloState_Deallocate                ( CMonteCarloState **self )
    cdef void              MonteCarloState_AdjustMoveSizes           ( CMonteCarloState  *self )
    cdef Status            MonteCarloState_MoveIsolate               ( CMonteCarloState  *self )
    cdef Status            MonteCarloState_MoveVolume                ( CMonteCarloState  *self )
    cdef void              MonteCarloState_StatisticsBlockAccumulate ( CMonteCarloState  *self )
    cdef void              MonteCarloState_StatisticsBlockStart      ( CMonteCarloState  *self )
    cdef void              MonteCarloState_StatisticsBlockStop       ( CMonteCarloState  *self )
    cdef void              MonteCarloState_StatisticsStop            ( CMonteCarloState  *self )
    cdef void              MonteCarloState_StatisticsStart           ( CMonteCarloState  *self )
