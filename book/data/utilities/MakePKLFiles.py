"""Make the pkl files for Examples 25 and 26."""

from pBabel           import MOLFile_ToSystem, XYZFile_ToSystem
from pCore            import Pickle
from pMolecule        import CrystalClassCubic, MMModelOPLS, NBModelMonteCarlo
from pMoleculeScripts import MergeByAtom

methane = MOLFile_ToSystem ( "../mol/methane.mol" )
water   = MOLFile_ToSystem ( "../mol/water.mol"   )

# . Systems to create.
_ToCreate = ( ( [ methane ] + 215 * [ water ], "Methane in Water.", "ch4_water215_cubicBox_mc.xyz", "ch4_water215_cubicBox_mc.pkl" ),
              (               216 * [ water ], "Water Box."       , "water216_cubicBox_mc.xyz"    , "water216_cubicBox_mc.pkl"     ) )

for ( molecules, label, xyzPath, pklPath ) in _ToCreate:

    system = MergeByAtom ( molecules )
    system.label = label

    xyzSystem  = XYZFile_ToSystem ( "../xyz/" + xyzPath )
    sideLength = float ( xyzSystem.label.split ( )[-1] )

    system.coordinates3 = xyzSystem.coordinates3

    system.DefineMMModel ( MMModelOPLS ( "bookSmallExamples" ) )
    system.DefineNBModel ( NBModelMonteCarlo ( ) )

    system.DefineSymmetry ( crystalClass = CrystalClassCubic ( ), a = sideLength )

    system.Summary ( )

    Pickle ( "../pkl/" + pklPath, system )
 
