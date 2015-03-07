"""Generate a MM model for a system."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions import dataPath, outPath
from pCore       import Pickle, Unpickle
from pMolecule   import MMModelOPLS

# . Set up a MM model.
mmModel = MMModelOPLS ( "pkaProtein", path = dataPath )

# . Get the system.
system = Unpickle ( os.path.join ( outPath, "step3.pkl" ) )
system.Summary ( )

# . Add the energy model.
system.DefineMMModel ( mmModel )
system.Summary ( )

# . Save the system.
mmModel.ClearModelBuildingData ( )
Pickle ( os.path.join ( outPath, "step4.pkl" ), system )
