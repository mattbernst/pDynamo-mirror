#-------------------------------------------------------------------------------
# . File      : Clone.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Clone functions."""

# . Clone function.
import copy

# . Generic.
Clone = copy.deepcopy

# . Specific.
ShallowClone = copy.copy
DeepClone    = copy.deepcopy

del copy
