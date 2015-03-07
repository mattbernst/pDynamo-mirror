"""Some definitions."""

import os, os.path

# . Output directory for large files.
scratchPath = os.getenv ( "PDYNAMO_SCRATCH" )
if scratchPath is None:
    outPath = ""
else:
    outPath = os.path.join ( scratchPath, "tutorials"      )
    if not os.path.exists  ( outPath ): os.mkdir ( outPath )
    outPath = os.path.join ( outPath    , "pdbFiles"       )
    if not os.path.exists  ( outPath ): os.mkdir ( outPath )
    outPath = os.path.join ( outPath    , "generatedFiles" )
    if not os.path.exists  ( outPath ): os.mkdir ( outPath )

# . Output directory for PDB components.
pdbDataPath = os.path.join ( outPath, "pdbData" )
if not os.path.exists ( pdbDataPath ): os.mkdir ( pdbDataPath )

# . Data path.
dataPath    = os.path.join ( os.getenv ( "PDYNAMO_ROOT" ), "tutorials", "pdbFiles", "data"  )
step2Path   = os.path.join ( os.getenv ( "PDYNAMO_ROOT" ), "tutorials", "pdbFiles", "step2" )
