#-------------------------------------------------------------------------------
# . File      : OOGLOffFileWriter.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for writing OOGL Off files."""

import math

from pCore        import logFile, LogFileActive, PolygonalSurface, TextFileWriter
from ExportImport import _Exporter

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class OOGLOffFileWriter ( TextFileWriter ):
    """Class for writing OOGL Off files."""

    def WritePolygonalSurface ( self, surface ):
        """Write a surface."""
        # . Initialization.
        self.Open ( )
        # . Header.
        ndimensions = surface.rank
        if ndimensions == 3: self.file.write ( "NOFF\n" )
        else:                self.file.write ( "NnOFF\n{:d}\n".format ( ndimensions ) )
        # . Number of vertices, faces and edges (latter not used).
        self.file.write ( "{:d} {:d} 0\n\n".format ( surface.NumberOfVertices ( ), surface.NumberOfPolygons ( ) ) )
        # . Vertices.
        format = 2 * ndimensions * " {:.10f}" + "\n"
        for items in surface.VertexIterator ( ):
            self.file.write ( format.format ( items ) )
        # . Polygons.
        self.file.write ( "\n" )
        format = "{:d}".format ( ndimensions ) + ndimensions * " {:d}" + "\n"
        for items in surface.PolygonIndexIterator ( ):
            self.file.write ( format.format ( items ) )
        # . Finish up.
        self.Close ( )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def OOGLOffFile_FromPolygonalSurface ( path, surface, log = logFile ):
    """Helper function that writes a polygonal surface to an OOGL Off file."""
    outfile = OOGLOffFileWriter ( path )
    outfile.WritePolygonalSurface ( surface )

# . Exporter definitions.
_Exporter.AddHandler ( { PolygonalSurface : OOGLOffFile_FromPolygonalSurface } , [ "off", "oogloff" ], "OOGL Off" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
