#-------------------------------------------------------------------------------
# . File      : OOGLOffFileReader.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for reading OOGL Off files."""

from pCore        import logFile, LogFileActive, PolygonalSurface, TextFileReader
from ExportImport import _Importer

#===================================================================================================================================
# . OOGLOff file reader class.
#===================================================================================================================================
class OOGLOffFileReader ( TextFileReader ):
    """OOGLOffFileReader is the class for OOGLOff files that are to be read."""

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.QPARSED:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Header line.
                line       = self.GetLine ( )
                hasNormals = ( line.find ( "N" ) > -1 )
                if line.find ( "n" ) > -1:
                    items = self.GetTokens ( converters = [ int ] )
                    ndimensions = items[0]
                else:
                    ndimensions = 3
                # . Counters line.
                items = self.GetTokens ( converters = [ int, int, int ] )
                nvertices = items[0]
                npolygons = items[1]
                # . Allocate the object.
                surface      = PolygonalSurface.WithSizes ( ndimensions, npolygons, nvertices, initialize = True )
                self.surface = surface
                # . Empty line.
                self.GetLine ( )
                # . Vertex lines.
                nfields = ndimensions
                if hasNormals: nfields *= 2
                converters = nfields * [ float ]
                for i in range ( nvertices ):
                    items = self.GetTokens ( converters = converters )
                    surface.SetVertex ( i, items[:ndimensions] )
                    if hasNormals: surface.SetNormal ( i, items[ndimensions:] )
                # . Empty line.
                self.GetLine ( )
                # . Face lines.
                nfields = ndimensions + 1
                converters = nfields * [ int ]
                for i in range ( npolygons ):
                    items = self.GetTokens ( converters = converters )
                    surface.SetPolygon ( i, items[1:] )
            except EOFError:
                pass
	    # . Close the file.
            self.WarningStop ( )
	    self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.QPARSED = True

    def ToPolygonalSurface ( self ):
        """Return a polygonal surface."""
        surface = None
        if self.QPARSED: surface = self.surface
        return surface

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def OOGLOffFile_ToPolygonalSurface ( path, log = logFile ):
    """Helper function that reads a surface from an OOGL Off file."""
    infile = OOGLOffFileReader ( path )
    infile.Parse ( )
    surface = infile.ToPolygonalSurface ( )
    return surface

# . Importer definitions.
_Importer.AddHandler ( { PolygonalSurface : OOGLOffFile_ToPolygonalSurface } , [ "off", "oogloff" ], "OOGL Off", defaultFunction = OOGLOffFile_ToPolygonalSurface )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
