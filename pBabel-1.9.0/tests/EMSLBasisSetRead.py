"""Test for reading EMSL basis set files in Gaussian 94 format."""

import glob, os, os.path

from pBabel    import EMSLG94File_ToGaussianBases
from pCore     import logFile, LogFileActive, TestCase, YAMLPickle, YAMLPickleFileExtension
from pMolecule import PeriodicTable

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The destination for results.
_Destination = "gaussianBasisSets"

# . The data path.
_Source = "emslG94"

# . The file extension.
_Extension = ".emslg94"

# . The data to try.
_EMSLG94Data = [ ( "321g" , False ), ( "631gs", False ), ( "ahlrichs", True ),
                 ( "demon", False ), ( "qzvp" , True  ), ( "sto3g"   , True ),
                 ( "svp"  , True  ), ( "tzvp" , True  ), ( "weigend" , True ) ]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class EMSLBasisSetReadTest ( TestCase ):
    """A test case for reading EMSL basis set files in Gaussian 94 format."""

    def runTest ( self ):
        """The test."""

        # . Initialization.
        isOK = True

        # . Output setup.
        dataPath = os.path.join ( os.getenv ( "PDYNAMO_PBABEL" ), "data", _Source )
        log = self.GetLog ( )

        # . Loop over the data.
        for ( label, isSpherical ) in _EMSLG94Data:

            # . Get the bases.
            bases = EMSLG94File_ToGaussianBases ( os.path.join ( dataPath, label + _Extension ), isSpherical = isSpherical )
            if log is not None:
                log.Paragraph ( "Processed bases for basis {:s} = {:d}.".format ( label, len ( bases ) ) )

            # . Check for an appropriate outPath.
            if self.resultPath is not None:
                outPath = os.path.join ( self.resultPath, _Destination )
                if not os.path.exists ( outPath ): os.mkdir ( outPath )
                outPath = os.path.join ( outPath, label.lower ( ) )
                if not os.path.exists ( outPath ): os.mkdir ( outPath )

                # . Save the bases.
                for basis in bases:
                    symbol = PeriodicTable.Symbol ( basis.atomicNumber )
                    YAMLPickle ( os.path.join ( outPath, symbol + YAMLPickleFileExtension ), basis )

        # . Success/failure.
        self.assertTrue ( isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    # . Run the test.
    test = EMSLBasisSetReadTest ( )
    test.run ( )
