#-------------------------------------------------------------------------------
# . File      : SystemSoftConstraintTrajectory.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for manipulating system soft-constraint trajectories."""

import glob, os, os.path, random

from pCore import Pickle                    , \
                  PickleFileExtension       , \
                  Real1DArray               , \
                  RegularHistogram          , \
                  RegularHistogramDimension , \
                  Unpickle

#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# . Remarks:
#
# * These trajectories have really been designed to work for potential of mean force calculations using WHAM and so they should
#   perhaps be renamed to SystemWHAMTrajectory?
#
# * The trajectories may be changed so that they store the unbiased energy - either directly or indirectly (the biased energy
#   along with the total softconstraint energy).
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Default block size.
_DefaultBlockSize = 1024

# . Block naming.
_BlockPrefix  = "block"
_BlockPostfix = PickleFileExtension

# . Header and footer naming.
_FooterName = "footer"
_HeaderName = "header"

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def SystemSoftConstraintTrajectory ( path, owner, mode = None ):
    """Function to return an appropriate system soft-constraint trajectory object."""
    # . Find the mode.
    if mode is None:
        if os.path.exists ( path ): mode == "r"
        else:                       mode == "w"
    # . Branch on the mode.
    if   mode == "a": return SystemSoftConstraintTrajectoryWriter ( path, owner, useAppendMode = True )
    elif mode == "r": return SystemSoftConstraintTrajectoryReader ( path, owner )
    elif mode == "w": return SystemSoftConstraintTrajectoryWriter ( path, owner )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SystemSoftConstraintTrajectoryDataHandler ( object ):
    """Handle the data in a series of system soft constraint trajectories."""

    defaultAttributes = { "data"   : None ,
                          "labels" : None ,
                          "paths"  : None ,
                          "system" : None }

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__.defaultAttributes.iteritems ( ): setattr ( self, key, value )

    def ExtractData ( self ):
        """Extract data from the trajectories."""
        if ( self.data is None ) and ( self.paths is not None ):
            # . Initialization.
            numberOfWindows = len ( self.paths )
            rank            = None
            data            = []
            energyModels    = []
            periods         = []
            pressures       = []
            temperatures    = []
            windowPoints    = []
            # . Loop over paths.
            for path in self.paths:
                # . Open the trajectory.
                trajectory = SystemSoftConstraintTrajectory ( path, self.system, mode = "r" )
                header     = trajectory.ReadHeader ( )
                # . First trajectory.
                if rank is None:
                    # . Basic information.
                    labels    = header["labels"]                      
                    rank      = len ( labels )                        
                    hasEnergy = header["hasPotentialEnergy"]          
                    hasVolume = header["hasVolume"]                   
                    # . Constraint information.
                    numberOfConstraints = rank
                    if hasEnergy: numberOfConstraints -= 1
                    if hasVolume: numberOfConstraints -= 1
                    models = header["softConstraintEnergyModels"]
                    for model in models:
                        if model.isPeriodic: periods.append ( model.period )
                        else:                periods.append ( None         )
                # . Subsequent trajectories.
                else:
                    # . Verify information.
                    newLabels = header["labels"]
                    isOK = ( len ( newLabels ) == rank ) and ( hasEnergy == header["hasPotentialEnergy"] ) and ( hasVolume == header["hasVolume"] )
                    if isOK:
                        for ( new, old ) in zip ( newLabels, labels ):
                            if new != old:
                                isOK = False
                                break
                    if not isOK: raise ValueError ( "Incompatible trajectories." )
                    # . Models.
                    models = header["softConstraintEnergyModels"]
                # . Save other data about the window.
                pressures.append    ( header.get ( "pressure",    None ) )
                temperatures.append ( header.get ( "temperature", None ) )
                # . Read data.
                frameData      = trajectory.ReturnAllFrameDataAsList ( )
                numberOfFrames = len ( frameData ) // rank
                windowPoints.append ( numberOfFrames )
                # . Unbias the energies (by removing the energies of the stored constraints).
                if hasEnergy and ( numberOfConstraints > 0 ):
                    d = 1
                    if hasVolume: d += 1
                    for i in range ( 0, numberOfFrames, rank ):
                        e = frameData[i]
                        for ( s, model ) in enumerate ( models ):
                            e -= model.Energy ( data[i+d+s] )[0]
                        frameData[i] = e
                # . Save the data.
                data.extend ( frameData )
                energyModels.append ( models )
                # . Finish up.
                trajectory.Close ( )
            # . Save all.
            self.data            = data
            self.energyModels    = energyModels
            self.hasEnergy       = hasEnergy
            self.hasVolume       = hasVolume
            self.labels          = labels
            self.numberOfWindows = numberOfWindows
            self.periods         = periods
            self.pressures       = pressures
            self.rank            = rank
            self.temperatures    = temperatures
            self.windowPoints    = windowPoints

    @classmethod
    def FromTrajectoryPaths ( selfClass, paths ):
        """Constructor given trajectory paths."""
        self       = selfClass ( )
        self.paths = paths
        self.ExtractData ( )
        return self

    def GetDataRanges ( self ):
        """Get the maximum and minimum values of each type of data."""
        maxima = None
        minima = None
        if self.data is not None:
            maxima = [ max ( data[i::self.rank] ) for i in range ( self.rank ) ]
            minima = [ min ( data[i::self.rank] ) for i in range ( self.rank ) ]
        return ( maxima, minima )

    def HistogramData ( self, bins ):
        """Histogram the data."""
        histogram = None
        if self.data is not None:
            # . Create the dimensions.
            dimensions = []
            for d in range ( self.rank ):
                dimensions.append ( RegularHistogramDimension.FromData ( self.data[d::self.rank], bins = bins[d], label = self.labels[d], period = self.periods[d] ) )
            # . Create the histogram.
            histogram = RegularHistogram.FromDimensions ( dimensions )
            histogram.BinData ( self.data )
        return histogram

    def Prune ( self, toKeep ):
        """Keep data corresponding to particular labels only."""
        if ( self.data is not None ) and ( self.labels is not None ):
            # . Find the label indices.
            indices = []
            for label in set ( toKeep ):
                try   : indices.append ( self.labels.index ( label ) )
                except: raise ValueError ( "Invalid prune label: {:s}".format ( label ) )
            if len ( indices ) <= 0: raise ValueError ( "All data cannot be pruned." )
            # . Processing.
            if len ( indices ) < len ( self.labels ):
                # . Energy and volume flags.
                hasEnergy = self.hasEnergy and ( "Potential Energy" in toKeep )
                hasVolume = self.hasVolume and ( "Volume"           in toKeep )
                # . Treat indices.
                indices.order ( )
                # . Indices of soft constraints.
                d = 0
                if hasEnergy: d += 1
                if hasVolume: d += 1
                scIndices = [ ( i - d ) for i in indices[d:] ]
                # . Old data.
                oldData         = self.data
                oldEnergyModels = self.energyModels
                oldLabels       = self.labels
                oldPeriods      = self.periods
                oldRank         = self.rank
                # . Data.
                numberOfFrames  = sum ( self.windowPoints )
                data            = []
                for c in range ( 0, numberOfFrames, oldRank ):
                    for r in indices: data.append ( oldData[c+r] )
                # . Energy models.
                energyModels = []
                for oldModels in oldEnergyModels:
                    energyModels.append ( [ oldModels[i] for i in scIndices ] ) 
                # . Save all.
                self.data         = data
                self.energyModels = energyModels
                self.hasEnergy    = hasEnergy
                self.hasVolume    = hasVolume
                self.labels       = [ oldLabels [i] for i in indices   ]
                self.periods      = [ oldPeriods[i] for i in scIndices ]
                self.rank         = len ( indices )

    def RehistogramData ( self, histogram ):
        """Histogram the data using an existing histogram."""
        if ( self.data is not None ) and ( histogram is not None ):
            histogram.InitializeCounts ( )
            histogram.BinData ( self.data )

    def ResampleData ( self ):
        """Resample the data with replacement."""
        if self.data is not None:
            oldData = getattr ( self, "originalData", None )
            if oldData is None:
                self.originalData = self.data
                oldData           = self.data
            newData = []
            rank    = self.rank
            start   = 0
            for points in self.windowPoints:
                end = start + points
                for i in range ( points ):
                    newStart = random.randrange ( start, end ) * rank
                    newData.extend ( oldData[newStart:newStart+rank] )
                start = end
# . To remove.
            if len ( newData ) != len ( oldData ): raise ValueError ( "Logic error in resampling." )
            self.data = newData

    def ResetData ( self ):
        """Reset to the original data."""
        originalData = getattr ( self, "originalData", None )
        if originalData is not None: self.data = originalData

    # . For E and V adjust the values by a constant amount to help numerical stability.
    # . This will not change the WHAM procedure as the changes are absorbed into the
    # . definition of the fs. This procedure is applicable to all items that enter
    # . linearly into the exponential.
    def TranslateEnergyVolumeData ( self ):
        """Translate energy and volume data to their minimum values."""
        if self.data is not None:
            n = 0
            if self.hasEnergy: n += 1
            if self.hasVolume: n += 1
            for d in range ( n ):
                minimum = min ( self.data[d::self.rank] )
                for i in range ( d, len ( self.data ), self.rank ): self.data[i] -= minimum

#===================================================================================================================================
# . Reader class.
#===================================================================================================================================
class SystemSoftConstraintTrajectoryReader ( object ):
    """Class for reading system soft constraint trajectories."""

    defaultAttributes = { "blocks"       :     0 ,
                          "blockSize"    :     0 ,
                          "currentBlock" :    -1 ,
                          "current"      :     0 ,
                          "data"         :  None ,
                          "owner"        :  None ,
                          "path"         :  None ,
                          "rank"         :     0 }

    def __init__ ( self, path, owner ):
        """Constructor."""
        self.__dict__.update ( self.__class__.defaultAttributes )
        # . Arguments.
        self.path  = path
        self.owner = owner
        # . Open the trajectory.
        self.Open ( )

    def Close ( self ):
        """Close the trajectory."""
        pass

    def Open ( self ):
        """Open the trajectory."""
        # . Check that the trajectory exists, is a directory and is readable.
        if os.access ( self.path, os.F_OK ) and os.access ( self.path, os.R_OK ) and os.path.isdir ( self.path ):
            # . Check for a valid header.
            if not os.path.exists ( os.path.join ( self.path, _HeaderName + _BlockPostfix ) ): raise IOError ( "Unable to find trajectory header." )
            # . Find the number of blocks.
            self.blocks = len ( glob.glob ( os.path.join ( self.path, _BlockPrefix + "*" + _BlockPostfix ) ) )
        # . Invalid trajectory.
        else: raise IOError ( "Invalid or non-existent trajectory." )

    def ReadBlock ( self ):
        """Read a block of data."""
        if self.currentBlock < self.blocks:
            self.data = Unpickle ( os.path.join ( self.path, _BlockPrefix + "{:d}".format ( self.currentBlock ) + _BlockPostfix ) )
            self.blockSize     = len ( self.data ) // self.rank
            self.current       = 0
            self.currentBlock += 1
        else: raise IndexError ( "Invalid block index." )

    def ReadFooter ( self ):
        """Read the footer."""
        pass

    def ReadHeader ( self ):
        """Read the trajectory header."""
        header = Unpickle ( os.path.join ( self.path, _HeaderName + _BlockPostfix ) )
        # . Set some object values.
        for ( key, value ) in header.iteritems ( ): setattr ( self, key, value )
        # . Set the number of dimensions.
        self.rank = len ( self.labels )
        # . Return the data.
        return header

    def RestoreOwnerData ( self ):
        """Restore data from a frame to the owner."""
        # . Data is not restored directly (as it would not make sense) but is put into an artificial state.
        if self.current >= self.blockSize:
            # . No more data.
            if self.currentBlock >= self.blocks:
                return False
            # . Read the next block.
            else:
                self.ReadBlock ( )
        # . Get the data.
        state = {}
        n     = self.current * self.rank
        for ( i, label ) in enumerate ( self.labels ):
            state[label] = self.data[n+i]
        self.current += 1
        # . Set the state.
        self.owner.configuration.SetTemporaryAttribute ( "scState", state )
        return True

    def ReturnAllFrameDataAsList ( self ):
        """Return all frame data as a list."""
        data = []
        self.currentBlock = 0
        for i in range ( self.blocks ):
            self.ReadBlock ( )
            data.extend ( list ( self.data ) )
        return data

#===================================================================================================================================
# . Writer class.
#===================================================================================================================================
class SystemSoftConstraintTrajectoryWriter ( object ):
    """Class for writing system soft constraint trajectories."""

    defaultAttributes = { "blocks"               :                 0 ,
                          "blockSize"            : _DefaultBlockSize ,
                          "current"              :                 0 ,
                          "data"                 :              None ,
                          "frames"               :                 0 ,
                          "owner"                :              None ,
                          "path"                 :              None ,
                          "hasPotentialEnergy"   :             False ,
                          "hasSoftConstraints"   :             False ,
                          "hasVolume"            :             False ,
                          "rank"                 :                 0 ,
                          "softConstraintLabels" :              None ,
                          "useAppendMode"        :             False }

    def __init__ ( self, path, owner, useAppendMode = False ):
        """Constructor."""
        self.__dict__.update ( self.__class__.defaultAttributes )
        # . Arguments.
        self.owner         = owner
        self.path          = path
        self.useAppendMode = useAppendMode
        # . Open the trajectory.
        self.Open ( )

    def Close ( self ):
        """Close the trajectory."""
        # . Flush out any remaining data.
        self.WriteBlock ( )
        # . Write the footer.
        #self.WriteFooter ( )

    def Open ( self ):
        """Open the trajectory."""
        # . Check to see if the path already exists and is a directory.
        pathExists = os.access ( self.path, os.F_OK )
        if pathExists:
            if not os.path.isdir ( self.path ): raise IOError ( "Trajectory exists that is not a directory." )
        else:
            os.mkdir ( self.path )
        # . Check for writeability.
        if not os.access ( self.path, os.W_OK ): raise IOError ( "Trajectory is not writeable." )
        # . Check the contents of an existing trajectory.
        if pathExists:
            # . Append mode - find the number of existing blocks.
            if self.useAppendMode:
                self.blocks = len ( glob.glob ( os.path.join ( self.path, _BlockPrefix + "*" + _BlockPostfix ) ) )
            #. Write mode - remove all existing files.
            else:
                for target in glob.glob ( os.path.join ( self.path, "*" ) ): os.remove ( target )

    def WriteBlock ( self ):
        """Write a block of data."""
        if self.current > 0:
            # . The following is a fudge until slicing is handled properly.
            if self.current == self.blockSize: outdata = self.data
            else:
                outdata = Real1DArray.WithExtent ( self.current * self.rank )
                for i in range ( len ( outdata ) ): outdata[i] = self.data[i]
            Pickle ( os.path.join ( self.path, _BlockPrefix + "{:d}".format ( self.blocks ) + _BlockPostfix ), outdata )
# . This is what it should be.
#            Pickle ( os.path.join ( self.path, _BlockPrefix + "{:d}".format ( self.blocks ) + _BlockPostfix, self.data[0:self.current*self.rank] )
            self.blocks += 1
            self.current = 0

    def WriteFooter ( self ):
        """Write a footer."""
        pass

    def WriteHeader ( self, pressure = None, temperature = None ):
        """Write the trajectory header."""
        # . Check for soft constraints.
        self.hasSoftConstraints = hasattr ( self.owner.energyModel, "softConstraints" ) and ( len ( self.owner.energyModel.softConstraints ) > 0 )
        # . Check the volume option.
        self.hasVolume          = self.hasVolume and hasattr ( self.owner, "symmetry" )
        # . Check that there will be data on the trajectory.
        if not ( self.hasPotentialEnergy or self.hasSoftConstraints or self.hasVolume ): raise ValueError ( "A soft-constraint trajectory contains no data." )
        # . Initialization.
        labels = []
        models = []
        # . Potential energy and volume.
        if self.hasPotentialEnergy: labels.append ( "Potential Energy" )
        if self.hasVolume:          labels.append ( "Volume"           )
        # . Soft constraints.
        if self.hasSoftConstraints:
            self.softConstraintLabels = sorted ( self.owner.energyModel.softConstraints )
            for label in self.softConstraintLabels: models.append ( self.owner.energyModel.softConstraints[label].energyModel )
            labels.extend ( self.softConstraintLabels )
        # . Set up the header.
        header   = { "labels" : labels, "hasPotentialEnergy" : self.hasPotentialEnergy, "hasVolume" : self.hasVolume, "softConstraintEnergyModels" : models }
        # . Optional information.
        if pressure    is not None: header["pressure"   ] = pressure
        if temperature is not None: header["temperature"] = temperature
        # . Write out.
        Pickle ( os.path.join ( self.path, _HeaderName + _BlockPostfix ), header )
        # . Setup the block for output.
        self.rank = len ( labels )
        self.data = Real1DArray.WithExtent ( self.rank * self.blockSize )

    def WriteOwnerData ( self ):
        """Write data from the owner to a frame."""
        if self.current >= self.blockSize: self.WriteBlock ( )
        n = self.current * self.rank
        if self.hasPotentialEnergy:
            self.data[n] = self.owner.configuration.energyTerms.potentialEnergy
            n += 1
        if self.hasVolume:
            self.data[n] = self.owner.symmetryParameters.volume
            n += 1
        if self.hasSoftConstraints:
            scState = self.owner.configuration.scState
            for label in self.softConstraintLabels:
                self.data[n] = scState[label][1]
                n += 1
        self.current += 1
        self.frames  += 1

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
