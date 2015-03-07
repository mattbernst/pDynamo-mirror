#-------------------------------------------------------------------------------
# . File      : SystemGeometryObjectiveFunction.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Defines the object function class for a system's geometrical variables."""

import math

from pCore     import Clone, CONSTANT_ATOMIC_MASS, CONSTANT_AVOGADRO_NUMBER, CONSTANT_BOLTZMANN, Coordinates3, logFile, LogFileActive, \
                      NormalDeviateGenerator, ObjectiveFunction, RandomNumberGenerator, Real1DArray, Real2DArray, UNITS_PRESSURE_ATMOSPHERES_TO_PASCALS
from pMolecule import CrystalClass

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The conversion factor from amu A^2 ps^-2 to kJ mol^-1 (equivalent to AMU_TO_KG * NAVOGADRO * 10^-3 / MS_TO_APS^2 ).
_AMUA2PS2_TO_KJMOL = 1.0e-2

# . The conversion factor from m s^-1 to A ps^-1.
_MS_TO_APS = 1.0e-2

# . The conversion factor from amu A^2 ps^-2 to Kelvin.
_AMUA2PS2_TO_K = CONSTANT_ATOMIC_MASS / ( CONSTANT_BOLTZMANN * _MS_TO_APS**2 )

# . The conversion factor from kJ mol^-1 to amu A^2 ps^-2.
_KJMOL_TO_AMUA2PS2 = 1.0e+2

# . The conversion factor from atm. Angstroms^3 to kJ mol^-1.
_PV_TO_KJMOL = UNITS_PRESSURE_ATMOSPHERES_TO_PASCALS * CONSTANT_AVOGADRO_NUMBER * 1.0e-33

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SystemGeometryObjectiveFunction ( ObjectiveFunction ):
    """The object function for geometry optimizations and dynamics of a single system."""

    # . Attributes.
    defaultAttributes = { "atomWeights"          :  None ,
                          "degreesOfFreedom"     :     0 ,
                          "fractional"           :  None ,
                          "freeAtoms"            :  None ,
                          "linearScalars"        :  None ,
                          "linearVectors"        :  None ,
                          "ncoordinatevariables" :     0 ,
                          "nsymmetryvariables"   :     0 ,
                          "nvariables"           :     0 ,
                          "QPRESSURE"            : False ,
                          "QSYMMETRY"            : False ,
                          "rtReference"          :  None ,
                          "startingHessian"      :  None ,
                          "system"               :  None ,
                          "variableWeights"      :  None }
    defaultAttributes.update ( ObjectiveFunction.defaultAttributes )

    def AccelerationConversionFactor ( self ):
        """Return the conversion factor from gradient units to acceleration units."""
        # . K to amu A^2 ps^-2.
        return ( - _KJMOL_TO_AMUA2PS2 )

    def Accelerations ( self, variables, accelerations ):
        """Evaluate the function and accelerations."""
        f = self.FunctionGradients ( variables, accelerations )
        accelerations.Scale ( - _KJMOL_TO_AMUA2PS2 )
        return f

    def AddLinearConstraint ( self, constraint ):
        """Add a linear constraint."""
        if len ( constraint ) != self.nvariables: raise ValueError ( "Invalid linear constraint length." )
        # . Orthogonalize to existing constraints.
        if self.linearVectors is not None:
            constraint = Clone ( constraint )
            self.linearVectors.ProjectOutOfArray ( constraint )
        # . Check to see if the constraint is valid.
        cnorm2 = constraint.Norm2 ( )
        if cnorm2 > 1.0e-10:
            constraint.Scale ( 1.0 / cnorm2 )
            # . Allocate space for new constraints.
            ncolumns = 1
            if self.linearVectors is not None: ncolumns += self.linearVectors.columns
            newconstraints = Real2DArray.WithExtents ( len ( constraint ), ncolumns )
            # . Copy over constraints.
            if self.linearVectors is not None:
                for r in range ( self.linearVectors.rows ):
                    for c in range ( self.linearVectors.columns ):
                        newconstraints[r,c] = self.linearVectors[r,c]
            for r in range ( len ( constraint ) ): newconstraints[r,ncolumns-1] = constraint[r]
            self.linearVectors = newconstraints
            # . Determine the linear scalars.
            self.linearScalars = Real1DArray.WithExtent ( self.linearVectors.columns )
            reference          = Real1DArray.WithExtent ( self.linearVectors.rows    )
            if self.rtReference is None: self.system.coordinates3.CopyToArray ( reference )
            else:                        self.rtReference.CopyToArray ( reference )
            self.linearVectors.VectorMultiply ( reference, self.linearScalars, 1.0, 0.0, transpose = True )
            # . Reset the number of degrees of freedom.
            self.degreesOfFreedom = self.nvariables - len ( self.linearScalars )

    def ApplyLinearConstraints ( self, vector ):
        """Apply linear constraints."""
        if self.linearVectors is not None: self.linearVectors.ProjectOutOfArray ( vector )

    def DefinePressure ( self ):
        """Set up pressure calculation."""
        try:
            cc = self.system.symmetry.crystalClass
            if not isinstance ( cc, CrystalClass ): raise
            self.QPRESSURE = True
        except:
            pass

    def DefineWeights ( self ):
        """Define atom weights - use masses by default."""
        self.atomWeights = self.system.atoms.GetItemAttributes ( "mass" )
        if self.nvariables > 0:
            self.variableWeights = Real1DArray.WithExtent ( self.nvariables )
            self.variableWeights.Set ( 0.0 )
            if self.freeAtoms is None: indices = range ( len ( self.system.atoms ) )
            else:                      indices = self.freeAtoms
            n = 0
            for iatom in indices:
                w = math.sqrt ( self.atomWeights[iatom] )
                for j in range ( 3 ):
                    self.variableWeights[n] = w
                    n += 1

    def DistanceSquared ( self, variables1, variables2, gradients ):
        """Calculate the distance squared and its gradients between two points."""
        # . This needs to be updated for crystal cases.
        variables1.CopyTo ( gradients )
        gradients.AddScaledArray ( -1.0, variables2 )
        d2 = gradients.Dot ( gradients )
        gradients.Scale ( 2.0 )
        return d2

    @classmethod
    def FromSystem ( selfClass, system ):
        """Constructor given a system."""
        # . Basic object.
        self = selfClass ( )
        # . Define the system.
        self.system = system
        # . Get the number of atoms.
        natoms = len ( system.atoms )
        # . Check for fixed atoms.
        hc = system.hardConstraints
        if hc is None:
            nfixed = 0
        else:
            nfixed         = hc.NumberOfFixedAtoms    ( )
            self.freeAtoms = hc.fixedAtoms.Complement ( upperBound = natoms )
        # . Set the number of variables.
        self.ncoordinatevariables = 3 * ( natoms - nfixed )
        self.nvariables           = self.ncoordinatevariables
        self.degreesOfFreedom     = self.nvariables
        # . Finish up.
        return self

    def Function ( self, variables ):
        """Evaluate the function."""
        self.VariablesPut ( variables )
        f = self.system.Energy ( log = None )
        return f

    def FunctionGradients ( self, variables, gradients ):
        """Evaluate the function and gradients."""
        self.VariablesPut ( variables )
        f = self.system.Energy ( doGradients = True, log = None )
        self.GradientsGet ( gradients )
        return f

    def GradientsGet ( self, gradients ):
        """Get the gradients."""
        if self.QSYMMETRY:
            spg = self.system.configuration.symmetryParameterGradients
            spg.MakeFractionalDerivatives ( self.system.configuration.symmetryParameters, self.system.coordinates3, self.system.configuration.gradients3 )
            spg.MakeCrystalDerivatives    ( self.system.configuration.symmetryParameters )
            self.system.configuration.gradients3.CopyToArray ( gradients, selection = self.freeAtoms )
            self.system.symmetry.crystalClass.EmptyGradientsToVector ( gradients, spg, self.ncoordinatevariables )
        else:
            self.system.configuration.gradients3.CopyToArray ( gradients, selection = self.freeAtoms )
        if self.variableWeights is not None: gradients.Divide ( self.variableWeights )
        if self.linearVectors is not None: self.linearVectors.ProjectOutOfArray ( gradients )

    def IncludeSymmetryParameters ( self ):
        """Flag the symmetry parameters as variables."""
        if ( self.system.symmetry is not None ) and ( self.system.symmetry.crystalClass is not None ):
            nvariables = len ( self.system.symmetry.crystalClass )
            if nvariables > 0:
                self.fractional         = Coordinates3.WithExtent ( self.ncoordinatevariables // 3 )
                self.nsymmetryvariables = nvariables
                self.nvariables        += nvariables
                self.QSYMMETRY          = True

    def NumberOfVariables ( self ):
        """Return the number of variables."""
        return self.nvariables

    def Pressure ( self, ke ):
        """Calculate the pressure and volume."""
        pressure = 0.0
        volume   = 0.0
        if self.QPRESSURE:
            cc       = self.system.symmetry.crystalClass
            sp       = self.system.configuration.symmetryParameters
            spg      = self.system.configuration.symmetryParameterGradients
            spg.MakeFractionalDerivatives ( sp, self.system.coordinates3, self.system.configuration.gradients3 )
            spg.MakeCrystalDerivatives    ( sp )
            volume   = sp.volume
            pressure = ( ( 2.0 * ke ) / ( 3.0 * volume ) - cc.GetVolumeDerivative ( sp, spg ) ) / _PV_TO_KJMOL
        return ( pressure, volume )

    def RemoveRotationTranslation ( self, reference = None ):
        """Remove rotation and translational degrees of freedom.

        |reference| is the reference coordinates.

        Nothing is done if there are fixed atoms.

        Weighting is done using |atomWeights| if it is has been defined.
        """
        if self.freeAtoms is None:
            # . Check that RT can be removed.
            #QNOTETHER    = ( SoftConstraints_Number_Of_Tethers ( system->constraints ) == 0 )
            #QROTATION    = QNOTETHER && ( ! Symmetry_Is_Periodic ( system->symmetry ) )
            #QTRANSLATION = QNOTETHER
            QROTATION    = ( self.system.symmetry is None ) or ( ( self.system.symmetry is not None ) and ( self.system.symmetry.crystalClass is None ) )
            QTRANSLATION = True
            #print "Rotation and Translation Flags = ", QROTATION, QTRANSLATION
            # . Remove RT.
            if QROTATION or QTRANSLATION:
                # . Assign the reference coordinates and modify the system coordinates if necessary.
                if reference is None:
                    self.rtReference = Clone ( self.system.coordinates3 )
                else:
                    self.rtReference = Clone ( reference )
                    self.system.coordinates3.Superimpose ( self.rtReference, weights = self.atomWeights )
                # . Get the RT vectors.
                ( self.linearVectors, self.linearScalars ) = self.rtReference.RotationTranslationVectors ( QROTATION, QTRANSLATION, dimension = self.nvariables, weights = self.atomWeights )
                self.degreesOfFreedom = self.nvariables - len ( self.linearScalars )

    def SetStartingHessian ( self, hessian ):
        """Set the starting hessian."""
        self.startingHessian = hessian

    def StartingHessian ( self, variables ):
        """Get a starting hessian."""
        if self.startingHessian is None: hessian = self.NumericalHessian ( variables )
        else:                            hessian = self.startingHessian
        if self.linearVectors is not None:
            hessian.ProjectOutVectors ( self.linearVectors          )
            hessian.Raise             ( self.linearVectors, 3.0e+04 )
        return hessian

    def Temperature ( self, velocities ):
        """Calculate the kinetic energy and temperature."""
        ke           = velocities.Dot ( velocities )
        temperature  = _AMUA2PS2_TO_K * ke / float ( self.degreesOfFreedom )
        ke          *= 0.5 * _AMUA2PS2_TO_KJMOL
        return ( ke, temperature )

    def TemperatureConversionFactor ( self ):
        """Return the conversion factor from temperature units to dynamics units."""
        # . K to amu A^2 ps^-2.
        return ( 1.0 / _AMUA2PS2_TO_K )

    def UndefineWeights ( self ):
        """Remove any weighting."""
        if ( self.atomWeights is not None ) or ( self.variableWeights is not None ):
            self.atomWeights     = None
            self.variableWeights = None
            if self.linearVectors is not None: self.RemoveRotationTranslation ( )

    def VariablesAllocate ( self ):
        """Return an object to hold the variables."""
        if self.NumberOfVariables ( ) > 0:
            variables = Real1DArray.WithExtent ( self.NumberOfVariables ( ) )
            variables.Set ( 0.0 )
            return variables
        else:
            return None

    def VariablesGet ( self, variables ):
        """Fill the variable array."""
        if self.QSYMMETRY:
            sp = self.system.configuration.symmetryParameters
            self.fractional.Gather ( self.system.coordinates3, selection = self.freeAtoms )
            self.fractional.Rotate ( sp.inverseM )
            self.fractional.CopyToArray ( variables )
            self.system.symmetry.crystalClass.EmptySymmetryParametersToVector ( variables, sp, self.ncoordinatevariables )
        else:
            self.system.coordinates3.CopyToArray ( variables, selection = self.freeAtoms )
        if self.variableWeights is not None: variables.Multiply ( self.variableWeights )

    def VariablesPut ( self, variables ):
        """Empty the variable array."""
        if self.variableWeights is not None: variables.Divide ( self.variableWeights )
        if self.QSYMMETRY:
            # . Get new box shape first.
            sp = self.system.configuration.symmetryParameters
            self.system.symmetry.crystalClass.FillSymmetryParametersFromVector ( variables, sp, self.ncoordinatevariables )
            # . Now get coordinates.
            self.fractional.CopyFromArray ( variables )
            self.fractional.Rotate  ( sp.M )
            self.fractional.Scatter ( self.system.coordinates3, selection = self.freeAtoms )
        else:
            self.system.coordinates3.CopyFromArray ( variables, selection = self.freeAtoms )
        if self.variableWeights is not None: variables.Multiply ( self.variableWeights )

    def VelocitiesAllocate ( self ):
        """Return an object with velocities.

        The velocities for the system must already exist in configuration.
        """
        return self.system.configuration.velocities

    def VelocitiesAssign ( self, temperature, normalDeviateGenerator = None ):
        """Set up the velocities for the system at a particular temperature.

        If |temperature| is None the velocities must already exist.
        """
        # . Check for an existing set of velocities of the correct size.
        velocities = self.system.configuration.__dict__.get ( "velocities", None )
        QASSIGN = ( velocities is None ) or ( ( velocities is not None ) and ( len ( velocities ) != self.nvariables ) )
        # . Assign velocities.
        if QASSIGN:
            if temperature is None:
                raise ValueError ( "Velocities need to be assigned for the system but a temperature has not been specified." )
            else:
                # . Get the generator.
                if normalDeviateGenerator is None:
                    normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithRandomSeed ( ) )
                # . Get the velocities.
                sigma      = _MS_TO_APS * math.sqrt ( CONSTANT_BOLTZMANN * temperature / CONSTANT_ATOMIC_MASS )
                velocities = Real1DArray.WithExtent ( self.NumberOfVariables ( ) )
                normalDeviateGenerator.NextStandardDeviates ( velocities )
                velocities.Scale ( ( sigma ) )
                self.system.configuration.velocities = velocities
        # . Project out linear constraints.
        if self.linearVectors is not None: self.linearVectors.ProjectOutOfArray ( velocities )
        # . Scale velocities if necessary (even for assigned velocities).
        if temperature is not None:
            ( ke, tactual ) = self.Temperature ( velocities )
            velocities.Scale ( math.sqrt ( temperature / tactual ) )

    def VolumeScale ( self, scale ):
        """Scale the volume of the system."""
        if self.QPRESSURE:
            self.system.symmetry.crystalClass.ScaleVolume ( self.system.configuration.symmetryParameters, scale )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
