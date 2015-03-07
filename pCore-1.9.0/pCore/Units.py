#-------------------------------------------------------------------------------
# . File      : Units.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Define physical and chemical units.

Currently only constants are provided.
"""

import math

# . Fundamental constants (NIST web site - October 2005).
CONSTANT_ATOMIC_MASS            = 1.66053886e-27  # . kg.
CONSTANT_AVOGADRO_NUMBER        = 6.0221415e23    # . Dimensionless.
CONSTANT_BOLTZMANN              = 1.3806505e-23   # . J K^-1.
CONSTANT_ELECTRON_CHARGE        = 1.60217653e-19  # . A s or C.
CONSTANT_MOLAR_GAS              = 8.314472        # . J mol^-1 K^-1.
CONSTANT_MOLAR_IDEAL_GAS_VOLUME = 22.413996e-3    # . m^3 mol^-1 (at 273.15 K and 101.325 kPa).
CONSTANT_PLANCK                 = 6.6260693e-34   # . J s.
CONSTANT_SPEED_OF_LIGHT         = 299792458.0     # . m s^-1.
CONSTANT_VACUUM_PERMITTIVITY    = 8.854187817e-12 # . A^2 s^4 / kg m^3 or F m^-1.

# . Angle.
UNITS_ANGLE_DEGREES_TO_RADIANS = math.pi / 180.0
UNITS_ANGLE_RADIANS_TO_DEGREES = 180.0 / math.pi

# . Dipole.
UNITS_DIPOLE_ATOMIC_UNITS_TO_DEBYES = 2.54176568

# . Energy.
UNITS_ENERGY_HARTREES_TO_ELECTRON_VOLTS                   = 27.2113845
UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE              = 2625.5
UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE = 4.184

UNITS_ENERGY_ELECTRON_VOLTS_TO_KILOJOULES_PER_MOLE        = UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE / UNITS_ENERGY_HARTREES_TO_ELECTRON_VOLTS

# . Conversion from e^2/Angstroms to kJ mol^-1.
UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE = ( ( 1.0e+7 * CONSTANT_AVOGADRO_NUMBER * CONSTANT_ELECTRON_CHARGE**2 ) / ( 4.0 * math.pi * CONSTANT_VACUUM_PERMITTIVITY ) )

# . Length.
UNITS_LENGTH_ANGSTROMS_TO_BOHRS = 1.0e-10 / 5.291772083e-11
UNITS_LENGTH_BOHRS_TO_ANGSTROMS = 5.291772083e-11 / 1.0e-10

# . Mass.
UNITS_MASS_AMU_TO_KG = CONSTANT_ATOMIC_MASS

# . Pressure.
UNITS_PRESSURE_ATMOSPHERES_TO_PASCALS             = 1.013250e+5
UNITS_PRESSURE_ATMOSPHERES_TO_KILOJOULES_PER_MOLE = UNITS_PRESSURE_ATMOSPHERES_TO_PASCALS * CONSTANT_AVOGADRO_NUMBER * 1.0e-33
