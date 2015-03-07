/*------------------------------------------------------------------------------
! . File      : Units.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==============================================================================
!=============================================================================*/
# ifndef _UNITS
# define _UNITS

# include <math.h>

/*------------------------------------------------------------------------------
! . Fundamental constants (NIST web site - October 2005).
!-----------------------------------------------------------------------------*/
/* . Avogadro's number - dimensionless. */
# define CONSTANT_AVOGADRO_NUMBER      6.0221415e+23

/* . Bohr radius - m. */
# define CONSTANT_BOHR_RADIUS          5.291772083e-11

/* . Boltzmann constant - J K^-1. */
# define CONSTANT_BOLTZMANN            1.3806505e-23

/* . Electron charge - A s or C. */
# define CONSTANT_ELECTRON_CHARGE      1.60217653e-19

/* . Vacuum permittivity - A^2 s^4 / kg m^3 or F m^-1. */
# define CONSTANT_VACUUM_PERMITTIVITY  8.854187817e-12

/*------------------------------------------------------------------------------
! . Various units.
!-----------------------------------------------------------------------------*/
/* . Angle. */
# define UNITS_ANGLE_DEGREES_TO_RADIANS ( M_PI / 180.0e+00 )
# define UNITS_ANGLE_RADIANS_TO_DEGREES ( 180.0e+00 / M_PI )

/* . Dipole. */
# define UNITS_DIPOLE_ATOMIC_UNITS_TO_DEBYES 2.54176568e+00

/* . Energy. */
# define UNITS_ENERGY_CALORIES_TO_JOULES 4.184e+00
# define UNITS_ENERGY_JOULES_TO_CALORIES ( 1.0e+00 / 4.184e+00 )


/* . Conversion from e^2/Angstroms to kJ mol^-1. */
# define UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE ( ( 1.0e+7 * CONSTANT_AVOGADRO_NUMBER * CONSTANT_ELECTRON_CHARGE * CONSTANT_ELECTRON_CHARGE ) / \
                                                           ( 4.0e+00 * M_PI * CONSTANT_VACUUM_PERMITTIVITY ) )

/* . Should be equivalent to UNITS_ENERGY_E2ANGSTROMS_TO_KILOJOULES_PER_MOLE * UNITS_LENGTH_ANGSTROMS_TO_BOHRS. */
# define UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE 2625.5e+00

/* . Length. */
# define UNITS_LENGTH_ANGSTROMS_TO_METRES 1.0e-10

# define UNITS_LENGTH_ANGSTROMS_TO_BOHRS ( UNITS_LENGTH_ANGSTROMS_TO_METRES / CONSTANT_BOHR_RADIUS )
# define UNITS_LENGTH_BOHRS_TO_ANGSTROMS ( CONSTANT_BOHR_RADIUS / UNITS_LENGTH_ANGSTROMS_TO_METRES )

/* . Pressure. */
# define UNITS_PRESSURE_ATMOSPHERES_TO_PASCALS               1.013250e+5
# define UNITS_PRESSURE_ATMOSPHERES_TO_KILOJOULES_PER_MOLE ( UNITS_PRESSURE_ATMOSPHERES_TO_PASCALS * CONSTANT_AVOGADRO_NUMBER * 1.0e-33 )

# endif
