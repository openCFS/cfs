// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SCFE_MYDEFS_2001
#define FILE_SCFE_MYDEFS_2001

#include <typeinfo>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

// includes for the C99 standard datatypes (e.g. uint32_t, long double)
#if defined(__GNUC__)
#include <stdint.h>
#else
#include "pstdint.h"
#endif
#include <cfloat>

#include "General/defs.hh"
#include "boost/shared_ptr.hpp"

#ifndef INTEGLIB
#include "olas.hh"
#endif

//! \file environment.hh
//! This file contains some global macro, class and enumeration data type
//! definitions for CFS++.


//! Maximal length of the trailing postfix of an auxilliary name,
//! i.e. the length of the extension after the basename
#define MAXPOSTFIX 15


namespace OutInfo {

  //! Global pointer to a stringstream used for generating error messages
  extern std::stringstream *error;

  //! Global pointer to a stringstream used for generating warning messages
  extern std::stringstream *warning;
}


namespace CoupledField {

  // Import Boost's namespace
  using namespace boost;


  // forward class declaration
  class Domain;
  class CFSMessenger;
  class LogConfigurator;

  //! redeclaration of types
  typedef int32_t Integer;
  typedef uint32_t UInt;
  typedef int16_t ShortInt;
  typedef float Float;
  typedef double Double;
  typedef std::complex<Double> Complex;
  typedef char Char;

  //! Define Enum types for numerical entries of vectors and matrices
  struct EntryType {
    typedef enum {NOENTRYTYPE, DOUBLE, FLOAT, COMPLEX, INTEGER, 
                  UINT, SHORTINT} ScalarType;
  };

  //! Define assignment of numerical types to their enum-representation
  template<class TYPE>
  struct EntryTypeMap {
    //! Associated enum representation of entry type
    static const EntryType::ScalarType S_TYPE = EntryType::NOENTRYTYPE;
  };
  
  // Explicit sepcialization for scalar entry types
#define DEFINE_SCALAR_TYPE(TYPE,TYPE_ENUM)                    \
 template<>                                                   \
  struct EntryTypeMap<TYPE> {                                 \
    static const EntryType::ScalarType S_TYPE = TYPE_ENUM;    \
  }
  DEFINE_SCALAR_TYPE( Double   , EntryType::DOUBLE);
  DEFINE_SCALAR_TYPE( Float    , EntryType::FLOAT);
  DEFINE_SCALAR_TYPE( Complex  , EntryType::COMPLEX);
  DEFINE_SCALAR_TYPE( Integer  , EntryType::INTEGER);
  DEFINE_SCALAR_TYPE( UInt     , EntryType::UINT);
  DEFINE_SCALAR_TYPE( ShortInt , EntryType::SHORTINT);
#undef DEFINE_SCALAR_TYPE
  
  // Type definitions for regions
  typedef int RegionIdType;
#define NO_REGION_ID -1
#define ALL_REGIONS -2

#define myEndl std::endl
#define myCout std::cout

  typedef Double (*pFncWith1Var)(const Double);

#define myEndl std::endl
#define myCout std::cout

#define myendl std::endl
#define mycout std::cout

  //! useful trick for testing problem
#define mark std::cout<<__FILE__<<__LINE__<<std::endl;

  const Double PI = acos(-1.0);
  const Double NORM_EPS = 1e-6;  // needed e.g. for lower bounds of norms in iteration loops
  const Double EPS = 1e-12;     // value for absolute precision (needed e.g. for lower bounds of norms in iteration loops)

  //! specifications of Lapack routines for different types of system matrices in 
  //! matrix.solveWithLapack
  //! Z - Complex valued matrix
  //! SY - symmetric matrix
  //! GE - general matrix
  //! HE - hermitian matrix
  typedef enum {ZSYSV, ZGESV, ZHESV} lapackSysMatType;
  

  //! logging configurator
  
  extern LogConfigurator * logConf;

#ifdef PROFILING
  //! Global memtrace pointer
  extern Profiler * profiler;
#endif

  //! Global pointer to domain object
  extern Domain *domain;

#ifdef USE_SCRIPTING
  //! Global pointer to messenger object
  extern CFSMessenger * messenger;
#endif

  // Definition of base datatypes

  typedef float            float32_t;
  typedef double           float64_t;
  typedef long double      float96_t;

  // enumeration with elements types.
  // enum ElementType{Line1, Triang1, Triang2, Quadrilateral1, Quadrilateral2};

  //! Type of finite element
  //! The enumeration contains the following values:
  //! - NOFETYPE
  //! - LINE
  //! - TRIA
  //! - QUAD
  //! - TET
  //! - HEX
  //! - PYR
  //! - WED
    //  typedef enum {NOFETYPE, LINE, TRIA, QUAD, TET, HEX, PYR, WED} FEType;

  // Definition of supported element types

  typedef enum  
  {
    ET_UNDEF      = 0,
    ET_POINT      = 1,
    ET_LINE2      = 2,
    ET_LINE3      = 3,
    ET_TRIA3      = 4,
    ET_TRIA6      = 5,
    ET_QUAD4      = 6,
    ET_QUAD8      = 7,
    ET_QUAD9      = 8,
    ET_TET4       = 9,
    ET_TET10      = 10,
    ET_HEXA8      = 11,
    ET_HEXA20     = 12,
    ET_HEXA27     = 13,
    ET_PYRA5      = 14,
    ET_PYRA13     = 15,
    ET_WEDGE6     = 16,
    ET_WEDGE15    = 17
  } FEType;

  const uint32_t NUM_ELEM_NODES[] = 
  {
    0,  // ET_UNDEF  
    1,  // ET_POINT
    2,  // ET_LINE2  
    3,  // ET_LINE3  
    3,  // ET_TRIA3  
    6,  // ET_TRIA6  
    4,  // ET_QUAD4  
    8,  // ET_QUAD8  
    9,  // ET_QUAD9  
    4,  // ET_TET4   
    10, // ET_TET10  
    8,  // ET_HEXA8  
    20, // ET_HEXA20 
    27, // ET_HEXA27 
    5,  // ET_PYRA5  
    13, // ET_PYRA13 
    6,  // ET_WEDGE6 
    15  // ET_WEDGE15
  };

  const uint32_t ELEM_DIM[] = 
  {
    0,  // ET_UNDEF  
    0,  // ET_POINT
    1,  // ET_LINE2  
    1,  // ET_LINE3  
    2,  // ET_TRIA3  
    2,  // ET_TRIA6  
    2,  // ET_QUAD4  
    2,  // ET_QUAD8  
    2,  // ET_QUAD9  
    3,  // ET_TET4   
    3,  // ET_TET10  
    3,  // ET_HEXA8  
    3,  // ET_HEXA20 
    3,  // ET_HEXA27 
    3,  // ET_PYRA5  
    3,  // ET_PYRA13 
    3,  // ET_WEDGE6 
    3   // ET_WEDGE15
  };

  const std::string ELEM_TYPE_NAMES[] = 
  {
    "ET_UNDEF",  
    "ET_POINT",         
    "ET_LINE2",         
    "ET_LINE3",         
    "ET_TRIA3",         
    "ET_TRIA6",         
    "ET_QUAD4",         
    "ET_QUAD8",         
    "ET_QUAD9",         
    "ET_TET4",          
    "ET_TET10",         
    "ET_HEXA8",         
    "ET_HEXA20",        
    "ET_HEXA27",        
    "ET_PYRA5",         
    "ET_PYRA13",        
    "ET_WEDGE6",        
    "ET_WEDGE15"       
  };

  const bool ELEM_QUADRATIC[] = 
  {
    false,  // ET_UNDEF  
    false,  // ET_POINT
    false,  // ET_LINE2  
    true,   // ET_LINE3  
    false,  // ET_TRIA3  
    true,   // ET_TRIA6  
    false,  // ET_QUAD4  
    true,   // ET_QUAD8  
    true,   // ET_QUAD9  
    false,  // ET_TET4   
    true,   // ET_TET10  
    false,  // ET_HEXA8  
    true,   // ET_HEXA20 
    true,   // ET_HEXA27 
    false,  // ET_PYRA5  
    true,   // ET_PYRA13 
    false,  // ET_WEDGE6 
    true    // ET_WEDGE15
  };


  /** enumeration with integration types - to be renamed into IntegrationType
   * ECONOMICAL are the "efficient" gaussian quadrature weights
   *            (->Solin, Segeth, Dolezel, Higher-Order Finite Element Methods) -> see also Elements/BaseFE 
   * CLASSICAL either the "original" implementation before ECONOMICAL refactoring or fixed product rule (wedge)
   * EXPERIMENTAL is only for developer/debugging use -> one does not need to modify "real" data
   * CARTESIAN is only for developer/debugging -> x1 =0-9, x2 = 1x-9x, x3 = 1xx-9xx
   * UNDEFINED only internal use!
   * LOBATTO and CHEBYSHEV are special line intgration methods (->Solin, Segeth, Dolezel, Higher-Order Finite Element Methods)*/
  enum IntegrationMethod {ECONOMICAL, CLASSICAL, LOBATTO, CHEBYSHEV, EXPERIMENTAL, CARTESIAN, SPECIAL, UNDEFINED}; 

  //! Damping type
  enum DampingType{NONE=0, RAYLEIGH=1, ABCDAMP=2, THERMOVISCOUS=3,
                   FRACTIONAL=4, FRACTIONAL_GL=5, FRACTIONAL_BLANK=6,
                   FRACTIONAL_GL_INT=7, FRACTIONAL_BLANK_INT=8,
                   PML, DAMPLAYER};

  //! Interpolation type used in fractional damping model
  enum InterpolType{NOTUSED=0, trueVAL=1, LIN1PT=2};

  //! Type of nonlinearity for certain pdes
  typedef enum { NO_NONLINEARITY, WESTERVELT, KUZNETSOV, VARIABLE_SOS_CN1, 
                 VARIABLE_SOS_CN2, VARIABLE_SOS_CN2Mean,
                 MATERIAL, GEOMETRIC, HYSTERESIS, PERMEABILITY  } NonLinType;


  //! Describes all possible solution types in a CFS simulation
  typedef enum{NO_SOLUTION_TYPE, MECH_DISPLACEMENT, MECH_ACCELERATION,
               MECH_VELOCITY, MECH_FORCE, MECH_STRESS, MECH_STRAIN, MECH_STRAIN_IRR,
               MECH_ENERGY, MECH_DEF_VOLUME, MECH_RHS_LOAD, MECH_PSEUDO_DENSITY,
               ELEC_POTENTIAL, ELEC_FIELD_INTENSITY, ELEC_FORCE_VWP, 
               ELEC_INTERFACE_FORCE, ELEC_CHARGE, ELEC_FLUX_DENSITY,
               ELEC_ENERGY, ELEC_POLARIZATION,ELEC_RHS_LOAD,ELEC_PSEUDO_POLARIZATION,
               SMOOTH_DISPLACEMENT, SMOOTH_VELOCITY, GRID_VELOCITY, SMOOTH_STRAIN,  
               ACOU_POTENTIAL, ACOU_PRESSURE, ACOU_PRESSURE_DERIV_1,
               ACOU_PRESSURE_DERIV_2,ACOU_VELOCITY,ACOU_FORCE, ACOU_POT_NRBC, 
               NRBC_PHI, ACOU_POWERDENSITY, ACOU_POWER, ACOU_INTENSITY,
               ACOU_POTENTIAL_DERIV_1, ACOU_POTENTIAL_DERIV_2, ACOU_RHS_LOAD,ACOU_RHSVAL, ACOUSURF_RHSVAL,
               ACOU_BUBBLE_RHS_VAL, ACOU_SOUND_SPEEED, ACOU_SURFINTENSITY,
               MAG_POTENTIAL, MAG_FLUX_DENSITY, MAG_HFIELD, MAG_EDDY_CURRENT, 
               MAG_POTENTIAL_DIV,
               MAG_FORCE_VWP, MAG_FORCE_LORENTZ, MAG_ENERGY,
               MAG_EDDY_POWER, MAG_RHS_LOAD,
               HEAT_TEMPERATURE, HEAT_RHS_LOAD,
               MPCCI, FLUID_FORCE,
               ACOU_PRESSUREXYZ,
               LAMBDA_K, FLUIDMECH_VELOCITY, FLUIDMECH_PRESSURE,
               FLUIDMECH_DENSITY, FLUIDMECH_FORCE, FLUIDMECH_TKE, 
               FLUIDMECH_STRESS, FLUIDMECH_STRAINRATE, FLUIDMECH_ENERGY, FLUIDMECH_STABILPARAM,
               BUBBLE_RADIUS, BUBBLE_RADIUS_DERIV_1, BUBBLE_VOLUME_FRAC,
               OPT_RESULT_1, OPT_RESULT_2, OPT_RESULT_3, LAGRANGE_MULT,
               THERMOMECH_FORCE, THERMOELEC_FORCE} SolutionType;


  //! describes the possible material types
  typedef enum{NO_MATERIAL, MAG_PERMEABILITY, MAG_RELUCTIVITY, MAG_CONDUCTIVITY, 
               MAG_PERMEABILITY_1, MAG_PERMEABILITY_2, MAG_PERMEABILITY_3,
               ELEC_PERMITTIVITY, MECH_STIFFNESS_TENSOR, MECH_EMODULUS,
               MECH_POISSON, MECH_KMODULUS, MECH_GMODULUS,
               MECH_LAME_MU, MECH_LAME_LAMBDA, COEFF_STRAIN_IRREVERSIBLE,
               MECH_EMODULUS_X, MECH_EMODULUS_Y, MECH_EMODULUS_Z,
               MECH_POISSON_XY, MECH_POISSON_YZ, MECH_POISSON_XZ,
               MECH_GMODULUS_YZ, MECH_GMODULUS_ZX, MECH_GMODULUS_XY,
               RAYLEIGH_ALPHA, RAYLEIGH_BETA, RAYLEIGH_FREQUENCY, RAYLEIGH_DELTA_FREQ, 
               BOVERA, LOSS_TANGENS_DELTA,
               DENSITY, ACOU_BULK_MODULUS, ACOU_SOUND_SPEED, 
               ACOU_ALPHA, FRACTIONAL_ALG, FRACTIONAL_MEMORY, FRACTIONAL_INTERPOL,
               FRACTIONAL_EXPONENT,
               HEAT_CONDUCTIVITY, HEAT_CAPACITY, PIEZO_TENSOR, HEAT_CONDUCTIVITY_TENSOR,
               X_SATURATION, Y_SATURATION, Y_REMANENCE, PREISACH_WEIGHTS, A_JILES, 
               ALPHA_JILES, K_JILES, C_JILES, P_DIRECTION, HYST_MODEL, 
               NONLIN_COEFFICIENT, NONLIN_DEPENDENCY, NONLIN_APPROXIMATION_TYPE,
               NONLIN_DATA_NAME, DYNAMIC_VISCOSITY, KINEMATIC_VISCOSITY,
               DATA_ACCURACY, MAX_APPROX_VAL,
               THERMAL_EXPANSION_TENSOR, PYROCOEFFICIENT_TENSOR} MaterialType;

  typedef enum{FULL, PLANE_STRAIN, PLANE_STRESS, PLANE, AXI} SubTensorType;

  typedef enum{ NO_CLASS, ELECTROMAGNETIC, ELECTROSTATIC, FLUID, FLOW,
                MECHANIC, PIEZO, THERMIC, PYROELECTRIC, THERMOELASTIC } MaterialClass;


  typedef enum{ noCurve, magBH } ApproxMaterialCurves;

  //! type of measured curve to be approximated
  //! GENERAL   = any curve
  //! BH        = magnetic BH curve
  typedef enum{ GENERAL, BH } ApproxCurveType;

  //! Enumberation for coupling method\n
  //! NO_COUPLING          = No coupling at all
  //! DIRECT_COUPLING      = Direct Coupling via matrix \n
  //! ITER_RHS_COUPLING    = Iterative via RHS \n
  //! ITER_MATRIX_COUPLING = Iterative via matrix
  typedef enum{NO_COUPLING, DIRECT_COUPLING, ITER_RHS_COUPLING,
                 ITER_MATRIX_COUPLING} CouplingMethod;

  //! Enumeration for Input Coupling types \n
  //! COORD = Coupling via coordinate displacement \n
  //! RHS   = Coupling via Right hand side \n
  //! ID_BC = Coupling via inhomogenous dirichlet bc \n
  //! MAT   = Coupling via material change \n
  typedef enum {COORD, RHS, ID_BC, MAT, GRID_VEL} CouplingInputType;

  //! Enumeration for Output Coupling types\n
  //! ELEM = Coupling via element quantities\n
  //! NODE = Coupling via node quantities\n
  typedef enum {NODE, ELEM} CouplingOutputType;

  //! Enumeration for types of coupling regions\n
  //! REGION = Coupling region is whole Subdomain\n
  //! NODES = Coupling region is specified as nodes in .conf file\n
  //! SURFACE = Coupling region is specified as 1D/2D surface elements
  typedef enum {REGION, NODES, SURFACE} CouplingRegionType;

  //! Enumeration for types of norms
  //! L2ABS = absolute L2-norm
  //! L2REL = relative L2 norm: (|val| - |oldval|) / |val|
  typedef enum {NO_NORM, L2ABS, L2REL} NormType;

  //! Enumeration for directions
  //! direction of various fields 
  //! "Rad" means readial, following two letters indicate the stress-plane
  enum  Directions {X=0, Y=1, Z=2, radXY=3, radXZ=4, radYZ=5};

  //! orientation of calculation plane in 2D
  //! (especially important for anisotropic simulations)
  enum orientation2D {xy, xz, yz};

  //! nonlinear method definition 
  enum NonLinMethod {FIXEDPOINT=1, NEWTON=2};

  //! output format for complex numbers
  typedef enum {REAL_IMAG, AMPLITUDE_PHASE} ComplexFormat;


  //! type of data
  typedef enum {INTEGER, REAL, IMAG, COMPLEX} DataType;

  //! Data type for specification of frequency sampling approach

  //! This enumeration data type is used for distinguishing the different
  //! approaches for sampling the frequency domain when performing a
  //! harmonic analysis with CFS++. The type currently allows the following
  //! values
  //! - NO_SAMPLING_TYPE
  //! - LINEAR_SAMPLING
  //! - LOG_SAMPLING
  //! - REVERSE_LOG_SAMPLING
  typedef enum { NO_SAMPLING_TYPE, LINEAR_SAMPLING, LOG_SAMPLING,
                 REVERSE_LOG_SAMPLING } FreqSamplingType;


  //--------------------- Stuff for handling different IO files -------------

  typedef enum {DEBUG_FILE, MEMTRACE_FILE, OLAS_FILE} AuxFileType;


  //------------------------ Stuff for bubble simulation --------------------

  typedef enum {NOBUBBLETYPE, KELLERMIKSIS, GILMORE} BubbleDynType;

  extern BubbleDynType bubbleDyn;
  

 

  //------------------------ Files for debug and information ---------

  // NOTE: OLAS uses the namespace 'OutInfo' for writing out data into the
  // different filestreams such as (*cla) etc. Therefore they are
  // explicitely imported into namespace CoupledField at this point
#ifndef INTEGLIB
  using OutInfo::debug;
  using OutInfo::cla;
  using OutInfo::memtrace;
#endif
  using OutInfo::error;
  using OutInfo::warning;

  // Forward declaration of class
  class WriteInfo;

  //! Global pointer to class performing logging to info file
  extern WriteInfo *Info;

  class ConfFile;
  //! pointer to class with methods for reading config-file. it is accessable
  //! from any place of the program
  //  extern ConfFile *conf; //name.conf

  class BaseFE;
  //! Pointers to derived classes of BaseElem. Initialized in grid.hh/grid.cc.
  //! They are used, when we read information about elements from mesh and
  //! create a pointer to the class containing the description of the Finite
  //! Element.
  extern BaseFE *ptQ1, *ptQ2, *ptQ9, *ptL1, *ptL2, *ptTet1, *ptTet2, *ptTr1, 
	*ptTr2, *ptHexa1, *ptHexa2, *ptHexa27, *ptPyra1, *ptPyra2, *ptWedge1, *ptWedge2;

  //! conversion from strings to enum types
  template <class TYPE>
  void String2Enum(const std::string &in, TYPE &out);

  //! conversion from enum types to strings
  template<class TYPE>
  void Enum2String(const TYPE &in, std::string &out);
 
  // Instantiation for all known enum types;
#if defined (__GNUC__)
#define DEFINE_ENUM_CONVERSION(TYPE)                                  \
  template<typename TYPE> void String2Enum(const std::string &in, TYPE &out); \
  template<typename TYPE> void Enum2String(const TYPE &in, std::string &out);

  DEFINE_ENUM_CONVERSION(AnalysisType)
  DEFINE_ENUM_CONVERSION(FreqSamplingType)
  DEFINE_ENUM_CONVERSION(CouplingInputType)
  DEFINE_ENUM_CONVERSION(CouplingOutputType)
  DEFINE_ENUM_CONVERSION(CouplingRegionType)
  DEFINE_ENUM_CONVERSION(NormType)
  DEFINE_ENUM_CONVERSION(ComplexFormat)
  DEFINE_ENUM_CONVERSION(EQNType)
  DEFINE_ENUM_CONVERSION(FEType)
  DEFINE_ENUM_CONVERSION(EntryType::ScalarType)
  DEFINE_ENUM_CONVERSION(DataType)
  DEFINE_ENUM_CONVERSION(MaterialClass)
  DEFINE_ENUM_CONVERSION(IntegrationMethod)
  DEFINE_ENUM_CONVERSION(NonLinType)
  DEFINE_ENUM_CONVERSION(DampingType)
#endif

  std::string MapSolTypeToUnit(SolutionType solType);

#ifdef INTEGLIB
  typedef enum {NOTYPE, SYSTEM, STIFFNESS, DAMPING, CONVECTION, MASS}
  FEMatrixType;
  
  //! class for flags of the FE matrix
	class FEMatrix_Flags {
		public:
			FEMatrix_Flags() {
				setCounterPart=false;
				setTransposeInt = true;
				setOnlyCounterPart = false;
			}
	
			//! Flag indicating assembling of the integrator
			//! in the counterpart of the pde location
			bool setCounterPart;
			//! Flag indicating the assembling of the integrator.
			//! in the counter part of the pde transposing it.
			//! Note: By default, we set the transpose of a matrix 
			//! true when assembling the counter part of the 
			//! matrix, i.e. the case of piezoelectric coupling.
			//! if we want set the counter part without transposing
			//! the integrator within the global FE-matrix, this 
			//! flag must be changed through 'SetTransposeInt()'
			bool setTransposeInt;
			//! Flag to set only the counter part of the element matrix.
			bool setOnlyCounterPart;
		
	};  
	//extern FEMatrix_Flags * matrix_flags;
	
#endif 
  
	
	
} // end of namespace

#endif // FILE_SCFE_MYDEFS
