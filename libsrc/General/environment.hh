#ifndef FILE_SCFE_MYDEFS_2001
#define FILE_SCFE_MYDEFS_2001

#include <typeinfo>
#include <iostream>
#include <vector>
//#include <math.h>
#include <cmath>
#include <complex>
#include "General/defs.hh"

#ifdef USE_OLAS
#include <olas.hh>
#else
#include <las_environment.hh>
#endif


//! \file environment.hh
//! This file contains some global macro, class and enumeration data type
//! definitions for CFS++.


//! Maximal length of the trailing postfix of an auxilliary name,
//! i.e. the length of the extension after the basename
#define MAXPOSTFIX 15


namespace CoupledField
{

  //! redeclaration of types
  typedef int Integer;
  typedef unsigned int UInt;
  typedef short ShortInt;
  typedef float Float;
  typedef double Double;
  typedef std::complex<Double> Complex;
  typedef char Char;
  typedef int Boolean;

#define FALSE 0
#define TRUE 1

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


  //! declaration sof functions. it is used in parsing functions from conf-file
  typedef Double (*pfn1var)(const Double);
  typedef Double (*pfn2var)(const Double, const Double);
  typedef Double (*pfn3var)(const Double, const Double, const Double);

  typedef enum {STATIC, TRANSIENT, HARMONIC, EIGENFREQUENCY, MULTI_SEQUENCE}
  AnalysisType;
  

  //! print grid only and then exit
  extern Boolean PrintGridOnly;


  //! enumeration with elements types.
  enum ElementType{Line1, Triang1, Triang2, Quadrilateral1, Quadrilateral2};

  //! enumeration with integration types. it is used in Elements classes
  enum IntegrationType {GaussOrder1, GaussOrder2, GaussOrder3, GaussOrder4,
			GaussOrder5, GaussOrder7, null};

  //! Damping type
  enum DampingType{NONE=0, RAYLEIGH=1, FRACTIONAL=2, ABCDAMP=3,
		   THERMOVISCOUS=4};

  //! Describes all possible solution types in a CFS simulation
  typedef enum{NO_SOLUTION_TYPE, 
		 MECH_DISPLACEMENT, MECH_ACCELERATION,
		 MECH_VELOCITY, MECH_FORCE, MECH_STRESS, MECH_STRAIN,
		 ELEC_POTENTIAL, ELEC_FIELD, ELEC_FORCE_VWP, 
		 ELEC_INTERFACE_FORCE, ELEC_CHARGE,  ELEC_FLUX_DENSITY,
		 SMOOTH_DISPLACEMENT, 
		 ACOU_POTENTIAL, ACOU_FORCE, 
		 ACOU_POTENTIAL_DERIV_1, ACOU_POTENTIAL_DERIV_2,
		 MAG_POTENTIAL, MAG_FIELD, MAG_EDDY_CURRENT, 
		 MAG_FORCE_VWP, MAG_FORCE_LORENTZ}
  SolutionType;

  //! Enumberation for coupling method\n
  //! NO_COUPLING          = No coupling at all
  //! DIRECT_COUPLING      = Direct Coupling via matrix\n
  //! ITER_RHS_COUPLING    = Iterative via RHS\n
  //! ITER_MATRIX_COUPLING = Iterative via matrix
  typedef enum{NO_COUPLING, DIRECT_COUPLING, 
		 ITER_RHS_COUPLING, ITER_MATRIX_COUPLING}
  CouplingMethod;

  //! Enumeration for Input Coupling types \n
  //! COORD = Coupling via coordinate displacement\n
  //! RHS   = Coupling via Right hand side\n
  //! ID_BC = Coupling via inhomogenous dirichlet bc\n
  //! MAT   = Coupling via material change\n
  typedef enum {COORD, RHS, ID_BC, MAT} CouplingInputType;

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

  //! Enumeration for type of equation numbering
  typedef enum {NODE_SCALAR, NODE_BLOCK, NODE_SCALAR_BLOCK,
		NODE_SUPER_BLOCK} EQNType;

  //! Enumeration for directions
  //! direction of various fields 
  //! "Rad" means readial, following two letters indicate the stress-plane
  enum  Directions {X, Y, Z, radXY, radXZ, radYZ};

  //! orientation of calculation plane in 2D
  //! (especially important for anisotropic simulations)
  enum orientation2D {xy, xz, yz};

  //! nonlinear method definition 
  enum NonLinMethod {FIXEDPOINT=1, NEWTON=2};

  //! output format for complex numbers
  typedef enum {REAL_IMAG, AMPLITUDE_PHASE} ComplexFormat;
  

//------------------------ Files for debug, trace and information ---------

  // NOTE: OLAS uses the namespace 'OutInfo' for writing out
  // data into the cla, trace, ... stream. Therefore they explicitely have to
  // be imported into namespace CoupledField.
#ifdef USE_OLAS
  using OutInfo::trace;
  using OutInfo::debug;
  using OutInfo::cla;
  using OutInfo::memtrace;
  using OutInfo::data;
#else
  extern std::ostream * trace; //name.trace
  extern std::ostream * debug; //name.deb
  extern std::ostream * cla; //name.cla
  extern std::ostream * memtrace; //name for name.mem 
  extern std::ostream * data; //name.data
#endif

  class WriteInfo;
  //! Global pointer to class performing logging to info file
  extern WriteInfo *Info;

  class ConfFile;
  //! pointer to class with methods for reading config-file. it is accessable
  //! from any place of the program
  extern ConfFile *conf; //name.conf

  class BaseParamHandler;
  //! Global pointer to class performing handling of steering parameters
  extern BaseParamHandler *params;

  //! class BaseElem;
  class BaseFE;

  //! Pointers to derived classes of BaseElem. Initialized in grid.hh/grid.cc.
  //! They are used, when we read information about elements from mesh and
  //! create a pointer to the class containing the description of the Finite
  //! Element.
  extern BaseFE *ptQ, *ptQ2, *ptL1, *ptL2, *ptTet, *ptTr1, *ptTr2, *ptHexa, *ptPyra, *ptWedge;

  //! class for flags of programm
  class Flags
  {
  public:
    Flags()
    { CalcErrorMap_=FALSE;
    adaptSpace_ = FALSE;
    }
 
    Boolean CalcErrorMap_;
    Boolean adaptSpace_;
  };
  
  extern Flags * flags;

  //! parameters necessary to describe coils

  enum COILTYPE {MEASUREMENT,CURRENT,VOLTAGE};
  
  struct coilDefStruct
  {
    int ID;
    double  current;
    double  voltage;
    double  coilArea;
    double  resistance;
    double  phase;
    std::string timefnc;
    COILTYPE type;
    std::string Lfile;
    std::string UIfile;
  };


  //! conversion from strings to enum types
  template <class TYPE>
  void String2Enum(const std::string &in, TYPE &out);

  //! conversion from enum types to strings
  template<class TYPE>
  void Enum2String(const TYPE &in, std::string &out);

  // Instantiation for all known enum types;
#ifdef __GNUC__
#define DEFINE_ENUM_CONVERSION(TYPE)                                  \
  template<class TYPE> void String2Enum(const std::string &in, TYPE &out); \
  template<class TYPE> void Enum2String(const TYPE &in, std::string &out);

  DEFINE_ENUM_CONVERSION(AnalysisType);
  DEFINE_ENUM_CONVERSION(CouplingInputType);
  DEFINE_ENUM_CONVERSION(CouplingOutputType);
  DEFINE_ENUM_CONVERSION(CouplingRegionType);
  DEFINE_ENUM_CONVERSION(NormType);
  DEFINE_ENUM_CONVERSION(ComplexFormat);
  DEFINE_ENUM_CONVERSION(EQNType);
}
#endif



#endif // FILE_SCFE_MYDEFS
