#ifndef FILE_SCFE_MYDEFS_2001
#define FILE_SCFE_MYDEFS_2001

#include <iostream>
#include <las_environment.hh>
#include <vector>
#include <math.h>

namespace CoupledField
{

  //! redeclaration of types
typedef int Integer;
typedef short ShortInt;
typedef float Float;
typedef double Double;
typedef char Char;
typedef int Boolean;

#define FALSE 0
#define TRUE 1
#define myEndl std::endl
#define myCout std::cout


  //! useful trick for testing problem
#define mark std::cout<<__FILE__<<__LINE__<<std::endl;

  const Double PI = acos(-1.0);
  

  //! declaration sof functions. it is used in parsing functions from conf-file
  typedef Double (*pfn1var)(const Double);
  typedef Double (*pfn2var)(const Double, const Double);
  typedef Double (*pfn3var)(const Double, const Double, const Double);

  enum AnalysisType {STATIC=1, TRANSIENT=2, HARMONIC=3, EIGENFREQUENCY=4};

  //! print grid only and then exit
  extern Boolean PrintGridOnly;

// #ifndef assert
// #define assert(ex) \
//        (void)((ex) ? 1 : \
//              (_error("Failed assertion " #ex " at line %d of `%s'.\n", \
//               __LINE__, __FILE__), 0))
// #endif

  //! enumeration with elements types.
  enum ElementType{Line1, Triang1, Triang2, Quadrilateral1, Quadrilateral2};
  //! enumeration with integration types. it is used in Elements classes
  enum IntegrationType {GaussOrder1, GaussOrder2, GaussOrder3, GaussOrder4, GaussOrder5, GaussOrder7, null};
  //! enumeration precondition's types. it is used in methods of LinAlg
  enum precond { non, Jacobi, SSOR, LU}; 

  std::ostream & operator << (std::ostream & out, const enum precond & type);
  
  std::ostream& operator << (std::ostream & outStr, std::vector<Double> xOut);
  std::ostream& operator << (std::ostream & outStr, std::vector<Integer> xOut);

  // Enumeration for Input Coupling types
  //   COORD = Coupling via coordinate displacement
  //   RHS   = Coupling via Right hand side
  //   ID_BC = Coupling via inhomogenous dirichlet bc
  //   MAT   = Coupling via material change
  enum CouplingInputType{COORD, RHS, ID_BC, MAT};

  // Enumeration for Output Coupling types
  //   ELEM = Coupling via element quantities
  //   NODE = Coupling via node quantities
  enum CouplingOutputType{NODE, ELEM};

  // Enumeration for types of coupling regions
  //   SUBDOMAIN = Coupling region is whole Subdomain
  //   NODES = Coupling region is specified as nodes in .conf file
  //   ELEMS1D = Coupling region is specified as 1D-Interface
  //   ELEMS2D = Coupling region is specified as 2D-interface
  enum CouplingRegionType{SUBDOMAIN, NODES, ELEMS1D, ELEMS2D};

  // Enumeration for types of norms
  //   L2ABS = absolute L2-norm
  //   L2REL = relative L2 norm: (|val| - |oldval|) / |val|
  enum NormType {L2ABS, L2REL};

  // Enumeration for directions
  // direction of various fields 
  // "Rad" means readial, following two letters indicate the stress-plane
  enum Directions {X, Y, Z, radXY, radXZ, radYZ};

  std::ostream & operator << (std::ostream & out, const enum precond & type);


// ------------------------ Files for debug, trace and information ---------
extern std::ostream * trace; //name.trace
extern std::ostream * debug; //name.deb
extern std::ostream * infofile; //name.info
extern std::ostream * cla; //name.cla
extern std::ostream * memtrace; //name for name.mem 
extern std::ostream * data; //name.data

  class WriteInfo;
  extern WriteInfo * Info; //class for log informations

class ConfFile;
  //! pointer to class with methods for reading config-file. it is accessable from any place of the program
extern ConfFile * conf; //name.conf

  //class BaseElem;
  class BaseFE;

  //! pointers to derived classes of BaseElem. it is initialized in grid.hh(grid.cc). 
  /// It is used, when we read information about elements from mesh and create pointer to class with description FE element.
  //extern BaseElem * ptQ, *ptTr, *ptTet, *ptL1, *ptHexa;

  extern BaseFE * ptQ, *ptL1, *ptTet, *ptTr1, *ptHexa;

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
}

#endif // FILE_SCFE_MYDEFS
