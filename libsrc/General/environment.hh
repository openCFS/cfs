#ifndef FILE_SCFE_MYDEFS_2001
#define FILE_SCFE_MYDEFS_2001

#include <las_environment.hh>

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

  //! useful trick for testing problem
#define mark std::cout<<__FILE__<<__LINE__<<std::endl;

  //! declaration sof functions. it is used in parsing functions from conf-file
  typedef Double (*pfn1var)(const Double);
  typedef Double (*pfn2var)(const Double, const Double);
  typedef Double (*pfn3var)(const Double, const Double, const Double);

  enum AnalysisType {STATIC=1, TRANSIENT=2, HARMONIC=3, EIGENFREQUENCY=4};

// #ifndef assert
// #define assert(ex) \
//        (void)((ex) ? 1 : \
//              (_error("Failed assertion " #ex " at line %d of `%s'.\n", \
//               __LINE__, __FILE__), 0))
// #endif

  //! enumeration with elements types.
  enum ElementType{Line1, Triang1, Triang2, Quadrilateral1, Quadrilateral2};
  //! enumeration with integration types. it is used in Elements classes
  enum IntegrationType {GaussOrder2, GaussOrder3, GaussOrder4, GaussOrder5, GaussOrder7, null};
  //! enumeration precondition's types. it is used in methods of LinAlg
 enum precond { non, Jacobi, SSOR, LU}; 

std::ostream & operator << (std::ostream & out, const enum precond & type);

// ------------------------ Files for debug, trace and information ---------
extern std::ostream * trace; //name.trace
extern std::ostream * debug; //name.deb
extern std::ostream * infofile; //name.info
extern std::ostream * cla; //name.cla

class ConfFile;
  //! pointer to class with methods for reading config-file. it is accessable from any place of the program
extern ConfFile * conf; //name.conf

  //! indicator for info-file. TRUE, if we need it, FALSE, otherwise
extern Boolean InfoPrint;

  //class BaseElem;
  class BaseFE;

  //! pointers to derived classes of BaseElem. it is initialized in grid.hh(grid.cc). it is used, when we read information about elements from mesh and create pointer to class with description FE element.
  //extern BaseElem * ptQ, *ptTr, *ptTet, *ptL1, *ptHexa;
  extern BaseFE * ptQ;

}

#endif // FILE_SCFE_MYDEFS
