#ifndef FILE_SCFE_MYDEFS_2001
#define FILE_SCFE_MYDEFS_2001

namespace CoupledField
{

typedef int Integer;
typedef short ShortInt;
typedef float Float;
typedef double Double;
typedef char Char;
typedef int Boolean;

#define FALSE 0
#define TRUE 1

#define mark std::cout<<__FILE__<<__LINE__<<std::endl;

  typedef Double (*pfn1var)(const Double);
  typedef Double (*pfn2var)(const Double, const Double);
  typedef Double (*pfn3var)(const Double, const Double, const Double);

// #ifndef assert
// #define assert(ex) \
//        (void)((ex) ? 1 : \
//              (_error("Failed assertion " #ex " at line %d of `%s'.\n", \
//               __LINE__, __FILE__), 0))
// #endif

  enum ElementType{Line1, Triang1, Triang2, Quadrilateral1, Quadrilateral2};
  enum IntegrationType {GaussOrder2, GaussOrder3, GaussOrder4, GaussOrder5, GaussOrder7, null};
 enum precond { non, Jacobi, SSOR, LU}; 

// ---------------solver specific enumerations: begin
enum GraphType   {NOGRAPH, NODEGRAPH, EDGEGRAPH, FACEGRAPH, VOLUMEGRAPH};
enum FEType      {NOFETYPE, TRIA, QUAD, TET, HEX};
enum MatrixType  {NOTYPE, SYSTEM, STIFFNESS, DAMPING, CONVECTION, MASS};
enum MatrixClass {NOCLASS, RSCALAR, CSCALAR, RBLOCK, CBLOCK, RFULL, CFULL, RSCALAREDGE, CSCALAREDGE, HBFEM};
enum PrecondType {NOPRECOND, ID, MG};
enum SolverType  {NOSOLVER, RealRichardson, ComplexRichardson, RealCG, ComplexCG, Lanczos, RealQMR, 
		  ComplexQMR};
enum CycleType   {NOCYCLE, VCYCLE, WCYCLE};

// ---------------solver specific enumerations: end

std::ostream & operator << (std::ostream & out, const enum precond & type);

// ------------------------ Files for debug, trace and information ---------
extern std::ostream * trace; //name.trace
extern std::ostream * debug; //name.deb
extern std::ostream * infofile; //name.info
extern std::ostream * cla; //name.cla

class ConfFile;
extern ConfFile * conf; //name.conf

extern Boolean InfoPrint;

class BaseElem;
extern BaseElem * ptQ, *ptTr, *ptTet, *ptL1, *ptHexa;

}

#endif // FILE_SCFE_MYDEFS
