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

//#ifndef assert
//#define assert(ex) \
//        (void)((ex) ? 1 : \
//              (_error("Failed assertion " #ex " at line %d of `%s'.\n", \
//               __LINE__, __FILE__), 0))
//#endif

#define TRIANGLE1 1
#define TRIANGLE2 2
#define QUADRILATERAL1 3
#define QUADRILATERAL2 4
#define TETRAHEDRAL1 11
#define TETRAHEDRAL2 12

//#define TETRAHEDRA1 10
//#define TETRAHEDRA2 11
//#define HEXAHEDRA1 12
//#define HEXAHEDRA2 13
enum ElementType {Triang1, Triang2, Quadrilateral1, Quadrilateral2};

enum IntegrationType {GaussOrder2, GaussOrder3, GaussOrder4, GaussOrder5, GaussOrder7};

/// for each level of remeshing we create own GridHierarchy with full information about grid
template <class Dim>
struct GridHierarchy
{
    Integer maxnumelem;
    Integer maxnumnode;
    Integer * Info; // Array of full info about element
    Dim * ptCoordinate;
    Integer * Connect;
    Integer * fp; //Array of first position in Info

};
//-------------------------  Used enumerations: ---------------------------- 
// enum for precondition matrix
 enum precond { non, Jacobi, SSOR, LU}; 

std::ostream & operator << (std::ostream & out, const enum precond & type);

// ------------------------ Files for debug, trace and information ---------
extern std::ostream * trace; //name.trace
extern std::ostream * debug; //name.deb
extern std::ostream * infofile; //name.info

// if InfoPrint=true, then we create file with information about methods, which are used 

extern Boolean InfoPrint;

}

#endif // FILE_SCFE_MYDEFS
