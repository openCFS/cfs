#ifndef FILE_INTERFACE_NetGen_2002
#define FILE_INTERFACE_NetGen_2002

#include "filetype.hh"
#include "grid.hh"

// from NetGen
#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace CoupledField
{
 
/// Interface to library NetGen
template<class Dim>
class InterfaceNetGen: public Grid<Dim>
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  InterfaceNetGen(FileType * const aptFileType);

  /// Deconstructor
  virtual ~InterfaceNetGen();
  
   /// Uniform subdivision of domain
  virtual void SubdivideUniform(const Integer level);

  /// Get coordinates of all nodes which belong to element
  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodes, Dim * ptCoordElem); 

   /// Get coordinates of node with global number inode
   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint);

   /// Get connection of element
  virtual void GetConnection(Integer * result, const Integer level, 
           const Integer numElem, const Integer numnodesPerElem);

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel);

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel);

  /// Return num of nodes per element i
  virtual Integer GetNumNodesPerElem(const Integer iElem, const Integer level);

  /// Print coordinates of grid in out
  virtual void PrintCoordinate(const Integer level, std::ostream * out) const;

  /// Put information about grid
  virtual void Read();

  //! Put global numbers of nodes with boundary condition in Vector
  virtual void GetNodesBoundaryCondition(Vector<Integer> & nodesDirBC, const Integer level);

  //! Get array of pointers to element type
  virtual BaseElem ** getptArrayElem() const
  { return ptArrayElem_;}
private:

  //!
  BaseElem * ptQ_, * ptTr_;  

  //! 
  FileType * ptFileType;
      
  //! 
  Mesh mesh;

  //! if we do subdivision, then this variable is TRUE
  Boolean DoesGridSubdivide;

  //!
  BaseElem ** ptArrayElem_;  

};

template<class Dim>
inline InterfaceNetGen<Dim>::InterfaceNetGen(FileType * aptFileType)
: Grid<Dim>(aptFileType)
{
#ifdef TRACE
 (*trace) << "Entering InterfaceNetGen<Dim>::InterfaceNetGen<Dim>" << std::endl;
#endif

  ptFileType=aptFileType;
  DoesGridSubdivide=FALSE;
  ptQ_=NULL;
  ptTr_=NULL;
  lastlevel_=0;
 
}

#ifdef __GNUC__
//template class InterfaceNetGen<Point3D>;
template class InterfaceNetGen<Point2D>;
#endif

} // end of namespace
#endif // 
