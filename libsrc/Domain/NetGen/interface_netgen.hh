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
class InterfaceNetGen: public Grid
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
 virtual void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level);

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel);

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel);

  /// Return num of nodes per element i
  Integer GetNumNodesPerElem(const Integer iElem, const Integer level);

  /// Put information about initial grid in mesh
  virtual void Read();

  ///
  virtual BaseElem * GetptElem(const Integer iElem);

  ///
  virtual Integer GetDim() { return dim_; } // mesh.GetDimension();}
   
  //! Here we mark elements for refinement: ei - number of elem
  virtual void SetRefinementFlag(const Integer ei);
  virtual void SetRefinementFlag(Vector<Integer> & ei);  

  //! Do refinement of elements, which we mark through function SetRefinementFlag
  virtual void Refine();

private:
  //! 
  FileType * ptFileType;
      
  //! 
  Mesh mesh;

  //!
  Integer dim_;

   //!
  void Init();  

  //! if we do subdivision, then this variable is TRUE
  Boolean DoesGridSubdivide;

  //! array of pointers to BaseElem
  std::vector<BaseElem*> allptElem;
};

template<class Dim>
inline InterfaceNetGen<Dim>::InterfaceNetGen(FileType * aptFileType)
: Grid(aptFileType)
{
#ifdef TRACE
 (*trace) << "Entering InterfaceNetGen<Dim>::InterfaceNetGen<Dim>" << std::endl;
#endif

  ptFileType=aptFileType;
  DoesGridSubdivide=FALSE;
  lastlevel_=0; 

  Init();
}

#ifdef __GNUC__
template class InterfaceNetGen<Point3D>;
template class InterfaceNetGen<Point2D>;
#endif

} // end of namespace
#endif // 
