#ifndef FILE_INTERFACE_AdaptGrid_2002
#define FILE_INTERFACE_AdaptGrid_2002

#include "filetype.hh"
#include "grid.hh"

// from Adapt Grid
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "Tetrahedron.h"
#include "Octahedron.h"
#include "MultilevelGrid.h"
#include "GeometrySensor.h"
#include "MeshReader.h"
#include "TetrahedronMeasure.h"
#include "MeshWriter.h"

namespace CoupledField
{
 
/// Interface to library Grid by Roberto G.
template<class Dim>
class InterfaceAdaptGrid: public Grid
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  InterfaceAdaptGrid(FileType * const aptFileType);

  /// Deconstructor
  virtual ~InterfaceAdaptGrid();
  
   /// Uniform subdivision of domain
  virtual void SubdivideUniform(const Integer level);

   /// Get coordinates of node with global number inode
   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint);

   /// Get connection of element
 virtual void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level);

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel);

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel);

  /// Put information about initial grid in mesh
  virtual void Read();

  ///
  virtual Integer GetDim() { return dim_; } // 
   
  //! Here we mark elements for refinement: ei - number of elem
  virtual void SetRefinementFlag(const Integer ei);
  virtual void SetRefinementFlag(const std::vector<Integer> ei);  

  //! Do refinement of elements, which we mark through function SetRefinementFlag
  virtual void Refine();

private:
  //! 
  FileType * ptFileType;
      
  //!
  grd::MultilevelGrid grid_;

  //!
  Integer dim_;

  //!
  std::vector<grd::Vertex*> vertex_;
  std::vector<grd::Element*> elems_;

};

template<class Dim>
inline InterfaceAdaptGrid<Dim>::InterfaceAdaptGrid(FileType * aptFileType)
: Grid(aptFileType)
{
#ifdef TRACE
 (*trace) << "Entering InterfaceAdaptGrid<Dim>::InterfaceAdaptGrid<Dim>" << std::endl;
#endif

  ptFileType=aptFileType;
  lastlevel_=0; 
}

#ifdef __GNUC__
template class InterfaceAdaptGrid<Point3D>;
template class InterfaceAdaptGrid<Point2D>;
#endif

} // end of namespace
#endif // 
