
#ifndef FILE_INTERFACE_GRIDLIB_2002
#define FILE_INTERFACE_GRIDLIB_2002

#include "filetype.hh"
#include "grid.hh"

// Include from Gridlib
#include <GoMesh.hh>
#include <GoTriangleMesh.hh>
#include <GoDefaultVertex.hh>
#include <GoTriangleElement.hh>

namespace CoupledField
{
 
/// Class for working with grid
template<class Dim> 
class InterfaceGridlib: public Grid<Dim>
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  InterfaceGridlib(FileType * const aptFileType);

  /// Deconstructor
  virtual ~InterfaceGridlib() { if (ptGoMesh ) delete ptGoMesh ;}
  
  /// Get coordinates of all nodes which belong to element
  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodes, Dim * ptCoordElem); 

   /// Get connection of element
  virtual void GetConnection(Integer * result, const Integer level, 
           const Integer numElem, const Integer numnodesPerElem);

  /// Return pointer to coordinates
  virtual Dim * GetptCoordinate(const Integer numlevel)
  { Error("Not implemented  GetptCoordinate"); return NULL;}

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel)
  { 
     return ptGoMesh->getNumVertices(numlevel);
  }

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel)
  { 
    return ptGoMesh->getNumElements(numlevel);
     }

  /// Return num of nodes per element i
  virtual Integer GetNumNodesPerElem(const Integer iElem, const Integer level)
  { 
//  { GoGeometryElement<float> * ptElem=ptGoMesh->getElement(iElem, level);
//      ptElem->getNumVertices();}
    return ptGoMesh->getElement(iElem, level)->getNumVertices();
  }

  /// Print coordinates of grid in out
  virtual void PrintCoordinate(const Integer level, std::ostream * out) const
  { Error("Not implemented yet  PrintCoordinate");}

  /// Put information about grid
  void Read();
private:

  //! 
  FileType * ptFileType;
      
  //! 
  GoMesh * ptGoMesh;
};

template<class Dim>
inline InterfaceGridlib<Dim>::InterfaceGridlib(FileType * aptFileType)
: Grid<Dim>(aptFileType)
{
#ifdef TRACE
 (*trace) <<
             "Entering InterfaceGridlib<Dim>::InterfaceCFS<Dim>" 
                << std::endl;
#endif
  ptFileType=aptFileType;
}

#ifdef __GNUC__
template class InterfaceGridlib<Point3D>;
template class InterfaceGridlib<Point2D>;
#endif

} // end of namespace
#endif // 
