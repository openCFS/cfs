#ifndef FILE_GRID_INTERFACE_CFS_2002
#define FILE_GRID_INTERFACE_CFS_2002

#include "grid.hh"
#include "grid_cfs.hh"

namespace CoupledField
{
 
class FileType;

/// Class for working with grid
template<class Dim> 
class GridInterfaceCFS: public Grid<Dim>
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  GridInterfaceCFS(FileType * const aptFileType);

  /// Deconstructor
  virtual ~GridInterfaceCFS() { if (ptGridCFS) delete ptGridCFS;}
  
  /// Uniform subdivision of domain
  virtual void Subdiv(const Integer level) {Error("Not implemented Subdiv()");} 

  //! Read of mesh 
  virtual void Read() 
  { ptGridCFS->Read();}

  /// Get coordinates of all nodes which belong to element
  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodes,  Dim * ptCoordElem) 
  { ptGridCFS->GetCoordOfNodesElem(numElem, numlevel, numnodes,ptCoordElem);}

   /// Get connection of element
  virtual void GetConnection(Integer * result, const Integer level, 
           const Integer numElem, const Integer numnodesPerElem)
  { ptGridCFS->GetConnection(result, level, numElem, numnodesPerElem);}

   /// Get coordinates of node with global number inode
  virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint) {ptGridCFS->GetCoordinateNode(inode,numlevel,rfPoint);}

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel)
  { return ptGridCFS->GetMaxnumnodes(numlevel); }

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel)
  { return ptGridCFS->GetMaxnumElem(numlevel);}

  /// Return num of nodes per element i
  virtual Integer GetNumNodesPerElem(const Integer iElem, const Integer level)
  { return ptGridCFS->GetNumNodesPerElem(iElem, level);}

  /// Print coordinates of grid in out
  virtual void PrintCoordinate(const Integer level, std::ostream * out) const
  { ptGridCFS->PrintCoordinate(level, out);}

protected:
private:
  GridCFS<Dim> * ptGridCFS;
  ///
};

template<class Dim>
inline GridInterfaceCFS<Dim>::GridInterfaceCFS(FileType * aptFileType)
: Grid<Dim>(aptFileType)
{
#ifdef TRACE
 (*trace) << "Entering GridInterfaceCFS<Dim>::GridInterfaceCFS<Dim>" << std::endl;
#endif
   ptGridCFS=new GridCFS<Dim>(ptFileType); 
}

#ifdef __GNUC__
template class GridInterfaceCFS<Point3D>;
template class GridInterfaceCFS<Point2D>;
#endif

} // end of namespace
#endif // FILE_GRID
