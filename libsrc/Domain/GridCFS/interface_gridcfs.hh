#ifndef FILE_GRID_INTERFACE_CFS_2002
#define FILE_GRID_INTERFACE_CFS_2002

#include "grid.hh"
#include "grid_cfs.hh"

namespace CoupledField
{
 
class FileType;

/// Class for working with grid
template<class Dim> 
class GridInterfaceCFS: public Grid
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  GridInterfaceCFS(FileType * const aptFileType);

  /// Deconstructor
  virtual ~GridInterfaceCFS() { if (ptGridCFS) delete ptGridCFS;}
  
  /// Uniform subdivision of domain
  virtual void SubdivideUniform(const Integer level) {Error("Not implemented Subdiv() in GridCFS",__FILE__,__LINE__);} 

  //! Read of mesh 
  virtual void Read() 
  { ptGridCFS->Read();}

  /// Get coordinates of all nodes which belong to element
  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodes,  Dim * ptCoordElem) 
  { ptGridCFS->GetCoordOfNodesElem(numElem, numlevel, numnodes,ptCoordElem);}

   /// Get connection of element
   virtual void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)
    { ptGridCFS->GetConnection(connect,iElem, level);}

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

  /// return pointer to BaseElem
    BaseElem * GetptElem(const Integer iElem)
 { return ptGridCFS->GetptElem(iElem);}

  /// return dimension of mesh
  virtual Integer GetDim(){ ptGridCFS->GetDim();}

protected:
private:
  GridCFS<Dim> * ptGridCFS;
  ///
};

template<class Dim>
inline GridInterfaceCFS<Dim>::GridInterfaceCFS(FileType * aptFileType)
: Grid(aptFileType)
{
#ifdef TRACE
 (*trace) << "Entering GridInterfaceCFS<Dim>::GridInterfaceCFS<Dim>" << std::endl;
#endif
 lastlevel_=0;
   ptGridCFS=new GridCFS<Dim>(ptFileType); 
}

#ifdef __GNUC__
template class GridInterfaceCFS<Point3D>;
template class GridInterfaceCFS<Point2D>;
#endif

} // end of namespace
#endif // FILE_GRID
