#ifndef FILE_SCFE_GRID_2001
#define FILE_SCFE_GRID_2001

#include "filetype.hh"
#include "baseelem.hh"

namespace CoupledField
{
 
/// Class for working with grid
class Grid
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  Grid(FileType * aptFileType) { ptFileType=aptFileType; }

  /// Deconstructor
  virtual ~Grid() { ;}

  //!
  virtual void Read()=0;

  /// 
  virtual void SubdivideUniform(const Integer level)=0;

  /// Get coordinates of all nodes which belong to element
  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevelGrid, const Integer numnodes, Point2D * ptCoordElem)
{ Error("Not implemented");}  
  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevelGrid, const Integer numnodes, Point3D * ptCoordElem)
{ Error(" Not implemented");}

   /// Get connection of element
  virtual void GetConnection(Integer * result, const Integer levelGrid, 
           const Integer numElem, const Integer numnodesPerElem)=0;

   /// Get coordinates of node with global number inode

   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point2D & rfPoint){ Error(" Not implemented");}
   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point3D & rfPoint){ Error(" Not implemented");}

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel)=0;

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel)=0;

  /// Return num of nodes per element i
  virtual Integer GetNumNodesPerElem(const Integer iElem, const Integer level)=0;

  //! Get array of nodes for boundary condition 
  virtual void GetNodesBoundaryCondition(Vector<Integer> & nodesDirBC, const Integer level)=0;

  //! Get array of pointers to element type
  virtual BaseElem ** getptArrayElem() const=0;

  //! Get last level of grid
  virtual Integer GetLastLevel() const { return lastlevel_;} 

  //! Get number of subdomains
  virtual Integer GetNumSubdomains() const=0; 
  
  //! Get pointer to array with nodes, that belongs to subdomain number num
  virtual Integer * GetElemSubdomain(const Integer num, const Integer level) const=0;
 
protected:

  FileType * ptFileType;

  Integer lastlevel_;

private:
  ///
};

} // end of namespace
#endif // FILE_GRID
