#ifndef FILE_SCFE_GRID_2001
#define FILE_SCFE_GRID_2001

#include "filetype.hh"

namespace CoupledField
{
 
/// Class for working with grid
template<class Dim> 
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
  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevelGrid, const Integer numnodes, Dim * ptCoordElem)=0;  

   /// Get connection of element
  virtual void GetConnection(Integer * result, const Integer levelGrid, 
           const Integer numElem, const Integer numnodesPerElem)=0;

   /// Get coordinates of node with global number inode
   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint)=0;

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel)=0;

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel)=0;

  /// Return num of nodes per element i
  virtual Integer GetNumNodesPerElem(const Integer iElem, const Integer level)=0;

  /// Print coordinates of grid in out
//  virtual void PrintCoordinate(const Integer level, ostream * out) const=0;

protected:
  FileType * ptFileType;
private:
  ///
};

#ifdef __GNUC__
template class Grid<Point3D>;
template class Grid<Point2D>;
#endif
} // end of namespace
#endif // FILE_GRID
