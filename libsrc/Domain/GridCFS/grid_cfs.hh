#ifndef FILE_SCFE_GRID_CFS_2001
#define FILE_SCFE_GRID_CFS_2001

#include "filetype.hh"
//#include "baseelem.hh"

namespace CoupledField
{

struct Elem;

/// Class for working with grid
template<class Dim> 
class GridCFS
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  GridCFS(FileType * const aptFileType);

  /// Deconstructor
  ~GridCFS();

  //! Read Grid Information
  void Read();
  
   /// Get connection of element
   void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level);

   /// Get coordinates of node with number inode
   void GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint);

  /// Return maximum number of nodes
  Integer GetMaxnumnodes(const Integer numlevel)
        { return maxnumnodes_;}

  /// Return maximum number of elements 
  Integer GetMaxnumElem(const Integer numlevel);

  /// Return num of nodes per element i
//  Integer GetNumNodesPerElem(const Integer iElem, const Integer level)
//{ return allelems[iElem].connect.size();}

  //!
  Integer GetDim() { return dim_;}

  //!
  void GetElemSD(std::vector<Elem> & els, const std::string sd, const Integer level);

  //!
  void GetCoordNodesElem(const Vector<Integer> connect, Dim * ptCoord);  

protected:
private:
  //!
  FileType *InFile;

  // 
  std::vector<Elem> * elems_;  

  std::vector<std::string> sd_;

  //
  Integer maxnumnodes_;

  //
  Dim * ptCoordinate_;

  //
  Integer dim_;

};

template class GridCFS<Point3D>;
template class GridCFS<Point2D>;

} // end of namespace
#endif // FILE_GRID_CFS
