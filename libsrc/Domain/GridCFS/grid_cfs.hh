#ifndef FILE_SCFE_GRID_CFS_2001
#define FILE_SCFE_GRID_CFS_2001

#include "filetype.hh"
#include "baseelem.hh"


namespace CoupledField
{

//! struct for Element
struct Elem
{
  BaseElem * ptElem;
  Vector<Integer> connect;
  std::string namesd;
};
 
/// Class for working with grid
template<class Dim> 
class GridCFS
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  GridCFS(FileType * const aptFileType);

  /// Deconstructor
  ~GridCFS(){ if (ptCoordinate_) delete [] ptCoordinate_;}

  //! Read Grid Information
  void Read();
  
  /// Get coordinates of all nodes which belong to element
  void GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodes,  Dim * ptCoordElem); 

   /// Get connection of element
   void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level);

   /// Get coordinates of node with number inode
   void GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint);

  /// Return maximum number of nodes
  Integer GetMaxnumnodes(const Integer numlevel)
        { return maxnumnodes_;}

  /// Return maximum number of elements 
  Integer GetMaxnumElem(const Integer numlevel)
        { return allelems.size();}

  /// Return num of nodes per element i
  Integer GetNumNodesPerElem(const Integer iElem, const Integer level)
{ return allelems[iElem].connect.size();}

  /// return pointer to pointer to BaseElem
  BaseElem * GetptElem(const Integer iElem)
 { return allelems[iElem].ptElem;}  

protected:
private:
  //!
  FileType *InFile;
  //
  std::vector<Elem> allelems;
  //
  Integer maxnumnodes_;
  //
  Dim * ptCoordinate_;

};

template class GridCFS<Point3D>;
template class GridCFS<Point2D>;

} // end of namespace
#endif // FILE_GRID_CFS
