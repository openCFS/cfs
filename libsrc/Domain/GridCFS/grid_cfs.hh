#ifndef FILE_SCFE_GRID_CFS_2001
#define FILE_SCFE_GRID_CFS_2001

#include "filetype.hh"
#include "baseelem.hh"

namespace CoupledField
{
 
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
  
  /// Print coordinates of grid in out
  void PrintCoordinate(const Integer level, std::ostream * out) const;

  /// Print information about each element in out
  void PrintInfoElem(const Integer level,const Integer i, std::ostream * out) const;

  /// Get coordinates of all nodes which belong to element
  void GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodes,  Dim * ptCoordElem); 

   /// Get connection of element
  void GetConnection(Integer * result, const Integer level, 
           const Integer numElem, const Integer numnodesPerElem);

   /// Get coordinates of node with number inode
   void GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint);

  /// Return pointer to coordinates
//  Dim * GetptCoordinate(const Integer numlevel)
//        { return gh[numlevel].ptCoordinate;}

  /// Return maximum number of nodes
  Integer GetMaxnumnodes(const Integer numlevel)
        { return gh[numlevel].maxnumnode;}

  /// Return maximum number of elements 
  Integer GetMaxnumElem(const Integer numlevel)
        { return gh[numlevel].maxnumelem;}

  /// Return num of nodes per element i
  Integer GetNumNodesPerElem(const Integer iElem, const Integer level);

  /// Return pointer to array of elements
  BaseElem ** getptArrayElem() const { return ptArrayElem_; }

   //! Get number of subdomains
  Integer getnumsubdomains() const {return maxnumsubdomain;}

   //! Get pointer to array with nodes, that belongs to subdomain number num
  Integer * getelemsubdomain(const Integer num, const Integer level) const
  { return pptelemsubdom[num];}  

protected:
private:
  //!
  FileType *InFile;
  ///
  Integer maxnumsubdomain;
  ///
  Integer numlevel;
  //
  GridHierarchy<Dim> gh[20];
  ///
  Integer sizeConnectElem;
  ///
  BaseElem ** ptArrayElem_;
  //
  BaseElem * ptQ_, * ptTr_;
  //
  Integer ** pptelemsubdom;

};


template class GridCFS<Point3D>;
template class GridCFS<Point2D>;

} // end of namespace
#endif // FILE_GRID_CFS
