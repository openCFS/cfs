#ifndef FILE_SCFE_GRID_CFS_2001
#define FILE_SCFE_GRID_CFS_2001

#include "filetype.hh"

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
  /// Return pointer to coordinates
  Dim * GetptCoordinate(const Integer numlevel)
        { return gh[numlevel].ptCoordinate;}

  /// Return maximum number of nodes
  Integer GetMaxnumnodes(const Integer numlevel)
        { return gh[numlevel].maxnumnode;}

  /// Return maximum number of elements 
  Integer GetMaxnumElem(const Integer numlevel)
        { return gh[numlevel].maxnumelem;}

  /// Return num of nodes per element i
  Integer GetNumNodesPerElem(const Integer iElem, const Integer level);

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
};


template class GridCFS<Point3D>;
template class GridCFS<Point2D>;

} // end of namespace
#endif // FILE_GRID_CFS
