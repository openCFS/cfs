#ifndef FILE_SCFE_GRID_2001
#define FILE_SCFE_GRID_2001

namespace CoupledField
{
 
class FileType;

/// Class for working with grid
template<class Dim> 
class Grid
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  Grid(FileType * const aptFileType);

  /// Deconstructor
  ~Grid();
  
  /// 
//  void SetCoordinate(Integer level, Integer nodenumber, Double * coord);

  /// Print coordinates of grid in out
  void PrintCoordinate(const Integer level, ostream * out) const;

  /// Print information about each element in out
  void PrintInfoElem(const Integer level,const Integer i, ostream * out) const;

  /// Get coordinates of all nodes which belong to element
  void GetCoordOfNodesElem(const Integer numElem, const Integer numlevelGrid,
                           Dim * ptCoordElem);

   /// Get connection of element
  void GetConnection(Integer * result, const Integer levelGrid, 
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


template class Grid<Point3D>;
template class Grid<Point2D>;

} // end of namespace
#endif // FILE_GRID
