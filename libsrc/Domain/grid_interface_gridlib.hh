#ifndef FILE_GRID_INTERFACE_GRIDLIB_2002
#define FILE_GRID_INTERFACE_GRIDLIB_2002

namespace CoupledField
{
 
class FileType;

/// Class for working with grid
template<class Dim> 
class GridInterfaceGridlib: public Grid<Dim>
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  GridInterfaceGridlib(FileType * const aptFileType);

  /// Deconstructor
  virtual ~GridInterfaceGridlib() { if (ptGridCFS) delete ptGridCFS;}
  
  /// Get coordinates of all nodes which belong to element
  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevel,  Dim * ptCoordElem) 
  { ptGridCFS->GetCoordOfNodesElem(numElem, numlevel,ptCoordElem);}

   /// Get connection of element
  virtual void GetConnection(Integer * result, const Integer level, 
           const Integer numElem, const Integer numnodesPerElem)
  { ptGridCFS->GetConnection(result, level, numElem, numnodesPerElem);}

  /// Return pointer to coordinates
//  virtual Dim * GetptCoordinate(const Integer numlevel)
//  { ptGridCFS->GetptCoordinate(numlevel);}

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel)
  { ptGridCFS->GetMaxnumnodes(numlevel); }

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel)
  { ptGridCFS->GetMaxnumElem(numlevel);}

  /// Return num of nodes per element i
  virtual Integer GetNumNodesPerElem(const Integer iElem, const Integer level)
  { ptGridCFS->GetNumNodesPerElem(iElem, level);}

  /// Print coordinates of grid in out
  virtual void PrintCoordinate(const Integer level, ostream * out) const
  { ptGridCFS->PrintCoordinate(level, out);}

private:

  GridCFS<Dim> * ptGridCFS;
};

template<class Dim>
inline GridInterfaceGridlib<Dim>::GridInterfaceGridlib(FileType * aptFileType)
: Grid<Dim>(aptFileType)
{
#ifdef TRACE
 (*trace) <<
             "Entering GridInterfaceGridlib<Dim>::GridInterfaceCFS<Dim>" 
                << endl;
#endif

   ptGridCFS=new GridCFS<Dim>(ptFileType); 

}

template class GridInterfaceGridlib<Point3D>;
template class GridInterfaceGridlib<Point2D>;

} // end of namespace
#endif // FILE_GRID
