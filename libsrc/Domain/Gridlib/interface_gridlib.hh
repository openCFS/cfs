
#ifndef FILE_INTERFACE_GRIDLIB_2002
#define FILE_INTERFACE_GRIDLIB_2002

namespace CoupledField
{
 
class FileType;
class GoMesh;

/// Class for working with grid
template<class Dim> 
class InterfaceGridlib: public Grid<Dim>
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  InterfaceGridlib(FileType * const aptFileType);

  /// Deconstructor
  virtual ~InterfaceGridlib() { if (ptGoMesh ) delete ptGoMesh ;}
  
  /// Get coordinates of all nodes which belong to element
  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevel,  Dim * ptCoordElem) 
  { Error("Not implemented yet");}

   /// Get connection of element
  virtual void GetConnection(Integer * result, const Integer level, 
           const Integer numElem, const Integer numnodesPerElem)
  { Error("Not implemented yet");}

  /// Return pointer to coordinates
//  virtual Dim * GetptCoordinate(const Integer numlevel)
//  { ptGridCFS->GetptCoordinate(numlevel);}

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel)
  { if (ptGoMesh) ptGoMesh->getNumVertices(numlevel);
     else Error("Mesh is undefined"); }

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel)
  { if (ptGoMesh) ptGoMesh->getNumElements(numlevel);
     else Error("Mesh is undefined");}

  /// Return num of nodes per element i
  virtual Integer GetNumNodesPerElem(const Integer iElem, const Integer level)
  { if (ptGoMesh)
  { GoGeometryElement<float> * ptElem=ptGoMesh->getElement(iElem, level);
      ptElem->getNumVertices();}
    else  Error("Mesh is undefined");
  }

  /// Print coordinates of grid in out
  virtual void PrintCoordinate(const Integer level, ostream * out) const
  { Error("Not implemented yet");}

private:

  //! Put information about grid
  void read();

  //! 
  FileType * ptFileType;
      
  //! 
  GoMesh * ptGoMesh;
};

template<class Dim>
inline InterfaceGridlib<Dim>::InterfaceGridlib(FileType * aptFileType)
: Grid<Dim>(aptFileType)
{
#ifdef TRACE
 (*trace) <<
             "Entering InterfaceGridlib<Dim>::InterfaceCFS<Dim>" 
                << std::endl;
#endif
  ptFileType=aptFileType;
  
  read();
}

#ifdef __GNUC__
template class InterfaceGridlib<Point3D>;
template class InterfaceGridlib<Point2D>;
#endif

} // end of namespace
#endif // 
