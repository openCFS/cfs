#ifndef FILE_SCFE_GRID_2001
#define FILE_SCFE_GRID_2001

#include "filetype.hh"
#include "baseelem.hh"

namespace CoupledField
{


struct Elem
{
  BaseElem * ptElem;
  Vector<Integer> connect;
  std::string namesd;
};


/// Class for working with grid
class Grid
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  Grid(FileType * aptFileType); 

  /// Deconstructor
  virtual ~Grid();

  //!
  virtual void Read()=0;

  /// 
  virtual void SubdivideUniform(const Integer level)=0;

  /// Get coordinates of all nodes which belong to element
//  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevelGrid, const Integer numnodes, Point2D * ptCoordElem)
//  { Error(" Not implemented",__FILE__,__LINE__);}

//  virtual void GetCoordOfNodesElem(const Integer numElem, const Integer numlevelGrid, const Integer numnodes, Point3D * ptCoordElem)
//  { Error(" Not implemented",__FILE__,__LINE__);}

   /// Get connection of element
   virtual void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)=0;

   /// Get coordinates of node with global number inode
   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point2D & rfPoint)
  { Error(" Not implemented",__FILE__,__LINE__);}

   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point3D & rfPoint)
  { Error(" Not implemented",__FILE__,__LINE__);}

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel)=0;

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel)=0;

  //! Get last level of grid
  virtual Integer GetLastLevel() const { return lastlevel_;} 

  //! return dimension of mesh
  virtual Integer GetDim()=0;

    //! Here we mark elements for refinement: ei - number of elem
  virtual void SetRefinementFlag(const Integer ei)
  { Error(" Not implemented",__FILE__,__LINE__);}

  virtual void SetRefinementFlag(Vector<Integer> & ei)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! Do refinement of elements, which we mark through function SetRefinementFlag
  virtual void Refine()
  { Error(" Not implemented",__FILE__,__LINE__);}

    //!
 virtual void GetElemSD(std::vector<Elem> &, const std::string sd, const Integer level)
   { Error(" Not implemented",__FILE__,__LINE__);}

  //!
  virtual void GetCoordNodesElem(const Vector<Integer> connect, Point2D * ptCoord)
  { Error(" Not implemented",__FILE__,__LINE__);}

  virtual void GetCoordNodesElem(const Vector<Integer> connect, Point3D * ptCoord)
  { Error(" Not implemented",__FILE__,__LINE__);}

protected:

  FileType * ptFileType;

  Integer lastlevel_;

  std::vector<std::string> listSD_;

private:
  ///
};

} // end of namespace
#endif // FILE_GRID
