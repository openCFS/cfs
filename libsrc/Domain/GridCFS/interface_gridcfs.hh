#ifndef FILE_GRID_INTERFACE_CFS_2002
#define FILE_GRID_INTERFACE_CFS_2002

#include "Domain/grid.hh"
#include "grid_cfs.hh"

namespace CoupledField
{
 
class FileType;

/// Class for working with grid
template<Integer dim> 
class GridInterfaceCFS: public Grid
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  GridInterfaceCFS(FileType * const aptFileType);

  /// Deconstructor
  virtual ~GridInterfaceCFS() { if (ptGridCFS) delete ptGridCFS;}
  
  //! Read of mesh 
  virtual void Read() 
  { ptGridCFS->Read();}

   /// Get connection of element
   virtual void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)
    { ptGridCFS->GetConnection(connect,iElem, level);}

   /// Get coordinates of node with global number inode
   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point<dim> & rfPoint) {ptGridCFS->GetCoordinateNode(inode,numlevel,rfPoint);}

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel)
  { return ptGridCFS->GetMaxnumnodes(numlevel); }

  //! Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel)
  { return ptGridCFS->GetMaxnumElem(numlevel);}

  //! Return maximum number of elements in subdomains
  virtual Integer GetMaxnumElem(const Integer numlevel, const std::vector<std::string> & subdoms)
  { return ptGridCFS->GetMaxnumElem(numlevel,subdoms);}

  /// return dimension of mesh
  virtual Integer GetDim(){ return ptGridCFS->GetDim();}

  //!
  void GetElemSD(std::vector<Elem*> & els, const std::string sd, const Integer level)
  { ptGridCFS->GetElemSD(els,sd,level);}

  //!
  std::vector<std::string>* GetAllSDs(){ return ptGridCFS->GetAllSDs();}

  //!
  virtual void GetCoordNodesElem(const Vector<Integer> connect, Point<dim> * ptCoord, const Integer level)
  { ptGridCFS->GetCoordNodesElem(connect, ptCoord, level);}

  //! gets the coordinates of the element nodes
  virtual void GetCoordNodesElemMat(const Vector<Integer> connect, Matrix<Double>& coordMat, const Integer level)
  { ptGridCFS->GetCoordNodesElemMat(connect, coordMat, level);}  

     //! return vector of element-neighbors for the element with number noOfElem
  virtual  std::vector<Elem*> *  GetNeighboursOfElem(const Integer noOfElem, std::string color)
  { return ptGridCFS->GetptNeighboursOfElem(noOfElem,color);}

  //! return vector of element-neighbors for the node with number noOfNode
  virtual void GetNeighboursOfNode(const Integer noOfNode, std::vector<Elem*> * neighbours)
  { ptGridCFS->GetNeighboursOfNode(noOfNode,neighbours);}
 
  //! in this function we calculate area of element
  virtual Double CalcAreaElem(const Elem* elem)
  { 
    return ptGridCFS->CalcAreaElem(elem);
  }  

  //!
  virtual void DefineBelonging4Elems(const std::vector<Elem*>& elemsSurf, const std::vector<Elem*>&elems, std::vector<Elem*> & belongingSE)
  { ptGridCFS->DefineBelonging4Elems(elemsSurf,elems,belongingSE);}

protected:
private:
  GridCFS<dim> * ptGridCFS;
  ///
};

template<Integer dim>
inline GridInterfaceCFS<dim>::GridInterfaceCFS(FileType * aptFileType)
: Grid(aptFileType)
{
#ifdef TRACE
 (*trace) << "Entering GridInterfaceCFS<Dim>::GridInterfaceCFS<Dim>" << std::endl;
#endif
 lastlevel_=0;
   ptGridCFS=new GridCFS<dim>(ptFileType); 
}

#ifdef __GNUC__
template class GridInterfaceCFS<3>;
template class GridInterfaceCFS<2>;
#endif

} // end of namespace
#endif // FILE_GRID
