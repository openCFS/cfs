#ifndef FILE_SCFE_GRID_2001
#define FILE_SCFE_GRID_2001

#include "filetype.hh"
#include "elem.hh"

namespace CoupledField
{

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

   /// Get connection of element
   virtual void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)=0;

   /// Get coordinates of node with global number inode
   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point<2> & rfPoint)
  { Error(" Not implemented",__FILE__,__LINE__);}

  virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point<3> & rfPoint)
  { Error(" Not implemented",__FILE__,__LINE__);}

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel)=0;

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel)=0;
  //! Return maximum number of elements, which belong to subdoms
  virtual Integer GetMaxnumElem(const Integer numlevel, const std::vector<std::string> & subdoms)
   { Error(" Not implemented",__FILE__,__LINE__);}  

  //! Get last level of grid
  virtual Integer GetLastLevel() const { return lastlevel_;} 

  //! return dimension of mesh
  virtual Integer GetDim()=0;
  
 //! prolongation of solution
  virtual void ProlongSol(const Vector<Double> sol_coarse, Vector<Double> &sol, const Integer alevel)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! update nodes for boundary conditions
  virtual void UpdateBCs(std::list<Integer> * bcs)
  { Error(" Not implemented",__FILE__,__LINE__);}

    //! return vector of element-neighbors for the element with number noOfElem
  virtual  std::vector<Elem*> *GetNeighboursOfElem(const Integer noOfElem, std::string color)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! return vector of element-neighbors for the node with number noOfNode
  virtual void GetNeighboursOfNode(const Integer noOfNode, std::vector<Elem*> * neighbours)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! Do refinement of elements, which we mark through function SetRefinementFlag
  virtual void Refine()
  { Error(" This fnc is valid only for adaptgrid. Change your config-file",__FILE__,__LINE__);}

  //! Do uniform refinement of elements, which we mark through function SetRefinementFlag
  virtual void RefineUniform()
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! restore initial coarse mesh
  virtual void RestoreInitialMesh()
  { Error(" Not implemented",__FILE__,__LINE__);}

  //!
 virtual void GetElemSD(std::vector<Elem*> &, const std::string sd, const Integer level)
   { Error(" Not implemented",__FILE__,__LINE__);}

  //!
  virtual std::vector<std::string>* GetAllSDs()
  { Error("Not implemented",__FILE__,__LINE__);}

  //!
  virtual void GetCoordNodesElem(const Vector<Integer> connect, Point<2> * ptCoord, const Integer level)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //!
  virtual void GetCoordNodesElem(const Vector<Integer> connect, Point<3> * ptCoord, const Integer level)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! in this function we calculate area of element
  virtual Double CalcAreaElem(const Elem* elem)
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
