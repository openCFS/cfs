#ifndef FILE_SCFE_GRID_2001
#define FILE_SCFE_GRID_2001

#include <DataInOut/filetype.hh>
#include "elem.hh"

namespace CoupledField
{

//! Class for working with grid
class Grid
{
public:
  //! Constructor 
  /*!
    \param aptFileType pointer to FileType for reading initial grid
  */
  Grid(FileType * aptFileType); 

  //! Deconstructor
  virtual ~Grid();

  //! reads the grid from input file
  virtual void Read()=0;

   //! Get connection of element
  /*!
    \param connect (output) contains global node numbers
    \param iElem (input) element level
    \param level (input) index for multilevel hierarchy
  */
   virtual void GetConnection(StdVector<Integer> & connect, const Integer iElem, const Integer level)=0;

   //! Get coordinates of node with global number inode
  /*!
    \param inode (input) node number
    \param numlevel (input) index for multilevel hierarchy
    \param rfPoint (output) coordinates of point 2D
  */
   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point<2> & rfPoint)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! gets coordinate of specified node
  /*!
    \param inode (input) node number
    \param numlevel (input) index for multilevel hierarchy
    \param rfPoint (output) coordinates of point 3D
  */
  virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point<3> & rfPoint)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! Return maximum number of nodes
  /*!
    \param numlevel (input) index for multilevel hierarchy
  */
  virtual Integer GetMaxnumnodes(const Integer numlevel)=0;
  
  //! Return maximum number of elements 
  /*!
    \param numlevel (input) index for multilevel hierarchy
  */
  virtual Integer GetMaxnumElem(const Integer numlevel)=0;

  //! Return maximum number of elements, which belong to subdoms
  /*!
    \param numlevel (input) index for multilevel hierarchy
    \param subdoms (input) contains the names of the subdomains
  */
  virtual Integer GetMaxnumElem(const Integer numlevel, const StdVector<std::string> & subdoms)
  { 
    Error(" Not implemented",__FILE__,__LINE__);
    return Dint;
  }  

  //! Get last level of grid
  virtual Integer GetLastLevel() const { return lastlevel_;} 

  //! return dimension of mesh
  virtual Integer GetDim()=0;
  
 //! prolongation of solution
  /*!
    \param sol_coarse (input) solution on coarse grid
    \param sol (output) contains the solution on the new grid (fine grid)
    \param alevel (input) index in multilevel hierarchy
  */
  virtual void ProlongSol(const Vector<Double> sol_coarse, Vector<Double> &sol, const Integer alevel)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! update nodes for boundary conditions
  /*!
    \param bcs list of boundary nodes
  */
  virtual void UpdateBCs(std::list<Integer> * bcs)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! return vector of element-neighbors for the element with number noOfElem
  /*!
    \param noOfElem (input) element level
    \param color (input) subdomain
  */ 
  virtual  StdVector<Elem*> *GetNeighboursOfElem(const Integer noOfElem, std::string color)
  { 
    Error(" Not implemented",__FILE__,__LINE__);
    return Evec;
  }

  //! return vector of element-neighbors for the node with number noOfNode
  /*!
    \param noOfNode (input) global number of node
    \param neighbours (output) list with neighbors
    \param subdomains (input) list of subdomains, on which neighbors are searched
  */
  virtual void GetNeighboursOfNode(const Integer noOfNode,
				   StdVector<Elem*> * neighbours)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! Do refinement of elements, which we mark through function SetRefinementFlag
  virtual void Refine(const Integer numLoops = 1)
  { Error(" This fnc is valid only for adaptgrid. Change your config-file",__FILE__,__LINE__);}

  //! Do uniform refinement of elements, which we mark through function SetRefinementFlag
  virtual void RefineUniform()
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! restore initial coarse mesh
//  virtual void ResetToCoarseGrid()
//  { Error(" Not implemented",__FILE__,__LINE__);}

  //!
  /*!
    \param els (output)
    \param sd (input) contains the name of the subdomain
    \param level (input) index for multilevel hierarchy
  */
  virtual void GetElemSD(StdVector<Elem*> & els, const std::string sd, const Integer level)
   { Error(" Not implemented",__FILE__,__LINE__);}

  //!
  virtual StdVector<std::string>* GetAllSDs()
  { 
    Error("Not implemented",__FILE__,__LINE__);
    return Dstr;
  }

  //! gets the coordinates of the element nodes
  /*!
    \param connect (input) global node numbers of element
    \param ptCoord (output) coordinates of the element nodes
    \param level (input) index for multilevel hierarchy
  */
  virtual void GetCoordNodesElem(const StdVector<Integer> connect, 
				 Point<2> * ptCoord, 
				 const Integer level)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! gets a matrix of the coordinates of the element nodes
  /*!
    \param connect (input) global node numbers of element
    \param ptCoord (output) coordinates of the element nodes (spaceDim \f$\times\f$ nrNodes);
    \param level (input) index for multilevel hierarchy
  */
  virtual void GetCoordNodesElemMat(const StdVector<Integer> connect, 
				    Matrix<Double>& coordMat, 
				    const Integer level)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! gets the coordinates of the element nodes
  /*!
    \param connect (input) global node numbers of element
    \param ptCoord (output) coordinates of the element nodes
    \param level (input) index for multilevel hierarchy
  */
  virtual void GetCoordNodesElem(const StdVector<Integer> connect, Point<3> * ptCoord, const Integer level)
  { Error(" Not implemented",__FILE__,__LINE__);}

  //! in this function we calculate area of element
  /*!
    \param elem (input) element object
  */
  virtual Double CalcAreaElem(const Elem* elem)
  { 
    Error(" Not implemented",__FILE__,__LINE__);
    return Ddummy;
  }

  //! Calculates the surface normal pointing OUT OF the neighbouring
  //! volume element
  //! \param n (output) Normal vector
  //! \param surfElem (input) Surface element
  //! \param surfElem (input) Surface element
  void CalcSurfNormalOutOfVol(Vector<Double> & n,
			      const Elem & surfElem,
			      const Elem & volElem);
  
  //! form list with interface-elements neighbours
  //! NOTE: an element is considered as neighbour, if both have 
  //! AT LEAST one common node
  /*!
    \param interfaceNodes (input) Nodes defining the interface between two domains
    \param subdoms (input) Subdomain adjacent to interface
    \param neighbours (output) Elements neighbouring (= have min. 1 node in common) to interface
    \param level (input) Refinement level
  */
  virtual void GetInterfaceNeighbours(StdVector<Integer> & interfaceNodes, 
				      StdVector<std::string> & subdoms, 
				      StdVector<Elem*> & neighbours,
				      Integer level) = 0;
  

  //! Find volume elems next to surface elems

  //! Get to a list of surface elements the neighbouring volume elements
  //! lying in one of the given regions.
  /*!
    \param surfElems (input) Vector of surface elems
    \param neighRegions (input) Region names, where the volume elems must lie
    \param volElems (output) Vector of surface elems.
    \param level (input) Refinement level
  */
  //!\note If not all surface elems were assigned to EXACT ONE volume
  //! element, an error is thrown. If the search was successfull, the
  //! i-the entry in the surfElems-vector corresponds to the i-th
  //! entry in the volElems-vector
  virtual void GetVolNeighboursForSurf(const StdVector<Elem*> & surfElems,
				       const StdVector<std::string> & neighRegions,
				       StdVector<Elem*> & volElems,
				       const Integer level) = 0;
    
  
   //! calculate number of nodes in patch of elements
  /*!
    \param patch (input)
    \param map (output)
  */
  virtual void CalcNumberOfNodesInPatch(const StdVector<Elem*> & patch, 
					StdVector<Integer> & map) = 0;


  /// resets the integration type of all known elements
  void SetIntTypeAllElems(IntegrationType aIntType);
  
  StdVector<std::string>& GetListSubDomains()
  { return listSD_;}
  
  
protected:

  FileType * ptFileType;   //!< pointer to input file
  Integer lastlevel_;      //!< last level in multilevel hierarchy
  StdVector<std::string> listSD_; //!< list of names of subdomains

  // dummies just for sun compiler
  Integer Dint;
  Double Ddummy;
  StdVector<Elem*> *Evec;
  StdVector<std::string>* Dstr;

  //! Auxiliary function: procedure for forming list with element-neighbors for nodes of patch of element
  /*!
    \param elems (input)
    \param nodeNeighbors (output)
    \param map (input)
  */
  virtual void FormNeighbors4NodesOfElements(const StdVector<Elem*> &elems, 
					     StdVector<StdVector<Elem*> > &nodeNeighbors, 
					     StdVector<Integer> & map) = 0;
private:
  ///
};

} // end of namespace
#endif // FILE_GRID
