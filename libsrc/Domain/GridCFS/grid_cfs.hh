#ifndef FILE_SCFE_GRID_CFS_2001
#define FILE_SCFE_GRID_CFS_2001

#include "DataInOut/filetype.hh"

#ifdef ADAPTGRID
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "Tetrahedron.h"
#include "Octahedron.h"
#include "MultilevelGrid.h"
#include "MeshReader.h"
#include "TetrahedronMeasure.h"
#include "MeshWriter.h"

#include "DataInOut/WriteInfo.hh"
#endif

namespace CoupledField
{

struct Elem;

/// Class for working with grid
template<Integer dim> 
class GridCFS
{
public:
  //! Constructor with parameter - pointer to FileType for reading initial grid
  GridCFS(FileType * const aptFileType);

  //! Deconstructor
  ~GridCFS();

  //! Read Grid Information
  void Read();
  
  //! Get connection of element
  /*!
    \param connect (output) contains global node numbers
    \param iElem (input) element level
    \param level (input) index for multilevel hierarchy
  */
  void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level);

  //! Get coordinates of node with number inode
  /*!
    \param inode (input) node number
    \param numlevel (input) index for multilevel hierarchy
    \param rfPoint (output) coordinates of point
  */
  void GetCoordinateNode(const Integer inode, const Integer numlevel, Point<dim> & rfPoint);

  //! Return maximum number of nodes
  /*!
    \param numlevel (input) index for multilevel hierarchy
  */
  Integer GetMaxnumnodes(const Integer numlevel)
        { return maxnumnodes_;}

  //! Return maximum number of elements 
  /*!
    \param numlevel (input) index for multilevel hierarchy
  */
  Integer GetMaxnumElem(const Integer numlevel);

  //! Return maximum number of elements in subdomains
  /*!
    \param numlevel (input) index for multilevel hierarchy
    \param subdoms (input) contains the names of the subdomains
  */
  Integer GetMaxnumElem( const Integer numlevel,
			 const std::vector<std::string> &subdoms );

  //! return dimension of grid
  Integer GetDim() { return dim_;}

  //!
  /*!
    \param  els (output)
    \param sd (input) contains the name of the subdomain
    \param level (input) index for multilevel hierarchy
  */
  void GetElemSD(std::vector<Elem*> & els, const std::string sd, const Integer level);

  //! returns all names of the subdomains contained in teh grid
  std::vector<std::string>* GetAllSDs(){ return &sd_;}

  void GetCoordNodesElem(const Vector<Integer> connect, Point<dim> * ptCoord, const Integer level); 


  //! gets the coordinates of the element nodes
  /*!
    \param connect (input) global node numbers of element
    \param ptCoord (output) coordinates of the element nodes (nrNodes \f$\times$ spaceDim);
    \param level (input) index for multilevel hierarchy
  */
  void GetCoordNodesElemMat( const Vector<Integer> connect,
			     Matrix<Double>& coordMat, const Integer level );
  


 //! return pointer to vector of element-neighbors for the element with number noOfElem
  /*!
    \param noOfElem (input) element level
    \param color (input) subdomain
  */ 
  std::vector<Elem*> * GetptNeighboursOfElem(const Integer noOfElem,
					     const std::string color);

  //! return vector of element-neighbors for the node with number noOfNode
  /*!
    \param  noOfNode (input) node number 
    \param neighbours (output) contains all elements, which contain the specified node
  */
  void GetNeighboursOfNode(const Integer noOfNode, std::vector<Elem*> * neighbours);

  //! procedure for forming list with element-neighbors for nodes of patch of element
  /*!
    \param elems (input)
    \param nodeNeighbors (output)
    \param map (input)
  */
  void FormNeighbors4NodesOfElements(const std::vector<Elem*> &elems, std::vector<std::vector<Elem*> > &nodeNeighbors, std::vector<Integer> & map);


  //! auxiliary function; to define belonging of one element to another from the list, for ex. surface element and boundary elements (min. 2 common nodes)
  /*!
    \param elemsSurf
    \param elems
    \param belongingSE
  */
  void DefineBelonging4Elems(const std::vector<Elem*>& elemsSurf, const std::vector<Elem*>&elems, std::vector<Elem*> & belongingSE);



  //! procedure for forming list with interface-elements neighbours
  /*!
    \param Interface (input) Nodes defining the interface between two domains
    \param Next2Surf (input) Subdomain adjacent to interface
    \param neighbours (output) Elements neighbouring (= have min. 1 node in common) to interface
  */
  void GetInterfaceNeighbours(std::vector<Integer> & interfaceNodes, 
			      std::vector<std::string> & subdoms, 
			      std::vector<Elem*> & neighbours,
			      Integer level);

  //! in this function we calculate area of element
  /*!
    \param elem (input) element object
  */
  Double CalcAreaElem(const Elem* elem);

  //! calculate number of nodes in patch of elements
  /*!
    \param patch (input)
    \param map (output)
  */
  void CalcNumberOfNodesInPatch(const std::vector<Elem*> & patch, std::vector<Integer>& map);

  //!
  void SetDim(const Integer dimension){dim_=dimension;}
  
#ifdef ADAPTGRID
  void putNodesFromGrid_RG(grd::MultilevelGrid * grid, const Integer level);

  void putElemsFromGrid_RG(grd::MultilevelGrid * grid, const Integer level);

  void Refine(grd::MultilevelGrid& grid);

  void ReRefine(grd::MultilevelGrid& grid);

  void RefineUniform(grd::MultilevelGrid& grid);

#endif

protected:
private:
  //<! file with initial mesh
  FileType * InFile;

  //<! list of elements for each subdomains
  std::vector<Elem*> * elems_;  

  //<! list of subdomains
  std::vector<std::string> sd_;

  //<! maximum number of nodes
  Integer maxnumnodes_;

  //<! pointer to array with coordinates of all nodes
  Point<dim> * ptCoordinate_;

  //<! dimension of grid
  Integer dim_;

  //! procedure for forming list with neighbors
  void FormNeighborsLists();

  //! list with neighbors for element
  std::vector<std::vector<Elem*> >** elNeighbors_;

  //! list with neighbors for nodes
  std::vector<std::vector<Elem*> > vtNeighbors_;

   //! class for mapping Grid_RG and GridCFS
  struct ElementMap {
    int sd;
    std::vector<Integer> map;
  };
  
  std::vector<ElementMap*> elemMap_; //!< mapping between GridRG and GridCFS

  //! only for test
  void SetRefinementFlag();

};

#ifdef __GNUC__
template class GridCFS<3>;
template class GridCFS<2>;
#endif

} // end of namespace
#endif // FILE_GRID_CFS
