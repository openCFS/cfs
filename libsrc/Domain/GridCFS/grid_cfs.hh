#ifndef FILE_SCFE_GRID_CFS_2001
#define FILE_SCFE_GRID_CFS_2001

#include "filetype.hh"
#include "acoustic2dPDE.hh"


#ifdef ADAPTGRID
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "Tetrahedron.h"
#include "Octahedron.h"
#include "MultilevelGrid.h"
#include "GeometrySensor.h"
#include "MeshReader.h"
#include "TetrahedronMeasure.h"
#include "MeshWriter.h"
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
  virtual Integer GetMaxnumElem(const Integer numlevel, const std::vector<std::string> & subdoms);

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

  //! gets the coordinates of the element nodes
  /*!
    \param connect (input) global node numbers of element
    \param ptCoord (output) coordinates of the element nodes
    \param level (input) index for multilevel hierarchy
  */
  void GetCoordNodesElem(const Vector<Integer> connect, Point<dim> * ptCoord, const Integer level); 

 //! return pointer to vector of element-neighbors for the element with number noOfElem
  /*!
    \param noOfElem (input) element level
    \param color (input) subdomain
  */ 
  std::vector<Elem*> * GetptNeighboursOfElem(const Integer noOfElem, const std::string color);

  //! return vector of element-neighbors for the node with number noOfNode
  /*!
    \param  noOfNode (input) node number 
    \param neighbours (output) contains all elements, which contain the specified node
  */
  void GetNeighboursOfNode(const Integer noOfNode, std::vector<Elem*> * neighbours);

  //! auxialary function; to define belonging of one element to another from the list, for ex. surface element and boundary elements
  /*!
    \param elemsSurf
    \param elems
    \param belongingSE
  */
  void DefineBelonging4Elems(const std::vector<Elem*>& elemsSurf, const std::vector<Elem*>&elems, std::vector<Elem*> & belongingSE);

  //! in this function we calculate area of element
  /*!
    \param elem (input) element object
  */
  Double CalcAreaElem(const Elem* elem);

#ifdef ADAPTGRID
  void putNodesFromGrid_RG(grd::MultilevelGrid * grid, const Integer level);

 void putElemsFromGrid_RG(grd::MultilevelGrid * grid, const Integer level);

  void Refine(grd::MultilevelGrid& grid);

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

  //! procedure for forming list with element-neighbors for nodes of patch of element
  /*!
    \param elems (input)
    \param nodeNeighbors (output)
    \param map (input)
  */
  void FormNeighbors4NodesOfElements(const std::vector<Elem*> &elems, std::vector<std::vector<Elem*> > &nodeNeighbors, std::vector<Integer> & map);

  //! calculate number of nodes in patch of elements
  /*!
    \param patch (input)
    \param map (output)
  */
  void CalcNumberOfNodesInPatch(const std::vector<Elem*> & patch, std::vector<Integer>& map);

  //! list with neighbors for element
  std::vector<std::vector<Elem*> >** elNeighbors_;

  //! list with neighbors for nodes
  std::vector<std::vector<Elem*> > vtNeighbors_;

   //! 
  struct ElementMap {
    int sd;
    std::vector<Integer> map;
  };
  
  //! useful for maping betw. elements of GridRG and GridCFS
  std::vector<ElementMap*> elemMap;

  //! only for test
  void SetRefinementFlag();

};

#ifdef __GNUC__
template class GridCFS<3>;
template class GridCFS<2>;
#endif

} // end of namespace
#endif // FILE_GRID_CFS
