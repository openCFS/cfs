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
  /// Constructor with parameter - pointer to FileType for reading initial grid
  GridCFS(FileType * const aptFileType);

  /// Deconstructor
  ~GridCFS();

  //! Read Grid Information
  void Read();
  
   /// Get connection of element
   void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level);

   /// Get coordinates of node with number inode
   void GetCoordinateNode(const Integer inode, const Integer numlevel, Point<dim> & rfPoint);

  /// Return maximum number of nodes
  Integer GetMaxnumnodes(const Integer numlevel)
        { return maxnumnodes_;}

  //! Return maximum number of elements 
  Integer GetMaxnumElem(const Integer numlevel);

  //! Return maximum number of elements in subdomains
 virtual Integer GetMaxnumElem(const Integer numlevel, const std::vector<std::string> & subdoms);

  //! return dimension of grid
  Integer GetDim() { return dim_;}

  //!
  void GetElemSD(std::vector<Elem*> & els, const std::string sd, const Integer level);

  ///
  std::vector<std::string>* GetAllSDs(){ return &sd_;}

  //!
  void GetCoordNodesElem(const Vector<Integer> connect, Point<dim> * ptCoord, const Integer level); 

 //! return pointer to vector of element-neighbors for the element with number noOfElem
  std::vector<Elem*> * GetptNeighboursOfElem(const Integer noOfElem, const std::string color);

  //! return vector of element-neighbors for the node with number noOfNode
  void GetNeighboursOfNode(const Integer noOfNode, std::vector<Elem*> * neighbours);

  //! auxialary function; to define belonging of one element to another from the list, for ex. surface element and boundary elements
  void DefineBelonging4Elems(const std::vector<Elem*>& elemsSurf, const std::vector<Elem*>&elems, std::vector<Elem*> & belongingSE);

  //! in this function we calculate area of element
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

  //<! dimention of grid
  Integer dim_;

  //! procedure for forming list with neighbors
  void FormNeighborsLists();

  //! procedure for forming list with element-neighbors for nodes of patch of element
  void FormNeighbors4NodesOfElements(const std::vector<Elem*> &elems, std::vector<std::vector<Elem*> > &nodeNeighbors, std::vector<Integer> & map);

  //! calculate number of nodes in patch of elements
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
  
  //! needful for maping betw. elements of GridRG and GridCFS
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
