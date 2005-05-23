#ifndef FILE_INTERFACE_AdaptGrid_2002
#define FILE_INTERFACE_AdaptGrid_2002

#include <DataInOut/filetype.hh>
#include <Domain/grid.hh>

#include <Domain/GridCFS/grid_cfs.hh>


// from Adapt Grid
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
#endif

namespace CoupledField
{

  class SetRefFlag;
  class SetRefFlagTest;

  /// Interface to library Grid by Roberto G.
  template<Integer dim>
  class InterfaceAdaptGrid: public Grid
  {
  public:
    //! Constructor with parameter - pointer to FileType for reading initial grid
    InterfaceAdaptGrid(FileType * const aptFileType);

    //! Deconstructor
    virtual ~InterfaceAdaptGrid();

    //! Get coordinates of node with global number inode
    virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Point<dim> & rfPoint);

    //! Get coordinate of all nodes that belong to elem ie
    virtual void GetCoordNodesElem(const Vector<Integer> connect, Point<dim> * ptCoord, const Integer level);

    //! 
    virtual void GetCoordNodesElemMat(const Vector<Integer> connect, Matrix<Double>& coordMat, const Integer level);

    //! Get connection of element
    virtual void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level);

    //! return vector of element-neighbors for the element with number noOfElem
    virtual std::vector<Elem*>* GetNeighboursOfElem(const Integer noOfElem, std::string color);

    //! Return maximum number of nodes
    virtual Integer GetMaxnumnodes(const Integer numlevel);
 
    //! Return maximum number of elements 
    virtual Integer GetMaxnumElem(const Integer numlevel);
    //! Return maximum number of elements in subdomains
    virtual Integer GetMaxnumElem(const Integer numlevel, const std::vector<std::string> & subdoms);

    //! Put information about initial grid in mesh
    virtual void Read();

    //! update nodes for boundary conditions
    virtual void UpdateBCs(std::list<Integer> * bcs);

    //! prolongation of solution
    virtual void ProlongSol(const Vector<Double> sol_coarse, Vector<Double> &sol, const Integer alevel);

    //! return dimension of grid
    virtual Integer GetDim() { return dim_; }  
  
    //! return vector of elements for this subdomain of grid
    virtual void GetElemSD(std::vector<Elem*> &, const std::string sd, const Integer level);

    //! return pointer to vector with all names of subdomains
    virtual std::vector<std::string> *GetAllSDs();

    //! restore initial coarse mesh
    //  virtual void ResetToCoarseGrid();

    //! in this function we calculate area of element
    virtual Double CalcAreaElem(const Elem* elem)
    { return ptgridcfs_->CalcAreaElem(elem);}

    //! Do refinement of elements, which we mark through function SetRefinementFlag
    virtual void Refine(const Integer numLoops = 1);
    virtual void RefineUniform();

    void FormNeighbors4NodesOfElements(const std::vector<Elem*> &elems, std::vector<std::vector<Elem*> > &nodeNeighbors, std::vector<Integer> & map)
    { ptgridcfs_->FormNeighbors4NodesOfElements(elems, nodeNeighbors,  map);}

    //!
    virtual void DefineBelonging4Elems(const std::vector<Elem*>& elemsSurf, const std::vector<Elem*>&elems, std::vector<Elem*> & belongingSE)
    { ptgridcfs_->DefineBelonging4Elems(elemsSurf,elems,belongingSE);}

    //!
    virtual void GetInterfaceNeighbours(std::vector<Integer> & Interface, 
                                        std::vector<std::string> & subdoms, 
                                        std::vector<Elem*> & Neighbours,
                                        Integer level)
    {  ptgridcfs_->GetInterfaceNeighbours(Interface, subdoms, Neighbours, level);}

    //!
    virtual void CalcNumberOfNodesInPatch(const std::vector<Elem*> & patch,
                                          std::vector<Integer> & map)
    { ptgridcfs_->CalcNumberOfNodesInPatch(patch,map);}

  private:

    //! 
    FileType * ptFileType;
       
    //!
    GridCFS<dim> * ptgridcfs_;

    //!
    Integer dim_;

#ifdef ADAPTGRID
    //!
    std::vector<grd::Vertex*> vertices_;
    std::vector<grd::Element*> elements_;

    //!
    grd::MultilevelGrid grid_;
    //!
    grd::ConformingClosure closure_;  
#endif

    //!
    void SetVertexNumbers();

    //! transformation from GridRG to CFS-Grid
    void Trans2CFSGrid(const Integer level=-1);

  };


  class SetRefFlag
  {
  public:

    SetRefFlag(){;}
    ~SetRefFlag(){;}

#ifdef ADAPTGRID
    void operator() (grd::Element * t)
    {
      t->markForRefinement();
    }
#endif

  private:
  };

  class SetRefFlagTest
  {
  public:

    SetRefFlagTest(){;}
    ~SetRefFlagTest(){;}

#ifdef ADAPTGRID
    void operator() (grd::Element * t)
    {
      t->markForRefinement();
    }
#endif
  };


#ifdef __GNUC__
  template class InterfaceAdaptGrid<2>;
  template class InterfaceAdaptGrid<3>;
#endif

} // end of namespace
#endif // 
