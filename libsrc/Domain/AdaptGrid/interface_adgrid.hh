#ifndef FILE_INTERFACE_AdaptGrid_2002
#define FILE_INTERFACE_AdaptGrid_2002

#include "filetype.hh"
#include "grid.hh"

#include "grid_cfs.hh"
#include "spaceerror.hh"

#include "bcs.hh"

// from Adapt Grid
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

class SetRefFlag;
class SetRefFlagTest;

/// Interface to library Grid by Roberto G.
template<class Dim>
class InterfaceAdaptGrid: public Grid
{
public:
  /// Constructor with parameter - pointer to FileType for reading initial grid
  InterfaceAdaptGrid(FileType * const aptFileType);

  /// Deconstructor
  virtual ~InterfaceAdaptGrid();

   /// Uniform subdivision of domain
  virtual void SubdivideUniform(const Integer level);

   /// Get coordinates of node with global number inode
   virtual void GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint);

   /// Get coordinate of all nodes that belong to elem ie
   virtual void GetCoordNodesElem(const Vector<Integer> connect, Dim * ptCoord, const Integer level);

   /// Get connection of element
 virtual void GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level);

  /// Return maximum number of nodes
  virtual Integer GetMaxnumnodes(const Integer numlevel);

  /// Return maximum number of elements 
  virtual Integer GetMaxnumElem(const Integer numlevel);

  /// Put information about initial grid in mesh
  virtual void Read();

  // update nodes for boundary conditions
  virtual void UpdateBCs(std::list<Integer> * bcs);

  // prolongation of solution
  virtual void ProlongSol(const Vector<Double> sol_coarse, Vector<Double> &sol, const Integer alevel);

  ///
  virtual Integer GetDim() { return dim_; } // 
  
  ///
  virtual void GetElemSD(std::vector<Elem*> &, const std::string sd, const Integer level);

  //! Here we mark elements for refinement: ei - number of elem
  virtual void SetRefinementFlag(const Integer ei);
  virtual void SetRefinementFlag(const std::vector<Integer> ei);  

  //! Do refinement of elements, which we mark through function SetRefinementFlag
  virtual void Refine();
  virtual void TestRefine();
  virtual void TestCoarse();

  //!
  void forEachElemSd(SetRefFlag & f,const std::string subdomain);
  void forEachElemSd(SetRefFlagTest & f,const std::string subdomain);

//   virtual void forEachElemSd(PutElemMatInAlgSys & f,const std::string subdomain);
//   virtual void forEachElemSd(PutElemMatAlgSysElst3d & f,const std::string subdomain);

  //!
 void  Trans2CFSGrid(const Integer level=-1);

private:
  BCs * ptBCs;
  //! 
  FileType * ptFileType;
       
  //!
  GridCFS<Dim> * ptgridcfs_;

  //!
  Integer dim_;

#ifdef ADAPTGRID
  //!
  std::vector<grd::Vertex*> vertex_;
  std::vector<grd::Element*> elems_;

  //!
  grd::MultilevelGrid grid_;
  //!
  grd::ConformingClosure closure_;  
#endif

  //!
  void SetVertexNumbers();
};

#ifdef ADAPTGRID
class SetRefFlag
{
public:

  SetRefFlag(SpaceErrorEstimator * apt){ ptError_=apt;}
  ~SetRefFlag(){;}

 void operator() (grd::Element * t)
 {
  if (ptError_->TestLocError(t)) t->markForRefinement();
 }

private:
 
  SpaceErrorEstimator * ptError_;

};

class SetRefFlagTest
{
public:

  SetRefFlagTest(){;}
  ~SetRefFlagTest(){;}

 void operator() (grd::Element * t)
 {
   t->markForRefinement();
 }

};
#endif

#ifdef __GNUC__
template class InterfaceAdaptGrid<Point3D>;
template class InterfaceAdaptGrid<Point2D>;
#endif

} // end of namespace
#endif // 
