#ifndef FILE_SCFE_GRID_CFS_2001
#define FILE_SCFE_GRID_CFS_2001

#include "filetype.hh"
#include "acoustic2dPDE.hh"
//#include "baseelem.hh"

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
template<class Dim> 
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
   void GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint);

  /// Return maximum number of nodes
  Integer GetMaxnumnodes(const Integer numlevel)
        { return maxnumnodes_;}

  /// Return maximum number of elements 
  Integer GetMaxnumElem(const Integer numlevel);

  //!
  Integer GetDim() { return dim_;}

  //!
  void GetElemSD(std::vector<Elem*> & els, const std::string sd, const Integer level);

  //!
  void GetCoordNodesElem(const Vector<Integer> connect, Dim * ptCoord, const Integer level); 


  std::vector<Integer*> * GetPtTestConnection(){ return &testconnect_;}

#ifdef ADAPTGRID
  void putNodesFromGrid_RG(grd::MultilevelGrid * grid, const Integer level);

 void putElemsFromGrid_RG(grd::MultilevelGrid * grid, const Integer level);
#endif

  //! Iterators
  // template<class T> void forEachElemSd(T & f,const std::string subdomain);
  //  void forEachElemSd(PutElemMatInAlgSys & f,const std::string subdomain);
  // void forEachElemSd(PutElemMatAlgSysElst3d & f,const std::string subdomain);

protected:
private:
  //!
  FileType *InFile;

  // 
  std::vector<Elem*> * elems_;  

  std::vector<std::string> sd_;

  std::vector<Integer*> testconnect_;

  //
  Integer maxnumnodes_;

  //
  Dim * ptCoordinate_;

  //
  Integer dim_;

};

template class GridCFS<Point3D>;
template class GridCFS<Point2D>;

} // end of namespace
#endif // FILE_GRID_CFS
