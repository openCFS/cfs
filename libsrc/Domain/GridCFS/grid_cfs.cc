#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "grid_cfs.hh"
#include "elements_header.hh"

namespace CoupledField
{

template<class Dim>
 GridCFS<Dim> :: GridCFS(FileType * const aptFileType)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GridCFS" << std::endl;
#endif

  InFile = aptFileType;
}

template<>
void GridCFS<Point2D> :: Read()
{
#ifdef TRACE
  (*trace) << "entering GridCFS::Read" << std::endl;
#endif
  
  InFile->ReadMaxnumnodes(maxnumnodes_);
  ptCoordinate_=new Point2D[maxnumnodes_];
  InFile->ReadCoordinate(ptCoordinate_, maxnumnodes_);

  InFile->ReadElems(allelems);

#ifdef TRACE
  (*trace) << "leaving GridCFS::Read" << std::endl;
#endif
}

template<>
void GridCFS<Point3D> :: Read()
{
#ifdef TRACE
  (*trace) << "entering GridCFS::Read" << std::endl;
#endif

  InFile->ReadMaxnumnodes(maxnumnodes_);
  ptCoordinate_=new Point3D[maxnumnodes_];
  InFile->ReadCoordinate(ptCoordinate_, maxnumnodes_);

  InFile->ReadElems(allelems);

#ifdef TRACE
  (*trace) << "leaving GridCFS::Read" << std::endl;
#endif
}

template<class Dim>
void GridCFS<Dim> :: GetCoordOfNodesElem(const Integer i, const Integer l, const Integer numnodes,  Dim * ptCoordElem)   
{
#ifdef TRACE
  (*trace) << "entering GridCFS:: GetCoordOfNodesElem" << std::endl;
#endif

  Integer k;
  for (k=0; k < numnodes; k++)
       ptCoordElem[k]=ptCoordinate_[allelems[i].connect[k]-1];
}

template< class Dim>
void GridCFS<Dim>::GetConnection(Vector<Integer> & connection, const Integer iElem, const Integer level)
{
  connection=allelems[iElem].connect;
}

template<class Dim>
void GridCFS<Dim> :: GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint) 
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GetCoordinateNode" << std::endl;
#endif
  rfPoint=ptCoordinate_[inode];
}

} // end namespace
