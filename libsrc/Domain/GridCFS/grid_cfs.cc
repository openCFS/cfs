#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "grid_cfs.hh"
#include "grid.hh"
#include "elements_header.hh"
#include "baseelem.hh"
#include "conffile.hh"

namespace CoupledField
{

template<class Dim>
 GridCFS<Dim> :: GridCFS(FileType * const aptFileType)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GridCFS" << std::endl;
#endif

  InFile = aptFileType;

  conf->getsubdom(sd_);
  elems_=new std::vector<Elem>[sd_.size()];  
}

template<>
void GridCFS<Point2D> :: Read()
{
#ifdef TRACE
  (*trace) << "entering GridCFS::Read 2D" << std::endl;
#endif
 
  dim_=InFile->ReadDim();
 
  InFile->ReadMaxnumnodes(maxnumnodes_);
  ptCoordinate_=new Point2D[maxnumnodes_];
  InFile->ReadCoordinate(ptCoordinate_, maxnumnodes_);

//  InFile->ReadElems(allelems);
  InFile->ReadEl(elems_,sd_);

#ifdef TRACE
  (*trace) << "leaving GridCFS::Read" << std::endl;
#endif
}

template<>
void GridCFS<Point3D> :: Read()
{
#ifdef TRACE
  (*trace) << "entering GridCFS::Read 3D" << std::endl;
#endif

  dim_=InFile->ReadDim();

  InFile->ReadMaxnumnodes(maxnumnodes_);
  ptCoordinate_=new Point3D[maxnumnodes_];
  InFile->ReadCoordinate(ptCoordinate_, maxnumnodes_);

//  InFile->ReadElems(allelems);
  InFile->ReadEl(elems_,sd_);

#ifdef TRACE
  (*trace) << "leaving GridCFS::Read" << std::endl;
#endif
}

/*
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
*/

template< class Dim>
void GridCFS<Dim>::GetConnection(Vector<Integer> & connection, const Integer iElem, const Integer level)
{
  Integer i;
  Integer sum=0,sum1=0;
  for (i=0; i<sd_.size(); i++)
  {
    sum1=sum;
    sum+=elems_[i].size();
    if (iElem < sum)
       { connection=elems_[i][iElem-sum1].connect;
         break;
       }
  }
}

template<class Dim>
void GridCFS<Dim> :: GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint) 
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GetCoordinateNode" << std::endl;
#endif
  rfPoint=ptCoordinate_[inode];
}

template<class Dim>
void GridCFS<Dim> :: GetElemSD(std::vector<Elem> & els, const std::string sd, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GetCoordinateNode" << std::endl;
#endif

 Integer i;
 for (i=0; i<sd_.size(); i++)
  if (sd_[i]==sd) els=elems_[i];

}

template<class Dim>
void  GridCFS<Dim> :: GetCoordNodesElem(const Vector<Integer> connect, Dim * ptCoord)
{
#ifdef TRACE
  (*trace) << "entering GridCFS :: GetCoordinateNode" << std::endl;
#endif

  Integer k;
  for (k=0; k < connect.size(); k++)
       ptCoord[k]=ptCoordinate_[connect[k]-1];
 
}

template<class Dim>
GridCFS<Dim> ::~GridCFS() 
{ if (ptCoordinate_) delete [] ptCoordinate_;
  if (elems_) delete [] elems_;
}

template<class Dim>
Integer GridCFS<Dim>::GetMaxnumElem(const Integer numlevel)
{
 Integer i;
 Integer res=0;
 for (i=0; i<sd_.size(); i++)
  res+=elems_[i].size();

 return res;
}

} // end namespace
