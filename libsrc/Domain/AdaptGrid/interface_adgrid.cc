#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "interface_adgrid.hh"
#include "elements_header.hh"

namespace CoupledField
{

template<>
void InterfaceAdaptGrid<Point2D>::Read() 
{
#ifdef TRACE
  (*trace)<< "Entering InterfaceAdaptGrid::Read 2D" << std::endl;
#endif
  
  dim_=ptFileType->ReadDim();

  Integer nnodes;
  ptFileType->ReadMaxnumnodes(nnodes); 

  Point2D * ptCoord=new Point2D[nnodes]; 
  ptFileType->ReadCoordinate(ptCoord, nnodes);

// read vertices  

  Integer inode;
  vertex_.resize(nnodes);  
  for (inode=0; inode<nnodes; inode++)
    {
      Double pos[3];
      pos[0]=ptCoord[inode].x;
      pos[1]=ptCoord[inode].y;
      pos[2]=0.0;

      grd::Vertex* tmpVt = new grd::Vertex(pos);
      tmpVt->setId(inode+1);

      vertex_[inode]= tmpVt;
    }

  if (ptCoord) delete [] ptCoord;

// read elements
 
  ptFileType->ReadElems4AdaptGrid(elems_);
 
 grid_.buildCoarseMesh(vertex_,elems_);
   
#ifdef TRACE
 (*trace) << "Leaving InterfaceAdaptGrid<Dim>::Read " << std::endl;
#endif
}

template<>
void InterfaceAdaptGrid<Point3D>::Read()
{
#ifdef TRACE
  (*trace)<< "Entering InterfaceAdaptGrid::Read 3D" << std::endl;
#endif

  dim_=ptFileType->ReadDim();

  Integer nnodes;
  ptFileType->ReadMaxnumnodes(nnodes);

  Point3D * ptCoord=new Point3D[nnodes];
  ptFileType->ReadCoordinate(ptCoord, nnodes);

// read vertices

  Integer inode;
  vertex_.resize(nnodes);
  for (inode=0; inode<nnodes; inode++)
    {
      Double pos[3];
      pos[0]=ptCoord[inode].x;
      pos[1]=ptCoord[inode].y;
      pos[2]=ptCoord[inode].z;

      grd::Vertex* tmpVt = new grd::Vertex(pos);
      tmpVt->setId(inode+1);

      vertex_[inode]= tmpVt;
    }

  if (ptCoord) delete [] ptCoord;

// read elements

  ptFileType->ReadElems4AdaptGrid(elems_);

  grid_.buildCoarseMesh(vertex_,elems_);

#ifdef TRACE
 (*trace) << "Leaving InterfaceAdaptGrid<Dim>::Read " << std::endl;
#endif
}

template<class Dim>
void InterfaceAdaptGrid<Dim>::SubdivideUniform(const Integer level)
{
#ifdef TRACE
 (*trace) << "entering InterfaceAdaptGrid<Dim>::SubdivideUniform 3D" << std::endl;
#endif

lastlevel_++;
;

#ifdef TRACE
 (*trace) << "leaving InterfaceAdaptGrid<Dim>::SubdivideUniform 3D" << std::endl;
#endif
}

template<class Dim>
void InterfaceAdaptGrid<Dim>::Refine()
{
 ;
}

template<>
void InterfaceAdaptGrid<Point3D>::SetRefinementFlag(const Integer ei)
{
 ;
}

template<>
void InterfaceAdaptGrid<Point2D>::SetRefinementFlag(const Integer ei)
{
;
}

template<>
void InterfaceAdaptGrid<Point3D>::SetRefinementFlag(const std::vector<Integer> ei)
{
 ;
}

template<>
void InterfaceAdaptGrid<Point2D>::SetRefinementFlag(const std::vector<Integer> ei)
{
 ;
}

template<>
void InterfaceAdaptGrid<Point2D>::GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)
{
;
}

template<>
void InterfaceAdaptGrid<Point3D>::GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)
{
 ;
}

template<>
void InterfaceAdaptGrid<Point2D>::GetCoordinateNode(const Integer inode, const Integer numlevel, Point2D & rfPoint) 
{   
; 
}

template<>
void InterfaceAdaptGrid<Point3D>::GetCoordinateNode(const Integer inode, const Integer numlevel, Point3D & rfPoint)
{
;
}

template<class Dim>
Integer InterfaceAdaptGrid<Dim>::GetMaxnumnodes(const Integer numlevel)
  { return 0; }

template<>
Integer InterfaceAdaptGrid<Point2D>::GetMaxnumElem(const Integer numlevel)
{
 return 0;
}

template<>
Integer InterfaceAdaptGrid<Point3D>::GetMaxnumElem(const Integer numlevel)
{
 return 0;
}

template<class Dim>
InterfaceAdaptGrid<Dim>::~InterfaceAdaptGrid()
{
;
}

} // end of namespace

