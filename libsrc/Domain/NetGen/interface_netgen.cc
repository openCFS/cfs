#include <iostream>
#include <fstream>
#include <string>
#include <vector.h>

#include "interface_netgen.hh"
#include "elements_header.hh"
#include  "nginterface.h"

#include <mystdlib.h>

#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>

#include <visual.hpp>

#include "nginterface.h"
#include <FlexLexer.h>

namespace CoupledField
{

//! struct for Element
struct Elem
{
  BaseElem * ptElem;
  Vector<Integer> connect;
  std::string namesd;
};

template<>
void InterfaceNetGen<Point2D>::Read() 
{
#ifdef TRACE
  (*trace)<< "Entering InterfaceNetGen::Read 2D" << std::endl;
#endif

  Integer nnodes;
  ptFileType->ReadMaxnumnodes(nnodes); 

  Point2D * ptCoord=new Point2D[nnodes]; 
  ptFileType->ReadCoordinate(ptCoord, nnodes);
  
  Integer inode;
  for (inode=0; inode<nnodes; inode++)
    {
      Point3d p;
      p.X()=ptCoord[inode].x;
      p.Y()=ptCoord[inode].y;
      p.Z()=0.0;
      mesh.AddPoint(p);
    }

  if (ptCoord) delete [] ptCoord;

  // put information about elements
  std::vector<Elem> allel;
  ptFileType->ReadElems(allel);

  Element2d el3(3);
  Element2d el4(4);

 allptElem.resize(allel.size());
 
 Integer i, elemsize;
 for (i=0; i < allel.size(); i++)
{ 

 elemsize=allel[i].connect.size();
 allptElem[i]=allel[i].ptElem;

 std::cout << " elemsize " << elemsize << std::endl;
 switch(elemsize)
{
 case 3:
      el3.SetIndex(1);

      el3.PNum(1)=allel[i].connect[0];
      el3.PNum(2)=allel[i].connect[1];
      el3.PNum(3)=allel[i].connect[2];

      mesh.AddSurfaceElement(el3);

      break;

  case 4:

      el4.SetIndex(1);

      el4.PNum(1)=allel[i].connect[0];
      el4.PNum(2)=allel[i].connect[1];
      el4.PNum(3)=allel[i].connect[2];
      el4.PNum(4)=allel[i].connect[3];

      mesh.AddSurfaceElement(el4);

      break;

    default:
	  Error("Unknown type of element");
	  break;
      
}
}

   mesh.ClearFaceDescriptors();
   mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));  
#ifdef TRACE
 (*trace) << "Leaving InterfaceNetGen<Dim>::Read " << std::endl;
#endif
}

template<>
void InterfaceNetGen<Point3D>::Read()
{
#ifdef TRACE
  (*trace)<< "Entering InterfaceNetGen::Read" << std::endl;
#endif

  Integer nnodes;
  ptFileType->ReadMaxnumnodes(nnodes);

  Point3D * ptCoord=new Point3D[nnodes];
  ptFileType->ReadCoordinate(ptCoord, nnodes);

  Integer inode;
  for (inode=0; inode<nnodes; inode++)
    {
      Point3d p;
      p.X()=ptCoord[inode].x;
      p.Y()=ptCoord[inode].y;
      p.Z()=ptCoord[inode].z;
      mesh.AddPoint(p);
    }

  if (ptCoord) delete [] ptCoord;

  // put information about elements
  std::vector<Elem> allel;
  ptFileType->ReadElems(allel);

  Element el(4);

  Integer i, elemsize;

for (i=0; i< allel.size(); i++)
{

elemsize=allel[i].connect.size();

switch (elemsize)
{
 case 4:
      el.SetIndex(1);

      el.PNum(1)=allel[i].connect[0];
      el.PNum(2)=allel[i].connect[1];
      el.PNum(3)=allel[i].connect[2];
      el.PNum(4)=allel[i].connect[3];

      mesh.AddVolumeElement(el);
      break;

 default:
       Error("Unknown type of element");
       break;

}
} // end of for


   mesh.ClearFaceDescriptors();
   mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
#ifdef TRACE
 (*trace) << "Leaving InterfaceNetGen<Dim>::Read " << std::endl;
#endif
}

template<class Dim>
void InterfaceNetGen<Dim>::SubdivideUniform(const Integer level)
{
#ifdef TRACE
 (*trace) << "entering InterfaceNetGen<Dim>::SubdivideUniform " << std::endl;
#endif

lastlevel_++;

Integer ei;
Integer flag=1;

Integer maxnumelem=GetMaxnumElem(level);
// mesh.SurfaceElement(1).SetRefinementFlag (1);

for(ei=1; ei<=maxnumelem; ei++)
  SetRefinementFlag(ei,flag);

Refine();
}

template<>
void InterfaceNetGen<Point3D>::SetRefinementFlag(const Integer ei, const Integer flag)
{
  mesh.VolumeElement(ei).SetRefinementFlag (flag != 0);
}

template<>
void InterfaceNetGen<Point2D>::SetRefinementFlag(const Integer ei, const Integer flag)
{
  mesh.SurfaceElement(ei).SetRefinementFlag (flag != 0);
}

template<class Dim>
void InterfaceNetGen<Dim>::Refine()
{
  BisectionOptions biopt;
  biopt.usemarkedelements = 1;

  std::cout << " 1 " << std::endl;

  Refinement ref;
  ref.Bisect (mesh, biopt);

  std::cout << " 2 " << std::endl;

  mesh.UpdateTopology();
  mesh.UpdateClusters();

 std::cout << " 3 " << std::endl;
}

void InterfaceNetGen<Point2D>::GetCoordOfNodesElem(const Integer iElem, const Integer numlevel, const Integer numnodes, Point2D * ptCoordElem)
{ 
  if (!ptCoordElem) Error("Allocate ptCoordElem before using in function InterfaceNetGen::GetCoordOfNodesElem");

  Element2d el=mesh.SurfaceElement(iElem+1);
  Integer result;
  Point3d point;
 
  Integer i;
  for (i=0; i<numnodes; i++)
  {
   result=el.PNum(i+1);
   point=mesh.Point(result);

   ptCoordElem[i].x=point.X();
   ptCoordElem[i].y=point.Y();
  }
}

void InterfaceNetGen<Point3D>::GetCoordOfNodesElem(const Integer iElem, const Integer numlevel, const Integer numnodes, Point3D * ptCoordElem)
{

 Error("Not implemented");
/*
  if (!ptCoordElem) Error("Allocate ptCoordElem before using in function InterfaceNetGen::GetCoordOfNodesElem");

  Element el=mesh.VolumeElement(iElem+1);
  Integer result;
  Point3d point;

  Integer i;
  for (i=0; i<numnodes; i++)
  {
   result=el.PNum(i+1);
   point=mesh.Point(result);

   ptCoordElem[i].x=point.X();
   ptCoordElem[i].y=point.Y();
   ptCoordElem[i].z=point.Z();
  }
*/
}

template<class Dim>
BaseElem * InterfaceNetGen<Dim>::GetptElem(const Integer iElem)
{
 return allptElem[iElem]; 
}

template<>
void InterfaceNetGen<Point2D>::GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)
{
  Element2d el=mesh.SurfaceElement(iElem+1);
  Integer elemsize=4;              //// !!!! Attention
  connect.Resize(elemsize);

  Integer i;
  for (i=0; i<elemsize; i++)
  {
   connect[i]=el.PNum(i+1);
  } 
}

template<>
void InterfaceNetGen<Point3D>::GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)
{
  Error("Not implemented",__FILE__,__LINE__);
  Element2d el=mesh.SurfaceElement(iElem+1);
  Integer elemsize=4;              //// !!!! Attention
  connect.Resize(elemsize);

  Integer i;
  for (i=0; i<elemsize; i++)
  {
   connect[i]=el.PNum(i+1);
  }
}

template<>
void InterfaceNetGen<Point2D>::GetCoordinateNode(const Integer inode, const Integer numlevel, Point2D & rfPoint) 
{   
    Point3d auxPoint=mesh.Point(inode+1);
    rfPoint.x=auxPoint.X();
    rfPoint.y=auxPoint.Y();
}

template<>
void InterfaceNetGen<Point3D>::GetCoordinateNode(const Integer inode, const Integer numlevel, Point3D & rfPoint)
{
    Point3d auxPoint=mesh.Point(inode+1);
    rfPoint.x=auxPoint.X();
    rfPoint.y=auxPoint.Y();
    rfPoint.z=auxPoint.Z();
}

template<class Dim>
Integer InterfaceNetGen<Dim>::GetMaxnumnodes(const Integer numlevel)
  { return mesh.GetNP(); }

template<>
Integer InterfaceNetGen<Point2D>::GetMaxnumElem(const Integer numlevel)
{
 return mesh.GetNSE();
}

template<>
Integer InterfaceNetGen<Point3D>::GetMaxnumElem(const Integer numlevel)
{
 return mesh.GetNE();
}

template<class Dim>
Integer InterfaceNetGen<Dim>::GetNumNodesPerElem(const Integer iElem, const Integer level)
{ 
  return mesh.SurfaceElement(iElem+1).GetNP();
}

template<>
Integer InterfaceNetGen<Point3D>::GetNumNodesPerElem(const Integer iElem, const Integer level)
{
  return mesh.VolumeElement(iElem+1).GetNP();
}

template<class Dim>
InterfaceNetGen<Dim>::~InterfaceNetGen()
{
;
}

} // end of namespace

