#include <iostream>
#include <fstream>
#include <string>
#include <vector.h>

#include "interface_netgen.hh"
#include "elements_header.hh"
#include "nginterface.h"

#include <mystdlib.h>

#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>

#include <visual.hpp>

#include "nginterface.h"
#include <FlexLexer.h>

namespace CoupledField
{

template<>
void InterfaceNetGen<Point2D>::Read() 
{
#ifdef TRACE
  (*trace)<< "Entering InterfaceNetGen::Read 2D" << std::endl;
#endif
  
  dim_=ptFileType->ReadDim();

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
//  ptFileType->ReadElems(allel);

  Element2d el3(3);
  Element2d el4(4);

 allptElem.resize(allel.size());
 
 Integer i, elemsize;
 for (i=0; i < allel.size(); i++)
{ 

 elemsize=allel[i].connect.size();
 allptElem[i]=allel[i].ptElem;

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
  mesh.AddFaceDescriptor (FaceDescriptor(1,0,1,1));  
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

  Integer dim_=ptFileType->ReadDim();

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
//  ptFileType->ReadElems(allel);

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
   mesh.AddFaceDescriptor (FaceDescriptor(1,0,1,1));
#ifdef TRACE
 (*trace) << "Leaving InterfaceNetGen<Dim>::Read " << std::endl;
#endif
}


template<class Dim>
void InterfaceNetGen<Dim>::SubdivideUniform(const Integer level)
{
#ifdef TRACE
 (*trace) << "entering InterfaceNetGen<Dim>::SubdivideUniform 3D" << std::endl;
#endif

lastlevel_++;
Refine();

#ifdef TRACE
 (*trace) << "leaving InterfaceNetGen<Dim>::SubdivideUniform 3D" << std::endl;
#endif
}

template<>
void InterfaceNetGen<Point3D>::SetRefinementFlag(const Integer ei)
{
  Integer flag=0;

  Integer i;
  Integer maxnumelem=mesh.GetNE();
  for (i=1; i<=maxnumelem; i++)
{
  if (i!=(ei+1))
    mesh.VolumeElement(i).SetRefinementFlag(flag);
}

}

template<>
void InterfaceNetGen<Point2D>::SetRefinementFlag(const Integer ei)
{
  Integer flag=0; 
 
  Integer i;
  Integer maxnumelem=mesh.GetNSE();
  for (i=1; i<=maxnumelem; i++)
{
  if (i!=(ei+1))
  mesh.SurfaceElement(i).SetRefinementFlag(flag);
}

}

template<>
void InterfaceNetGen<Point3D>::SetRefinementFlag(Vector<Integer> & ei)
{
 Integer noref=0;
 sort(ei.get(),ei.size());

 Integer maxnumelem=mesh.GetNE(); 
 Integer i,j=0;
 for (i=0; i<maxnumelem; i++)
  {
    if (i==ei[j]) {
                     if (j<(ei.size()-1)) j++;
                  }
    else
        mesh.VolumeElement(i+1).SetRefinementFlag(noref);
  }
}

template<>
void InterfaceNetGen<Point2D>::SetRefinementFlag(Vector<Integer> & ei)
{
 Integer noref=0;
 sort(ei.get(),ei.size());

 Integer maxnumelem=mesh.GetNSE();
 std::cout << maxnumelem << std::endl;

 Integer i,j=0;
 for (i=0; i<maxnumelem; i++)
  {
    if (i==ei[j]) {
                     if (j<(ei.size()-1)) j++;
                  }
    else
        mesh.SurfaceElement(i+1).SetRefinementFlag(noref);
  }
}

template<class Dim>
void InterfaceNetGen<Dim>::Refine()
{
  BisectionOptions biopt;
  biopt.usemarkedelements = 1;

  Refinement ref;
  ref.Bisect (mesh, biopt);

  mesh.UpdateTopology();
  mesh.UpdateClusters();
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
  Integer elemsize=el.GetNP();           
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
  Element el=mesh.VolumeElement(iElem+1);
  Integer elemsize=el.GetNP();         
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

template<class Dim>
void InterfaceNetGen<Dim>::Init()
{
 mycout=&cout;
 myerr=&cerr;
 testout=new ofstream("test.out");
}

} // end of namespace

