#include <iostream>
#include <fstream>
#include <string>
#include <vector.h>

#include "interface_netgen.hh"
//#include "filetype.hh"
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

template<>
void InterfaceNetGen<Point2D>::Read() 
{
#ifdef TRACE
  (*trace)<< "Entering InterfaceNetGen::Read" << std::endl;
#endif

  Integer data[1];
//  ptFileType-> ReadGeneralAnalChoice(data, FileType::numnode, FileType::endGAnal);
 
  Integer nnodes;
  ptFileType->ReadMaxnumnodes(nnodes); 

//  ptFileType->ReadGeneralAnalChoice(data,FileType::numgroup,FileType::endGAnal); 
   maxnumsubdomain_=1;
  pptelemsubdom_=new Integer*[maxnumsubdomain_];

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
  
  // ######################### number of groupes

//  Integer data1[2];
//  ptFileType->ReadGeneralElemChoice(0,data1, FileType::numelem, FileType::maxnode, FileType::endGElem);
  Integer nelems;
  ptFileType->ReadMaxnumelem(nelems);

  Integer nelemNodes;
  ptFileType->ReadNumberNodesPerElem(nelemNodes);

  Integer * Connect=new Integer[nelems*nelemNodes];

  ptArrayElem_=new BaseElem*[nelems+1];
  // ########################## number of groupes

  Integer i;

//  for (i=0; i<maxnumsubdomain_; i++)
  ptFileType->ReadElemConnectionGH(nelems, Connect, nelemNodes, 0,0);

  Element2d el3(3);
  Element2d el4(4);

  ptQ_=new Quad1();
  ptTr_=new Triangle1();

  switch (nelemNodes)
    {
      case 3:

  for (i=0; i<nelems; i++)
    {
      el3.SetIndex(1);

      el3.PNum(1)=Connect[i*nelemNodes];
      el3.PNum(2)=Connect[i*nelemNodes+1];
      el3.PNum(3)=Connect[i*nelemNodes+2];

      mesh.AddSurfaceElement(el3);

      ptArrayElem_[i]=ptTr_;
   }

      break;

      case 4:

  for (i=0; i<nelems; i++)
    {
      el4.SetIndex(1);

      el4.PNum(1)=Connect[i*nelemNodes];
      el4.PNum(2)=Connect[i*nelemNodes+1];
      el4.PNum(3)=Connect[i*nelemNodes+2];
      el4.PNum(4)=Connect[i*nelemNodes+3];

      mesh.AddSurfaceElement(el4);

      ptArrayElem_[i]=ptQ_;
   }

      break;

    default:
	  Error("Unknown type of element");
	  break;
      
   }

   ptArrayElem_[nelems]=NULL;

   delete [] Connect;

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

  Integer data[1];
//  ptFileType-> ReadGeneralAnalChoice(data, FileType::numnode, FileType::endGAnal);

  Integer nnodes;
  ptFileType->ReadMaxnumnodes(nnodes);

//  ptFileType->ReadGeneralAnalChoice(data,FileType::numgroup,FileType::endGAnal);
  maxnumsubdomain_= 1;
  pptelemsubdom_=new Integer*[maxnumsubdomain_];

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

  // ######################### number of groupes

  Integer nelems;
  ptFileType->ReadMaxnumelem(nelems);

  Integer nelemNodes;
  ptFileType->ReadNumberNodesPerElem(nelemNodes);

  Integer * Connect=new Integer[nelems*nelemNodes];

  ptArrayElem_=new BaseElem*[nelems+1];
  // ########################## number of groupes

  Integer i;

//  for (i=0; i<maxnumsubdomain_; i++)
  ptFileType->ReadElemConnectionGH(nelems, Connect, nelemNodes, 0,0);

  Element el (4);

  ptTet_=new Tetrahedral1();

  switch (nelemNodes)
    {
      case 4:

  for (i=0; i<nelems; i++)
    {
      Element el (4);
      el.SetIndex(1);

      el.PNum(1)=Connect[i*nelemNodes];
      el.PNum(2)=Connect[i*nelemNodes+1];
      el.PNum(3)=Connect[i*nelemNodes+2];
      el.PNum(4)=Connect[i*nelemNodes+3];

      std::cout << Connect[i*nelemNodes] << " " << Connect[i*nelemNodes+1] << " " << Connect[i*nelemNodes+2] << " " << Connect[i*nelemNodes+3] << std::endl;
 
      mesh.AddVolumeElement(el);

      ptArrayElem_[i]=ptTet_;
   }

      break;

    default:
          Error("Unknown type of element");
          break;

   }

   ptArrayElem_[nelems]=NULL;

   delete [] Connect;

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

template<class Dim>
void InterfaceNetGen<Dim>::GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodes, Dim * ptCoordElem)
{ 
  if (!ptCoordElem) Error("Allocate ptCoordElem before using in function InterfaceNetGen::GetCoordOfNodesElem");

  Element2d el=mesh.SurfaceElement(numElem+1);
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

template<>
void InterfaceNetGen<Point2D>::GetConnection(Integer * result, const Integer level,   const Integer numElem, const Integer numnodesPerElem)
{
 Element2d el=mesh.SurfaceElement(numElem+1);
 Integer i;
 for (i=0; i < el.GetNP(); i++)
  {
   result[i]=el.PNum(i+1);
  }
}

template<>
void InterfaceNetGen<Point3D>::GetConnection(Integer * result, const Integer level,   const Integer numElem, const Integer numnodesPerElem)
{
 Element el=mesh.VolumeElement(numElem+1);
 Integer i;
 for (i=0; i < el.GetNP(); i++)
  {
   result[i]=el.PNum(i+1);
   std::cout << " el num " << numElem << " " << result[i] << std::endl;
  }

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
void InterfaceNetGen<Dim>::PrintCoordinate(const Integer level, std::ostream * out) const
  { Error("Not implemented yet",__FILE__,__LINE__); }

template<class Dim>
InterfaceNetGen<Dim>::~InterfaceNetGen()
{
  if (ptQ_) delete ptQ_;
  if (ptTr_) delete ptTr_;
  if (ptTet_) delete ptTet_;
}

} // end of namespace

