#include <iostream>
#include <fstream>
#include <string>
#include <vector.h>

#include "interface_netgen.hh"
//#include "filetype.hh"
#include "elements_header.hh"

namespace CoupledField
{

template<>
void InterfaceNetGen<Point2D>::Read() 
{
#ifdef TRACE
  (*trace)<< "Entering InterfaceNetGen::Read" << std::endl;
#endif

  Integer data[1];
  ptFileType-> ReadGeneralAnalChoice(data, FileType::numnode, FileType::endGAnal);
 
  Integer nnodes=data[0]; 

  ptFileType->ReadGeneralAnalChoice(data,FileType::numgroup,FileType::endGAnal);  maxnumsubdomain_=data[0];
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

  Integer data1[2];
  ptFileType->ReadGeneralElemChoice(0,data1, FileType::numelem, FileType::maxnode, FileType::endGElem);
  Integer nelems=data1[0];
  Integer nelemNodes=data1[1];

  Integer * Connect=new Integer[nelems*nelemNodes];

  ptArrayElem_=new BaseElem*[nelems+1];
  // ########################## number of groupes

  Integer i;

//  for (i=0; i<maxnumsubdomain_; i++)
  ptFileType->ReadElemConnectionGH(nelems, Connect, nelemNodes, 0,0);

  Element2d el3(3);
  Element2d el4(4);

  ptQ_=new Quad1(GaussOrder2);
  ptTr_=new Triangle1(GaussOrder3);

  switch (nelemNodes)
    {
      case 3:

  for (i=0; i<nelems; i++)
    {
      el3.SetIndex(i+1);

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
      el4.SetIndex(i+1);

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

template<class Dim>
void InterfaceNetGen<Dim>::SubdivideUniform(const Integer level)
{ Error("Not implemented yet",__FILE__,__LINE__); }

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

template<class Dim>
void InterfaceNetGen<Dim>::GetConnection(Integer * result, const Integer level,   const Integer numElem, const Integer numnodesPerElem)
{
 Element2d el=mesh.SurfaceElement(numElem+1);
 Integer i;
 for (i=0; i < el.GetNP(); i++)
  {
   result[i]=el.PNum(i+1);
  }
}

template<class Dim>
Integer InterfaceNetGen<Dim>::GetMaxnumnodes(const Integer numlevel)
  { return mesh.GetNP(); }

template<class Dim>
Integer InterfaceNetGen<Dim>::GetMaxnumElem(const Integer numlevel)
{
 return mesh.GetNSE();
}

template<class Dim>
Integer InterfaceNetGen<Dim>::GetNumNodesPerElem(const Integer iElem, const Integer level)
{ 
  return mesh.SurfaceElement(iElem+1).GetNP();
}

template<class Dim>
void InterfaceNetGen<Dim>::PrintCoordinate(const Integer level, std::ostream * out) const
  { Error("Not implemented yet",__FILE__,__LINE__); }

template<class Dim>
void InterfaceNetGen<Dim>::GetNodesBoundaryCondition(Vector<Integer> & nodesDirBC, const Integer level)
{
 if (level==0) ptFileType->ReadDirichletBC(nodesDirBC);
}

template<class Dim>
InterfaceNetGen<Dim>::~InterfaceNetGen()
{
  if (ptQ_) delete ptQ_;
  if (ptTr_) delete ptTr_;
}

template<class Dim>
void InterfaceNetGen<Dim>::GetNodesBoundaryCondition(Vector<Integer> & nodesDirBC, const Integer level)
{
 if (level==0) ptFileType->ReadDirichletBC(nodesDirBC);
}

} // end of namespace

