#include <iostream>
#include <fstream>
#include <string>
#include <vector.h>


#include "interface_gridlib.hh"
#include <GoTriangleMesh.hh>
#include <GoVolumeMesh.hh>
#include <GoDefaultVertex.hh>
#include <GoTriangleElement.hh>
#include <GoTetrahedronElement.hh> 
#include <GoHexahedronElement.hh>
#include <GoFlowVertex.hh>

namespace CoupledField
{

template<class Dim>
void InterfaceGridlib<Dim>::GetConnection(Integer * result, const Integer level,
           const Integer numElem, const Integer numnodesPerElem)
{
 Integer i;
 for (i=0; i<numnodesPerElem; i++)
 result[i]=ptGoMesh->getElement(numElem,level)->getVertex(i)->getId()+1;
}

template<>
void InterfaceGridlib<Point2D>::GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodesPerElem, Point2D * ptCoordElem)
{
 float x,y,z;
 Integer i;
 for (i=0; i<numnodesPerElem; i++)
 { 
 ptGoMesh->getElement(numElem, numlevel)->getVertex(i)->getPosition(x,y,z);
 ptCoordElem[i].x= double(x);
 ptCoordElem[i].y= double(y);
 }
}

template<>
void InterfaceGridlib<Point3D>::GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodesPerElem, Point3D * ptCoordElem)
{
 float x,y,z;
 Integer i;
 for (i=0; i<numnodesPerElem; i++)
 {
 ptGoMesh->getElement(numElem, numlevel)->getVertex(i)->getPosition(x,y,z);
 ptCoordElem[i].x= double(x);
 ptCoordElem[i].y= double(y);
 ptCoordElem[i].z= double(z);
 }
}

template<>
void InterfaceGridlib<Point2D>::GetCoordinateNode(const Integer inode, const Integer numlevel, Point2D & rfPoint)
{
 float x,y,z;
 ptGoMesh->getVertex(inode,numlevel)->getPosition(x,y,z);
 rfPoint.x= double(x);
 rfPoint.y= double(y);
}

template<>
void InterfaceGridlib<Point3D>::GetCoordinateNode(const Integer inode, const Integer numlevel, Point3D & rfPoint)
{
 float x,y,z;
 ptGoMesh->getVertex(inode,numlevel)->getPosition(x,y,z);
 rfPoint.x= double(x);
 rfPoint.y= double(y);
 rfPoint.z= double(z);
}

template<>
void InterfaceGridlib<Point2D>::Read() 
{
#ifdef TRACE
  (*trace)<< "Entering InterfaceGridlib::Read" << std::endl;
#endif

  Integer data[1];
  ptFileType-> ReadGeneralAnalChoice(data, FileType::numnode, FileType::endGAnal);
 
  Integer nnodes=data[0]; 

  ptGoMesh=new GoTriangleMesh;

  vector<GoDefaultVertex*> vv;


  GoDefaultVertex * v=NULL;
  vv.reserve(nnodes); 

  Point2D * ptCoord=new Point2D[nnodes]; 
  ptFileType->ReadCoordinate(ptCoord, nnodes);
   
  Integer inode;
  for (inode=0; inode<nnodes; inode++)
    {
      v=new GoDefaultVertex;
      v->setPosition(ptCoord[inode].x, ptCoord[inode].y, 0.0);
      v->setId(inode);
      ptGoMesh->addVertex(v);
      vv.push_back(v);
    }

  if (ptCoord) delete [] ptCoord;

  // put information about elements
  
  // ######################### number of groupes

  Integer data1[2];
  ptFileType->ReadGeneralElemChoice(0,data1, FileType::numelem, FileType::maxnode, FileType::endGElem);
  Integer nelems=data1[0];
  Integer nelemNodes=data1[1];

  Integer * Connect=new Integer[nelems*nelemNodes];
  // ########################## number of groupes
  ptFileType->ReadElemConnectionGH(nelems, Connect, nelemNodes, 0);

  GoGeometryElement<float > * t=NULL;
  Integer i;
  for (i=0; i<nelems; i++)
    {
      switch (nelemNodes)
	{
	case 3:
	  t=new GoTriangleElement;//######################;

	  t->setVertex(0,vv[Connect[i*nelemNodes]-1]);
	  vv[Connect[i*nelemNodes]-1]->setElement(t);
	  vv[Connect[i*nelemNodes]-1]->setValence(vv[Connect[i*nelemNodes]-1]->getValence()+1);
    
	  t->setVertex(1,vv[Connect[i*nelemNodes+1]-1]);
	  vv[Connect[i*nelemNodes+1]-1]->setElement(t);
	  vv[Connect[i*nelemNodes+1]-1]->setValence(vv[Connect[i*nelemNodes+1]-1]->getValence()+1);
 
	  t->setVertex(2,vv[Connect[i*nelemNodes+2]-1]);
	  vv[Connect[i*nelemNodes+2]-1]->setElement(t);
	  vv[Connect[i*nelemNodes+2]-1]->setValence(vv[Connect[i*nelemNodes+2]-1]->getValence()+1);

	  break;

	default:
	  Error("Unknown type of element");
	  t=NULL;
	  break;
	}
  
      if (t) {
	t->setId(i-1);
	ptGoMesh->addElement(t);
      }
    }
}

template<>
void InterfaceGridlib<Point3D>::Read()
{
#ifdef TRACE
  (*trace)<< "Entering InterfaceGridlib::Read - Point3D" << std::endl;
#endif

  Integer data[1];
  ptFileType-> ReadGeneralAnalChoice(data, FileType::numnode, FileType::endGAnal);

  Integer nnodes=data[0];

  ptGoMesh=new GoVolumeMesh;

  vector<GoFlowVertex*> vv;

  GoFlowVertex * v=NULL;
  vv.reserve(nnodes);

  Point3D * ptCoord=new Point3D[nnodes];
  ptFileType->ReadCoordinate(ptCoord, nnodes);

  Integer inode;
  for (inode=0; inode<nnodes; inode++)
    {
      v=new GoFlowVertex;
      v->setPosition(ptCoord[inode].x, ptCoord[inode].y, ptCoord[inode].z);
      v->setId(inode);
      ptGoMesh->addVertex(v);
      vv.push_back(v);
    }

  if (ptCoord) delete [] ptCoord;

  // put information about elements

  // ######################### number of groupes

  Integer data1[2];
  ptFileType->ReadGeneralElemChoice(0,data1, FileType::numelem, FileType::maxnode, FileType::endGElem);
  Integer nelems=data1[0];
  Integer nelemNodes=data1[1];
  std::cout << data1[0] << " num element nodes " << data1[1] << std::endl;

  Integer * Connect=new Integer[nelems*nelemNodes];
  // ########################## number of groupes
  ptFileType->ReadElemConnectionGH(nelems, Connect, nelemNodes, 0);

  GoGeometryElement<float > * t=NULL;
  Integer i;
  for (i=0; i<nelems; i++)
    {
      switch (nelemNodes)
        {
        case 4:
          t=new GoTetrahedronElement;//######################;

          t->setVertex(0,vv[Connect[i*nelemNodes]-1]);
          vv[Connect[i*nelemNodes]-1]->setElement(t);
          vv[Connect[i*nelemNodes]-1]->setValence(vv[Connect[i*nelemNodes]-1]->getValence()+1);

          t->setVertex(1,vv[Connect[i*nelemNodes+1]-1]);
          vv[Connect[i*nelemNodes+1]-1]->setElement(t);
          vv[Connect[i*nelemNodes+1]-1]->setValence(vv[Connect[i*nelemNodes+1]-1]->getValence()+1);

          t->setVertex(2,vv[Connect[i*nelemNodes+2]-1]);
          vv[Connect[i*nelemNodes+2]-1]->setElement(t);
          vv[Connect[i*nelemNodes+2]-1]->setValence(vv[Connect[i*nelemNodes+2]-1]->getValence()+1);

          t->setVertex(3,vv[Connect[i*nelemNodes+3]-1]);
          vv[Connect[i*nelemNodes+3]-1]->setElement(t);
          vv[Connect[i*nelemNodes+3]-1]->setValence(vv[Connect[i*nelemNodes+3]-1]->getValence()+1);

          break;

        case 8:
          t=new GoHexahedronElement;
          
          t->setVertex(0,vv[Connect[i*nelemNodes]-1]);
          vv[Connect[i*nelemNodes]-1]->setElement(t);
          vv[Connect[i*nelemNodes]-1]->setValence(vv[Connect[i*nelemNodes]-1]->getValence()+1);

          t->setVertex(1,vv[Connect[i*nelemNodes+1]-1]);
          vv[Connect[i*nelemNodes+1]-1]->setElement(t);
          vv[Connect[i*nelemNodes+1]-1]->setValence(vv[Connect[i*nelemNodes+1]-1]->getValence()+1);

          t->setVertex(2,vv[Connect[i*nelemNodes+2]-1]);
          vv[Connect[i*nelemNodes+2]-1]->setElement(t);
          vv[Connect[i*nelemNodes+2]-1]->setValence(vv[Connect[i*nelemNodes+2]-1]->getValence()+1);

          t->setVertex(3,vv[Connect[i*nelemNodes+3]-1]);
          vv[Connect[i*nelemNodes+3]-1]->setElement(t);
          vv[Connect[i*nelemNodes+3]-1]->setValence(vv[Connect[i*nelemNodes+3]-1]->getValence()+1);

          t->setVertex(4,vv[Connect[i*nelemNodes+4]-1]);
          vv[Connect[i*nelemNodes+4]-1]->setElement(t);
          vv[Connect[i*nelemNodes+4]-1]->setValence(vv[Connect[i*nelemNodes+4]-1]->getValence()+1);

          t->setVertex(5,vv[Connect[i*nelemNodes+5]-1]);
          vv[Connect[i*nelemNodes+5]-1]->setElement(t);
          vv[Connect[i*nelemNodes+5]-1]->setValence(vv[Connect[i*nelemNodes+5]-1]->getValence()+1);

          t->setVertex(6,vv[Connect[i*nelemNodes+6]-1]);
          vv[Connect[i*nelemNodes+6]-1]->setElement(t);
          vv[Connect[i*nelemNodes+6]-1]->setValence(vv[Connect[i*nelemNodes+6]-1]->getValence()+1);

          t->setVertex(7,vv[Connect[i*nelemNodes+7]-1]);
          vv[Connect[i*nelemNodes+7]-1]->setElement(t);
          vv[Connect[i*nelemNodes+7]-1]->setValence(vv[Connect[i*nelemNodes+7]-1]->getValence()+1); 

         break;

        default:
          Error("Unknown type of element",__FILE__,__LINE__);
          t=NULL;
          break;
        }

      if (t) {
        t->setId(i-1);
        ptGoMesh->addElement(t);
      }
    }
}
} // end of namespace
