#include <iostream>
#include <fstream>
#include <string>
#include <vector.h>


#include "interface_gridlib.hh"

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

template<class Dim>
void InterfaceGridlib<Dim>::GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodesPerElem, Dim * ptCoordElem)
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

template<class Dim>
void InterfaceGridlib<Dim>::GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint)
{
 float x,y,z;
 ptGoMesh->getVertex(inode,numlevel)->getPosition(x,y,z);
 rfPoint.x= double(x);
 rfPoint.y= double(y);
}

template<class Dim>
void InterfaceGridlib<Dim>::Read() 
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

  Dim * ptCoord=new Dim[nnodes]; 
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
} // end of namespace
