#include <iostream>
#include <fstream>
#include <string>
#include <vector.h>

// from NetGen
#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

#include "interface_netgen.hh"
#include "filetype.hh"

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
  // ########################## number of groupes
  ptFileType->ReadElemConnectionGH(nelems, Connect, nelemNodes, 0);

  Integer i;
  for (i=0; i<nelems; i++)
    {
      Element2d el;
      el.SetIndex(i+1);

      switch (nelemNodes)
	{
	case 3:

            el.PNum(1)=Connect[i*nelemNodes];
            el.PNum(2)=Connect[i*nelemNodes+1];
            el.PNum(2)=Connect[i*nelemNodes+2];

	  break;

	default:
	  Error("Unknown type of element");
	  break;
	}
      
      mesh.AddSurfaceElement (el);
   }

   delete [] Connect;

      mesh.ClearFaceDescriptors();
      mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));  
}

  template<class Dim>
  void InterfaceNetGen<Dim>::SubdivideUniform(const Integer level)
  { Error("Not implemented yet",__FILE__,__LINE__); }

  template<class Dim>
  void InterfaceNetGen<Dim>::GetCoordOfNodesElem(const Integer numElem, const Integer numlevel, const Integer numnodes, Dim * ptCoordElem)
   { Error("Not implemented yet",__FILE__,__LINE__); }

  template<class Dim>
   void InterfaceNetGen<Dim>::GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint) 
  { Error("Not implemented yet",__FILE__,__LINE__); }

  template<class Dim>
     void InterfaceNetGen<Dim>::GetConnection(Integer * result, const Integer level,
           const Integer numElem, const Integer numnodesPerElem)
   { Error("Not implemented yet",__FILE__,__LINE__); }

  template<class Dim>
     Integer InterfaceNetGen<Dim>::GetMaxnumnodes(const Integer numlevel)
  { Error("Not implemented yet",__FILE__,__LINE__); }

  template<class Dim>
     Integer InterfaceNetGen<Dim>::GetMaxnumElem(const Integer numlevel)
      { Error("Not implemented yet",__FILE__,__LINE__); }

  template<class Dim>
    Integer InterfaceNetGen<Dim>::GetNumNodesPerElem(const Integer iElem, const Integer level)
      { Error("Not implemented yet",__FILE__,__LINE__); }

  template<class Dim>
   void InterfaceNetGen<Dim>::PrintCoordinate(const Integer level, std::ostream * out) const
    { Error("Not implemented yet",__FILE__,__LINE__); }

  template<class Dim>
   void InterfaceNetGen<Dim>::GetNodesBoundaryCondition(Vector<Integer> & nodesDirBC, const Integer level)
   { Error("Not implemented yet",__FILE__,__LINE__); }

} // end of namespace

