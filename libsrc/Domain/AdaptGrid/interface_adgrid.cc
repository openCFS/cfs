#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "conffile.hh"
#include "interface_adgrid.hh"
#include "elements_header.hh"

#include "acoustic2dPDE.hh"

namespace CoupledField
{

template<class Dim>
InterfaceAdaptGrid<Dim>::InterfaceAdaptGrid(FileType * aptFileType)
: Grid(aptFileType)
{
#ifdef TRACE
 (*trace) << "Entering InterfaceAdaptGrid<Dim>::InterfaceAdaptGrid<Dim>" << std::endl;
#endif

  ptFileType=aptFileType;
  conf->getsubdom(listSD_);
  lastlevel_=0;
}

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
 
 ptFileType->ReadGrid_RG(elems_,&vertex_,listSD_);

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

  ptFileType->ReadGrid_RG(elems_,&vertex_,listSD_);

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

// regular and unregular elements
template<>
void InterfaceAdaptGrid<Point2D>::GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)
{
 Integer i, nnodes, count=0;
 typedef std::list<grd::Element*>::iterator ElmI;
 ElmI p;

 std::list<grd::Element*> *le=grid_.getGridLevel(level)->getTriangleList();
 
 for ( p=le->begin() ; p!= le->end(); ++p, count++) 
{
 if (count==iElem)
{
 if ((*p)->isRegular())
  {
   Integer nnodes=(*p)->getNoOfVertices();
   connect.Resize(nnodes);
   for (i=0; i<nnodes; i++)
    {
     connect[i]=(*p)->getVertex(i)->getId();
    }    
  }
 break;
}
}

} // end of getConnection

template<>
void InterfaceAdaptGrid<Point3D>::GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)
{
 Error("Not implemented", __FILE__,__LINE__);
 ;
}

template<>
void InterfaceAdaptGrid<Point2D>::GetCoordinateNode(const Integer inode, const Integer level, Point2D & rfPoint) 
{ 
  typedef std::list<grd::Vertex*>::iterator VerI;  
  VerI p;
  double * ps;
  std::list<grd::Vertex*> *le=grid_.getGridLevel(level)->getVertexList();

  for (p=le->begin(); p!=le->end(); ++p) 
  {
    if ((*p)->getId() == (inode+1) )
    {
     ps=(*p)->getPosition();
     rfPoint.x=ps[0];
     rfPoint.y=ps[1];
     break; 
    }
  }
; 
}

template<>
void InterfaceAdaptGrid<Point3D>::GetCoordinateNode(const Integer inode, const Integer numlevel, Point3D & rfPoint)
{
 Error("Not implemented", __FILE__,__LINE__);
}

template<class Dim>
void  InterfaceAdaptGrid<Dim>::GetCoordNodesElem(const Vector<Integer> connect, Dim * ptCoord, const Integer level)
{
  Integer k;
  for (k=0; k < connect.size(); k++)
     GetCoordinateNode(connect[k]-1,level,ptCoord[k]);
}

template<class Dim>
Integer InterfaceAdaptGrid<Dim>::GetMaxnumnodes(const Integer numlevel)
  { return grid_.getNoOfVertices(numlevel); }

template<>
Integer InterfaceAdaptGrid<Point2D>::GetMaxnumElem(const Integer numlevel)
{
 return grid_.getNoOfElements();
}

template<>
Integer InterfaceAdaptGrid<Point3D>::GetMaxnumElem(const Integer numlevel)
{
 return grid_.getNoOfElements();
}

template<class Dim>
InterfaceAdaptGrid<Dim>::~InterfaceAdaptGrid()
{
;
}

template<class Dim>
void InterfaceAdaptGrid<Dim>::forEachElemSd(PutElemMatInAlgSys & f,const std::string subdomain)
{
#ifdef TRACE
  (*trace) << " entering InterfaceAdaptGrid<Dim>::forEachElemSd " << std::endl;
#endif

  Integer numsd;

 Integer i;
 for (i=0; i<listSD_.size(); i++)
 {
   if (listSD_[i]==subdomain) numsd=i;
 }

 
 Elem t;

 Integer toplevel=grid_.getTopLevel();
 grd::GridLevel * gridlv;
 grd::ConformingClosure * closure=grid_.getClosure();

 typedef std::list<grd::Element*>::iterator ElmI;
 typedef grd::ConformingClosure::triangleIterator  TriI;
  for (int j = 0; j <= toplevel; j++)
 {
    gridlv=grid_.getGridLevel(j);
    list<grd::Element*>** lt = gridlv->getElementList();
    int type = 0;
    while (lt[type]) {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->getValue() == numsd)
	  {
        if ((*p)->isRegular()) {
          Integer tp=(*p)->type();
	  if (tp == GRD_TRIANGLE) t.ptElem=ptTr;
          else
      Error("this type of elements is not implemented yet",__FILE__,__LINE__);

          Integer j;
          Integer nnodes=(*p)->getNoOfVertices();
          t.connect.Resize(nnodes);
          for(j=0; j<nnodes; j++)
	    t.connect[j]=(*p)->getVertex(j)->getId();

          f(t);
        }

	//        else if ((*p)->isIrregular()) {
        //  (*p)->close(*closure);
        //  for (TriI tri =(*closure).beginTriangle(); tri != (*closure).endTriangle(); ++tri)
        //    f(*tri);
       
	  } // check that from the subdomain
      } // loop elems
      // next element type
      type++;
 } // type loop
} // level loop
}

template<class Dim>
void InterfaceAdaptGrid<Dim>::forEachElemSd(PutElemMatAlgSysElst3d & f,const std::string subdomain)
{
 Integer numsd;

 Integer i;
 for (i=0; i<listSD_.size(); i++)
 {
   if (listSD_[i]==subdomain) numsd=i;
 }

  // typedef std::list<grd::Element*>::iterator ElmI;
//   typedef grd::ConformingClosure::tetrahedralIterator  TetI;
//   for (int j = 0; j <= topLevel; j++) {
//     list<Element*>** lt = grid[j]->getElementList();
//     int type = 0;
//     while (lt[type]) {
//       for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
//       {
//         if ((*p)->getValue() == color)
// 	  {
//         if ((*p)->isRegular()) {
//           f(*p);
//         }
//         else if ((*p)->isIrregular()) {
//           (*p)->close(closure);
//           for (TetI tri = closure.beginTriangle(); tri != closure.endTriangle(); ++tri)
//             f(*tri);
//         }
//           }
//       }
//       // next element type
//       type++;
//     }
//   }
//  //  grid_.forEachElementSd(f,numsd);
}

} // end of namespace

