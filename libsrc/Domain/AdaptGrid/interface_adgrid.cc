#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "conffile.hh"
#include "interface_adgrid.hh"
#include "elements_header.hh"

#include "acoustic2dPDE.hh"

#include "outUnverg.hh"

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

  ptgridcfs_=NULL;

 ptBCs=new BCs(ptFileType);
 ptBCs->ReadBCs();
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
 // set nodes at level 0 as processed
 std::list<grd::Vertex*>* lt = (grid_.getGridLevel(0))->getVertexList();
 std::list<grd::Vertex*>::iterator p;
 for (p = lt->begin(); p != lt->end(); ++p) {
   (*p)->setProcessed();
 }

 TestRefine();
 

 // Refine();
 Trans2CFSGrid();

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
#ifdef TRACE
  (*trace) << " entering InterfaceAdaptGrid<Dim>::Refine() " << std::endl;
#endif

 grid_.refine();
 grid_.setVertexNumbers();
 lastlevel_++;

 Trans2CFSGrid();
}

template<class Dim>
void InterfaceAdaptGrid<Dim>::TestRefine()
{
#ifdef TRACE
  (*trace) << " entering InterfaceAdaptGrid<Dim>::Refine() " << std::endl;
#endif

 SetRefFlagTest f;
 forEachElemSd(f,listSD_[0]);
 grid_.refine();
 SetVertexNumbers();
 lastlevel_++;

 forEachElemSd(f,listSD_[0]);
 grid_.refine();
 SetVertexNumbers();

 forEachElemSd(f,listSD_[0]);
 grid_.refine();
 SetVertexNumbers();

 forEachElemSd(f,listSD_[0]);
 grid_.refine();
 SetVertexNumbers();

 Trans2CFSGrid();

//  ptBCs->Update(this);
//  ptBCs->printBCs(); 

}

template<class Dim>
void InterfaceAdaptGrid<Dim>::TestCoarse()
{
#ifdef TRACE
  (*trace) << " entering InterfaceAdaptGrid<Dim>::Refine() " << std::endl;
#endif

 SetRefFlagTest f;
 forEachElemSd(f,listSD_[0]);
 grid_.refine();
 SetVertexNumbers();
 lastlevel_++;

 forEachElemSd(f,listSD_[0]);
 grid_.refine();
 SetVertexNumbers();

 Trans2CFSGrid();

}

template<class Dim>
void InterfaceAdaptGrid<Dim>::UpdateBCs(std::list<Integer> * data, std::list<Integer> * result, std::vector<std::string> levels)
{
  Integer i;
  for (i=0; i<levels.size(); i++)
    {
      ;
    }
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

  ptgridcfs_->GetConnection(connect,iElem,level);

//  Integer i, j, nnodes, count=0,startId;
//  if (level==0) startId=0;
//  else startId=1;


//  typedef std::list<grd::Element*>::iterator ElmI;
//  ElmI p;

//  typedef grd::ConformingClosure::triangleIterator TriI;
//  std::list<grd::Element*> *le;
//  //=grid_.getGridLevel(level)->getTriangleList();
 
//  // for (j=0; j<=level; j++)
//  //  {
//  le=grid_.getGridLevel(level)->getTriangleList();

//  for ( p=le->begin() ; p!= le->end(); ++p, count++) 
// {
//  if (count==iElem)
// {
//   //  if ((*p)->isRefined()) std::cout << " refined " << std::endl;
//   //  if ((*p)->isRegular() || (*p)->isRefined()) 
//   //  {
//   //     if ( (*p)->isIrregular())
//     {
//     std::cout << " iireg " << std::endl;
//    Integer nnodes=(*p)->getNoOfVertices();
//    connect.Resize(nnodes);
//    for (i=0; i<nnodes; i++)
//     {
//      connect[i]=(*p)->getVertex(i)->getId()+startId;
//     }    
//    //  }
//    //  else if ((*p)->isIrregular())
//    // {
//    //    std::cout << " iireg " << std::endl;
// //     (*p)->close(closure_);
// //       for (TriI t=closure_.beginTriangle(); t != closure_.endTriangle(); ++t)
// //         {
// //           std::cout << " iireg " << std::endl;
// //              Integer nnodes=(*p)->getNoOfVertices();
// //    connect.Resize(nnodes);
// //    for (i=0; i<nnodes; i++)
// //     {
// //      connect[i]=(*p)->getVertex(i)->getId()+startId;
// //     } 
	
//   }

//  break;
// }
// } // end for elems
// // } // end for levels

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
  ptgridcfs_->GetCoordinateNode(inode,level,rfPoint);

//   typedef std::list<grd::Vertex*>::iterator VerI;  
  
//   Double * ps;
//   Integer ilev, startId;

//   if (level == 0) startId=1;
//   else startId=0;

//   for (ilev=0; ilev<=level; ilev++)
//     {
//   std::list<grd::Vertex*> *le=grid_.getGridLevel(ilev)->getVertexList();

//   for (VerI p=le->begin(); p!=le->end(); ++p) 
//   {
//     if ((*p)->getId() == (inode+startId) )
//     {
//      ps=(*p)->getPosition();
//      rfPoint.x=ps[0];
//      rfPoint.y=ps[1];
//      break; 
//     }
//   }
//     }
// ; 
}

template<>
void InterfaceAdaptGrid<Point3D>::GetCoordinateNode(const Integer inode, const Integer numlevel, Point3D & rfPoint)
{
 Error("Not implemented", __FILE__,__LINE__);
}

template<class Dim>
void  InterfaceAdaptGrid<Dim>::GetCoordNodesElem(const Vector<Integer> connect, Dim * ptCoord, const Integer level)
{
  ptgridcfs_->GetCoordNodesElem(connect,ptCoord,level);
//   Integer k;
//   for (k=0; k < connect.size(); k++)
//      GetCoordinateNode(connect[k]-1,level,ptCoord[k]);
}

template<class Dim>
Integer InterfaceAdaptGrid<Dim>::GetMaxnumnodes(const Integer numlevel)
  {
    return ptgridcfs_->GetMaxnumnodes(numlevel);
   
    //  return grid_.getNoOfVertices(numlevel);
 }

template<>
Integer InterfaceAdaptGrid<Point2D>::GetMaxnumElem(const Integer numlevel)
{
    return ptgridcfs_->GetMaxnumElem(numlevel);
  // return grid_.getNoOfElements();
}

template<>
Integer InterfaceAdaptGrid<Point3D>::GetMaxnumElem(const Integer numlevel)
{
 return ptgridcfs_->GetMaxnumnodes(numlevel);
//  return grid_.getNoOfElements();
}

 template<class Dim>
 void InterfaceAdaptGrid<Dim>::GetElemSD(std::vector<Elem> & elems, const std::string sd, const Integer level)
{
  return ptgridcfs_->GetElemSD(elems,sd,level);
}

template<class Dim>
InterfaceAdaptGrid<Dim>::~InterfaceAdaptGrid()
{
 if (ptgridcfs_) delete ptgridcfs_;
}


template<class Dim>
void InterfaceAdaptGrid<Dim>::forEachElemSd(SetRefFlag & f,const std::string subdomain)
{
#ifdef TRACE
 (*trace) << " entering InterfaceAdaptGrid<Dim>::forEachElemSd Refinement " << std::endl;
#endif
 
 Integer numsd;
 Integer i;
 for (i=0; i<listSD_.size(); i++)
 {
   if (listSD_[i]==subdomain) numsd=i;
 }

 Integer toplevel=grid_.getTopLevel();
 grd::GridLevel * gridlv;
 grd::ConformingClosure * closure=grid_.getClosure();

 Integer counter=0;
 typedef std::list<grd::Element*>::iterator ElmI;
 typedef grd::ConformingClosure::triangleIterator  TriI;
  for (int j = 0; j <= toplevel; j++)
 {
    gridlv=grid_.getGridLevel(j);
    list<grd::Element*>** lt = gridlv->getElementList();
    int type = 0;
    while (lt[type]) {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p, counter++)
      {
        if ((*p)->getValue() == numsd)
          {
        if ((*p)->isRegular()) {
           if (counter<6) f(*p);
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
void InterfaceAdaptGrid<Dim>::forEachElemSd(SetRefFlagTest & f,const std::string subdomain)
{
#ifdef TRACE
 (*trace) << " entering InterfaceAdaptGrid<Dim>::forEachElemSd Refinement " << std::endl;
#endif
 
 Integer numsd;
 Integer i;
 for (i=0; i<listSD_.size(); i++)
 {
   if (listSD_[i]==subdomain) numsd=i;
 }

 Integer toplevel=grid_.getTopLevel();
 grd::GridLevel * gridlv;
 grd::ConformingClosure * closure=grid_.getClosure();

//  Integer counter=0;
//  typedef std::list<grd::Element*>::iterator ElmI;
//  typedef grd::ConformingClosure::triangleIterator  TriI;
//   for (int j = 0; j <= toplevel; j++)
//  {
//     gridlv=grid_.getGridLevel(j);
//     list<grd::Element*>** lt = gridlv->getElementList();
//     int type = 0;
//     while (lt[type]) {
//       for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p, counter++)
//       {
//         if ((*p)->getValue() == numsd)
//           {
//         if ((*p)->isRegular()) {
//            if (counter<6) f(*p);
//         }
//           } // check that from the subdomain
//       } // loop elems
//       // next element type
//       type++;
//  } // type loop
// } // level loop

  Integer topLevel = grid_.getTopLevel();
  typedef std::list<grd::Element*>::iterator ElmI;
  for (int j = 0; j <= toplevel; j++)
  {
    Integer counter=0;
    gridlv=grid_.getGridLevel(j);
    list<grd::Element*>* lt = gridlv->getTriangleList();
    for (ElmI p = lt->begin(); p != lt->end(); ++p, counter++)
    {
      if ((*p)->getValue() == numsd)
      {
	if (topLevel < 3) {
	  if (!(*p)->isRefined()) {
	    if (counter<1) f(*p);
	  }
	}
	else {
	  if (j == 3) {
	    std::cerr << "Mark element for coarsening at level " << j << "\n";
	    (*p)->markForCoarsening();
	  }
	}
	  
      }
    } // for p list of triangle
  } // for j topLevel
 
} // end of functions

template<class Dim>
void  InterfaceAdaptGrid<Dim>::SetVertexNumbers()
 {
  Integer level;
  Integer counter;
  Integer noOfVertices = 0;
  Integer noOfLevels = grid_.getNoOfLevels();
  std::cerr << "No. of levels: " << noOfLevels << '\n';
  for (level = 0; level < noOfLevels; level++) {
    typedef std::list<grd::Vertex*>::iterator VI;
    std::list<grd::Vertex*>* lt = (grid_.getGridLevel(level))->getVertexList();
    for (VI p = lt->begin(); p != lt->end(); ++p) {
      grd::Vertex* vt = *p;
      if (vt->getId() != GRD_NEWVERTEX) {
        noOfVertices++;
      }
    }
  }

  std::cerr << "No. of vertices: " << noOfVertices << '\n';

  counter = 1;
  for (level = 1; level < noOfLevels; level++) {
    typedef std::list<grd::Vertex*>::iterator VI;
    std::list<grd::Vertex*>* lt = (grid_.getGridLevel(level))->getVertexList();
    for (VI p = lt->begin(); p != lt->end(); ++p) {
      grd::Vertex* vt = *p;
      if (vt->getId() == GRD_NEWVERTEX) {
	std::cerr << "Setting node number: " << (noOfVertices+counter) << "  at level: " << level << '\n';
	vt->setId(noOfVertices+counter);
        counter++;
      }
    }
  }
  std::cerr << "No. of vertices: " << (noOfVertices+counter-1) << '\n';
//   Integer toplevel=grid_.getTopLevel();
//   grd::GridLevel * pgl=grid_.getGridLevel(toplevel);
//   Integer novergrdprev=grid_.getGridLevel(toplevel-1)->noOfVertices();
//   std::cout << " nodes " << novergrdprev << std::endl;
//   novergrdprev++;
//   pgl->setVertexNumbers(novergrdprev);
 }

template<class Dim>
void InterfaceAdaptGrid<Dim>::Trans2CFSGrid(const Integer alevel)
{
#ifdef TRACE
  (*trace) << " entering InterfaceAdaptGrid<Dim>::Transform2CFSGrid()" << std::endl;
#endif

  if (ptgridcfs_) delete ptgridcfs_;  

  ptgridcfs_=new GridCFS<Dim>(ptFileType);

  Integer level=alevel;
  if (level==-1) level=grid_.getTopLevel();
  ptgridcfs_->putNodesFromGrid_RG(&grid_,level); 
  ptgridcfs_->putElemsFromGrid_RG(&grid_,level);
}

template<class Dim>
void InterfaceAdaptGrid<Dim>::ProlongSol(Vector<Double>& sol,const Integer alevel)
{
#ifdef TRACE
  (*trace) << " entering InterfaceAdaptGrid<Dim>::Transform2CFSGrid()" << std::endl;
#endif
  Integer level;
  Integer noOfLevels = grid_.getNoOfLevels();
  typedef std::list<grd::Edge*>::iterator EdgI;
  std::list<grd::Edge*>* edglt;
  typedef std::list<grd::Element*>::iterator ElmI;
  std::list<grd::Element*>* elmlt;

  // Process first the edges
  for (level = 1; level < (noOfLevels-1); level++) {
    edglt = (grid_.getGridLevel(level))->getEdgeList();
    for (EdgI e = edglt->begin(); e != edglt->end(); ++e) {
      grd::Edge* edge = *e;
      if (edge->isRefined()) {
	grd::Vertex* v = edge->getMidPoint();
	if (!v->isProcessed()) {
	  grd::Vertex* v0 = edge->getVertex(0);
	  grd::Vertex* v1 = edge->getVertex(1);
	  Integer id0 = v0->getId();
	  Integer id1 = v1->getId();
	  Integer id  = v->getId();

	  // compute new function values
	    
	  // set nodes as processed
	  v->setProcessed();
        } // if isProcessed 
      } // if isRefined
    } // for edge
  } // for level
 
  // Process the elements
  for (level = 1; level < noOfLevels; level++) {
    elmlt = (grid_.getGridLevel(level))->getTriangleList();
    for (ElmI p = elmlt->begin(); p != elmlt->end(); ++p) {
      grd::Element* elm = *p;
      Integer noOfVertices = elm->getNoOfVertices();
      for (Integer vertex = 0; vertex < noOfVertices; vertex++) {
	grd::Vertex* v = elm->getVertex(vertex);
	if (!v->isProcessed()) {
	  grd::Element* parentElm = elm->getParent();
	  
	  // compute new function values

	  // set node as processed
	  v->setProcessed();
	}
      }
    }
  } //for level
}

} // end of namespace

//   std::list<grd::Element *> ** pe=pgl->getElementList();
//     typedef std::list<grd::Element*>::iterator EllI;

//    Integer counter=0,type=0;
//    while(pe[type])
//      { for (EllI p=(*pe[type]).begin(); p!=(*pe[type]).end(); ++p)
//        {
// 	 std::cout << counter << " " << (*p)->getNoOfVertices() << std::endl;
//          counter++;
//          grd::Vertex * h;
//          std::cout << " connect ";
//          for (Integer k=0; k<(*p)->getNoOfVertices(); k++)
//        {
//          h=(*p)->getVertex(k);
// 	 std::cout << h->getId() <<" "; 
// 	     }
// 	 std::cout << std::endl;
//        } 
//      type++; 
//      }

// template<class Dim>
// void InterfaceAdaptGrid<Dim>::forEachElemSd(PutElemMatInAlgSys & f,const std::string subdomain)
// {
// #ifdef TRACE
//   (*trace) << " entering InterfaceAdaptGrid<Dim>::forEachElemSd " << std::endl;
// #endif

//   Integer numsd;

//  Integer i;
//  for (i=0; i<listSD_.size(); i++)
//  {
//    if (listSD_[i]==subdomain) numsd=i;
//  }
 
//  Elem t;

//  Integer toplevel=grid_.getTopLevel();
//  grd::GridLevel * gridlv;
//  grd::ConformingClosure * closure=grid_.getClosure();

//  typedef std::list<grd::Element*>::iterator ElmI;
//  typedef grd::ConformingClosure::triangleIterator  TriI;
//   for (int j = 0; j <= toplevel; j++)
//  {
//     gridlv=grid_.getGridLevel(j);
//     list<grd::Element*>** lt = gridlv->getElementList();
//     int type = 0;
//     while (lt[type]) {
//       for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
//       {
//         if ((*p)->getValue() == numsd)
// 	  {
//         if ((*p)->isRegular()) {
//           Integer tp=(*p)->type();
// 	  if (tp == GRD_TRIANGLE) t.ptElem=ptTr;
//           else
//       Error("this type of elements is not implemented yet",__FILE__,__LINE__);

//           Integer j;
//           Integer nnodes=(*p)->getNoOfVertices();
//           t.connect.Resize(nnodes);
//           for(j=0; j<nnodes; j++)
// 	    t.connect[j]=(*p)->getVertex(j)->getId();

//           f(t);
//         }

// 	//        else if ((*p)->isIrregular()) {
//         //  (*p)->close(*closure);
//         //  for (TriI tri =(*closure).beginTriangle(); tri != (*closure).endTriangle(); ++tri)
//         //    f(*tri);
       
// 	  } // check that from the subdomain
//       } // loop elems
//       // next element type
//       type++;
//  } // type loop
// } // level loop
// }

// template<class Dim>
// void InterfaceAdaptGrid<Dim>::forEachElemSd(PutElemMatAlgSysElst3d & f,const std::string subdomain)
// {
//  Integer numsd;

//  Integer i;
//  for (i=0; i<listSD_.size(); i++)
//  {
//    if (listSD_[i]==subdomain) numsd=i;
//  }

//  Error("Not implemented yet", __FILE__, __LINE__);

// }
