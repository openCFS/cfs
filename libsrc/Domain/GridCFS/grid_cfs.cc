#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "grid_cfs.hh"
#include "grid.hh"
#include "elements_header.hh"
#include "baseelem.hh"
#include "conffile.hh"

namespace CoupledField
{

template<class Dim>
 GridCFS<Dim> :: GridCFS(FileType * const aptFileType)
{
#ifdef TRACES
  (*trace) << "entering GridCFS::GridCFS" << std::endl;
#endif

  InFile = aptFileType;

  conf->getsubdom(sd_);
  elems_=new std::vector<Elem*>[sd_.size()];  
  dim_ = 0;
}

template<>
void GridCFS<Point2D> :: Read()
{
#ifdef TRACE
  (*trace) << "entering GridCFS::Read 2D" << std::endl;
#endif
 
  dim_=InFile->ReadDim();
 
  InFile->ReadMaxnumnodes(maxnumnodes_);
  ptCoordinate_=new Point2D[maxnumnodes_];
  InFile->ReadCoordinate(ptCoordinate_, maxnumnodes_);
  InFile->ReadEl(elems_,sd_);
//   Integer i,j;
//   for(i=0; i<elems_[0].size(); i++){
//     Integer elemsize=(elems_[0][i]->connect).size();
//     testconnect_[i]=new Integer[elemsize];
//     for (j=0; j<elemsize; j++) {
//       testconnect_[i][j]=(elems_[0][i]->connect)[j];
//     }
//   }
      
#ifdef TRACE
  (*trace) << "leaving GridCFS::Read" << std::endl;
#endif
}

template<>
void GridCFS<Point3D> :: Read()
{
#ifdef TRACE
  (*trace) << "entering GridCFS::Read 3D" << std::endl;
#endif

  dim_=InFile->ReadDim();

  InFile->ReadMaxnumnodes(maxnumnodes_);
  ptCoordinate_=new Point3D[maxnumnodes_];
  InFile->ReadCoordinate(ptCoordinate_, maxnumnodes_);

  InFile->ReadEl(elems_,sd_);

#ifdef TRACE
  (*trace) << "leaving GridCFS::Read" << std::endl;
#endif
}

#ifdef ADAPTGRID
template<>
void GridCFS<Point2D>::putNodesFromGrid_RG(grd::MultilevelGrid * grid, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::putNodesFromGrid_R" << std::endl;
#endif

 Integer maxnumnodes= (*grid).getNoOfVertices();
 maxnumnodes_=maxnumnodes;
 ptCoordinate_=new Point2D[maxnumnodes];

  typedef std::list<grd::Vertex*>::iterator VerI;  
  
  Double * ps;
  Integer ilev, i=0;
  std::cout << " Total no. of vertices " << (*grid).getNoOfVertices() << std::endl;
  Integer topLevel = grid->getTopLevel();
  for (ilev=0; ilev<=topLevel; ilev++) {
    std::list<grd::Vertex*> *le=(*grid).getGridLevel(ilev)->getVertexList();
    for (VerI p=le->begin(); p!=le->end(); ++p) {
      Integer index = (*p)->getId();
      index--;
      if (index >= maxnumnodes) {
	std::cerr << " ERROR: catastrophic error, index overflow\n "
	          << " The index: " << index << "  the maxnumnodes: " << maxnumnodes << '\n';
      }
      ps=(*p)->getPosition();
      ptCoordinate_[index].x=ps[0];
      ptCoordinate_[index].y=ps[1];
      i++;
    }
  }
} // end of function putNodesFromGrid_RG

template<>
void GridCFS<Point3D>:: putNodesFromGrid_RG(grd::MultilevelGrid * grid, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::putNodesFromGrid_R" << std::endl;
#endif
  ;
}

template<>
void GridCFS<Point2D>::putElemsFromGrid_RG(grd::MultilevelGrid * grid, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::putElemsFromGrid_R" << std::endl;
#endif

  typedef std::list<grd::Element*>::iterator ElmI;
  ElmI p;

 Integer i, j, nnodes,type;

 std::list<grd::Element*> * le;
 std::list<grd::Element*> ** lt;

 Integer noOfLevels=grid->getNoOfLevels();
 for (j=0; j< noOfLevels; j++)
 {
   lt=grid->getGridLevel(j)->getElementList();
   type=0;

   while(lt[type]) {
     for ( p=lt[type]->begin() ; p!= lt[type]->end(); ++p) {
          if ((*p)->isRegular()) {
        Elem * el=new Elem();
	 
       Integer nnodes=(*p)->getNoOfVertices();
       Integer etype = (*p)->type();

       switch(etype) {
       case GRD_TRIANGLE:
        el->ptElem=ptTr;
        break;
       case GRD_QUADRANGLE:
        el->ptElem=ptQ;
        break;
       default:
	 Error("Unknown type of element in GridRG", __FILE__,__LINE__);
	 break;
       }

       (*el).connect.Resize(nnodes);
 
       for (i=0; i<nnodes; i++)
       {	
         (*el).connect[i]=((*p)->getVertex(i))->getId();	 
       }
 
       // test connection array
       Integer * testconnection=new Integer[nnodes];
       for (i=0; i<nnodes; i++)
	 testconnection[i]=((*p)->getVertex(i))->getId();
	 testconnect_.push_back(testconnection);

       Integer sd=(*p)->getValue();
       if (sd >= sd_.size()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);

       (*el).namesd=sd_[sd];

       elems_[sd].push_back(el);
     } // if isRegular
     else if ((*p)->isIrregular()) {
       grd::ConformingClosure closure;
       typedef grd::ConformingClosure::triangleIterator TriI;
       typedef grd::ConformingClosure::quadrangleIterator QuadI;
       (*p)->close(closure);
       
       // Process closing triangles
       for (TriI tri = closure.beginTriangle(); tri != closure.endTriangle(); ++tri) {
	 Elem * el=new Elem();
         Integer nnodes=(*tri)->getNoOfVertices();
         el->ptElem=ptTr;
         (*el).connect.Resize(nnodes);
         for (i=0; i<nnodes; i++) {
	   (*el).connect[i]=((*tri)->getVertex(i))->getId();
	 }
	 Integer sd=(*tri)->getValue();
	 if (sd >= sd_.size()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);

	 (*el).namesd=sd_[sd];
	 elems_[sd].push_back(el);
       } // for tri
       // Process now the quads
       for (QuadI quad = closure.beginQuadrangle(); quad != closure.endQuadrangle(); ++quad) {
	 Elem * el=new Elem();
         Integer nnodes=(*quad)->getNoOfVertices();
         el->ptElem=ptQ;
         (*el).connect.Resize(nnodes);
         for (i=0; i<nnodes; i++) {
	   (*el).connect[i]=((*quad)->getVertex(i))->getId();
	 }
	 Integer sd=(*quad)->getValue();
	 if (sd >= sd_.size()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);

	 (*el).namesd=sd_[sd];
	 elems_[sd].push_back(el);
       } // for quad
     } // else if isIrregular
   } // for element
   type++;
   } // end while(); list of elements types
 } // for level

} // end of function 

template<>
void GridCFS<Point3D>:: putElemsFromGrid_RG(grd::MultilevelGrid * grid, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::putElemsNodesFromGrid_R" << std::endl;
#endif

  ;
}


#endif // end of ADAPTGRID

template< class Dim>
void GridCFS<Dim>::GetConnection(Vector<Integer> & connection, const Integer iElem, const Integer level)
{
  Integer i,j;
  Integer sum=0,sum1;

  for (i=0; i<sd_.size(); i++) {
    sum1=sum;
    sum+=elems_[i].size();
    if (iElem < sum) {
      Integer elemsize=(elems_[i][iElem-sum1]->connect).size();
       connection.Resize(elemsize);
       for ( j=0; j < elemsize; j++) {
         connection[j]=(elems_[i][iElem-sum1]->connect)[j];
 	      } 
      break;
    }
  }
}

template<class Dim>
void GridCFS<Dim> :: GetCoordinateNode(const Integer inode, const Integer numlevel, Dim & rfPoint) 
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GetCoordinateNode" << std::endl;
#endif
  rfPoint=ptCoordinate_[inode];
}

template<class Dim>
void GridCFS<Dim> :: GetElemSD(std::vector<Elem*> & els, const std::string sd, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GetElemSD" << std::endl;
#endif

 Integer i;
 for (i=0; i<sd_.size(); i++)
  if (sd_[i]==sd) els=elems_[i];

}

template<class Dim>
void  GridCFS<Dim> :: GetCoordNodesElem(const Vector<Integer> connect, Dim * ptCoord, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS :: GetCoordinateNode" << std::endl;
#endif

  Integer k;
  for (k=0; k < connect.size(); k++)
       ptCoord[k]=ptCoordinate_[connect[k]-1];
 
}

template<class Dim>
GridCFS<Dim> ::~GridCFS() 
{ 
#ifdef TRACE
  (*trace) << " entering GridCFS<Dim> ::~GridCFS() \n";
#endif

if (ptCoordinate_) delete [] ptCoordinate_;

 Integer i,j; 
 if (elems_) {
   for (i=0; i<sd_.size(); i++) {
     for (j=0; j<elems_[i].size(); j++) 
       delete elems_[i][j];     
   }
   delete [] elems_;
 }

 for (i=0; i<testconnect_.size(); i++)
     if (testconnect_[i]) delete [] testconnect_[i];
}

template<class Dim>
Integer GridCFS<Dim>::GetMaxnumElem(const Integer numlevel)
{
 Integer i;
 Integer res=0;
 for (i=0; i<sd_.size(); i++)
  res+=elems_[i].size();

 return res;
}

} // end namespace

/*
template<class Dim>
void GridCFS<Dim>::forEachElemSd(PutElemMatInAlgSys & f,const std::string subdomain)
{
#ifdef TRACE
  (*trace) << "entering iterator GridCFS::forEachElemSD" << std::endl;
#endif

 Integer numsd, nesd; 
 std::vector<Elem> * els;

 Integer i;
 for (i=0; i<sd_.size(); i++)
 {
   if (sd_[i]==subdomain) numsd=i;
 }

 nesd=elems_[numsd].size();

 for (i=0; i< nesd; i++) 
  {
   f(elems_[numsd][i]);
  }

#ifdef TRACE
  (*trace) << "leaving iterator GridCFS::forEachElemSD" << std::endl;
#endif
}

template<class Dim>
void GridCFS<Dim>::forEachElemSd(PutElemMatAlgSysElst3d & f,const std::string subdomain)
{
#ifdef TRACE
  (*trace) << "entering iterator GridCFS::forEachElemSD" << std::endl;
#endif

 Integer numsd, nesd; 
 std::vector<Elem> * els;

 Integer i;
 for (i=0; i<sd_.size(); i++)
 {
   if (sd_[i]==subdomain) numsd=i;
 }

 nesd=elems_[numsd].size();

 for (i=0; i< nesd; i++) 
  {
   f(elems_[numsd][i]);
  }

#ifdef TRACE
  (*trace) << "leaving iterator GridCFS::forEachElemSD" << std::endl;
#endif
}
*/
