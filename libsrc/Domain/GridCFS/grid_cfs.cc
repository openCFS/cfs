#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

#include "grid_cfs.hh"
#include <Domain/grid.hh>
#include <Elements/elements_header.hh>
//#include <Elements/basefe.hh>
#include <DataInOut/conffile.hh>

namespace CoupledField
{

template<Integer dim>
 GridCFS<dim> :: GridCFS(FileType * const aptFileType)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GridCFS" << std::endl;
#endif

  InFile = aptFileType;

  conf->getsubdom(sd_);
  elems_=new std::vector<Elem*>[sd_.size()];  

  elNeighbors_=NULL;
} 

template<Integer dim>
void GridCFS<dim> :: Read()
{
#ifdef TRACE
  (*trace) << "entering GridCFS::Read " << std::endl;
#endif
 
  dim_=InFile->ReadDim();
 
  InFile->ReadMaxnumnodes(maxnumnodes_);
  ptCoordinate_=new Point<dim>[maxnumnodes_];
  InFile->ReadCoordinate(ptCoordinate_, maxnumnodes_);

  InFile->ReadEl(elems_,sd_);
      
  FormNeighborsLists();

#ifdef TRACE
  (*trace) << "leaving GridCFS::Read" << std::endl;
#endif
}

#ifdef ADAPTGRID
template<Integer dim>
void GridCFS<dim>::putNodesFromGrid_RG(grd::MultilevelGrid * grid, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::putNodesFromGrid_R" << std::endl;
#endif

 Integer maxnumnodes = (*grid).getNoOfVertices();
 maxnumnodes_=maxnumnodes;
 ptCoordinate_=new Point<dim>[maxnumnodes];

  typedef std::list<grd::Vertex*>::iterator VerI;  
  
  Double * ps;
  Integer ilev, i=0;
  std::cout << "\t\033[32m no. of vertices: \033[0m " << (*grid).getNoOfVertices() << std::endl;

  if (InfoPrint)
    (*infofile) <<  " total no. of vertices: " << (*grid).getNoOfVertices() << std::endl;

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
      ptCoordinate_[index][0]=ps[0];
      ptCoordinate_[index][1]=ps[1];
   if (dim==3)
   ptCoordinate_[index][2]=ps[2];
      i++;
    }
  }
} // end of function putNodesFromGrid_RG

template<>
void GridCFS<2>::putElemsFromGrid_RG(grd::MultilevelGrid * grid, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::putElemsFromGrid_RG" << std::endl;
#endif

  typedef std::list<grd::Element*>::iterator ElmI;
  ElmI p;

 Integer i, j, nnodes,type;

 std::list<grd::Element*> * le;
 std::list<grd::Element*> ** lt;

 // Element Map
 if (!elemMap_.empty()) {
   for (i = 0; i < elemMap_.size(); i++) {
     ElementMap* tmp = elemMap_[i];
     delete tmp;
   }
   elemMap_.clear();
 } 
 
 Integer noOfLevels=grid->getNoOfLevels();
 for (j=0; j< noOfLevels; j++)
 {
   lt=grid->getGridLevel(j)->getElementList();
   type=0;

   while(lt[type]) {
     for ( p=lt[type]->begin() ; p!= lt[type]->end(); ++p)
       {
	 if ((*p)->isRegular())
	   {
	     Elem * el=new Elem();
	     // Element maping
	     ElementMap* tmpMap = new ElementMap;
	     
	     Integer nnodes=(*p)->getNoOfVertices();
	     Integer etype = (*p)->type();
	     
	     switch(etype)
	       {
	       case GRD_TRIANGLE:
		 el->ptElem=ptTr1;
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
	     
	     Integer sd=(*p)->getValue();
	     if (sd >= sd_.size()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);
	     
	     (*el).namesd=sd_[sd];

	     elems_[sd].push_back(el);
	     
	     // put info in elemMap
	     Integer position = elems_[sd].size()-1;
	     tmpMap->map.resize(1);
	     tmpMap->map[0] = position;
	     tmpMap->sd = sd;
	     elemMap_.push_back(tmpMap);
	     
	   } // if isRegular
	 else if ((*p)->isIrregular()) {
	   ElementMap* tmpMap = new ElementMap;
	   tmpMap->sd = (*p)->getValue();
	   grd::ConformingClosure closure;
	   typedef grd::ConformingClosure::triangleIterator TriI;
	   typedef grd::ConformingClosure::quadrangleIterator QuadI;
	   (*p)->close(closure);
	   
	   // Process closing triangles
	   for (TriI tri = closure.beginTriangle(); tri != closure.endTriangle(); ++tri) {
	     Elem * el=new Elem();
	     Integer nnodes=(*tri)->getNoOfVertices();
	     el->ptElem=ptTr1;
	     (*el).connect.Resize(nnodes);
	     for (i=0; i<nnodes; i++)
	       {
		 (*el).connect[i]=((*tri)->getVertex(i))->getId();
	       }
	     Integer sd=(*tri)->getValue();
	     if (sd >= sd_.size()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);
	     
	     (*el).namesd=sd_[sd];
	     elems_[sd].push_back(el);
	     
	     // maping
	     Integer position = elems_[sd].size()-1;
	     tmpMap->map.push_back(position);
	   } // for tri
	   
	   // Process now the quads
	   for (QuadI quad = closure.beginQuadrangle(); quad != closure.endQuadrangle(); ++quad) {
	     Elem * el=new Elem();
	     Integer nnodes=(*quad)->getNoOfVertices();
	     el->ptElem=ptQ;
	     (*el).connect.Resize(nnodes);
	     for (i=0; i<nnodes; i++)
	       {
		 (*el).connect[i]=((*quad)->getVertex(i))->getId();
	       }
	     Integer sd=(*quad)->getValue();
	     if (sd >= sd_.size()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);
	     
	     (*el).namesd=sd_[sd];
	     elems_[sd].push_back(el);

	     // maping
	     Integer position = elems_[sd].size() - 1;
	     tmpMap->map.push_back(position);
	   } // for quad
	   // update element map list
	   elemMap_.push_back(tmpMap);
	 } // else if isIrregular
       } // for element
     type++;
   } // end while(); list of elements types
 } // for level
 
 if (InfoPrint)
   (*infofile) << " total number of elements (only for first subdomains): " << elems_[0].size() << std::endl;
     std::cerr << "\t\033[32m no. of elements: \033[0m " << elems_[0].size() << std::endl;
     
     FormNeighborsLists();
     
} // end of function 

template<>
void GridCFS<3>::putElemsFromGrid_RG(grd::MultilevelGrid * grid, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::putElemsNodesFromGrid_R" << std::endl;
#endif

  typedef std::list<grd::Element*>::iterator ElmI;
  ElmI p;

  std::vector<grd::Element*> tets;

  Integer i, j, nnodes,type;

  std::list<grd::Element*> * le;
  std::list<grd::Element*> ** lt;

  // Init Tets buffer
  tets.resize(4);
  for (i = 0; i < 4; i++)
    {
      tets[i] = new grd::Tetrahedron;
    }

 // Element Map
 if (!elemMap_.empty()) {
   for (i = 0; i < elemMap_.size(); i++) {
     ElementMap* tmp = elemMap_[i];
     delete tmp;
   }
   elemMap_.clear();
 } 
 
 Integer i1,i2,lev;
 Integer noOfLevels=grid->getNoOfLevels();
  for (j=0; j< noOfLevels; j++)
  {
    lt=grid->getGridLevel(j)->getElementList();
    type=2; // take only volume elementes

    while(lt[type])
    {
      for ( p=lt[type]->begin() ; p!= lt[type]->end(); ++p)
      {
        if ((*p)->isRegular())
        {
          Elem * el=new Elem();
          // Element maping
          ElementMap* tmpMap = new ElementMap;

          Integer nnodes=(*p)->getNoOfVertices();
          Integer etype = (*p)->type();

          Integer sd,position,nnds;

          switch(etype)
          {
            case GRD_TRIANGLE:
              break;
            case GRD_QUADRANGLE:
              break;
            case GRD_TETRAHEDRON:
              el->ptElem=ptTet;
              (*el).connect.Resize(nnodes);

              for (i=0; i<nnodes; i++) {
                (*el).connect[i]=((*p)->getVertex(i))->getId();
              }

              sd=(*p)->getValue();

              if (sd >= sd_.size())
                Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);

              (*el).namesd=sd_[sd];

              elems_[sd].push_back(el);

              // put info in elemMap
              position = elems_[sd].size()-1;
              tmpMap->map.resize(1);
              tmpMap->map[0] = position;
              tmpMap->sd = sd;

              elemMap_.push_back(tmpMap);
              break;
            case GRD_OCTAHEDRON:
              (*p)->getTetras(tets);
              Elem * elT[4];

              Integer it;
              nnds=4;

              // loop over tetrahedrals of octahedron
              for (it=0; it<4; it++)
              {
                elT[it]=new Elem();
                elT[it]->ptElem=ptTet;
                elT[it]->connect.Resize(nnds);
                // copy of connection array
                for (i=0; i<nnds; i++)
                {
                  elT[it]->connect[i]=(tets[it]->getVertex(i))->getId();
                }

                sd=(*p)->getValue();

                if (sd >= sd_.size()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);
                elT[it]->namesd=sd_[sd];

                elems_[sd].push_back(elT[it]);
                // mapping
                Integer position = elems_[sd].size()-1;
                tmpMap->map.push_back(position);
		tmpMap->sd=sd;
              } // end: loop over tetrahedrals of octahedron

              // update element map list
              elemMap_.push_back(tmpMap);
              break;
            case GRD_HEXAHEDRON:
           el->ptElem=ptHexa;
           break;
         default:
           Error("Unknown type of element in GridRG", __FILE__,__LINE__);
           break;
         } // end of switch
          }
        else if ((*p)->isIrregular()) {
          ElementMap* tmpMap = new ElementMap;
         tmpMap->sd = (*p)->getValue();
         grd::ConformingClosure closure;
         typedef grd::ConformingClosure::tetrahedronIterator TetI;
         (*p)->close(closure);
       
         // Process closing tetrahedrons
         for (TetI tet = closure.beginTetrahedron(); tet != closure.endTetrahedron(); ++tet)
           {
             Elem * el=new Elem();
             Integer nnodes=(*tet)->getNoOfVertices();
             el->ptElem=ptTet;
             (*el).connect.Resize(nnodes);
             for (i=0; i<nnodes; i++)
          {
            (*el).connect[i]=((*tet)->getVertex(i))->getId();
          }
             Integer sd=(*tet)->getValue();
             if (sd >= sd_.size()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);
             (*el).namesd=sd_[sd];
             elems_[sd].push_back(el);

             // maping
             Integer position = elems_[sd].size()-1;
             tmpMap->map.push_back(position);
           } // for tet
         // update element map list
         elemMap_.push_back(tmpMap);
       } // else if isIrregular
          } // for element
        type++;
      } // end while(); list of elements types
       } // for level

     // clean buffer of tets
     if (!tets.empty()) 
       {
    for (i = 0; i < 4; i++)
      {
        delete tets[i];
      }
    tets.clear();
       }

     if (InfoPrint)
       (*infofile) << " total number of elements: " << elems_[0].size() << std::endl;
     std::cerr << "\t\033[32m no. of elements: \033[0m " << elems_[0].size() << std::endl;
     FormNeighborsLists();
}

template<Integer dim>
void GridCFS<dim>::Refine(grd::MultilevelGrid& grid)
{
#ifdef TRACE
  (*trace) << " entering  GridCFS<Dim>::Refine " << std::endl;
#endif

    Integer i,j;
     Integer counter = 0;

     // Mesh refinement
     Integer noOfLevels = grid.getNoOfLevels();
     typedef std::list<grd::Element*>::iterator ElmI;

     Integer k;
     for (Integer j = 0; j < noOfLevels; j++)
       {
	 grd::GridLevel* gridlv = grid.getGridLevel(j);
	 list<grd::Element*>** lt = gridlv->getElementList();
	 Integer type;
	 if (dim==3) type=2;
	 else type=0;
	 while (lt[type])
	   {
	     for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
	       {
		 if (!(*p)->isRefined())
		   {
		     Integer  sde = (*p)->getValue();
		     Integer flag = 0;
		     ElementMap* map = elemMap_[counter];
		     
		     Integer sdm = map->sd;
		     if (sde != sdm) 
		       Error("Wrong number of subdomain",__FILE__,__LINE__);
		     
		     for (i = 0; i < map->map.size(); i++)
		       {
			 Integer elmId = map->map[i];
			 flag = elems_[sdm][elmId]->refinementFlag;
			 Integer numRefs=elems_[sdm][elmId]->refinementNumber;
			 if (flag == 1) 
			   {
			    
			     (*p)->markForRefinement(numRefs-1);
		    
			     break;
			   }

			 if (flag == -1)
			   {
			     (*p)->markForCoarsening();
			     break;
			   }
			 
		       }
		     
		     // update counter
		     counter++;
		   }
	       } // for loop elems
	     // next element type
	     type++;
	   } // type loop
       } // level loop

     grid.refine();      
}

template<Integer dim>
void GridCFS<dim>::ReRefine(grd::MultilevelGrid& grid)
{
#ifdef TRACE
     (*trace) << " entering  GridCFS<Dim>::ReRefine " << std::endl;
#endif
    
     Integer i,j;
     Integer counter = 0;

     // Mesh refinement
     Integer topLevel = grid.getTopLevel();
     typedef std::list<grd::Element*>::iterator ElmI;

     Integer k;
     grd::GridLevel* gridlv = grid.getGridLevel(topLevel);
     list<grd::Element*>** lt = gridlv->getElementList();
     Integer type;
     if (dim==3) type=2;
     else type=0;
     while (lt[type])
       {
	 for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
	   {
	     if (!(*p)->isRefined())
	       {
		 Integer  sde = (*p)->getValue();
		 grd::Element *parent = (*p)->getParent();
		 
		 Integer numRefs = parent->getNumOfRefinements(); 
		 if (numRefs)
		   (*p)->markForRefinement(numRefs-1);
	       } 
	   } // for loop elems
	     // next element type
	 type++;
       } // type loop
     grid.refine();     
}

template<Integer dim>
void GridCFS<dim>::RefineUniform(grd::MultilevelGrid& grid)
{
 Integer i,j;

 // Mesh refinement
 Integer noOfLevels = grid.getNoOfLevels();
 typedef std::list<grd::Element*>::iterator ElmI;

 for (Integer j = 0; j < noOfLevels; j++) {
   grd::GridLevel* gridlv = grid.getGridLevel(j);
   list<grd::Element*>** lt = gridlv->getElementList();
   Integer type = 0;
   if (dim == 3)
     type = 2;
   while (lt[type]) {
     for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p) {
       if (!(*p)->isRefined())
	 (*p)->markForRefinement();
     } // for loop elems
     // next element type
     type++;
   } // type loop
 } // level loop

 grid.refine();
}

template<Integer dim>
void GridCFS<dim>::SetRefinementFlag()
{
  Integer i,j;
  for (i=0; i<sd_.size(); i++) {
    for (j=0; j<elems_[i].size(); j++) {
      elems_[i][j]->refinementFlag=TRUE;
    }
  }
}

#endif // end of ADAPTGRID

template<Integer dim>
void GridCFS<dim>::GetConnection(Vector<Integer> & connection, const Integer iElem, const Integer level)
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

template<Integer dim>
void GridCFS<dim> :: GetCoordinateNode(const Integer inode, const Integer numlevel, Point<dim> & rfPoint) 
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GetCoordinateNode" << std::endl;
#endif
  rfPoint=ptCoordinate_[inode];
}

template<Integer dim>
void GridCFS<dim> :: GetElemSD(std::vector<Elem*> & els, const std::string sd, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS::GetElemSD" << std::endl;
#endif

 Integer i;
 for (i=0; i<sd_.size(); i++)
  if (sd_[i]==sd) els=elems_[i];

}

template<Integer dim>
void  GridCFS<dim> :: GetCoordNodesElem(const Vector<Integer> connect, Point<dim> * ptCoord, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS :: GetCoordinateNodesElem" << std::endl;
#endif

  Integer k;
  for (k=0; k < connect.size(); k++) 
       ptCoord[k]=ptCoordinate_[connect[k]-1];

}


template<Integer dim>
void  GridCFS<dim> :: GetCoordNodesElemMat(const Vector<Integer> connect, Matrix<Double>& coordMat, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS :: GetCoordinateNodesElemMat" << std::endl;
#endif

  coordMat.Resize(dim, connect.size());
  
  Integer k;
  for (k=0; k < connect.size(); k++)    
    for (int actDim=0; actDim < dim; actDim++)
      coordMat[actDim][k] = ptCoordinate_[connect[k]-1][actDim];
  
}


template<Integer dim>
GridCFS<dim> ::~GridCFS() 
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

 //! list with neighbors for element
 if (elNeighbors_) {    
   for (i = 0; i < sd_.size(); i++) {
     delete elNeighbors_[i];    
   }
   delete [] elNeighbors_;
 }

 // deconstructor for elemMap
 if (elemMap_.size() > 0) {
   for (i = 0; i < elemMap_.size(); i++) {
     ElementMap* tmp = elemMap_[i];
     delete tmp;
   }
   elemMap_.clear();
 }
}

template<Integer dim>
Integer GridCFS<dim>::GetMaxnumElem(const Integer numlevel)
{
 Integer i;
 Integer res=0;
 for (i=0; i<sd_.size(); i++)
  res+=elems_[i].size();

 return res;
}

template<Integer dim>
Integer GridCFS<dim>::GetMaxnumElem(const Integer numlevel, const std::vector<std::string> & subdoms)
{
  Integer res=0;
  Integer i,j;
  for (j=0; j<subdoms.size(); j++)
    for (i=0; i<sd_.size(); i++)
      if (sd_[i]==subdoms[j]) res+=elems_[i].size();

  return res;
}

template<Integer dim>
void GridCFS<dim>::FormNeighborsLists()
{
#ifdef TRACE
  (*trace) << " entering  GridCFS<Dim>::FormNeighborsLists " << std::endl;
#endif

  Integer i,j,k,n,m;

  Integer size=sd_.size(); // number of subdomains in domain
  elNeighbors_ = new  std::vector<std::vector<Elem*> >*[size];
  // initialize                                           
  // for each subdomain of grid we create vector( noOfElems) with vector of neighbors of each element
  for (i = 0; i < sd_.size(); i++) {
    size=elems_[i].size();
    elNeighbors_[i] = new std::vector<std::vector<Elem*> >(size);    
  }

  vtNeighbors_.resize(maxnumnodes_);     // vector with vector of neighbors-elements for each node

  // first main loop: form list with element-neighbors for each node
  for (i = 0; i < sd_.size(); i++) {
    for (j = 0; j < elems_[i].size(); j++) {
      Integer noOfVertices = elems_[i][j]->connect.size();
      for (k = 0; k < noOfVertices; k++) {
	Integer id = elems_[i][j]->connect[k];
	vtNeighbors_[id-1].push_back(elems_[i][j]);
      } // for loop over vertices, index k
    } // for loop over elements of the subdomain i, index j
  } // for loop over sd, index i

//   // check of nodes list
//   for (i=0; i<maxnumnodes_; i++) {
//     for (j=0; j<vtNeighbors_[i].size(); j++) {
//       cout << " size " << vtNeighbors_[i].size() << endl;
//       Elem* ptE=vtNeighbors_[i][j];
//       cout << " no of nodes " << i+1 << endl;
//       cout << ptE->connect << endl;
//     }
//   }
      
  // second main loop: form list with element-neighbors for each element
  for (i = 0; i < sd_.size(); i++) {   // do loop over subdomains
    for (j = 0; j < elems_[i].size(); j++) {   // do loop over elements of subdomain
      Integer noOfVertices = elems_[i][j]->connect.size(); 
      for (k = 0; k < noOfVertices; k++) {   // do loop over vertices of elements of the subdomain
	Integer id = elems_[i][j]->connect[k]; 
	for (n= 0; n < vtNeighbors_[id-1].size(); n++) {  // do loop over list of neighbors for node
	  Elem* ptel = vtNeighbors_[id-1][n];           

	  if (ptel != elems_[i][j]) {           // check that this element is not the same element for which we are looking for neighbors
	    Boolean flag = false;
	    for (m = 0; m < (*elNeighbors_[i])[j].size(); m++) { // check that this element is new in list of neighbors
	      Elem* ptel_tmp = (*elNeighbors_[i])[j][m];
	      if (ptel_tmp == ptel)
		flag = true;
	    } // for loop over neighbors, index m
	    if (!flag)
	      (*elNeighbors_[i])[j].push_back(ptel);
	  } // if (ptle != elem)
	} // for loop over vtNeighbors, index n
      } // for loop over vertices, index k
    } // for loop over elements of the subdomain, index j
  } // for loop over sd, index i

  //check of element-neighbors
//   for (i=0; i<sd_.size(); i++) {
//       for (j = 0; j < elems_[i].size(); j++) {  
// 	cout << " element: " << j << " with connect " << elems_[i][j]->connect << endl;	

// 	for (m = 0; m < (*elNeighbors_[i])[j].size(); m++) { // loop over neigbors for elem j
// 	  cout << " neighbors: " << m << endl;   
// 	  Elem * ptElm=(*elNeighbors_[i])[j][m];
// 	  cout << ptElm->connect;
// 	}
//       }
//   }

} // end of function FormNeighborsLists

//! return pointer to vector of element-neighbors for the element with number noOfElem
template<Integer dim>
std::vector<Elem*> * GridCFS<dim>::GetptNeighboursOfElem(const Integer noOfElem, const std::string color)
{
  
  std::vector<Elem*> * result;

  if (elNeighbors_==NULL) Error("You can't use function GetptNeighboursOfElem, since list of neighbors is not formed. You should use function FormNeighborsList before.");

  Integer i;
  for (i=0; i<sd_.size(); i++) {
    if (sd_[i]==color) {          
        result=&(*elNeighbors_[i])[noOfElem];
    }
  }

  return result;
}

 //! return vector of element-neighbors for the node with number noOfNode
template<Integer dim>
void  GridCFS<dim>::GetNeighboursOfNode(const Integer noOfNode, std::vector<Elem*> * neighbours)
{
  neighbours=&vtNeighbors_[noOfNode];
}

template<Integer dim>
void GridCFS<dim>::FormNeighbors4NodesOfElements(const std::vector<Elem*>& elems, std::vector<std::vector<Elem*> > &nodeNeighbors, std::vector<Integer> & map)
{
#ifdef TRACE
  (*trace) << " entering  GridCFS<Dim>::FormNeighbors4NodeOfElements " << std::endl;
#endif

  Integer iel,k;
  CalcNumberOfNodesInPatch(elems,map);
  Integer maxnumnodes=map.size();
 // calculation number of nodes in patch of elements

   // vector with vector of neighbors-elements for each node
  nodeNeighbors.resize(maxnumnodes);    

  // first main loop: form list with element-neighbors for each node
    for (iel = 0; iel < elems.size(); iel++) {
      Integer noOfVertices = elems[iel]->connect.size();
      for (k = 0; k < noOfVertices; k++) {
	Integer id = elems[iel]->connect[k];
	Integer imp;
	for (imp=0;imp<map.size();imp++) 
	  {
	    if (id==map[imp]) 
	      {
		break;
	      }	    
	  }	
	nodeNeighbors[imp].push_back(elems[iel]);
      } // for loop over vertices, index k
    } // for loop over elements iel

} // end of function FormNeighbors4NodesOfElements

template<Integer dim>
void GridCFS<dim>::DefineBelonging4Elems(const std::vector<Elem*>& elemsSurf, const std::vector<Elem*>&elems, std::vector<Elem*> & belongingSE)
{
#ifdef TRACE
  (*trace) << " entering  GridCFS<Dim>:DefineBelonging4Elems " << std::endl;
#endif

  Integer noOfSurfElems=elemsSurf.size();
  belongingSE.resize(noOfSurfElems);

  // TEST Routine for checking elements
  // std::cerr << "elems.size() = " << elems.size() << std::endl;
  // for (Integer i=0; i<elems.size(); i++)
  //   std::cerr << "elems[" << i <<"].connect = " << elems[i]->connect << std::endl;

  // form list with neighbors for each nodes in patch of boundary elements
  std::vector<Integer> map;
  std::vector<std::vector<Elem*> > listNeighbors;
  FormNeighbors4NodesOfElements(elems,listNeighbors,map);

  Integer ise,je;
  for (ise=0; ise<noOfSurfElems; ise++) { // loop over surface elements

    Boolean FoundNd=FALSE;
    Elem * ptAuxElem;

    Vector<Integer> &connectSE=elemsSurf[ise]->connect;

    // get list of neighbors for first node of the surface element
    Integer imp;   // get local number for this node
    for (imp=0; imp<map.size(); imp++) 
      {
	if (connectSE[0]==map[imp]) 
	  break;
      }
    
    std::vector<Elem*> &listNeigh4Elem=listNeighbors[imp];
       
    // loop over list of neighbors
    Integer ine;
    for (ine=0;ine<listNeigh4Elem.size();ine++) {
      ptAuxElem=listNeigh4Elem[ine];

      // check is there other vertices of element
      // loop over other nodes of surf element
      for (je=1;je<connectSE.size();je++) {
	Integer verSE=connectSE[je];

	//loop over vertices of the element
	FoundNd=FALSE;
	Vector<Integer> &vertices=ptAuxElem->connect;
	Integer ivt;
	for(ivt=0;ivt<vertices.size();ivt++) {
	  if (verSE==vertices[ivt]) {
	    FoundNd=TRUE;
	    break;
	  }
	} // end of loop over vertices of neigh-element
	
	if (!FoundNd) break;
      } // end of loop over nodes of surf element

      if (FoundNd) {
	belongingSE[ise]=ptAuxElem;
	break;
      }
    } // end loop over neighbors 

  } // loop over Surf element
	     
}


template<Integer dim>
void GridCFS<dim>::CalcNumberOfNodesInPatch(const std::vector<Elem*> & patch, std::vector<Integer> & map)
{
#ifdef TRACE
  (*trace) << "entering GridCFS<Dim>::CalcNumberOfNodesInPatch" << std::endl;
#endif

  Integer iels,ivc,imp;
  Vector<Integer> vec_connect;
  Boolean NewNode;

  for (iels=0; iels<patch.size(); iels++) // loop over elements in patch
    {
     
      vec_connect=patch[iels]->connect;
      
      for (ivc=0; ivc<vec_connect.size(); ivc++) { 
	 NewNode=TRUE;
	// loop over vector with global nodes for previous elements
	for (imp=0; imp<map.size(); imp++) {
	  // check that this node is not new
	  if (map[imp] == vec_connect[ivc]) { 
	    NewNode=FALSE;
	  }	 
	}

      if (NewNode) {
	map.push_back(vec_connect[ivc]);
      }

      } // end of loop over nodes in element

    } // end of loop over elements in patch   
}

template<Integer dim>
void GridCFS<dim>::GetInterfaceNeighbours(std::vector<Integer> & interfaceNodes, std::vector<std::string> & subdoms, std::vector<Elem*> & neighbours, Integer level)
{
#ifdef TRACE
  (*trace) << "entering GridCFS<Dim>::GetInterfaceNeighbours" << std::endl;
#endif
  
  Boolean belongs2Interface;
  std::vector<Elem*> elems;
  std::vector<Integer> map;
  
  //std::cerr << "In GetInterfaceNeighbourS" << std::endl;
  //std::cerr << "InterfaceNodes.size = " << interfaceNodes.size() << std::endl;

  // loop over all subdomains
  for (Integer isd=0; isd<subdoms.size(); isd++)
    {
      //std::cerr << "-----------------------------------" << std::endl;
      //std::cerr << "examing Subdomain " << subdoms[isd] << std::endl;
      GetElemSD(elems, subdoms[isd], level);

      // loop over all elements in subdomain
      for (Integer iNS=0; iNS < elems.size(); iNS++)
	{

	  Elem *aux = elems[iNS];
	  Vector<Integer>  aux_connect = aux->connect;
	  //std::cerr << " ====" << std::endl;
	  // std::cerr << "examing element Nr " << iNS << std::endl;
	  //std::cerr << "connect = " << aux_connect << std::endl;
	  
	  belongs2Interface = false;
	  
	  // check if any node is common in Interface
	  for (Integer inode=0; inode<aux_connect.size(); inode++) {

	    for (Integer nnode=0; nnode<interfaceNodes.size(); nnode++) {

	      if (interfaceNodes[nnode] == aux_connect[inode]) {
		belongs2Interface = true;
		break;
	      }
	    }
	  }
	  
	  if (belongs2Interface)
	    neighbours.push_back(elems[iNS]);
	}
    }
}
  
template<>
Double GridCFS<2>::CalcAreaElem(const Elem* elem)
{
#ifdef TRACE
  (*trace) << "entering GridCFS<Dim>::CalcAreaElem" << std::endl;
#endif

  Double res;
  Double a,b,c,s;
  Point<2> A,B,C,D;

  Integer type=elem->ptElem->feType();
  switch(type) {
  case 0: // line
    A=ptCoordinate_[elem->connect[0]-1];
    B=ptCoordinate_[elem->connect[1]-1];
    res=dist(A,B);
    break;
  case 1:     // triangle
    A=ptCoordinate_[elem->connect[0]-1];
    B=ptCoordinate_[elem->connect[1]-1];
    C=ptCoordinate_[elem->connect[2]-1];
  
    a=dist(A,B);
    b=dist(B,C);
    c=dist(A,C);
  
    s=(a+b+c)/2;
    res=sqrt(s*(s-a)*(s-b)*(s-c));

    break; 
  case 2:     // quadrilateral
    Double res1, res2;
    A=ptCoordinate_[elem->connect[0]-1];
    B=ptCoordinate_[elem->connect[1]-1];
    C=ptCoordinate_[elem->connect[2]-1];
    D=ptCoordinate_[elem->connect[3]-1];
   
    a=dist(A,B);
    b=dist(B,C);
    c=dist(A,C);
    s=(a+b+c)/2;
    res1=sqrt(s*(s-a)*(s-b)*(s-c));
    a=dist(C,D);
    b=dist(A,D);
    s=(a+b+c)/2;
    res2=sqrt(s*(s-a)*(s-b)*(s-c));
    res=res1+res2;
    break;

  default:
    Error("There isn't such kind of element's type",__FILE__,__LINE__);
  }

  return res;
}




template<>
Double GridCFS<3>::CalcAreaElem(const Elem* elem)
{
  Error("Not implemented yet",__FILE__,__LINE__);

  //just for SUN compiler
  Double dummy;
  return dummy;
}

} // end namespace


