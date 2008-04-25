// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <DataInOut/ParamHandling/ConfFile.hh>
#include "interface_adgrid.hh"
#include <Elements/elements_header.hh>

namespace CoupledField
{

  template<Integer Dim>
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

  }

  template<Integer dim>
  void InterfaceAdaptGrid<dim>::Read() 
  {
#ifdef TRACE
    (*trace)<< "Entering InterfaceAdaptGrid::Read 2D" << std::endl;
#endif
  
    dim_=ptFileType->ReadDim();

    Integer nnodes;
    ptFileType->ReadMaxnumnodes(nnodes); 

    Point<dim> * ptCoord=new Point<dim>[nnodes]; 
    ptFileType->ReadCoordinate(ptCoord, nnodes);

    // read id of bnd nodes
    std::vector<Integer> idBCs;
    std::vector<Integer> colorBCs;
    ptFileType->ReadBCs_GridRG(idBCs,colorBCs);

    // read vertices  
    Integer inode;
    Integer ibnd=0;
    vertices_.resize(nnodes);  
    for (inode=0; inode<nnodes; inode++)
      {
        Double pos[3];
        pos[0]=ptCoord[inode][0];
        pos[1]=ptCoord[inode][1];
        if (dim==2)
          pos[2]=0.0;
        else
          pos[2]=ptCoord[inode][2];

        grd::Vertex* tmpVt = new grd::Vertex(pos);
        tmpVt->setId(inode+1);

        //check is that vertice boundary or not
        std::vector<Integer>::const_iterator p;
        ibnd = 0;
        for ( p=idBCs.begin(); p!=idBCs.end(); p++,ibnd++) {
          if ((inode+1) == *p) {      
            tmpVt->setBoundaryNode();
            tmpVt->setColorBndNode(colorBCs[ibnd]);
          }
        }
        vertices_[inode]= tmpVt;
      }

    if (ptCoord) delete [] ptCoord;

    // read elements
    ptFileType->ReadGrid_RG(elements_,&vertices_,listSD_);

    grid_.buildCoarseMesh(vertices_,elements_);

    // set nodes at level 0 as processed: 
    // needed for transient problem in restoring solution
    std::list<grd::Vertex*>* lt = (grid_.getGridLevel(0))->getVertexList();
    std::list<grd::Vertex*>::iterator p;
    for (p = lt->begin(); p != lt->end(); ++p) {
      (*p)->setProcessed();
    }
  
    // set anisotropic refinement. only for Quadrangle
    typedef std::list<grd::Element*>::iterator Eit;
    std::list<grd::Element*>* elt = (grid_.getGridLevel(0))->getQuadrangleList();
    for (Eit p = elt->begin(); p != elt->end(); ++p) {
      const double dir[3] = {0.0,1.0,0.0};
      (*p)->setAnisotropy(dir);
    }
  
    // transfer this grid to GridCFS, since we do computation only with GridCFS 
 
    Trans2CFSGrid();

#ifdef TRACE
    (*trace) << "Leaving InterfaceAdaptGrid<Dim>::Read " << std::endl;
#endif
  }

  template<Integer dim>
  void InterfaceAdaptGrid<dim>::Refine(const Integer numLoops)
  {
#ifdef TRACE
    (*trace) << " entering InterfaceAdaptGrid<Dim>::Refine() " << std::endl;
#endif

    ptgridcfs_->Refine(grid_);
  
    Integer i;
    for(i=0; i<numLoops-1; i++)
      ptgridcfs_->ReRefine(grid_);

    SetVertexNumbers();
    Trans2CFSGrid();

  }

  template<Integer dim>
  void InterfaceAdaptGrid<dim>::RefineUniform()
  {
#ifdef TRACE
    (*trace) << " entering InterfaceAdaptGrid<Dim>::RefineUniform() " << std::endl;
#endif

    ptgridcfs_->RefineUniform(grid_);
    SetVertexNumbers();
    Trans2CFSGrid();
  }

  template<Integer dim>
  void InterfaceAdaptGrid<dim>::UpdateBCs(std::list<Integer> * bcs)
  {
#ifdef TRACE
    (*trace) << " entering InterfaceAdaptGrid<Dim>::UpdateBCs " << std::endl;
#endif

    Integer level;
    Integer noOfLevels = grid_.getNoOfLevels();

    std::list<grd::Vertex*> * lv;
    typedef std::list<grd::Vertex*>::iterator VerI;
    VerI p;

    for (level=0; level<noOfLevels; level++)
      {
        lv=(grid_.getGridLevel(level))->getVertexList();
        for (p=lv->begin(); p!=lv->end(); ++p)
          {
            if ((*p)->isBoundaryNode())
              { 
                bcs[(*p)->getColorBndNode()].push_back((*p)->getId());
              }
          }
      }

    std::cout << " size of BCs: " << (*bcs).size() << std::endl;

  }

  // regular and unregular elements
  template<Integer dim>
  void InterfaceAdaptGrid<dim>::GetConnection(Vector<Integer> & connect, const Integer iElem, const Integer level)
  {
    ptgridcfs_->GetConnection(connect,iElem,level);
  }

  template<Integer dim>
  void InterfaceAdaptGrid<dim>::GetCoordinateNode(const Integer inode, const Integer level, Point<dim> & rfPoint) 
  { 
    ptgridcfs_->GetCoordinateNode(inode,level,rfPoint);
  }

  template<Integer dim>
  void  InterfaceAdaptGrid<dim>::GetCoordNodesElem(const Vector<Integer> connect, Point<dim> * ptCoord, const Integer level)
  {
    ptgridcfs_->GetCoordNodesElem(connect,ptCoord,level);
  }

  template<Integer dim>
  void  InterfaceAdaptGrid<dim>::GetCoordNodesElemMat(const Vector<Integer> connect,
                                                        Matrix<Double>& coordMat, const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering InterfaceAdaptGrid::GetCoordinateNodesElemMat" << std::endl;
#endif

    ptgridcfs_->GetCoordNodesElemMat(connect, coordMat, level);

  }

  template<Integer dim>
  Integer InterfaceAdaptGrid<dim>::GetMaxnumnodes(const Integer numlevel)
  {
    return ptgridcfs_->GetMaxnumnodes(numlevel);
  }

  template<Integer dim>
  Integer InterfaceAdaptGrid<dim>::GetMaxnumElem(const Integer numlevel)
  {
    return ptgridcfs_->GetMaxnumElem(numlevel);
  }

  template<Integer dim>
  Integer InterfaceAdaptGrid<dim>::GetMaxnumElem(const Integer numlevel, const std::vector<std::string> & subdoms)
  {
    return ptgridcfs_->GetMaxnumElem(numlevel,subdoms); 
  }

  template<Integer dim>
  void InterfaceAdaptGrid<dim>::GetElemSD(std::vector<Elem*> & elems, const std::string sd, const Integer level)
  {
    return ptgridcfs_->GetElemSD(elems,sd,level);
  }

  template<Integer dim>
  std::vector<std::string>* InterfaceAdaptGrid<dim>::GetAllSDs()
  {
    return ptgridcfs_->GetAllSDs();
  }

  template<Integer dim>
  InterfaceAdaptGrid<dim>::~InterfaceAdaptGrid()
  {
    delete ptBCs;
    if (ptgridcfs_) delete ptgridcfs_;
  }

  template<Integer dim>
  void  InterfaceAdaptGrid<dim>::SetVertexNumbers()
  {
    Integer level;
    Integer counter;
    Integer noOfVertices = 0;
    Integer noOfLevels = grid_.getNoOfLevels();

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

    counter = 1;
    for (level = 1; level < noOfLevels; level++) {
      typedef std::list<grd::Vertex*>::iterator VI;
      std::list<grd::Vertex*>* lt = (grid_.getGridLevel(level))->getVertexList();
      for (VI p = lt->begin(); p != lt->end(); ++p) {
        grd::Vertex* vt = *p;
        if (vt->getId() == GRD_NEWVERTEX) {
          vt->setId(noOfVertices+counter);
          counter++;
        }
      }
    }
  }

  template<Integer dim>
  void InterfaceAdaptGrid<dim>::Trans2CFSGrid(const Integer alevel)
  {
#ifdef TRACE
    (*trace) << " entering InterfaceAdaptGrid<Dim>::Transform2CFSGrid()" << std::endl;
#endif

    if (ptgridcfs_) { delete ptgridcfs_;} 

    ptgridcfs_=new GridCFS<dim>(ptFileType);

    Integer level=alevel;
    if (level==-1) level=grid_.getTopLevel();
    ptgridcfs_->SetDim(dim_);
    ptgridcfs_->putNodesFromGrid_RG(&grid_,level); 
    ptgridcfs_->putElemsFromGrid_RG(&grid_,level);

  }

  /*
    template<Integer dim>
    void InterfaceAdaptGrid<dim>::ResetToCoarseGrid()
    {
    #ifdef TRACE
    (*trace) << " entering InterfaceAdaptGrid<Dim>::ResetToCoarseGrid()" << std::endl;
    #endif

    if (ptgridcfs_) { delete ptgridcfs_;} 
    ptgridcfs_=NULL;
  
    grid_.reset();

    Read();

    std::cout << " no. of elements: " << grid_.getNoOfTriangles() << " level: " << grid_.getTopLevel() << std::endl;
    std::cout <<  " no. of nodes: " << grid_.getNoOfVertices() << " level: " << grid_.getTopLevel() << std::endl;
    }
  */

  template<Integer dim>
  void InterfaceAdaptGrid<dim>::ProlongSol(const Vector<Double> sol_coarse, Vector<Double> &sol, const Integer alevel)
  {
#ifdef TRACE
    (*trace) << " entering InterfaceAdaptGrid<Dim>::ProlongSol()" << std::endl;
#endif
  
    Integer numnodes=grid_.getNoOfVertices();
    sol.Resize(numnodes);
    sol.Init();
    Integer i;
    for (i=0; i < sol_coarse.size(); i++) {
      sol[i]=sol_coarse[i];
    }    

    Integer level;
    Integer noOfLevels = grid_.getNoOfLevels();
    typedef std::list<grd::Edge*>::iterator EdgI;
    std::list<grd::Edge*>* edglt;
    typedef std::list<grd::Element*>::iterator ElmI;
    std::list<grd::Element*>* elmlt;

    // Process first the edges
    for (level = 0; level < noOfLevels; level++) {
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

          
            sol[id-1]=(sol[id0-1]+sol[id1-1])/2; // compute new function values
            
            // set nodes as processed
            v->setProcessed();
          } // if isProcessed 
        } // if isRefined
      } // for edge
    } // for level
 
    // Process the elements
    for (level = 1; level < noOfLevels; level++) {
      elmlt = (grid_.getGridLevel(level))->getQuadrangleList();
      for (ElmI p = elmlt->begin(); p != elmlt->end(); ++p) {
        grd::Element* elm = *p;
        Integer noOfVertices = elm->getNoOfVertices();
        for (Integer vertex = 0; vertex < noOfVertices; vertex++) {
          grd::Vertex* v = elm->getVertex(vertex);
          if (!v->isProcessed()) {        
            grd::Element* parentElm = elm->getParent();
          
            grd::Vertex* v0=elm->getVertex(0);
            grd::Vertex* v1=elm->getVertex(1);
            grd::Vertex* v2=elm->getVertex(2);
            grd::Vertex* v3=elm->getVertex(3);

            Double * p0=v0->getPosition();
            Double * p1=v1->getPosition();
            Double * p2=v2->getPosition();
            Double * p3=v3->getPosition();
            Double * pmid=v->getPosition();

            Integer id0=v0->getId();        
            Integer id1=v1->getId();
            Integer id2=v2->getId();
            Integer id3=v3->getId();
            Integer id=v->getId();
          
            Double d0=grd::distance(p0,pmid);
            Double d1=grd::distance(p1,pmid);
            Double d2=grd::distance(p2,pmid);
            Double d3=grd::distance(p3,pmid);

            sol[id-1]=(d0*sol[id0-1]+d1*sol[id1-1]+d2*sol[id2-1]+d3*sol[id3-1])/(d0+d1+d2+d3);

            // set node as processed
            v->setProcessed();
          }
        }
      }
    } //for level
  }  

  //! return vector of element-neighbors for the element with number noOfElem
  template <Integer dim>
  std::vector<Elem*>* InterfaceAdaptGrid<dim>::GetNeighboursOfElem(const Integer noOfElem, std::string color)
  {
#ifdef TRACE
    (*trace) << " entering InterfaceAdaptGrid<Dim>::GetNeighboursOfElem " << std::endl;
#endif
    std::vector<Elem*>*  neighbors=ptgridcfs_->GetptNeighboursOfElem(noOfElem,color);
    return neighbors;
  }


  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate InterfaceAdaptGrid<2>
#pragma instantiate InterfaceAdaptGrid<3>
#endif
} // end of namespace
