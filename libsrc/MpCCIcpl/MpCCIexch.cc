#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "MpCCIexch.hh"
#include <DataInOut/ParamHandling/ConfFile.hh>
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

 //#ifdef MpCCI
#include <cci.h>
 //#endif

namespace CoupledField
{

MpCCIexch::MpCCIexch(Grid *aptgrid, Integer nNodesSD)
{
#ifdef TRACE
  (*trace) << "entering MpCCIexch::MpCCIexch " << std::endl;
#endif
    std::cout<<"Nodes SD: "<<nNodesSD<<std::endl;
  ptgrid_ = aptgrid;
  // MpCCInodes_ = ptgrid_->GetMaxnumnodes(0);
  MpCCInodes_ = nNodesSD;

  Dim_ = ptgrid_->GetDim();


  //Get general specific coupling description for CFS++ side
  params->Get("meshId", meshId_, "MpCCI-flownoise");
  params->Get("partId",partId_, "MpCCI-flownoise");
  params->Get("nNodeIds", nNodeIds_,"MpCCI-flownoise");
  params->Get("GlobalDim", GlobalDim_, "MpCCI-flownoise");

  nodeIds_   = new Integer[MpCCInodes_]; 
  nodeIds_[0] = 0;

}

MpCCIexch::~MpCCIexch()
{
if (nodeIds_)  delete [] nodeIds_;
}

void MpCCIexch::PutExchangeGrid2MpCCI(StdVector<std::string> coupledsubdoms)
{
#ifdef TRACE
  (*trace) << "entering MpCCIexch::PutExchangeGrid2MpCCI" << std::endl;
#endif


  // Here starts part for giving mixed mesh to MpCCI
  
  CCI_Def_partition(meshId_, partId_);
 
  actlevel_= 0;

  Matrix<Double> ptCoordNodes;
  StdVector<Integer> connecth;
  Integer i,j,ii;
  Integer elsize = -1;

#define REALTYPE CCI_DOUBLE
  typedef double Realtype;

  std::cout<<"Nodes: "<<MpCCInodes_<<std::endl;
  Realtype * NODEDATA=new Realtype[3*MpCCInodes_];
  int ** TOPOLOGYDATA;
  TOPOLOGYDATA=new int*[coupledsubdoms.GetSize()];
  StdVector<Elem*> elemssd;     

    for (i=0; i<(coupledsubdoms.GetSize()); i++)
      {
      ptgrid_->GetElemSD(elemssd,coupledsubdoms[i],actlevel_);

      elsize=(elemssd[0]->connect).GetSize();
      TOPOLOGYDATA[i]=new int[elsize*elemssd.GetSize()];	  
      int k=0;
      for (j=0; j< elemssd.GetSize(); j++)
	{

	  connecth=elemssd[j]->connect;
	  // 	  ptCoordNodes=new Point<2>[connecth.size()];
	  //	  ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,actlevel_);
	  ptgrid_->GetCoordNodesElemMat(connecth,ptCoordNodes,actlevel_);
	  for (ii=0; ii<elsize; ii++, k++)
	    {  
	      TOPOLOGYDATA[i][k]=connecth[ii];
	      NODEDATA[3*(connecth[ii]-1)]=ptCoordNodes[0][ii];		      
	      NODEDATA[3*(connecth[ii]-1)+1]=ptCoordNodes[1][ii];
	      if(Dim_==3)
		NODEDATA[3*(connecth[ii]-1)+2]= ptCoordNodes[2][ii]; // z-component for the 2d case is zero
	      else
		NODEDATA[3*(connecth[ii]-1)+2]= 0.0; // z-component for the 2d case is zero		      
	    }	
	}


      // 	      std::cout<<"NODEDATA 1d Array:"<<std::endl;
      // 	      for (ii=0; ii<(3*size_); ii++)
      // 		std::cout<<NODEDATA[ii]<<"\t";
      // 	      std::cout<<std::endl;
      // 	      std::cout<<"TOPOLOGYDATA 1d Array:"<<std::endl;
      // 	      for (ii=0; ii<(elsize*maxnumelem_); ii++)
      // 		std::cout<<TOPOLOGYDATA[ii]<<"\t";		    
      // 	      std::cout<<std::endl;

    }

  //define the nodes
  CCI_Def_nodes(meshId_, partId_, GlobalDim_, MpCCInodes_, nNodeIds_, nodeIds_, REALTYPE, NODEDATA);

   for (i=0; i<coupledsubdoms.GetSize(); i++)
     {
      ptgrid_->GetElemSD(elemssd,coupledsubdoms[i],actlevel_);
      int k=0;
      nElemIds_=0;
      Integer nElemSD=elemssd.GetSize();
      elemIds_ = new Integer[nElemSD];
      elemIds_[0] = 0;
      nElemTypes_ = 1;
      nNodesPerElem_ = new Integer[nElemSD];
      elemTypes_ = new Integer[nElemSD];
      elsize=(elemssd[0]->connect).GetSize();// each subdomain must have only one elementtype
      switch(elsize)
	{
	case 2:
	  {
	    nNodesPerElem_[0] = 2;      
	    elemTypes_[0] = CCI_ELEM_LINE;
	    break;
	  }
	case 3:
	  {
	    nNodesPerElem_[0] = 3;      
	    elemTypes_[0] = CCI_ELEM_TRIANGLE;
	    break;
	  }
	case 4:
	  {
	    nNodesPerElem_[0] = 4;  
	    if(Dim_==3)
	      elemTypes_[0] = CCI_ELEM_TETRAHEDRON;
	    else
	      elemTypes_[0] = CCI_ELEM_QUAD;
	    break;
	  }
	case 5:
	  {
	    nNodesPerElem_[0] = 5;      
	    elemTypes_[0] = CCI_ELEM_PYRAMID;
	    break;
	  }
	case 6:
	  {
	    nNodesPerElem_[0] = 6;
	    elemTypes_[0] = CCI_ELEM_PRISM;
	    break;
	  }
	case 8:
	  {
	    nNodesPerElem_[0] = 8;      
	    elemTypes_[0] = CCI_ELEM_HEXAHEDRON;
	    break;
	  }
	}
      
      //define the elements
      CCI_Def_elems(meshId_, partId_, nElemSD, nElemIds_, elemIds_, nElemTypes_, elemTypes_, nNodesPerElem_, TOPOLOGYDATA[i]);
          }
  

  //Close the definition phase; contact detection.
  //i take part in the coupling
  MpCCIprocess_ = 1;
  CCI_Close_setup(MpCCIprocess_);

if (NODEDATA)  delete [] NODEDATA;
if (TOPOLOGYDATA)  delete [] TOPOLOGYDATA; 

if (elemIds_)  delete [] elemIds_;
if ( nNodesPerElem_)  delete [] nNodesPerElem_;
if (elemTypes_)  delete [] elemTypes_;

}

void MpCCIexch::CouplCompPhase(Matrix<Double> & flowdata, Integer timestep)
{
#ifdef TRACE
  (*trace) << "entering MpCCIexch::CouplCompPhase" << std::endl;
#endif

  //Coupling parameters for MpCCI
  Integer nQuantityIds;
  Integer *quantityIds = new Integer[CCI_MAX_NQUANTITIES];
  //    quantityIds[0] = 0;
  Integer nLocalMeshIds;
  Integer *localMeshIds = new Integer[CCI_MAX_NQUANTITIES];
  localMeshIds[0] = 0;
  Integer remoteCodeId;


  Integer quantityId1; // Id for pressure
  Integer quantityId2; // Id for vel
  Integer quantityDim1; //pressure (scalar)
  Integer quantityDim2; //velx, vely, velz (vector)
  Integer maxnEmptyNodes;

  params->Get("nQuantityIds",nQuantityIds, "MpCCI-flownoise");
  params->Get("nLocalMeshIds",nLocalMeshIds, "MpCCI-flownoise");
  params->Get("remoteCodeId",remoteCodeId, "MpCCI-flownoise");
  params->Get("quantityId1",quantityId1, "MpCCI-flownoise");
  params->Get("quantityId2",quantityId2, "MpCCI-flownoise");
  params->Get("quantityDim1",quantityDim1, "MpCCI-flownoise");
  params->Get("quantityDim2",quantityDim2, "MpCCI-flownoise");
  params->Get("maxnEmptyNodes",maxnEmptyNodes, "MpCCI-flownoise");

  CCI_Status status;
  Integer comm = CCI_COMM_RCODE[remoteCodeId];
  quantityIds[0] = quantityId1;
  quantityIds[1] = quantityId2;

  Integer *emptyNodes = new Integer[MpCCInodes_];
  Integer nEmptyNodes;

  Double *value_Press = new Double[MpCCInodes_]; // pressure
  Double *value_VxVy = new Double[MpCCInodes_*3]; //*3 due to velx, vely, velz (for 2D=0)

    Integer globalConvergence = CCI_CONTINUE;
    Integer myConvergence     = CCI_CONTINUE;

    //Retrieve valuedata via MpCCI
    //   while ( globalConvergence != CCI_STOP ) 
    //     {

    CCI_Recv(nQuantityIds, quantityIds, nLocalMeshIds, localMeshIds, comm, &status);

    //Get_nodes is a local operation
    //PRESSURE
    CCI_Get_nodes( meshId_, partId_, quantityId1, quantityDim1, MpCCInodes_, nNodeIds_, nodeIds_,
		   CCI_DOUBLE, value_Press, maxnEmptyNodes, emptyNodes, &nEmptyNodes);
    //VELOCITY
    CCI_Get_nodes( meshId_, partId_, quantityId2, quantityDim2, MpCCInodes_, nNodeIds_, nodeIds_,
		   CCI_DOUBLE, value_VxVy, maxnEmptyNodes, emptyNodes, &nEmptyNodes);

    // Putting values in our matrix flowdata
    Integer k = 0;
    for (Integer inode=0; inode<MpCCInodes_; inode++)
      {
	flowdata[0][inode] = value_Press[inode];
	for(Integer i=0;i<Dim_;i++)
	  flowdata[i+1][inode] = value_VxVy[k+i];
	
	k = k+3;
      }

 
    //check covergence: do another time step or not!
    CCI_Check_convergence(myConvergence,&globalConvergence,CCI_ANY_CODE);

    //     }

if (quantityIds)  delete [] quantityIds;
if (value_Press)  delete [] value_Press;
if (value_VxVy)  delete []  value_VxVy;
    std::cout<<"Leaving CplCompPhase"<<std::endl;
}
      
} // end of namespace
