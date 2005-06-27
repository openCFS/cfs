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
    ENTER_FCN("entering MpCCIexch::MpCCIexch");

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

  void MpCCIexch::PutExchangeGrid2MpCCI(StdVector<RegionIdType> coupledsubdoms)
  {
    ENTER_FCN("entering MpCCIexch::PutExchangeGrid2MpCCI");

    // Here starts part for giving mixed mesh to MpCCI
  
    CCI_Def_partition(meshId_, partId_);
 
    Matrix<Double> ptCoordNodes;
    StdVector<UInt> connecth;
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
	ptgrid_->GetVolElems(elemssd,coupledsubdoms[i]);
	elsize=(elemssd[0]->connect).GetSize();
	//workaround for computing with quadratic hexas 
	if  ((elsize==20)&&(Dim_==3))
	  {
	    elsize=8;
	  }
	else if ((elsize==10)&&(Dim_==3))
	  {
	    elsize=4;
	  }
	else
	  elsize=(elemssd[0]->connect).GetSize();
      
	TOPOLOGYDATA[i]=new int[elsize*elemssd.GetSize()];	  



	int k=0;
	for (j=0; j< elemssd.GetSize(); j++)
	  {
	    connecth=elemssd[j]->connect;
	    // 	  ptCoordNodes=new Point<2>[connecth.size()];
	    //	  ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,actlevel_);
	    ptgrid_->GetElemNodesCoord(ptCoordNodes,connecth);
	    for (ii=0; ii<elsize; ii++, k++)
	      {  
		if ((elsize==8)&&(Dim_==2))
		  {
		    // A bit hard coded node order transformation
		    if ((ii==0)||(ii==7))
		      TOPOLOGYDATA[i][k]=connecth[ii];
		    else
		      if ((ii==2)||(ii==4)||(ii==6))
			TOPOLOGYDATA[i][k]=connecth[ii/2];
		      else
			{
			  switch(ii) 
			    {
			    case 1:
			      {
				TOPOLOGYDATA[i][k]=connecth[ii+3];
				break;
			      }
			    case 3:
			      {
				TOPOLOGYDATA[i][k]=connecth[ii+2];
				break;
			      }
			    case 5:
			      {
				TOPOLOGYDATA[i][k]=connecth[ii+1];
				break;
			      }
			    }
			}
		  }
		else
		  TOPOLOGYDATA[i][k]=connecth[ii];
	      
		NODEDATA[3*(connecth[ii]-1)]=ptCoordNodes[0][ii];		      
		NODEDATA[3*(connecth[ii]-1)+1]=ptCoordNodes[1][ii];
		if(Dim_==3)
		  NODEDATA[3*(connecth[ii]-1)+2]= ptCoordNodes[2][ii];
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
	ptgrid_->GetVolElems(elemssd,coupledsubdoms[i]);
	int k=0;
	nElemIds_=0;
	Integer nElemSD=elemssd.GetSize();
	elemIds_ = new Integer[nElemSD];
	elemIds_[0] = 0;
	nElemTypes_ = 1;
	nNodesPerElem_ = new Integer[nElemSD];
	elemTypes_ = new Integer[nElemSD];
	// each subdomain must have only one elementtype
	//workaround for computing with quadratic hexas 
	elsize=(elemssd[0]->connect).GetSize();
	if  ((elsize==20)&&(Dim_==3))
	  {
	    elsize=8;
	  }
	else if ((elsize==10)&&(Dim_==3))
	  {
	    elsize=4;
	  }
	else
	  elsize=(elemssd[0]->connect).GetSize();
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
	      if (Dim_==3)
		elemTypes_[0] = CCI_ELEM_HEXAHEDRON;
	      else
		elemTypes_[0] = CCI_ELEM_QUAD8;
	      break;
	      
	    }
	  }
      
	//define the elements
	CCI_Def_elems(meshId_, partId_, nElemSD, nElemIds_, elemIds_, 
		      nElemTypes_, elemTypes_, nNodesPerElem_, TOPOLOGYDATA[i]);
	// =======
	//         elsize=(elemssd[0]->connect).GetSize();
	//         TOPOLOGYDATA[i]=new int[elsize*elemssd.GetSize()];          
	//         int k=0;
	//         for (j=0; j< elemssd.GetSize(); j++)
	//           {
	//             connecth=elemssd[j]->connect;
	//             //      ptCoordNodes=new Point<2>[connecth.size()];
	//             //      ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,actlevel_);
	//             ptgrid_->GetCoordNodesElemMat(connecth,ptCoordNodes,actlevel_);
	//             for (ii=0; ii<elsize; ii++, k++)
	//               {  
	//                 if ((elsize==8)&&(Dim_==2))
	//                   {
	//                     // A bit hard coded node order transformation
	//                     if ((ii==0)||(ii==7))
	//                       TOPOLOGYDATA[i][k]=connecth[ii];
	//                     else
	//                       if ((ii==2)||(ii==4)||(ii==6))
	//                         TOPOLOGYDATA[i][k]=connecth[ii/2];
	//                       else
	//                         {
	//                           switch(ii) 
	//                             {
	//                             case 1:
	//                               {
	//                                 TOPOLOGYDATA[i][k]=connecth[ii+3];
	//                                 break;
	//                               }
	//                             case 3:
	//                               {
	//                                 TOPOLOGYDATA[i][k]=connecth[ii+2];
	//                                 break;
	//                               }
	//                             case 5:
	//                               {
	//                                 TOPOLOGYDATA[i][k]=connecth[ii+1];
	//                                 break;
	//                               }
	//                             }
	//                         }
	//                   }
	//                 else
	//                   TOPOLOGYDATA[i][k]=connecth[ii];
              
	//                 NODEDATA[3*(connecth[ii]-1)]=ptCoordNodes[0][ii];               
	//                 NODEDATA[3*(connecth[ii]-1)+1]=ptCoordNodes[1][ii];
	//                 if(Dim_==3)
	//                   NODEDATA[3*(connecth[ii]-1)+2]= ptCoordNodes[2][ii];
	//                 else
	//                   NODEDATA[3*(connecth[ii]-1)+2]= 0.0; // z-component for the 2d case is zero                   
	//               }   
	// >>>>>>> 1.11

    //              std::cout<<"NODEDATA 1d Array:"<<std::endl;
    //              for (ii=0; ii<(3*size_); ii++)
    //                std::cout<<NODEDATA[ii]<<"\t";
    //              std::cout<<std::endl;
    //              std::cout<<"TOPOLOGYDATA 1d Array:"<<std::endl;
    //              for (ii=0; ii<(elsize*maxnumelem_); ii++)
    //                std::cout<<TOPOLOGYDATA[ii]<<"\t";                  
    //              std::cout<<std::endl;

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
    ENTER_FCN("entering MpCCIexch::CouplCompPhase");

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


    Boolean nodalSrc;
    //check type of flow data
    if( params->HasValue( "type", "nodalSrc", "acoustic", "flowData" ) ) {
      nodalSrc = TRUE;
      Info->PrintF("acoustic", "In MpCCIexch Using FlowData as RHS nodal source\n" );
    }
    



    if (nodalSrc == TRUE)
      for (Integer inode=0; inode<MpCCInodes_; inode++)
	{
	  flowdata[0][inode] = value_Press[inode]; // Getting first column as INT(dTdxdw)
	}
    else
      {
	for (Integer inode=0; inode<MpCCInodes_; inode++)
	  {
	    flowdata[0][inode] = value_Press[inode];
	    for(Integer i=0;i<Dim_;i++)
	      flowdata[i+1][inode] = value_VxVy[k+i];
              
	    k = k+3;
	  }
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
