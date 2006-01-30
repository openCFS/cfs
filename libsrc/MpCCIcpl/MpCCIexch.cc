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

MpCCIexch::MpCCIexch(Grid *aptgrid, StdVector<UInt> & mapSD)
{
  ENTER_FCN("entering MpCCIexch::MpCCIexch");

  //std::cout<<"Nodes SD: "<<nNodesSD<<std::endl;
  ptgrid_ = aptgrid;
  mapSD_ = mapSD;
  MpCCInodes_ = mapSD_.GetSize();

  Dim_ = ptgrid_->GetDim();


  //Get general specific coupling description for CFS++ side
  params->Get("meshId", meshId_, "MpCCI-flownoise");
  params->Get("partId",partId_, "MpCCI-flownoise");
  params->Get("nNodeIds", nNodeIds_,"MpCCI-flownoise");
  params->Get("GlobalDim", GlobalDim_, "MpCCI-flownoise");

  nodeIds_   = new Integer[MpCCInodes_]; 
  nodeIds_[0] = 0;

  //To activate storage of node numbers in an array
  nNodeIds_=MpCCInodes_;

  // For writing the grid definition and the interpolated 
  // source values in files
      writeGridFile_ = FALSE;
      writeSrcFileperTS_ = FALSE;
      writeSrcFileperNode_ = FALSE;
    if( params->IsSet( "writeCoupledGrid","pdeList" ,"acoustic" ) ) {
      writeGridFile_ = TRUE;
      Info->PrintF("acoustic","Writing grid def. file of coupled domain\n");
    }
    if( params->IsSet( "writeSrcFileperTS","pdeList" ,"acoustic" ) ) {
      writeSrcFileperTS_ = TRUE;
      Info->PrintF("acoustic","Writing coarse sources of coupled domain in time files (NrFiles=NrTimeSteps)\n");
    }
    if( params->IsSet( "writeSrcFileperNode","pdeList" ,"acoustic" ) ) {
      writeSrcFileperNode_ = TRUE;
      Info->PrintF("acoustic","Writing coarse transient sources in single nodal files (NrFiles=NrNodes)\n");
    }
    if (writeGridFile_)
      {
        std::string filename;
        filename = "elemfile";
        filename.append( ".elem" );
        outelemfile_ = new std::ofstream(filename.c_str());
        filename = "nodefile";
        filename.append( ".node" );
        outnodefile_ = new std::ofstream(filename.c_str());
      }

}
  //constructer used for FSI
  MpCCIexch::MpCCIexch(Grid *aptgrid, StdVector<RegionIdType>  subdoms)
  {
    ENTER_FCN( "entering MpCCIexch::MpCCIexch" );

    ptgrid_ = aptgrid;
    Dim_ = ptgrid_->GetDim();
    subdoms_ = subdoms;
    
    TOPOLOGYDATA_=new int*[subdoms_.GetSize()];
    NODEDATA_ =new double*[subdoms_.GetSize()];
    //Get general specific coupling description for CFS++ side
    params->Get("GlobalDim", GlobalDim_, "mpcciFSI");
  }

MpCCIexch::~MpCCIexch()
{
  if (nodeIds_)  delete [] nodeIds_;
 if (outelemfile_)
   delete outelemfile_;
 if (outnodefile_)
   delete outnodefile_;
 if (outsrcfile_) 
      delete outsrcfile_;
}

void MpCCIexch::PutExchangeGrid2MpCCI(StdVector<RegionIdType> coupledsubdoms)
{
    ENTER_FCN("entering MpCCIexch::PutExchangeGrid2MpCCI");

    // Here starts part for giving mixed mesh to MpCCI
  
  CCI_Def_partition(meshId_, partId_);

  Matrix<Double> ptCoordNodes;
  StdVector<UInt> connecth;
  UInt i,j,ii;
  UInt elsize = 0;

#define REALTYPE CCI_DOUBLE
  typedef double Realtype;

  std::cout<<"Nodes: "<<MpCCInodes_<<std::endl;
  Realtype * NODEDATA=new Realtype[3*MpCCInodes_];
  int ** TOPOLOGYDATA;
  TOPOLOGYDATA=new int*[coupledsubdoms.GetSize()];
  StdVector<Elem*> elemssd;     

  Point<2> ptNodalCoord2D;      
  Point<3> ptNodalCoord3D;      

  if (writeGridFile_)
    {
      (*outnodefile_) <<MpCCInodes_<<std::endl;
    }

  for (i=0; i<(coupledsubdoms.GetSize()); i++)
      {
        // This is part of workaround for 3D quadratic elements
	for (UInt n=0; n<(mapSD_.GetSize()); n++)
	  {
	    nodeIds_[n] = (int)mapSD_[n];
//   	    std::cout<<"nodeIds["<<n<<"]= "<<nodeIds_[n]<<std::endl;
//   	    std::cout<<"mapSD_.GetSize: "<<mapSD_.GetSize()<<std::endl;
            if (Dim_==2)
              {
                ptgrid_->GetNodeCoordinate(ptNodalCoord2D, nodeIds_[n]);
                NODEDATA[3*n]=ptNodalCoord2D[0];		      
                NODEDATA[3*n+1]=ptNodalCoord2D[1];  
                NODEDATA[3*n+2]= 0.0; // z-component for the 2d case is zero
                if (writeGridFile_)
                {
                  (*outnodefile_) << NODEDATA[3*n]<<" "<< NODEDATA[3*n+1]<<std::endl;
                }	      
              }
            else
              {
                ptgrid_->GetNodeCoordinate(ptNodalCoord3D, nodeIds_[n]);
                NODEDATA[3*n]=ptNodalCoord3D[0];		      
                NODEDATA[3*n+1]=ptNodalCoord3D[1];
                NODEDATA[3*n+2]=ptNodalCoord3D[2];
              }
	  }
    
        ptgrid_->GetElems(elemssd,coupledsubdoms[i]);
        //ptgrid_->GetVolElems(elemssd,coupledsubdoms[i]);
        //ptgrid_->GetElemSD(elemssd,coupledsubdoms[i],actlevel_);
      
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

      TOPOLOGYDATA[i]=new int[elsize*elemssd.GetSize()];	  

      if (writeGridFile_)
        {
          (*outelemfile_) << elemssd.GetSize()<<" "<<elsize<< std::endl;
        }

      int k=0;
      for (j=0; j< elemssd.GetSize(); j++)
	{
	  connecth=elemssd[j]->connect;
	  //	  ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,actlevel_);
	  // ptgrid_->GetCoordNodesElemMat(connecth,ptCoordNodes,actlevel_);
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
              
              if (writeGridFile_ && (ii<(elsize-1)))
                {
                  (*outelemfile_) << TOPOLOGYDATA[i][k]<<" ";
                }	      
              else
                {
                  if (writeGridFile_)
                  (*outelemfile_) << TOPOLOGYDATA[i][k];
                }
            }
          if (writeGridFile_)
               {
                 (*outelemfile_) << std::endl;
               }	     
        } // loop over all elements
      }

  //define the nodes
          CCI_Def_nodes((int)meshId_, (int)partId_, (int)GlobalDim_, MpCCInodes_, (int)nNodeIds_, &nodeIds_[0], REALTYPE, NODEDATA);

   for (i=0; i<coupledsubdoms.GetSize(); i++)
     {
       ptgrid_->GetElems(elemssd,coupledsubdoms[i]);
       //ptgrid_->GetVolElems(elemssd,coupledsubdoms[i]);
      nElemIds_=0;
      UInt nElemSD=elemssd.GetSize();
      elemIds_ = new Integer[nElemSD];
      elemIds_[0] = 0;
      nElemTypes_ = 1;
      nNodesPerElem_ = new Integer[nElemSD];
      elemTypes_ = new Integer[nElemSD];
      // each subdomain must have only one elementtype
      elsize=(elemssd[0]->connect).GetSize();
      //workaround for computing with quadratic hexas and tetras
      if  ((elsize==20)&&(Dim_==3))
	{
	  elsize=8;
	}
      else if ((elsize==10)&&(Dim_==3))
	{
	  elsize=4;
	}

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
	    if (Dim_==3)
	      elemTypes_[0] = CCI_ELEM_PRISM;
	    else
	      elemTypes_[0] = CCI_ELEM_TRIANGLE6;
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
      CCI_Def_elems((int)meshId_, (int)partId_, (int)nElemSD, (int)nElemIds_, elemIds_, 
		    (int)nElemTypes_, elemTypes_, (int*)nNodesPerElem_, TOPOLOGYDATA[i]);
          }
  

  //Close the definition phase; contact detection.
  //i take part in the coupling
  MpCCIprocess_ = 1;
    SETPROFILE("Before CCI_Close_setup()");
  CCI_Close_setup(MpCCIprocess_);
    SETPROFILE("After CCI_Close_setup()");

if (NODEDATA)  delete [] NODEDATA;
if (TOPOLOGYDATA)  delete [] TOPOLOGYDATA; 

if (elemIds_)  delete [] elemIds_;
if ( nNodesPerElem_)  delete [] nNodesPerElem_;
if (elemTypes_)  delete [] elemTypes_;

//  if (outelemfile_)
//    delete [] outelemfile_;
//  if (outnodefile_)
//    delete [] outnodefile_;
}
    
void MpCCIexch::DefMpcciPartition(UInt meshId, UInt partId)
{
  ENTER_FCN( "MpCCIexch::DefMpcciPartition" );
  CCI_Def_partition(meshId, partId);
}

void MpCCIexch::DefMpcciNodes(UInt meshId, UInt partId, UInt nrNodesSD, 
                              UInt* nodeIds,  NodeEQN & eqnData)
{
  ENTER_FCN( "MpCCIexch::DefMpcciNodes" );

  UInt i;
  UInt globalNode;

  Point<3> ptPoint;

#define REALTYPE CCI_DOUBLE


  NODEDATA_[partId-1]=new double[3*nrNodesSD];

  for (i=0; i<nrNodesSD ; i++)
    {
      globalNode = eqnData.PDE2MeshNode(nodeIds[i]);
      ptgrid_->GetNodeCoordinate(ptPoint, globalNode);
      
      NODEDATA_[partId-1][3*i]  =ptPoint[0];		      
      NODEDATA_[partId-1][3*i+1]=ptPoint[1];
      
      if(Dim_==3)
	NODEDATA_[partId-1][3*i+2]= ptPoint[2]; // z-component for the 2d case is zero
      else
	NODEDATA_[partId-1][3*i+2]= 0.0; // z-component for the 2d case is zero	

//       Info->PrintF("","local:%d;\t global:%d; \t%f\t%f\t%f\n",nodeIds[i],eqnData.PDE2MeshNode(nodeIds[i]), 
// 		   NODEDATA_[partId-1][3*i], NODEDATA_[partId-1][3*i+1], NODEDATA_[partId-1][3*i+2]);

    }

  //define the nodes
  CCI_Def_nodes( (int) meshId, (int) partId, (int) GlobalDim_, (int) nrNodesSD, 
                 (int) nrNodesSD, (int*) nodeIds, REALTYPE, NODEDATA_[partId-1]);
}

void MpCCIexch::DefMpcciElements(UInt meshId, UInt partId, NodeEQN & eqnData)
{
  UInt i,j;
  UInt NrOfElemsInSD, NodesPerElem;
  UInt localPDENode;
  StdVector<Elem*> elemssd;     
  StdVector<UInt> connecth;

  ptgrid_->GetElems(elemssd,subdoms_[partId-1]);
//   ptgrid_->GetVolElems(elemssd,subdoms_[partId-1]);

  NrOfElemsInSD = elemssd.GetSize(); // evtl check: if (NrOfElemsInSD==-1) dann exit

  NodesPerElem  =(elemssd[0]->connect).GetSize(); // evtl check: if (NrOfElemsInSD==-1) dann exit

  TOPOLOGYDATA_[partId-1]=new int[NrOfElemsInSD*NodesPerElem];	  

  for (i=0; i<NrOfElemsInSD ; i++) // loop over all elements belonging to subdoms_[partId-1]
    {
      //get all element nodes
      connecth=elemssd[i]->connect;
//       Info->PrintF("","local elem:%d\t",i);
      for (j=0; j<NodesPerElem; j++)
	{
	  localPDENode = eqnData.Mesh2PDENode(connecth[j]);
	  TOPOLOGYDATA_[partId-1][NodesPerElem*i+j]=localPDENode;

// 	  Info->PrintF("","Top[%d]=%d, %d\t",NodesPerElem*i+j, TOPOLOGYDATA_[partId-1][NodesPerElem*i+j], eqnData.PDE2MeshNode(localPDENode));
	}  
//       Info->PrintF("","\n");
    }

  nElemIds_=0;
  nElemTypes_ = 1;
  nNodesPerElem_ = new Integer[NrOfElemsInSD];
  elemTypes_ = new Integer[NrOfElemsInSD];
  elemIds_ = new Integer[NrOfElemsInSD];
  elemIds_[0] = 0;
      
  switch(NodesPerElem)
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
	  {
	    elemTypes_[0] = CCI_ELEM_QUAD;
	    //elemTypes_[0] = CCI_ELEM_TETRAHEDRON;
	  }
	else
	  {
	    elemTypes_[0] = CCI_ELEM_QUAD;
	  }
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
  CCI_Def_elems(meshId, partId, NrOfElemsInSD, nElemIds_, elemIds_, nElemTypes_, elemTypes_, nNodesPerElem_, TOPOLOGYDATA_[partId-1]);
}

void MpCCIexch::FinishMpcciSetup(std::string couplingType)
{
  Integer *quantityIds = new Integer[CCI_MAX_NQUANTITIES];
  Integer *localMeshIds = new Integer[CCI_MAX_NQUANTITIES];
  
  if(couplingType=="shell")
    {
      quantityIds[0] = 31;
      quantityIds[1] = 33;
      quantityIds[2] = 31;
      quantityIds[3] = 33;
  
      localMeshIds[0] = 102;
      localMeshIds[1] = 102;
      localMeshIds[2] = 103;
      localMeshIds[3] = 103;
      
      CCI_Def_comm(211,1,111,4,quantityIds,4,localMeshIds);
      
      quantityIds[0] = 15;
      quantityIds[1] = 16;
      
      localMeshIds[0] = 102;
      
      CCI_Def_comm(212,1,112,2,quantityIds,1,localMeshIds);
    }
  else if(couplingType=="solid")
    {
      quantityIds[0] = 31; //pressure
      quantityIds[1] = 33; //tau
  
      localMeshIds[0] = 102;

      //Define communicator for transfer of data from Fastest to this code
      CCI_Def_comm(211,1,111,2,quantityIds,1,localMeshIds);
      
      quantityIds[0] = 15; // displacement
      quantityIds[1] = 16; // iter
      
      localMeshIds[0] = 102;

      //Define communicator for transfer of data from this code to Fastest
      CCI_Def_comm(212,1,112,2,quantityIds,1,localMeshIds);
    }
  else
    {
      std::string errmsg = "MpCCIType '" + couplingType;
      errmsg += "does not exist";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

  //Close the definition phase; contact detection.
  //i take part in the coupling
  MpCCIprocess_ = 1;
  CCI_Close_setup(MpCCIprocess_);

  if (NODEDATA_) delete [] NODEDATA_;
  if (TOPOLOGYDATA_) delete [] TOPOLOGYDATA_; 

  if (elemIds_) delete [] elemIds_;
  if ( nNodesPerElem_) delete [] nNodesPerElem_;
  if (elemTypes_) delete [] elemTypes_;
}

void MpCCIexch::CouplCompPhase(Matrix<Double> & flowdata, Double acttime)
{
  ENTER_FCN("entering MpCCIexch::CouplCompPhase");


  //Coupling parameters for MpCCI
  Integer nQuantityIds;
  Integer *quantityIds = new Integer[CCI_MAX_NQUANTITIES];
  //    quantityIds[0] = 0;
  UInt nLocalMeshIds;
  Integer *localMeshIds = new Integer[CCI_MAX_NQUANTITIES];
  localMeshIds[0] = 0;
  UInt remoteCodeId;


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
//     std::cout<<"nodeIds size: "<<nodeIds_[MpCCInodes_-1]<<std::endl;
//     std::cout<<"nNodeIds= "<<nNodeIds_<<std::endl;
    //Get_nodes is a local operation
    //PRESSURE
    CCI_Get_nodes((int)meshId_, (int)partId_, quantityId1, quantityDim1, (int)MpCCInodes_, nNodeIds_,  &nodeIds_[0],
		   CCI_DOUBLE, value_Press, maxnEmptyNodes, emptyNodes, &nEmptyNodes);
    //VELOCITY
    CCI_Get_nodes((int) meshId_, (int)partId_, quantityId2, quantityDim2, MpCCInodes_, (int)nNodeIds_,  &nodeIds_[0],
		   CCI_DOUBLE, value_VxVy, maxnEmptyNodes, emptyNodes, &nEmptyNodes);

    // Putting values in our matrix flowdata
    Integer k = 0;


    Boolean nodalSrc=FALSE;
    //check type of flow data
    if( params->HasValue( "type", "nodalSrc", "acoustic", "flowData" ) ) {
      nodalSrc = TRUE;
      Info->PrintF("acoustic", "In MpCCIexch Using FlowData as RHS nodal source\n" );
    }
    

    std::cout<<"flowdata length= "<<flowdata.GetSizeCol()<<std::endl;
    if (nodalSrc == TRUE)
      {
        if (writeSrcFileperTS_)
          {
            // static variable for suffix of output src file
            static Integer filenum=1;
            std::string filename;
            filename = "srcfile";
            filename.append( ".ts" );
            if ( filenum < 10 ) filename.append( "000" );
            else if ( filenum < 100 ) filename.append( "00" );
            else if ( filenum < 1000 ) filename.append( "0" );
            else if ( filenum > 10000 ) {
              Info->Error( "Number of src file exceeds 9999!",
                           __FILE__, __LINE__ );
            }
            filename.append( Info->GenStr( filenum ) );
            filenum++;
            outsrcfile_ = new std::ofstream(filename.c_str());
          }
        
        static Integer TimeStepCtr=1;
        for (UInt inode=0; inode<MpCCInodes_; inode++)
          {
            // Getting first column (no node identifier)
            // flowdata[0][inode] = value_Press[inode];
            
            // Getting first column as INT(dTdxdw) with node identifier
            flowdata[0][nodeIds_[inode]-1] = value_Press[inode];

            // std::cout<<"flowdata[0]["<<nodeIds_[inode]-1<<"]= "<<value_Press[inode]<<std::endl;
            // std::cout<<"value_Press["<<inode<<"]= "<<value_Press[inode]<<std::endl;
            if (writeSrcFileperTS_)
              {
                (*outsrcfile_) <<  value_Press[inode]<<std::endl;
                // !!!!THIS IS ONLY TO ALLOW Flemming (Stuttgart) test his reading of the data
                // we assign the nodal index as output to the nodal value just for checking!!!
                //(*outsrcfile_) << nodeIds_[inode] <<std::endl;
              }
            if (writeSrcFileperNode_)
              {
                std::string filename;
                std::ofstream outsrcnodalfile;
                filename = "nodalSrcs/timesrcfile";
                filename.append( ".node" );
                filename.append( Info->GenStr( nodeIds_[inode] ) );
                //create the file if it doesn't exist yet
                
                if (TimeStepCtr==1)
                  {
                    // Create directory to save in the source files
                    std::string S="mkdir -p nodalSrcs";
                    system(S.c_str());
                    std::cout<<"In first time step"<<std::endl;
                    outsrcnodalfile.open(filename.c_str(), std::ios::out | std::ios::trunc);
                  }
                else
                  {
                    outsrcnodalfile.open(filename.c_str(), std::ios::app);
                  }

                if (!outsrcnodalfile) 
                  {
                    std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
                      ") In MpCCIexch: Can't open src nodal file for output:" << filename.c_str() << std::endl;
                    exit(1);
                  }

                outsrcnodalfile<< std::setiosflags(std::ios::uppercase | std::ios::scientific)
                               << acttime << " " << value_Press[inode] << std::endl;
                // !!!!THIS IS ONLY TO test reading of the data
                // we assign the nodal index as output to the nodal value just for checking!!!
                //outsrcnodalfile << nodeIds_[inode] <<std::endl;

                outsrcnodalfile.close();
              }
          }
    TimeStepCtr++;
      }
    else
      {
        for (UInt inode=0; inode<MpCCInodes_; inode++)
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

// if (outsrcfile_) 
//      delete outsrcfile_;


    std::cout<<"Leaving CplCompPhase"<<std::endl;
}

void MpCCIexch::RecvAllPartitions(std::string couplingType)
{
  ENTER_FCN( "MpCCIexch::RecvAllPartitions" );

  //Coupling parameters for MpCCI
  Integer nQuantityIds;
  Integer *quantityIds = new Integer[CCI_MAX_NQUANTITIES];
  
  Integer nLocalMeshIds;
  Integer *localMeshIds = new Integer[CCI_MAX_NQUANTITIES];

  Integer remoteCodeId;

  // courrently the FSI is done on two sides. Side 'a' and side 'b'.
  // For each side the force due to fluid pressure (p) and due to 
  // molecular momentum exchange (tau) has to be received.
  Integer quantityId1; // Id for force from p 'a'
  Integer quantityId2; // Id for force from tau 'a'

  Integer quantityDim1; // fx, fy, fz (p_a   force vector)
  Integer quantityDim2; // fx, fy, fz (tau_a   force vector)

  CCI_Status status;
  //Integer comm = CCI_COMM_RCODE[remoteCodeId];


  if(couplingType=="shell")
    {
      nQuantityIds = 4;
      nLocalMeshIds = 4;

      quantityIds[0] = 31;
      quantityIds[1] = 33;
      quantityIds[2] = 31;
      quantityIds[3] = 33;

      //  Integer *localMeshIds = new Integer[CCI_MAX_NQUANTITIES];
      localMeshIds[0] = 102;
      localMeshIds[1] = 102;
      localMeshIds[2] = 103;
      localMeshIds[3] = 103;
      
      CCI_Recv(nQuantityIds, quantityIds, nLocalMeshIds, localMeshIds, 211, &status);
    }
  else if(couplingType=="solid")
    {
      nQuantityIds = 2;
      nLocalMeshIds = 1;

      quantityIds[0] = 31;
      quantityIds[1] = 33;

      //  Integer *localMeshIds = new Integer[CCI_MAX_NQUANTITIES];
      localMeshIds[0] = 102;
      
      CCI_Recv(nQuantityIds, quantityIds, nLocalMeshIds, localMeshIds, 211, &status);
    }
  else
    {
      std::string errmsg = "MpCCIType '" + couplingType;
      errmsg += "does not exist";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }
  
}

void MpCCIexch::GetNodalValOfOnePartition(UInt partId, Vector<Double> & forceData, UInt nrNodesSD, 
                                          UInt* nodeIds, std::string couplingType)
{
  ENTER_FCN( "MpCCIexch::GetNodalValOfOnePartition" );

  Integer quantityId1; // Id for force from p 'a'
  Integer quantityId2; // Id for force from tau 'a'
  Integer quantityId3; // Id for force from p 'b'
  Integer quantityId4; // Id for force from tau 'b'

  Integer quantityDim1; // fx, fy, fz (p_a   force vector)
  Integer quantityDim2; // fx, fy, fz (tau_a   force vector)
  Integer quantityDim3; // fx, fy, fz (p_b   force vector)
  Integer quantityDim4; // fx, fy, fz (tau_b   force vector)
  Integer maxnEmptyNodes;

  Integer *emptyNodes = new Integer[nrNodesSD];
  Integer nEmptyNodes;

  if(couplingType=="shell")
    {
      Double *value_force_pa   = new Double[nrNodesSD*3]; //*3 due to fx, fy, fz (In 2D 3rd comp: fz=0)
      Double *value_force_taua = new Double[nrNodesSD*3]; //    -- || --
      Double *value_force_pb   = new Double[nrNodesSD*3]; //    -- || --
      Double *value_force_taub = new Double[nrNodesSD*3]; //    -- || --

      quantityId1   = 31;
      quantityId2   = 33;
      quantityId3   = 31;
      quantityId4   = 33;

      quantityDim1  = 3;
      quantityDim2  = 3;
      quantityDim3  = 3;
      quantityDim4  = 3;

      maxnEmptyNodes= 1;

      CCI_Get_nodes( 102, (int) partId, quantityId1, quantityDim1, (int) nrNodesSD, (int) nrNodesSD, 
                     (int*) nodeIds, CCI_DOUBLE, value_force_pa, maxnEmptyNodes, emptyNodes, &nEmptyNodes);

      CCI_Get_nodes( 102, (int) partId, quantityId2, quantityDim2, (int) nrNodesSD, (int) nrNodesSD, 
                     (int*) nodeIds, CCI_DOUBLE, value_force_taua, maxnEmptyNodes, emptyNodes, &nEmptyNodes);

      CCI_Get_nodes( 103, (int) partId, quantityId3, quantityDim3, (int) nrNodesSD, (int) nrNodesSD, 
                     (int*) nodeIds, CCI_DOUBLE, value_force_pb, maxnEmptyNodes, emptyNodes, &nEmptyNodes);

      CCI_Get_nodes( 103, (int) partId, quantityId4, quantityDim4, (int) nrNodesSD, (int) nrNodesSD, 
                     (int*) nodeIds,CCI_DOUBLE, value_force_taub, maxnEmptyNodes, emptyNodes, &nEmptyNodes);

      for (UInt k=0; k<3*nrNodesSD; k++)
	{
	  forceData[k] = value_force_pa[k] + value_force_taua[k] + value_force_pb[k] + value_force_taub[k];
	}
      if (value_force_pa)    delete []  value_force_pa;
      if (value_force_taua)  delete []  value_force_taua;
      if (value_force_pb)    delete []  value_force_pb;
      if (value_force_taub)  delete []  value_force_taub;
    }


  else if(couplingType=="solid")
    {
      Double *value_force_pa   = new Double[nrNodesSD*3]; //*3 due to fx, fy, fz (In 2D 3rd comp: fz=0)
      Double *value_force_taua = new Double[nrNodesSD*3]; //    -- || --

      quantityId1   = 31;
      quantityId2   = 33;

      quantityDim1  = 3;
      quantityDim2  = 3;
      maxnEmptyNodes= 1;

      CCI_Get_nodes( 102, (int) partId, quantityId1, quantityDim1, (int) nrNodesSD, (int) nrNodesSD, 
                     (int*) nodeIds, CCI_DOUBLE, value_force_pa, maxnEmptyNodes, emptyNodes, &nEmptyNodes);

      CCI_Get_nodes( 102, (int) partId, quantityId2, quantityDim2, (int) nrNodesSD, (int) nrNodesSD, 
                     (int*) nodeIds, CCI_DOUBLE, value_force_taua, maxnEmptyNodes, emptyNodes, &nEmptyNodes);

      for (UInt k=0; k<3*nrNodesSD; k++)
	{
	  forceData[k] = value_force_pa[k] + value_force_taua[k];
	}
      if (value_force_pa)    delete []  value_force_pa;
      if (value_force_taua)  delete []  value_force_taua;

    }
  else
    {
      std::string errmsg = "MpCCIType '" + couplingType;
      errmsg += "does not exist";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

  if (emptyNodes)  delete []  emptyNodes;
}

void MpCCIexch::PutPartition(UInt partId, const Vector<Double> & displData, UInt nrNodesSD, 
                             UInt* nodeIds, Boolean conv)
{
  ENTER_FCN("entering MpCCIexch::PutPartition");

  //Coupling parameters for MpCCI
  Integer nQuantityIds;
  Integer *quantityIds = new Integer[CCI_MAX_NQUANTITIES];
  //    quantityIds[0] = 0;
  Integer nLocalMeshIds;
  Integer *localMeshIds = new Integer[CCI_MAX_NQUANTITIES];
  localMeshIds[0] = 0;
  Integer remoteCodeId;

  Integer quantityId1; // Id for displacement dx, dy, dz
  Integer quantityId2; // Id for convergence check; this should be done cheaper!!!

  Integer quantityDim1; // for coupled quantity (vector)
  Integer quantityDim2; // for coupled quantity (vector)

  Integer maxnEmptyNodes;

//   nodeIds_   = new int[nrNodesSD]; 
//   nodeIds_[0] = 0;

  Integer globalConvergence = CCI_CONTINUE;
  Integer myConvergence     = CCI_CONTINUE;

  double *value_disp = new double[nrNodesSD*3]; //*3 due three dims (for 2D, 3rdcomp.=0)
  double *value_conv = new double[nrNodesSD]; 

  //Retrieve valuedata via MpCCI
  //   while ( globalConvergence != CCI_STOP ) 
  //     {

  quantityId1 = 15; //should be read out of name.cci file
  quantityId2 = 16; //should be read out of name.cci file

  quantityDim1 = 3; 
  quantityDim2 = 1; 

  for (UInt k=0; k<3*nrNodesSD; k++)
    {
      value_disp[k]   = displData[k];
    }

  if (conv == TRUE) 
    {
      for (UInt k=0; k<nrNodesSD; k++)
	{
	  value_conv[k]    = 0.0; // the mechanic PDE converged
	}
    }
  else 
    {
      for (UInt k=0; k<nrNodesSD; k++)
	{
	  value_conv[k]    = 1.0; // the mechanic PDE is not converged
	}
    }

  CCI_Put_nodes( 102, (int) partId, quantityId1, quantityDim1, (int) nrNodesSD, (int) nrNodesSD, 
                 (int*) nodeIds, CCI_DOUBLE, value_disp);

  CCI_Put_nodes( 102, (int) partId, quantityId2, quantityDim2, (int) nrNodesSD, (int) nrNodesSD, 
                 (int*) nodeIds, CCI_DOUBLE, value_conv);

  std::cout<<"Leaving PutPartition"<<std::endl;

  if (value_disp)  delete [] value_disp;
  if (value_conv)  delete [] value_conv;

}

void MpCCIexch::SendAllPartitions()
{
  ENTER_FCN("entering MpCCIexch::SendAllPartitions");

  //Coupling parameters for MpCCI
  Integer nQuantityIds;
  Integer *quantityIds = new Integer[CCI_MAX_NQUANTITIES];
  //    quantityIds[0] = 0;
  Integer nLocalMeshIds;
  Integer *localMeshIds = new Integer[CCI_MAX_NQUANTITIES];
  localMeshIds[0] = 0;
  Integer remoteCodeId;


  Integer globalConvergence = CCI_CONTINUE;
  Integer myConvergence     = CCI_CONTINUE;

  nQuantityIds=2;
  quantityIds[0] = 15;
  quantityIds[1] = 16;
  nLocalMeshIds=1;
  localMeshIds[0] = 102;
  remoteCodeId = 1;
      
  Integer comm = CCI_COMM_RCODE[remoteCodeId];

  CCI_Send(nQuantityIds, quantityIds, nLocalMeshIds, localMeshIds, 212);

  CCI_Check_convergence(myConvergence,&globalConvergence,CCI_ANY_CODE);

  std::cout<<"CCI_CONTINUE_VALUE"<<CCI_CONTINUE<<std::endl;
  std::cout<<"Leaving SendAllPartitions"<<std::endl;
}      
} // end of namespace
