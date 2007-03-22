// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "mpcciPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/solveStepMpCCI.hh"
#include "Forms/forms_header.hh"
#include "CoupledPDE/pdecoupling.hh"

#ifdef MpCCI
#include <MpCCIcpl/MpCCIexch.hh>
#endif

//#include "nodeEQN.hh"


namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  MpcciPDE::MpcciPDE( Grid * aptgrid, ParamNode* paramNode )
    :SinglePDE( aptgrid, paramNode ) 
  {

    ENTER_FCN( "MpcciPDE::MpcciPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================
    //dofspernode_ = 1;
    //solTypes_ = MPCCI;
    //solDofs_ = 1;
    pdename_          = "mpcci";
    pdematerialclass_ = MECHANIC;
    flagFirstTimeStep_= true;
    converged_=false;

    // get firs region node
    StdVector<ParamNode*> regionNodes = 
      myParam_->Get("regionList")->GetList("region");
    if( regionNodes.GetSize() == 0 ) {
      EXCEPTION( "No region defined in MpcciPDE!" );
    }
    // take 'type' parameter of first region node
    regionNodes[0]->Get( "type", MpCCIType_ );

    // get node 'mpcciFSI'
    ParamNode * mpcciNode = param->Get( "mpcciFSI" );
    
    
    if (MpCCIType_ == "shell")
      {
	NumMeshIds_=2;
	meshId_= new UInt[NumMeshIds_];
        Info->PrintF( pdename_, 
                      "thin structure is assumed (forces on both sides!\n" );
        mpcciNode->Get( "meshIdA", meshId_[0] );
        mpcciNode->Get( "meshIdB", meshId_[0] );        
      }
    else if (MpCCIType_ == "solid")
      {
        Info->PrintF( pdename_, 
                      "thick structure is assumed (forces on one side!\n" );
	NumMeshIds_=1;
	meshId_ = new UInt[NumMeshIds_];
        mpcciNode->Get( "meshIdA", meshId_[0] );
      }
    else
      {
        EXCEPTION( "MpCCIType '" << MpCCIType_ << "does not exist" );
      }

    // Create dummy result dof object
    shared_ptr<ResultInfo> res1(new ResultInfo);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = MECH_FORCE;
    res1->dofNames = "";
    //res1->dofNames = "mpcci";
    res1->definedOn = ResultInfo::NODE;
    res1->fctType = fct;
    results_.Push_back( res1 );

  }

  void MpcciPDE::Init(UInt bcSequenceIndex )
  {

    ENTER_FCN( "MpCCI::Init()" );
    
    sequenceStep_ = bcSequenceIndex;
    
    // =====================================================================
    // get regions/subdomains for PDE
    // =====================================================================
    // Obtain regions the pde is defined on
    StdVector<ParamNode*> regionNodes = 
      myParam_->Get("regionList")->GetList("region");
    
    bool usePenalty = false;
    eqnMap_ = shared_ptr<EqnMap>(new EqnMap( ptgrid_, pdeId_, !usePenalty ));

    Info->PrintF( pdename_, " %s lives on regions:\n", pdename_.c_str());
    for ( UInt k = 0; k < regionNodes.GetSize(); k++ ) {
      std::string actRegionName = regionNodes[k]->Get("name")->AsString();
      RegionIdType actRegionId = ptgrid_->RegionNameToId( actRegionName );
      subdoms_.Push_back( actRegionId );
      Info->PrintF( pdename_, " %s\n", actRegionName.c_str() );
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[k] );
      
      //std::cerr << "results_->size():" << results_->size() << std::endl;
      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );
        
      }

    eqnMap_->Finalize();
 
    //eqnData_ = new ScalarNodeEQN( ptgrid_, subdoms_, dofspernode_, false );
    //eqnData_->CalcMpcciMapping();
    numPDENodes_ = eqnMap_->GetNumLocalNodes();
    numElems_ = eqnMap_->GetNumLocalElems();

    solVec_ = new Vector<Double>;


    PreparePDE4Computation();
    DefineSolveStep();
  }



  void MpcciPDE::PreparePDE4Computation()
  {
    ENTER_FCN( "MpcciPDE::PreparePDE4Computation" );
#ifdef MpCCI

  StdVector<Elem*> elemssd;

  MpCCInodes_ = 0;
  for (UInt i=0; i<subdoms_.GetSize(); i++)
    {
//       ptgrid_->GetVolElems(elemssd,subdoms_[i]);
      ptgrid_->GetElems(elemssd,subdoms_[i]);
      mapSD_.Clear(); // needed to call CalcNumberOfNodesInPatch(.,.) correctly
      ptgrid_->GetNodesOfElemList(mapSD_, elemssd);
      MpCCInodes_ += mapSD_.GetSize();
    }

  GetNodesOfSubdomain();
  SetupNodesSubdomainsMapping();
  ptMpCCIexch_ = new MpCCIexch(ptgrid_, subdoms_);
  for (UInt i=0; i<NumMeshIds_; i++)
    {
      for (UInt j=0; j<subdoms_.GetSize(); j++)
	{
	  ptMpCCIexch_->DefMpcciPartition(meshId_[i], j+1);
	  ptMpCCIexch_->DefMpcciNodes(meshId_[i], j+1, numOfNodesInSD_[j], localNodes_[j], eqnMap_);
	  ptMpCCIexch_->DefMpcciElements(meshId_[i], j+1, eqnMap_);
	  //possible alternative:
	  //ptMpCCIexch_->DefMpcciPartNodeElem(i+1, MpCCInodes_,*eqnData_);
	}
    }
  ptMpCCIexch_->FinishMpcciSetup(MpCCIType_);
#endif
  }

  void MpcciPDE::DefineAlgSys() {
    
    ENTER_FCN( "MpcciPDE::DefineAlgSys" );
  }

  void MpcciPDE::DefineIntegrators()
  {
    ENTER_FCN( "MpcciPDE::DefineIntegerators" );
  }


  void MpcciPDE::DefineSolveStep()
  {
    ENTER_FCN( "MpcciPDE::DefineSolveStep" );
  
    solveStep_ = new SolveStepMpCCI(*this);
  }


  
  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void MpcciPDE::InitCoupling(PDECoupling * Coupling)
  {
    ENTER_FCN( "MpcciPDE::InitCoupling" );
  
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;

    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();
  

    nonLin_ = false;

    // Initialization of coupling helper arrays
    std::string quantity;
    StdVector<UInt> * couplingnodes = NULL;
    StdVector<Elem*> interface_tmp;
    StdVector<StdVector<ShortInt> > isBoundaryNode_tmp;
//     StdVector<std::string> * neighRegions = NULL;
    StdVector<StdVector<UInt> > elemNodeToCouplingNode_tmp;
    F_Interface_.Resize(numCouplings);
    F_Interface_.Init();
    isBoundaryNode_.Resize(numCouplings);
    isBoundaryNode_.Init();
    elemNodeToCouplingNode_.Resize(numCouplings);
    elemNodeToCouplingNode_.Init();

    for (UInt actCoupling=0; actCoupling<numCouplings; actCoupling++)
      {
        // Initialize arrays for coupling surface elements
        if (ptCoupling_->GetOutputQuantity(actCoupling) == FLUID_FORCE)
          {

          
            ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
            if (couplingnodes == 0)
              std::cerr << "Couplingnodes = 0!!!!" << std::endl;
          
            F_Interface_[actCoupling] = interface_tmp;

            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);

            isBoundaryNode_tmp.Resize(interface_tmp.GetSize());
            isBoundaryNode_tmp.Init();
            elemNodeToCouplingNode_tmp.Resize(interface_tmp.GetSize());
            elemNodeToCouplingNode_tmp.Init();
          
            for (UInt ielem=0; ielem<interface_tmp.GetSize(); ielem++)
              {
                isBoundaryNode_tmp[ielem].Resize(interface_tmp[ielem]->connect.GetSize());
                isBoundaryNode_tmp[ielem].Init();
                elemNodeToCouplingNode_tmp[ielem].Resize(interface_tmp[ielem]->connect.GetSize());
                elemNodeToCouplingNode_tmp[ielem].Init();

                // Determine Boundary Nodes
                for (UInt ielemnode=0; ielemnode<isBoundaryNode_tmp[ielem].GetSize(); ielemnode++)
                  for (UInt inodes=0; inodes<(*couplingnodes).GetSize(); inodes++)
                    if (interface_tmp[ielem]->connect[ielemnode] == (*couplingnodes)[inodes] )
                      {
                        isBoundaryNode_tmp[ielem][ielemnode] = 1;
                        elemNodeToCouplingNode_tmp[ielem][ielemnode] = inodes;
                        break;
                      } // end if

              } // end for (ielems)

            isBoundaryNode_[actCoupling] = isBoundaryNode_tmp;
            elemNodeToCouplingNode_[actCoupling]  = elemNodeToCouplingNode_tmp;
          } // end if
            
      } // end for (actNode)
  }
  


  void MpcciPDE::CalcInputCoupling()
  {

    ENTER_FCN( "MpcciPDE::CalcInputCoupling" );

    std::string errMsg;
    StdVector<UInt> * nodes;
    CFSVector * val;
    UInt pdeNode;
    UInt k, couplingDof;
    UInt * nodeIds;

    // Reset counter for boundary conditions
    couplingBCsCounter_ = 0;
  
    // Outer loop over all INPUT coupling terms
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
      {

        //    ptCoupling_ = &ptCoupling_[i];
        ptCoupling_->GetInputValues(i, val);
        couplingDof = ptCoupling_->GetInputDof(i);
    
        // Up to now, Coupling is only possible with
        // Real valued solutions
        Vector<Double> const & help = dynamic_cast<Vector<Double>&>(*val);
        Vector<Double> displ;
        switch(ptCoupling_->GetInputType(i))
          {
            // -------------------
            // COORDINATE COUPLING
            // -------------------
          case COORD:
          
            ptCoupling_->GetInputNodes(i, nodes);

            for (UInt ii=0; ii<subdoms_.GetSize(); ii++)
              {
                displ.Resize( numOfNodesInSD_[ii] * ptCoupling_->GetInputDof(i) );
                displ.Init(-1);
                nodeIds=new UInt[numOfNodesInSD_[ii]];
                k=0;

                for (UInt j=0; j<nodes->GetSize(); j++)
                  {
                    pdeNode = eqnMap_->Mesh2PdeNode((*nodes)[j]);

                    // begin FASTEST NETZ CHECK
                    Point ptPoint;
                    ptgrid_->GetNodeCoordinate(ptPoint, (*nodes)[j]);
                    //Info->PrintF( pdename_, "pdeNode=%d \t (*nodes)[%d]=%d \t x=%e \t y=%e \t z=%e \n", pdeNode, j, (*nodes)[j], ptPoint[0], ptPoint[1], ptPoint[2] );
                    // end FASTEST NETZ CHECK

                    //Info->PrintF( pdename_, "pdeNode=%d \t (*nodes)[%d]=%d\n", pdeNode, j, (*nodes)[j] );
                  
                    if(NodeBelongsToSD_(pdeNode,ii)==true)
                      {
                        nodeIds[k]=pdeNode;
                        for (UInt dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
                          {
                            displ[dof + k*dim_] = help[dof + j*dim_];
                            //Test der Netzadaption in FASTEST x-Displ=yCoordinate
                            //if(dof==0) displ[dof + k*dim_] = ptPoint[1];

                            //Info->PrintF( pdename_, "%e \t ", displ[dof + k*dim_]);
                          }
                        //Info->PrintF( pdename_, "\n" );
                        k++;
                      }
                  }
                //Info->PrintF(pdename_, "k =%d \t  numOfNodesInSD_[%d]=%d \n", k, ii, numOfNodesInSD_[ii]);
      
                if (flagFirstTimeStep_== true) {;} // do nothing
                else
#ifdef MpCCI
                  ptMpCCIexch_->PutPartition(ii+1, displ, numOfNodesInSD_[ii], nodeIds, converged_);
#endif
                if(nodeIds) delete [] nodeIds;
              } // end for ii

            if (flagFirstTimeStep_== true)
              {
                //for the first iteration of the first time step 
                //do not call ptMpCCIexch_->SendAllPartitions();
                flagFirstTimeStep_=false;
              }
#ifdef MpCCI
            // Send the displacements to FASTEST-3D
            else ptMpCCIexch_->SendAllPartitions();
#endif
            break;
          case RHS:
            EXCEPTION(" No use for RHS coupling!" );
            break;
          case ID_BC:
            EXCEPTION(" No use for ID_BC coupling!" );
            break;
          case MAT:
            EXCEPTION(" No use for MAT coupling!" );
            break;
          }  // end switch
      } // end for
  }


  void MpcciPDE::CalcOutputCoupling()
  {
    ENTER_FCN( "MpcciPDE::CalcOutputCoupling" );

    SolutionType quantity;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * values = 0;
    UInt forcesCount = 0;

    std::string errMsg;
    UInt pdeNode, k; //, eqnNr,eqnDof;
    UInt * nodeIds;


    // loop over all output coupling quantities
    for (UInt actCoupling=0; actCoupling<ptCoupling_->GetNumOutputCouplings(); actCoupling++)
      {
        quantity = ptCoupling_->GetOutputQuantity(actCoupling);
        ptCoupling_->GetOutputValues(actCoupling, values);

        Vector<Double> * temp = dynamic_cast<Vector<Double> *>(values);
        Vector<Double> subdomForces;

        temp->Init(0.0);

        switch(ptCoupling_->GetOutputType(actCoupling))
          {
	  
          case NODE:	  
            
            ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
            
            if (quantity == FLUID_FORCE)
              {
#ifdef MpCCI
                if (MpCCIType_=="shell")
                  {
                    Info->PrintF( pdename_, 
                        "thin structure is assumed (forces on both sides!\n" );
                  }
                else if (MpCCIType_=="solid")
                  {
                    Info->PrintF( pdename_, 
                        "thick structure is assumed (forces on oned side!\n" );
                  }
                else
                  {
                    std::string errmsg = "MpCCIType '" + MpCCIType_;
                  }

                ptMpCCIexch_->RecvAllPartitions(MpCCIType_);
                //ptMpCCIexch_->RecvAllPartitions("shell");
#endif
                forcesCount++;
                for (UInt ii=0; ii<subdoms_.GetSize(); ii++)
                  {
                    subdomForces.Resize(numOfNodesInSD_[ii] * ptCoupling_->GetInputDof(actCoupling));
                    subdomForces.Init(0.0);
                    nodeIds=new UInt[numOfNodesInSD_[ii]];
#ifdef MpCCI
                    ptMpCCIexch_->GetNodalValOfOnePartition(ii+1 , subdomForces, numOfNodesInSD_[ii], localNodes_[ii], MpCCIType_);
#endif
                    for (k=0; k<numOfNodesInSD_[ii]; k++)
                      {
                        for (UInt j=0; j<couplingNodes->GetSize(); j++)
                          {
                            pdeNode = eqnMap_->Mesh2PdeNode((*couplingNodes)[j]);
                         
                            if(pdeNode==localNodes_[ii][k])
                              {
                                for (UInt dof=0; dof<ptCoupling_->GetInputDof(actCoupling); dof++)
                                  {
                                    temp[actCoupling][dof + j*dim_] += subdomForces[dof + k*dim_];
                                  }
                              }
                          }
                      }
                  }
              }
            break;
          case ELEM:
            EXCEPTION("No Element coupling output" );
          }
      }
    
  }


  bool MpcciPDE::HasOutput(SolutionType output)
  {
    ENTER_FCN( "MpcciPDE::HasOutput" );
  
    switch (output)
      {
      case FLUID_FORCE:
        return true;
        break;
      default:
        return false;
        break;
      }
    return false;
  }



  void MpcciPDE::GetNodesOfSubdomain()
  {
    ENTER_FCN( "MpcciPDE::GetNodesOfSubdomain" );
    UInt i, j;
    UInt localPDENode;
    UInt numOfSubdom=subdoms_.GetSize();
    StdVector<UInt> globalNodes;
    StdVector<Elem*> elemsInSD;

    localNodes_=new UInt*[numOfSubdom];

    numOfNodesInSD_=new UInt[numOfSubdom];
    for (i=0; i<numOfSubdom; i++)
      {
	globalNodes.Resize(0);
// 	ptgrid_->GetVolElems(elemsInSD,subdoms_[i]);
	ptgrid_->GetElems(elemsInSD,subdoms_[i]);
	ptgrid_->GetNodesOfElemList(globalNodes, elemsInSD);
	numOfNodesInSD_[i] = globalNodes.GetSize();
	localNodes_[i]  =new UInt[numOfNodesInSD_[i]];
	for (j=0; j<numOfNodesInSD_[i]; j++)
	  {
	    localPDENode = eqnMap_->Mesh2PdeNode(globalNodes[j]);
	    localNodes_[i][j]=localPDENode;
	  }
      }
  }

  void MpcciPDE::SetupNodesSubdomainsMapping()
  {
    ENTER_FCN( "MpcciPDE::SetupNodesSubdomainsMapping" );

    UInt i, j, k;
    UInt numOfSubdom=subdoms_.GetSize();

    NodeBelongsToSD_.Resize(numPDENodes_+1, numOfSubdom);
    NodeBelongsToSD_.Init(false);

    for (i=1; i<=numPDENodes_; i++)
      {
	for (j=0; j<numOfSubdom; j++)
	  {
	    for (k=0; k<numOfNodesInSD_[j]; k++)
	      {
		if(i==localNodes_[j][k])
		  {
		    NodeBelongsToSD_(i,j)=true;
		    Info->PrintF( pdename_, "pdeNode:%d is in SD %d \t", i, j );
		    break;
		  }
	      }
	  }
	Info->PrintF( pdename_, "\n");
      }
  }

} // end of namespace

