// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
//#include <stdio.h>
#include <complex>

#include "acouFlowNoise.hh"

#include "DataInOut/freqfunc.hh"
#include "Forms/forms_header.hh"
#include "MpCCIcpl/MpCCIexch.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "Driver/stdSolveStep.hh"

#ifdef MpCCI

#if (MpCCI_RELEASE == 305)
#include <mpcci.h>
#else
#include <cci.h>
#endif

#endif

namespace CoupledField
{

  AcouFlowNoise::AcouFlowNoise(Grid *aptgrid, ParamNode* paramNode )
    :AcousticPDE(aptgrid, paramNode)
  {
    ENTER_FCN( "AcouFlowNoise::AcouFlowNoise" );

    // pdename_ is also acoustic for this case
    pdename_          = "acoustic";
    pdematerialclass_ = FLUID;
    nodalSrc_ = false;
    vortexSrc_ = false;
    plotRHS_ = false; 
    plotRHSVel_= false;
    isHarmonic_=false;
    vortexFlag_=0;
    pressFormul_=false;

    myParam_->Get( "pressFormul", pressFormul_, false );
    if( pressFormul_ ) {
      Info->PrintF("acoustic","Using sources computed from fluid pressure\n");
    }    

    myParam_->Get( "valRHS", plotRHS_, false );
    if( plotRHS_ ) {
      Info->PrintF(pdename_, "Writing RHS as solution of problem\n" );
    }

    StdVector<Elem*> elemssd;
    StdVector<std::string> regionNames, coupledRegionNames;

    // -----------------------------------------------------------------------
    // @Max: Please note, that the values "writeVortexFineGrid"
    // and "writeVortexSrcsperTS" are not defined in the xsd-schema definition,
    // so, in fact, they are NEVER set!
    // -----------------------------------------------------------------------
    //For writing fine vortex flow field in files for later MpCCI exchange
      writeGridFile_ = false;
      writeSrcFileperTS_ = false;
      myParam_->Get( "writeVortexFineGrid", writeGridFile_, false );
    if( writeGridFile_ ) {
      Info->PrintF("acoustic","Writing fine grid def. file of coupled domain\n");
    }
    myParam_->Get( "writeVortexSrcsperTS", writeSrcFileperTS_, false );
    if( writeSrcFileperTS_ ) {
      // Create directory to save in the source files
      std::string S="mkdir -p fineSrcs_perTS";
      system(S.c_str());
      Info->PrintF("acoustic","Writing fine sources in time files (NrFiles=NrTimeSteps)\n");
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

#ifdef MpCCI
    //check type of flow data
    ParamNode * flowNode = myParam_->Get("flowData", false );
    if (flowNode) 
      {
        nodalSrc_ = myParam_->Has("flowData", "type", "nodalSrc" );
      }
    
    if( nodalSrc_ ) {
      Info->PrintF(pdename_, " Using FlowData as RHS nodal source\n" );
    }

    StdVector<ParamNode*> regionNodes = 
      myParam_->Get("regionList")->GetList( "region" );

    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
      regionNames.Push_back( regionNodes[i]->Get("name")->AsString() );
    }

    StdVector<ParamNode*> cplRegionNodes = 
      param->Get( "MpCCI-flownoise" )->GetList( "coupledregion" );
    std::cerr << "cplRegions = " << cplRegionNodes.Serialize() << std::endl;
    for( UInt i = 0; i < cplRegionNodes.GetSize(); i++ ) {
      coupledRegionNames.Push_back( cplRegionNodes[i]->Get("name")->AsString() );
    }

    StdVector<RegionIdType> regionIds;
    ptgrid_->RegionNameToId( regionIds, regionNames );
    ptgrid_->RegionNameToId( couplSubDomId_, coupledRegionNames );
    
    bool OnlyLinNodes=true;
    bool Find=false;
    //    std::cerr << "\nName of all regions:\n" << regionNames << std::endl;
    //    std::cerr << "Coupled regions:\n" << coupledRegionNames  << std::endl;
    

    for (UInt i=0; i<regionIds.GetSize(); i++)
      {
        for (UInt j=0; j<couplSubDomId_.GetSize(); j++)
          {
            if (coupledRegionNames[j] == regionNames[i])
              {
                std::cout<<"\ncoupledRegionNames[j]: "<<coupledRegionNames[j]<<std::endl;
                std::cout<<"\nregionNames[i]: "<<regionNames[i]<<std::endl;
                ptgrid_->GetVolElems(elemssd,regionIds[i]); //couplSubDomId_[j]);
		
                if (dim_ == 3)
                  { //TO SEND ONLY LINEAR 3D ELEMS TO MPCCI SINCE NO QUADRATIC 
                    //ARE ALLOWED
                    ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, false);
                    ptgrid_->GetNodesOfElemList(mapSD_onlyLinNodes_, elemssd, OnlyLinNodes);
                    ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_onlyLinNodes_);
                  }
                else
                  {
                    ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, false);
                    //ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_.GetSize() );
                    ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_allNodes_);
                    std::cout<<"\nmapSD_allNodes_.GetSize(): "<<mapSD_allNodes_.GetSize()<<std::endl;
                  }
		
                Find=true;

              }
          }
      }
    if (!Find) 
      {
        EXCEPTION( "Subdom to be coupled with MpCCI is not in list of PDE subdoms. "
                   << "Please, check .xml-file" );
      }
    // On CFS side we always use the complete coupled region and later give a zero src to the
    // quadratic (or mid) nodes in case of 3D quadratric acoustic grid.
    MpCCInodes_=mapSD_allNodes_.GetSize();
    ptMpCCIexch_->PutExchangeGrid2MpCCI(couplSubDomId_);

#else
    
    ParamNode * flowDataNode = myParam_->Get( "flowData" );
    if( flowDataNode ) 
      vortexSrc_ = flowDataNode->Get("type")->AsString() == "vortexSrc";
    StdVector<RegionIdType> regionIds;
    if( vortexSrc_  ) {
      StdVector<ParamNode*> regionNodes = 
        myParam_->Get("regionList")->GetList( "region" );
      for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
        regionNames.Push_back( regionNodes[i]->Get("name")->AsString() );
      }
      StdVector<ParamNode*> cplRegionNodes = 
        param->Get( "MpCCI-flownoise" )->GetList( "coupledregion" );
      for( UInt i = 0; i < cplRegionNodes.GetSize(); i++ ) {
        coupledRegionNames.Push_back( cplRegionNodes[i]->Get("name")->AsString() );
      }
      ptgrid_->RegionNameToId( regionIds, regionNames );
      ptgrid_->RegionNameToId( couplSubDomId_, coupledRegionNames );
      std::cerr << "\nSubdoms:\n" << regionNames << std::endl;
      std::cerr << "CoupledSubdoms:\n" << coupledRegionNames << std::endl;
      Info->PrintF(pdename_, " Using Vortex function as RHS nodal source\n" );
    }
    else 
      {
        nodalSrc_ = flowDataNode->Get("type")->AsString() == "nodalSrc";
        if( nodalSrc_ ) {
          // Now verify that the type of analysis is HARMONIC
          // Determine type of analysis
          if ( analysistype_==HARMONIC )
            {
              isHarmonic_=true;
              ComputeRHSforHarm_=true;    
              Info->PrintF(pdename_, "Using FlowData from dataset as RHS nodal source\n" );
              Info->PrintF(pdename_, "Computing using nodal frequency files (No MpCCI used)\n" );



              // Now getting the correct size of coupled domain. It must have the same
              // numbering as the nodal frequency file, if not freqfunc will complain!!
              StdVector<ParamNode*> regionNodes = 
                myParam_->Get("regionList")->GetList( "region" );
              for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
                regionNames.Push_back( regionNodes[i]->Get("name")->AsString() );
              }
              StdVector<ParamNode*> cplRegionNodes = 
                param->Get( "MpCCI-flownoise" )->GetList( "coupledregion" );
              for( UInt i = 0; i < cplRegionNodes.GetSize(); i++ ) {
                coupledRegionNames.Push_back( cplRegionNodes[i]->Get("name")->AsString() );
              }
              ptgrid_->RegionNameToId( regionIds, regionNames );
              ptgrid_->RegionNameToId( couplSubDomId_, coupledRegionNames );

              bool Find=false;
              //              std::cerr << "Name of all regions:\n" << regionNames << std::endl;
              //              std::cerr << "Coupled regions:\n" << coupledRegionNames  << std::endl;
    

              for (UInt i=0; i<regionIds.GetSize(); i++)
                {
                  for (UInt j=0; j<couplSubDomId_.GetSize(); j++)
                    {
                      if (coupledRegionNames[j] == regionNames[i])
                        {
                          ptgrid_->GetVolElems(elemssd,couplSubDomId_[j]);
                          ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, false);
                          //true returns only the corner nodes (for 3D MpCCI linear coupling)
                          if (dim_==3)
                            ptgrid_->GetNodesOfElemList(mapSD_onlyLinNodes_, elemssd, true);
                          Find=true;
                        }
                    }
                }
              if (!Find) 
                {
                  EXCEPTION( "Subdom to be coupled with NodalFreqFilesSet is not "
                             << "in list of PDE subdoms. Please, check .xml-file" );
                }
              UInt CoupledNodes;
              if (dim_==3)
                CoupledNodes=mapSD_onlyLinNodes_.GetSize();
              else 
                CoupledNodes=mapSD_allNodes_.GetSize();
              flowdata_.Resize(1, CoupledNodes);
              flowdata_.Init();
            }
          else
            {
              EXCEPTION( "Only Nodal Frequency Sources or Vortex Source are allowed "
                         << "if not using MpCCI" );
            }
        }
        
      }

#endif

  }

  AcouFlowNoise::~AcouFlowNoise()
  {
    ENTER_FCN( "AcouFlowNoise::AcouFlowNoise" );
#ifdef MpCCI
    delete ptMpCCIexch_;
#endif
  }

  void AcouFlowNoise::DefineSolveStep()
  {
    ENTER_FCN( "AcouFlowNoise::DefineSolveStep" );

    solveStep_ = new StdSolveStep(*this);
  }

  void AcouFlowNoise::ComputeRHS(const Double atime)
  {
    ENTER_FCN( "AcouFlowNoise::ComputeRHS" );

  
    Vector<Double> elemvec, nodalval;
    UInt i;
    Matrix<Double> ptCoordNodes;
    StdVector<UInt> connecth;
    // For changing connecth to PDE
    StdVector<Integer> connect_PDE; 
    BaseFE * ptEl;
    Double actFreq;   
    UInt j;
    UInt elsize = 0;
    StdVector<Elem*> elemssd;     
    Double valmult=1.0;
    Vector<Double> valVec(1);

    static UInt timestep=1;

#ifdef MpCCI
    double starttime, endtime;
    //std::cout<<"MpCCInodes_: "<< MpCCInodes_ << " dimension: " 
    //             << dim_ << std::endl;
    if (nodalSrc_ == true) 
      {
      //we get already the integrated acoustic source term
//      if ( variableSpeedOfSoundCN_ )  {
         flowdata_.Resize(1, MpCCInodes_);
//      }
 //     else 
 //        flowdata_.Resize(1, MpCCInodes_);
      }
    else
      flowdata_.Resize(1+dim_, MpCCInodes_);

    flowdata_.Init();

    //starttime = CCI_Wtime();
 
    //Integer timestep = 0;
    ptMpCCIexch_->CouplCompPhase(flowdata_, atime);
 

    //endtime = CCI_Wtime();
      

    //    std::cout<<"Transfer of Data CouplCompPhase() for 1 time step took: "
    //         <<(endtime-starttime)<<" seconds"<<std::endl;

    // Quadrupole computation
    std::cout << "Receiving "<<MpCCInodes_<<" acoustic nodal sources for time t = "<< atime<<"s"<< std::endl;


#else
    // If data from fluid file, call to get fluid flow data in flowdata_  
    std::string flowdatafile;
    //    Warning( "Using no MpCCI, vortexSrc is assumed.", __FILE__, __LINE__ );
      
    //     conf->get("acousrc_file",flowdata);
    //     ReadFlowData(flowdatafile.c_str(), timestep, flowdata_); 
#endif 


    ///////For ramping
      Vector<Double> fmin, fmax, facRampMin, facRampMax, bndoffsetMin, bndoffsetMax;
      fmin.Resize(dim_);
      fmax.Resize(dim_);
      facRampMin.Resize(dim_);
      facRampMax.Resize(dim_);
      bndoffsetMin.Resize(dim_);
      bndoffsetMax.Resize(dim_);
      //For damping also in harmonic case in the outflow region of the fluid domain
      if (nodalSrc_)
      {

        // Fetch MpCCI-Node and query data
        ParamNode * mpcciNode = param->Get( "MpCCI-flownoise" );
        mpcciNode->Get("xfmin",fmin[0]);
        mpcciNode->Get("yfmin",fmin[1]);
        mpcciNode->Get("xfmax",fmax[0]);
        mpcciNode->Get("yfmax",fmax[1]);
        mpcciNode->Get("facrampXmin",facRampMin[0]);
        mpcciNode->Get("facrampYmin",facRampMin[1]);
        mpcciNode->Get("facrampXmax",facRampMax[0]);
       mpcciNode ->Get("facrampYmax",facRampMax[1]);

        if (dim_==3)
          {
            mpcciNode->Get("zfmin",fmin[2], "MpCCI-flownoise");
            mpcciNode->Get("zfmax",fmax[2], "MpCCI-flownoise");
            mpcciNode->Get("facrampZmin",facRampMin[2], "MpCCI-flownoise");
            mpcciNode->Get("facrampZmax",facRampMax[2], "MpCCI-flownoise");
            
          } 

        for (UInt actDim=0; actDim<dim_; actDim++)
          {
            bndoffsetMin[actDim]=facRampMin[actDim]*fmin[actDim];
            bndoffsetMax[actDim]=facRampMax[actDim]*fmax[actDim];
          }
      }
    
#ifdef MpCCI
    //    starttime = CCI_Wtime();
#endif 
    if (nodalSrc_ == false) 
      {
      
        if (vortexSrc_ == true)
          {
            std::cout<<"Calling ComputeRHSwithVortexSource(atime)"<<std::endl;
            ComputeRHSwithVortexSource(atime);
          }
        else
          {
            //Using interpolated velocity field and computing source in CFS
            // Correct valmult value is -1.0, 
            // if plugging in source (ddTij/dxidxj) directly in weak form then 1.0
            valmult=-1.0;
      
            for (i=0; i<couplSubDomId_.GetSize(); i++)
              {
                ptgrid_->GetVolElems(elemssd,couplSubDomId_[i]);
              
                for (j=0; j< elemssd.GetSize(); j++)
                  {
                    ptEl=elemssd[j]->ptElem;

                    LinearFlowNoiseInt * linear_load = 
                      new LinearFlowNoiseInt(ptEl);
                  
                    UInt ii;
                    elsize=(elemssd[j]->connect).GetSize();
                    connecth.Resize(elsize);
                    for (ii=0; ii<elsize; ii++)
                      connecth[ii]=(elemssd[j]->connect)[ii];
                  
                    Matrix<Double> ptCoordNodes;
                    ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);        
                    linear_load->CalcElemVector4Quad(ptCoordNodes, connecth,
                                                     flowdata_, elemvec);
                    elemvec*=valmult;
                  
                    // Ramping before adding elemrhs to global vector 
                    // to avoid spurious effect at bnd. of fluid dom.
                  
                    //               for (ii=0; ii<elsize; ii++)
                    //                 {
                    //                   if (ptCoordNodes[0][ii]<bndoffsetMin[0])
                    //                     {
                    //                       elemvec[ii]-=elemvec[ii]*
                    //                   (ptCoordNodes[0][ii]-bndoffsetMin[0])/(fmin[0]-bndoffsetMin[0]);
                    //                     }
                  
                    //                   else
                    //                     if (ptCoordNodes[0][ii]>bndoffsetMax[0])
                    //                       elemvec[ii]-=elemvec[ii]*
                    //                   (ptCoordNodes[0][ii]-bndoffsetMax[0])/(fmax[0]-bndoffsetMax[0]);
                    //                   if (ptCoordNodes[1][ii]<bndoffsetMin[1])
                    //                     elemvec[ii]-=elemvec[ii]*
                    //                   (ptCoordNodes[1][ii]-bndoffsetMin[1])/(fmin[1]-bndoffsetMin[1]);
                    //                   else    
                    //                     if (ptCoordNodes[1][ii]>bndoffsetMax[1])
                    //                       elemvec[ii]-=elemvec[ii]*
                    //                   (ptCoordNodes[1][ii]-bndoffsetMax[1])/(fmax[1]-bndoffsetMax[1]);
                    //                 }
                  
                    //end ramping
                  
                    // CHANGE connecth
                    //Mesh2PDENode(connect_PDE,connecth,mesh2PDENode_);
                    eqnMap_->GetNodeEqn(connecth, connect_PDE);
                    // for setting with homogeneous rhs       
                    //linear_load->CalcElemVector(ptCoordNodes, elemvec); 
                  
                    //    std::cout<<"elemvect QUADRUPOLE: "<<elemvec<<std::endl;
                    // Quadrupole activated!!   
                    algsys_->SetElementRHS(&elemvec[0], pdeId_, 
                                           connect_PDE.GetPointer(), connect_PDE.GetSize());
                  
                    delete linear_load;
                  }
              }
          }
      }
    
    else
      {     //for assigning acoustic nodal sources
        UInt node, dof;
        Integer eqnNr;
        StdVector<UInt> connect(1);
        Double srcVal;
        Double pointX, pointY, pointZ;
        bool inBox; 

        valmult = -1.0;    
        dof     = 1;

        if (!isHarmonic_)
          {
            for (UInt idx=0; idx<flowdata_.GetSizeCol() ; idx++)
              {
                srcVal = flowdata_[0][idx];
                node = idx + 1;
              
                // Ramping before adding to RHS vector
                Matrix<Double> ptCoordNodes;
                connecth.Resize(1);
                connecth[0] = node;
                ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);         

		pointX = ptCoordNodes[0][0];
		pointY = ptCoordNodes[1][0];
		if(dim_==3) {
		  pointZ = ptCoordNodes[2][0];
		}
		//check, if point is in specified (xml-file) box 
		inBox = false;
		if ( pointX > fmin[0] &&  pointX < fmax[0] ) {
                  if ( pointY > fmin[1] &&  pointY < fmax[1] ) {
                    inBox = true;
		    if ( dim_ == 3 ) {
		      inBox = false;
		      if ( pointZ > fmin[2] &&  pointZ < fmax[2] ) {
			inBox = true;
		      }
		    }
		  }
		}
		if ( !inBox ) {
                  srcVal = 0.0;
		}
		else {
                  for (UInt actDim=0; actDim<dim_; actDim++)
                    {
                      if (ptCoordNodes[actDim][0]<bndoffsetMin[actDim])
                        srcVal -= srcVal*(ptCoordNodes[actDim][0]-bndoffsetMin[actDim])/
                          (fmin[actDim]-bndoffsetMin[actDim]);
                      else
                        if (ptCoordNodes[actDim][0]>bndoffsetMax[actDim])
                          srcVal -= srcVal*(ptCoordNodes[actDim][0]-bndoffsetMax[actDim])/
                            (fmax[actDim]-bndoffsetMax[actDim]);
                    }
                }
                //End of ramping
            		
		srcVal *= valmult;//Due to weak formulation, factor for RHS sources
		
		//add to RHS
		eqnNr = eqnMap_->GetNodeEqn(node,dof );
		algsys_->SetNodeRHS(srcVal, pdeId_, eqnNr, 1);   
		
// 		if ( saveNodalSourcesRHS_ ) {
// 		  valVec[0] = srcVal;
// 		  rhsNodalSrc_.SetNodalResult(node, valVec);
// 		}
              }
          } 
        else
          {
            std::cout << "Using frequency source files..."<< std::endl;
            actFreq = atime;
            Double omega = actFreq * 2 * PI;
            FreqFunc * ptFreqFunc = new FreqFunc();
            StdVector<Double> Ampl_Phase;
            Ampl_Phase.Resize(2);
            Ampl_Phase.Init();
            //Getting freq. file name without node number extension from xml file
            std::string nameFreqFile;
            param->Get("sequenceStep", "index", GenStr(sequenceStep_) )
              ->Get("analysis")->Get("harmonic")->Get("freqDataFile", nameFreqFile );
            for (UInt idx=0; idx<flowdata_.GetSizeCol() ; idx++) {
              if (dim_==3)
                node = mapSD_onlyLinNodes_[idx];
              else
                node = mapSD_allNodes_[idx];    
          
              Ampl_Phase=ptFreqFunc->NodalFreqFuncAtFreq(actFreq,nameFreqFile.c_str(),
                                                         node);              
              Double valAmpl = Ampl_Phase[0];
              Double valPhase = Ampl_Phase[1];

              if (pressFormul_)
                {
                  valAmpl = -omega*omega*valAmpl;
                }
            
              valAmpl*=valmult; //Due to weak formulation, factor for RHS sources

              Complex complexValue( valAmpl * cos( valPhase / 180 * PI ),
                                    valAmpl * sin( valPhase / 180 * PI ) );

              // Ramping before adding to RHS vector (NOW IT MAKES ZERO THE RHS ENTRY!)
              Matrix<Double> ptCoordNodes;
              connecth.Resize(1);
              connecth[0] = node;
              ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);         
            
              for (UInt actDim=0; actDim<dim_; actDim++)
                {
                  if (ptCoordNodes[actDim][0]<bndoffsetMin[actDim])
                    //valAmpl -= val*(ptCoordNodes[actDim][0]-bndoffsetMin[actDim])/(fmin[actDim]-bndoffsetMin[actDim]);
                    complexValue = 0;
                  else
                    if (ptCoordNodes[actDim][0]>bndoffsetMax[actDim])
                      //valAmpl -= val*(ptCoordNodes[actDim][0]-bndoffsetMax[actDim])/(fmax[actDim]-bndoffsetMax[actDim]);
                      complexValue = 0;
                }
              //end ramping
            
              //add to RHS
              //Since it can be that there is no sequential order of nodes
              //we use a map with index idx
              if (dim_ == 3)
                {
                  eqnNr = eqnMap_->GetNodeEqn(mapSD_onlyLinNodes_[idx],dof );
                  //std::cout<<"EqnDOF= "<<eqnDof<<" EqnNr= "<<eqnNr<<" dof= "<<dof
                  //         <<" NodeNumber= "<<mapSD_onlyLinNodes_[idx]
                  //         <<" Ampl= "<<Ampl_Phase[0]<<std::endl;
                }
              else
                eqnNr = eqnMap_->GetNodeEqn(mapSD_allNodes_[idx],dof );
            
              algsys_->SetNodeRHS(complexValue, pdeId_, eqnNr, 1);
            }
          }//end else for frequency analysis
      }//end else in case nodalSrc is true
    
#ifdef MpCCI
    //    endtime = CCI_Wtime();
#endif    
  
    if (plotRHS_ && !isHarmonic_ && vortexFlag_!=6 && vortexFlag_!=7){
      // For plotting the RHS as solution for analysing it
      rhs_.SetNumSolutions(1);
      rhs_.SetNumNodes(numPDENodes_);
      rhs_.SetSolutionType(ACOU_RHS_LOAD);
      rhs_.SetResult( results_[0] );
      rhs_.SetNumDofs(1);
      rhs_.SetPtrEQNData( eqnMap_.get(), ptgrid_);
      rhs_.SetRegions( subdoms_ );
      rhs_.Init(0.0);
      
      Double *ptRHS;
      algsys_->GetRHSVal( ptRHS );
      rhs_.CopyFromAlgSysDataPointer(ptRHS);

      //For writing out the topology and RHSval in files (only for Vortex simul.)
      if (writeGridFile_ && writeSrcFileperTS_ && vortexSrc_)
        {
          WriteOutVortexFineSources(timestep);
        }
    }
    
    timestep=timestep+1;

  } 

  void  AcouFlowNoise::WriteOutVortexFineSources(UInt timestep)
  {
    ENTER_FCN( "AcouFlowNoise::WriteOutVortexFineSources" );

    // static variable for suffix of output src file
    static Integer filenum=1;
    std::string filename;
    StdVector<UInt> connecth;
    filename = "fineSrcs_perTS/srcfile_1p0mill_TS0p25_";
    if ( filenum < 10 ) filename.append( "000" );
    else if ( filenum < 100 ) filename.append( "00" );
    else if ( filenum < 1000 ) filename.append( "0" );
    else if ( filenum > 10000 ) {
      EXCEPTION( "Number of src file exceeds 9999!" );
    }
    filename.append( GenStr( filenum ) );
    filename.append( ".dat" );
    filenum++;
    outsrcfile_ = new std::ofstream(filename.c_str());

    StdVector<Elem*> elemssd;  
    Vector<Double> ptNodalCoord2D;   
    for (UInt i=0; i<couplSubDomId_.GetSize(); i++)
      {
        UInt node, dof;
        Integer eqnNr;
        dof    = 1;
        UInt elsize = 0;
        //First we get the elem definitions
        ptgrid_->GetVolElems(elemssd,couplSubDomId_[i]);
        if(timestep==1)
          {
                  
            for (UInt j=0; j< elemssd.GetSize(); j++)
              {
                elsize=(elemssd[j]->connect).GetSize();
                connecth.Resize(elsize);
                connecth.Init();
                for (UInt ii=0; ii<elsize; ii++)
                  {
                    connecth[ii]=(elemssd[j]->connect)[ii];
                    if ((ii<(elsize-1)))
                      {
                        (*outelemfile_) << connecth[ii]<<" ";
                      }	      
                    else
                      {
                        (*outelemfile_) << connecth[ii];
                      }
                    //std::cout<<connecth[ii]<<" ";
                  }
                (*outelemfile_) << std::endl;
                //std::cout<<std::endl;

                //Matrix<Double> ptCoordNodes;
                //ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);
              }
          }
              
        //Now we get the node numbers and corresponding RHS values

        ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, false);
                               
        for (UInt idx=0; idx<mapSD_allNodes_.GetSize() ; idx++)
          {
            node = mapSD_allNodes_[idx];
            //std::cout<<std::endl;
            Double RHSnodalVal=0;
            if (timestep==1)
              {
                ptgrid_->GetNodeCoordinate(ptNodalCoord2D, node);
                (*outnodefile_) << ptNodalCoord2D[0]<<" "<< ptNodalCoord2D[1]<<" "<<0.0<<std::endl;
                //    std::cout<<"x: "<<ptNodalCoord2D[0]<<", y: "<<ptNodalCoord2D[1]<<std::endl; 
              }
            eqnNr = eqnMap_->GetNodeEqn(node,dof );
            rhs_.Get(idx,1,RHSnodalVal);
            (*outsrcfile_) <<  RHSnodalVal<<std::endl;
            //std::cout<<"node: "<<node<<" eqnNr: "<<eqnNr<<" RHSnodalVal: "<<RHSnodalVal<<std::endl;
          }
      }    
    
  }
  

 void AcouFlowNoise::ComputeRHSwithVortexSource(const Double atime)
 {
   ENTER_FCN( "AcouFlowNoise::ComputeRHSwithVortexSource" );

   Matrix<Double> Src_Matrix;
   Vector<Double> elemvec, nodalval;
   UInt j;
   UInt elsize = 0;
   StdVector<Elem*> elemssd;  
   StdVector<UInt> connecth;
   // For changing connecth to PDE
   StdVector<Integer> connect_PDE; 
   BaseFE * ptEl; 
    
   // vortexFlag_: 1: Pak, 2: -dPdt, 3: dTijdi vector, 4: Only vel. vals
   myParam_->Get( "flowData" )->Get( "vortexSrc" )->Get("outType", vortexFlag_ );
   for (UInt i=0; i<couplSubDomId_.GetSize(); i++)
     {
       std::cerr << "Current time: "<<atime;
       std::cerr <<std::endl;
       ptgrid_->GetVolElems(elemssd,couplSubDomId_[i]);
       if (vortexFlag_==1 || vortexFlag_==5 || vortexFlag_==6 || vortexFlag_==7)
         //Getting Analytical solution (P_ak) or Tangential velocity (arg) as RHSval
         {
           std::cout<<"Getting Analytical solution (P_ak) or Tangential velocity (arg) as RHSval"<<std::endl;
           UInt node, dof;
           Integer eqnNr;
           StdVector<UInt> connect(1);
           dof    = 1;
           if (vortexFlag_==6 || vortexFlag_==7)
             {
               plotRHSVel_= true;
               std::cout<<"Init rhs_ for putting vel vector field in RHSval for visualization"<<std::endl;
               std::cout<<"numPDENodes_: "<<numPDENodes_<<std::endl;
               std::cout<<"elemssd.GetSize(): "<<elemssd.GetSize()<<std::endl;
               std::cout<<"couplSubDomId_[i]: "<<couplSubDomId_[i]<<std::endl;
               rhs_.SetNumSolutions(1);
               rhs_.SetNumNodes(numPDENodes_);
               //             rhs_.SetNumNodes(elemssd.GetSize());
               rhs_.SetSolutionType(ACOU_RHS_LOAD);
               rhs_.SetResult( results_[0] );
               rhs_.SetNumDofs(1);
               rhs_.SetResult( results_[0] );
               rhs_.SetPtrEQNData( eqnMap_.get(), ptgrid_);
               rhs_.SetRegions( subdoms_ );
               rhs_.Init(0.0);
               
               rhs2_.SetNumSolutions(1);
               rhs2_.SetNumNodes(numPDENodes_);
               //               rhs2_.SetNumNodes(elemssd.GetSize());
               rhs2_.SetSolutionType(ACOU_POT_NRBC);
               rhs2_.SetResult( results_[0] );
               rhs2_.SetNumDofs(1);
               rhs2_.SetResult( results_[0] );
               rhs2_.SetPtrEQNData(eqnMap_.get(),ptgrid_);
               rhs2_.SetRegions( subdoms_ );
               rhs2_.Init(0.0);

             }

           ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, false);
           std::cout<<"mapSD_allNodes_SIZE: "<<mapSD_allNodes_.GetSize()<<std::endl;
           for (UInt idx=0; idx<mapSD_allNodes_.GetSize() ; idx++)
             {
               Double val = mapSD_allNodes_[idx];
               //node = mapSD_allNodes_[idx];
               node = idx + 1;

               Matrix<Double> ptCoordNodes;
               connecth.Resize(1);
               connecth[0] = node;
               ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);         
               
               Double P_ak;
               Vector<Double> vectorVal;
               VortexAnalytical(P_ak, vectorVal, ptCoordNodes[0][0],ptCoordNodes[1][0] , atime, vortexFlag_);
               //std::cout<<"After getting P_ak in ComputeRHSwithVortexSource(), P_ak= "<<P_ak<<std::endl;
               val=P_ak;
               Vector<Double>tempVelX, tempVelY;
               tempVelX.Resize(1);
               tempVelY.Resize(1);
               
               //add to RHS
               eqnNr = eqnMap_->GetNodeEqn(node,dof );
               algsys_->SetNodeRHS(val, pdeId_, eqnNr, 1 ); 
                             
               if (vortexFlag_==6 || vortexFlag_==7)
                 {
                   eqnMap_->GetNodeEqn(connecth, connect_PDE);
                   tempVelX[0]=vectorVal[0];
                   tempVelY[0]=vectorVal[1];
                   rhs_.SetNodalResult(node,tempVelX);
                   rhs2_.SetNodalResult(node,tempVelY);
                 }
             }
         }
       else
         { //Using Vortex Analytical with potential function
           
           for (j=0; j< elemssd.GetSize(); j++)
             {
               UInt ii;
               elsize=(elemssd[j]->connect).GetSize();
               elemvec.Resize(elsize);
               elemvec.Init();
               nodalval.Resize(elsize);
               nodalval.Init();
               connecth.Resize(elsize);
               for (ii=0; ii<elsize; ii++)
                 connecth[ii]=(elemssd[j]->connect)[ii];
               Matrix<Double> ptCoordNodes;
               ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);

//                if (vortexFlag_==3 || vortexFlag_==4)
//                  {
                   Src_Matrix.Resize(dim_,elsize);
                   for (UInt i=0;i<elsize;i++)
                     for (UInt j=0;j<dim_;j++)
                       Src_Matrix[j][i]=0.0;
//                 }
                  
                   Double r_sqr;
                   
               for (ii=0; ii<elsize; ii++)
                 {
                   r_sqr=((ptCoordNodes[0][ii])*(ptCoordNodes[0][ii])+
                                 (ptCoordNodes[1][ii])*(ptCoordNodes[1][ii]));

                   //Original core (about 2.5m)  
                   if (r_sqr<=(200./81.)*(200./81.) && (vortexFlag_==2))
                     //Smaller core source (about 1.5m core)  
                   //  if (r_sqr<=((140./81.)*(140./81.))) 
                   //Bigger core (about 12m)
                   //if (r_sqr<=((1000./81.)*(1000./81.))) //Give shifted results (above zero)
                   //Full core   
                     //      if (r_sqr<=-1)
                   //Very small core just to avoid central singularity (center of fluid domain)
                            //    if (r_sqr<=0.025)
                     {
                       nodalval[ii]=0; 
                       Src_Matrix[0][ii]=0;
                       Src_Matrix[1][ii]=0;
                     }
                   else
                     {
                       Double srcVal=0;
                       Vector<Double> vectorVal;
                       //std::cout << "In else before calling VortexAnalytical"<<std::endl;
                       VortexAnalytical(srcVal, vectorVal, ptCoordNodes[0][ii],
                                        ptCoordNodes[1][ii], atime, vortexFlag_);



                       if (vortexFlag_==3 || vortexFlag_==4)
                         {
                           for (UInt j=0;j<dim_;j++)
                             {
                               Src_Matrix[j][ii] = vectorVal[j];
                               //std::cout<<"Src_Matrix["<<ii<<"]["<<j<<"]= "<<Src_Matrix[ii][j]<<std::endl;
                             }
                         }
                       else
                         nodalval[ii]=srcVal;
                     }
                 }

               if (vortexFlag_==3 || vortexFlag_==4) //Computing element vector using dTijdi from VortexAnalytical
                 {
                   ptEl=elemssd[j]->ptElem;
                   LinearFlowNoiseInt * linear_load = new LinearFlowNoiseInt(ptEl);
               
                   //To test implementation of CalcElemVec_withVortexVel()
                   //                Matrix<Double>  ptCoordtotest;
                   //                Matrix<Double>  NodalVel;
                   //                ptCoordtotest.Resize(2,4);
                   //                NodalVel.Resize(2,4);
                   //                ptCoordtotest[0][0] =  -1;
                   //                ptCoordtotest[1][0] =  -1;
                   //                ptCoordtotest[0][1] = 1;
                   //                ptCoordtotest[1][1] =  -1;
                   //                ptCoordtotest[0][2] = 1;
                   //                ptCoordtotest[1][2] = 1;
                   //                ptCoordtotest[0][3] = -1;
                   //                ptCoordtotest[1][3] = 1;
   
                   //                NodalVel[0][0] = 10.0;
                   //                NodalVel[0][1] = 8.5;
                   //                NodalVel[0][2] = 9.0;
                   //                NodalVel[0][3] = 8.0;
                   //                NodalVel[1][0] = 0.0;
                   //                NodalVel[1][1] = 1.0;
                   //                NodalVel[1][2] = 2.0;
                   //                NodalVel[1][3] = 1.5;
      
                   if (vortexFlag_==3)
                     linear_load->CalcElemVec_withdTijdi(ptCoordNodes,Src_Matrix, elemvec);
                   else
                     {
                       //To test with test vels:
                       //linear_load->CalcElemVec_withVortexVel(ptCoordtotest,NodalVel, elemvec);
                       // Original
                       linear_load->CalcElemVec_withVortexVel(ptCoordNodes,Src_Matrix, elemvec);
                     }
      
                   //                for (ii=0; ii<elsize; ii++)
                   //                  {
                   //                    if(((ptCoordNodes[0][ii])*(ptCoordNodes[0][ii])+
                   //                        (ptCoordNodes[1][ii])*(ptCoordNodes[1][ii]))<=((400./81.)*(400./81.)))
                   //                      {
                   //                        elemvec[ii]=0;
                   //                        std::cout<<"Making elemvec["<<ii<<"] equal 0"<<std::endl;
                   //                      }
                   //                  }

                   eqnMap_->GetNodeEqn(connecth, connect_PDE);

                   //std::cout<<"Elem vector obtained using Src_Matrix from VortexAnalytical: "<<std::endl;
                   //std::cout<<elemvec<<std::endl;
                  
                   algsys_->SetElementRHS(&elemvec[0], pdeId_, 
                                          connect_PDE.GetPointer(), connect_PDE.GetSize());
                  
                   delete linear_load;
                 }
               else
                 {
                   // Here the RHS integration of the nodal forces takes place
                   Vector<Double>  Sf;
                   ptEl=elemssd[j]->ptElem;
                   const UInt nrIntPts= ptEl->GetNumIntPoints();
                   const Vector<Double> & intWeights = ptEl->GetIntWeights();
                   Vector<Double> valueAtIP;
                   valueAtIP.Resize(nrIntPts);
                      
                   for (UInt actIntPt=1; actIntPt<=nrIntPts; actIntPt++)
                     {
                       ptEl->GetShFncAtIp(Sf, actIntPt, elemssd[j]);
                       valueAtIP[actIntPt-1]= nodalval*Sf;
                     }
                   Double jacDet;
                   for (UInt actIntPt=1; actIntPt<=nrIntPts;  actIntPt++) {
                        
                     jacDet = ptEl->CalcJacobianDetAtIp(actIntPt, ptCoordNodes,
                                                        elemssd[j]);
                     ptEl -> GetShFncAtIp(Sf, actIntPt, elemssd[j]);
                        
                     if (isaxi_) {
                       Vector<Double> coordAtIP;
                       coordAtIP = ptCoordNodes * Sf;
                       elemvec += Sf * (coordAtIP[0] * valueAtIP[actIntPt-1]* intWeights[actIntPt-1]
                                        * 2 * PI * jacDet);
                     }
                     else {
                       elemvec += Sf * intWeights[actIntPt-1] 
                         * valueAtIP[actIntPt-1]*jacDet;
                     }
                   }
                   eqnMap_->GetNodeEqn(connecth, connect_PDE);
                   algsys_->SetElementRHS(&elemvec[0], pdeId_, 
                                          connect_PDE.GetPointer(), connect_PDE.GetSize());
                 }
             }
           //ALSO TO FINISH THE OUTPUT OF CONTROL DATA 
           std::cerr <<std::endl;
         }
     }
}

  void AcouFlowNoise::ReadFlowData(const Char * aname, UInt timestep,
                                   Matrix<Double> &flowdata_)
  {
    ENTER_FCN( "AcouFlowNoise::ReadFlowData" );

    std::ifstream flowdatafile;
    UInt maxnumnodes;
    UInt maxnumqtts;
    std::ofstream testflowf;
    Char * aux=new Char[2];
    Char * anameloc=new Char[30];

    strcpy(anameloc,aname);     

    sprintf(aux,"%i",timestep);

    //   if (timestep/10 < 1) strcat(anameloc,"0");
     
    strcat(anameloc,aux);
    strcat(anameloc,".dat");
    std::cout<<"name of currentfile: "<<anameloc<<std::endl;
  
    /* Open as a binary file */
    //flowdatafile.open(aname, std::ios::in | std::ios::binary); 

    flowdatafile.open(anameloc,std::ios::in);
    if (!flowdatafile) 
      {
        std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
          ") Can't open Time-FlowSrc-File:" << anameloc << std::endl;
        exit(1);
      }

    std::string buffer, buffer2;

    // for maxnumqtts, currently for 2D we are working with 3, p and vx, vy

    UInt ibuf;

    /* Set pointer to beginning of file: */
    std::string::size_type pos=0;
    flowdatafile.seekg(pos, std::ios::beg);

    flowdatafile >> ibuf >> maxnumqtts >> maxnumnodes >> ibuf;

    flowdata_.Resize(maxnumqtts,maxnumnodes);
    flowdata_.Init();
    UInt i,j;

    UInt buffernodenum = 0;
    for (i=0; i < maxnumnodes; i++)
      {
        flowdatafile >> buffernodenum;
        for (j=0; j < maxnumqtts; j++)
          flowdatafile >> flowdata_[j][buffernodenum-1];
      }

    flowdatafile.close();


  }

  void AcouFlowNoise::VortexAnalytical(Double &srcVal, Vector<Double>& dTij_di, const Double x,
                                       const Double y, const Double t, 
                                       const UInt outType)
  {
    ENTER_FCN( "AcouFlowNoise::VortexAnalytical" );

    // From Vortex Implementation of HydSol by username iagschwa

    // local variable declaration

    Double      tPhys, XPhys, YPhys, arg;
    Double      p0Phys, rho0Phys, r_0Phys, GammaZirkPhys;
    Double      c0Phys, MaRot, iMaRot, qPhys, omegaPhys, akkPhys;
    Double      Equationgamma; // before it was EQN%gamma
     
    Double      U_inc;                         // Variablen, Stroemungs-
    Double      V_inc;                         // groessen aus Potential-
    Double      P_inc;                         // theorie
    Double      P_ak;                          // Akust. Druck aus MAE
    Double      Phi_t, Phi_tt, Phi_ttt;        // und seine Ableitungen

    Double      U_x, V_x, P_x;
    Double      U_y, V_y, P_y;
    Double      U_t, V_t, P_t;
    
    Double      U_xx, V_xx, P_xx;
    Double      U_xy, V_xy;
    Double      U_yy, V_yy, P_yy;
    Double      U_tt, V_tt, P_tt;
    
    Double      U_xt, V_xt, P_xt;
    Double      U_yt, V_yt, P_yt;
    
    Double      r, theta;                      // Polarkoordinaten r, theta
    //std::complex <double>    z, b;                          // Kompl.Koordinate z, 
    // Wirbelbahn b
    //std::complex<double>    zq_bq;                         // z**2-b**2
    //std::complex<double>    Im;                            // i defined below
    //std::complex<double>    w;                             // Komplexe Potential-
    // funktion
    //std::complex<double>    w_z, w_t, w_zz, w_zt;          // und ihre Ableitungen
    //std::complex<double>    w_zzz, w_tt, w_ttt, w_ztt;
    //std::complex<double>    w_zzt;
    //std::complex<double>    H2_2;
    // dimensionsbehaftete Var. die noch nicht in common definiert sind
    Double     KappaPhys, T_umlaufPhys, lambdaPhys;                    

    //--------------------------------------------------------------------------
    // ...Phys kennzeichnet die dimensionsbehafteten Groessen
    // ...     (ohne Endung) sind entdimensionalisierte Groessen
    //--------------------------------------------------------------------------
    //
    // Gegeben sind die folgenden Groessen :         
    //                                               
    p0Phys        = 0.714285714;                     // Umgeb.druck  [kg/(m*s^2)]
    rho0Phys      = 1.;                              // Umgeb.dichte     [kg/m^3]
    r_0Phys       = 1.;                              // halb. Wirbelabst. [  m  ]
    GammaZirkPhys = 1.00531;                         // Zirkulation       [m^2/s]
    Equationgamma = 1.4;                             // Need to be defined here as well


                                                                             //  Value
    c0Phys        = sqrt(Equationgamma*p0Phys/rho0Phys);// Schallgeschw.[ m/s ]    1.0 
    MaRot         = GammaZirkPhys/4./PI /r_0Phys/c0Phys;// Rotations-Mach-Zahl     0.08
    iMaRot        = 1. / MaRot;                         // 1/Ma_rot               12.5 
    
    qPhys         = MaRot * c0Phys;                  // induz. Geschw.    [ m/s ]  0.08
    KappaPhys     = GammaZirkPhys / (2.*PI);         // Wirbelstaerke pro          0.16
                                                     // Wirbel            [m^2/s]  
    T_umlaufPhys  = 2. * PI * r_0Phys / qPhys;       // Zeit fuer einen           78.53
                                                     // Wirbelumlauf      [  s  ]
    omegaPhys     = 2. * PI / T_umlaufPhys;          // Winkelgeschw.des           0.08
                                                     // Wirbel            [ 1/s ] 
    akkPhys       = 2. * omegaPhys / c0Phys;         // Argument,Berechnung        0.16
                                                     // akust.Druck,MAE   [ 1/m ]
    lambdaPhys    = c0Phys * T_umlaufPhys;           // Wellenlaenge      [  m  ] 78.53 

    std::complex<Double> Im(0., 1.);
      
    tPhys = t;
    XPhys = x;
    YPhys = y;

    //--------------------------------------------------------------------------
    // Vortex computation with dimensions
    //--------------------------------------------------------------------------
    //
    r     = sqrt(XPhys*XPhys + YPhys*YPhys);             // Komplexe Koordinate
    std::complex<Double> z( XPhys       , YPhys      );  // des Wirbelpaars

    arg   = omegaPhys*tPhys;
    std::complex<Double> b( r_0Phys * cos(arg) ,           
              r_0Phys * sin(arg)  );        // Komplexe Bahn des Wirbelpaars

    std::complex<double> zq_bq(pow(z, 2) - pow(b, 2));    // z**2-b**2
    
    // Complex potential function
                                  
    std::complex<double> w(GammaZirkPhys/(2.*PI*Im) * log(zq_bq));
    
    // Ableitungen nach z, 1. Ordnung
    //std::complex<double> w_z;
    Double r_c=0.10;
                         
    std::complex<double> w_z(GammaZirkPhys*z/(PI*Im*zq_bq));

    Double r_limit;
    r_limit=0.15;
    Double x_vort, y_vort;

    //Vortex core model only in a region with rad=r_limit around the two votices
    if ( ( abs(z-b) <=r_limit ) || ( abs(z+b) <=r_limit ) )
      {
        if  ( abs(z-b) <=r_limit ) 
          {
            y_vort = YPhys-imag(b);
            x_vort = XPhys-real(b);
          }
        if ( abs(z+b) <=r_limit ) 
          {
            y_vort = YPhys+imag(b);
            x_vort = XPhys+real(b);          
          }
        U_inc= 1.35* -GammaZirkPhys/(2*PI)* y_vort / (r_c*r_c  + y_vort*y_vort + x_vort*x_vort);
        V_inc= 1.35* GammaZirkPhys/(2*PI)* x_vort / (r_c*r_c  + y_vort*y_vort + x_vort*x_vort);
      }
    else
      {
        //Original
        U_inc  =    std::real(w_z);                    
        V_inc  = - std::imag(w_z);
      }
    
        
    // Ableitungen nach z, 2. Ordnung                              
    std::complex<double> w_zz(-GammaZirkPhys*(pow(z,2) + pow(b,2)) / (PI*Im*pow(zq_bq,2)));

    U_x   =   real(w_zz);
    V_x   = -imag(w_zz);
    U_y   =  V_x;
    V_y   = -U_x;
    
    // Ableitungen nach z, 3. Ordnung
                                  
    std::complex<double> w_zzz(-GammaZirkPhys*(-2.*pow(z,3) - 6.*z*pow(b,2))/(PI*Im*pow(zq_bq,3)));
    U_xx =  real( w_zzz);
    U_xy = -imag(w_zzz);
    
    V_yy =  imag(w_zzz);
    V_xy = -real( w_zzz);
    
    U_yy =  V_xy;
    V_xx =  U_xy;
    
    // Ableitungen nach z, nach t, 2. Ordnung
                                  
    std::complex<double> w_zt(2.* omegaPhys* GammaZirkPhys* z* pow((b),2) / (PI*pow((zq_bq),2)));
    U_t  =  real( w_zt);
    V_t  = -imag(w_zt);
    
    // Ableitungen nach t, 1. und 2. und 3. Ordnung
                                  
    std::complex<double> w_t(-omegaPhys*GammaZirkPhys * pow(b,2) / ( PI*zq_bq ));
    Phi_t   = real(w_t);
    
    std::complex<double> w_tt(-(2.*Im*pow(omegaPhys,2)*GammaZirkPhys/PI)*pow(b,2)*pow(z,2)/pow(zq_bq,2));
    Phi_tt  = real(w_tt);
    
    std::complex<double> w_ttt( (4*pow(omegaPhys,3)*GammaZirkPhys/PI)*pow(b,2)*pow(z,2)*(pow(z,2)+pow(b,2))/pow(zq_bq,3));
    Phi_ttt = real(w_ttt);
    
    // Ableitungen nach z, nach t, Gesamtordnung 3
                                  
    std::complex<double> w_ztt( (4.*Im*pow(omegaPhys,2)*GammaZirkPhys/PI)*(pow(b,2)*pow(z,3) + z*pow(b,4))/pow(zq_bq,3));
    U_tt =  real( w_ztt);
    V_tt = -imag(w_ztt);
    
    std::complex<double> w_zzt( -(2.*omegaPhys*GammaZirkPhys/PI)*(3.*(pow(b,2))*(pow(z,2)) + pow(b,4))/pow(zq_bq,3));
    U_xt =  real(w_zzt);
    U_yt =  real(w_zzt*Im);
    V_xt = -imag(w_zzt);
    V_yt = -U_xt;
    
    P_inc = p0Phys + rho0Phys*( -Phi_t  -(pow(U_inc,2)+pow(V_inc,2))/2. );
    P_t   = -rho0Phys*( Phi_tt      + U_inc*U_t + V_inc*V_t );
    P_x   = -rho0Phys*( U_t    + U_inc*U_x + V_inc*V_x );
    P_y   = -rho0Phys*( V_t    + U_inc*U_y + V_inc*V_y );
    
    P_tt  = -rho0Phys*( Phi_ttt + pow(U_t,2) + U_inc*U_tt + pow(V_t,2) + V_inc*V_tt );
    
    P_xt  = -rho0Phys*( U_tt + U_t*U_x + U_inc*U_xt + V_t*V_x + V_inc*V_xt );
    
    P_yt  = -rho0Phys*( V_tt + U_t*U_y + U_inc*V_xt + V_t*V_y - V_inc*U_xt );
    P_xx  = -rho0Phys*( U_xt + pow(U_x,2)  + U_inc*U_xx + pow(V_x,2)  + V_inc*V_xx );
    P_yy  = -rho0Phys*( V_yt + pow(U_y,2)  + U_inc*U_yy + pow(V_y,2)  + V_inc*V_yy );
    

    switch( outType ) 
      {
        //-----------------Berechne den akkustischen Druck--------------------------
      case 1: 
        {
          // Wenn outAkk gleich 1, dann soll der
          // akkustische Druck ausgegeben werden
        
              // Akustischer Druck
              if (r > 0.5)
                { 
                  if (abs(x)<1.e-18) 
                    {
                      std::cout<<"In if for converting Theta in analytical MAE"<<std::endl;
                      std::cout<<"x= "<<x<<std::endl;
                      std::cout<<"abs(x)= "<<abs(x)<<std::endl;
                      
                      if (y>0.0)
                        {
                          theta = PI/2.;

                        }
                      else
                        {
                          theta = 3.*PI/2.;
                        }
                    }
                  else
                    {
                      theta = atan(y/x);   // Polarkoordinaten r, theta
                    }
            
                  // Hankel Funktion 2nd order and 2 kind H2_2
                  // with Bessel function 2nd order, 1st kind JJ2 and 2nd kind YY2
            
                  std::complex<Double> H2_2( JJ2(akkPhys*r), -YY2(akkPhys*r) ) ;
            
                  P_ak = std::real( (-Im * rho0Phys * (GammaZirkPhys*GammaZirkPhys*GammaZirkPhys
                                                       *GammaZirkPhys)    * H2_2
                                     * exp(Im * 2. * (omegaPhys * tPhys-theta)))   
                                    /(64. * (PI*PI*PI) * (r_0Phys*r_0Phys*
                                                          r_0Phys*r_0Phys) * (c0Phys*c0Phys))          );
                }
              else
                {
                  P_ak = 0.;
                }
              srcVal=P_ak;
          break;
        }
      case 2:
        {
          P_ak = 0.;
          srcVal = -P_tt;
          break;
        }
      case 3:
        {
          //Calculation of Lighthill's RHS term
          //Using weak formulation of Lighthill's quadrupole term
          //
          //This is only a 2D implementation!!
          // Sending back the complete divergence of T
          dTij_di.Resize(2);
          dTij_di[0] = 2.0 * U_inc * U_x + V_inc * U_y + U_inc * V_y; 
          dTij_di[1] = 2.0 * V_inc * V_y + U_inc * V_x + V_inc * U_x; 
          Double density=1.0;
          dTij_di*= density;
          break;
        }
      case 4:
        {
          // Sending back only the velocity components for evaluation of
          // dTijdi in linearForm.cc
          dTij_di.Resize(2);   
          dTij_di[0] = U_inc; 
          dTij_di[1] = V_inc;
          break;
        }
      case 5:
        {// Sending back the tangential velocity as srcVal
          P_ak = 0.;
          //srcVal=omegaPhys*r;
          srcVal=sqrt(U_inc*U_inc+V_inc*V_inc);
          break;
        }
      case 6:
        {// Sending back the velocity field for plotting it as RHS
          srcVal=sqrt(U_inc*U_inc+V_inc*V_inc);
          P_ak = 0.;
          //srcVal=omegaPhys*r;
          dTij_di.Resize(2);   
          dTij_di[0] = U_inc; 
          dTij_di[1] = V_inc; 
          break;
        }
      case 7:
        {// Sending back the gradient field for plotting it as RHS
          srcVal=sqrt(U_inc*U_inc+V_inc*V_inc);
          P_ak = 0.;
          //srcVal=omegaPhys*r;
          dTij_di.Resize(2);   
          if (((x<(cos(2*arg)+0.25)) && (y<(sin(2*arg)+0.25)) && (x>(cos(2*arg)-0.25)) && (y>(sin(2*arg)-0.25))) ||
              ((x<(-cos(2*arg)+0.25)) && (y<(-sin(2*arg)+0.25)) && (x>(-sin(2*arg)-0.25)) && (y>(-sin(2*arg)-0.25))))
            {
              std::cout<<"aqui primera excepcion, "<<"Real b: "<<real(b)+0.25<<"Imag b: "<<imag(b)+0.25<<std::endl;
              dTij_di[0] = 0.; 
              dTij_di[1] = 0.; 
            }
          else
            {
              dTij_di[0] = 2.0 * U_inc * U_x + V_inc * U_y + U_inc * V_y; 
              dTij_di[1] = 2.0 * V_inc * V_y + U_inc * V_x + V_inc * U_x;           
            }
          break;
        }
      default:
        EXCEPTION( "Need to give an outType value between 1 and 7 "
                   << "Current outType value: " << outType );
      }
    

  }
  
  //___________________________________
  //
  // Bessel function 2nd order and 1st kind
  //___________________________________
  
  Double AcouFlowNoise::JJ2(const Double x)
  {
    
    //--------------------------------------------------------------------------
    
    Double JJ2=0;
    Double RESULT=0;
    //--------------------------------------------------------------------------
    Double PI=3.1415926535;
    //--------------------------------------------------------------------------
    
    if (x<2.)
      {
        RESULT = 0.5*(0.5*x)*(0.5*x);
      }
    
    if (x>=2.) 
      {
        RESULT = sqrt(2./PI/x)*cos(x-1.25*PI);
      }
    
    JJ2 =RESULT;
    
    return JJ2;
    
      }
  
  //___________________________________
  //
  // Bessel function 2nd order and 2nd kind 
  //___________________________________

  Double AcouFlowNoise::YY2(const Double x)
  {
    //--------------------------------------------------------------------------
    Double YY2=0;
    Double RESULT=0;
    Double PI=3.1415926535;
    
    if (x<2.) 
      {
        RESULT = -1./PI*sqrt(0.5*x);
      }
    
    if (x>=2.)
      {
        RESULT = sqrt(2./PI/x)*sin(x-1.25*PI);
      }
    
    YY2 =RESULT;
    
    return YY2;
    
  }
  
} // end of namespace
