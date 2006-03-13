#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
//#include <stdio.h>
#include <complex>

#include "acouFlowNoise.hh"

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/freqfunc.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "MpCCIcpl/MpCCIexch.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Driver/solveStepAcouFlowNoise.hh"

#ifdef MpCCI
#include <cci.h>
#endif

namespace CoupledField
{

  AcouFlowNoise::AcouFlowNoise(Grid *aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut)
    :AcousticPDE(aptgrid, aptTimeFunc, aptOut)
  {
    ENTER_FCN( "AcouFlowNoise::AcouFlowNoise" );

    // pdename_ is also acoustic for this case
    pdename_          = "acoustic";
    pdematerialclass_ = "fluid";
    nodalSrc_ = FALSE;
    vortexSrc_ = FALSE;
    plotRHS_ = FALSE; 
    plotRHSVel_= FALSE;
    isHarmonic_=FALSE;
    vortexFlag_=0;

    if( params->IsSet( "valRHS","pdeList" ,"acoustic" ) ) {
      plotRHS_ = TRUE;
      Info->PrintF(pdename_, "Writing RHS as solution of problem\n" );
    }

    StdVector<Elem*> elemssd;
    StdVector<std::string> regionNames, coupledRegionNames;

    //For writing fine vortex flow field in files for later MpCCI exchange
      writeGridFile_ = FALSE;
      writeSrcFileperTS_ = FALSE;
    if( params->IsSet( "writeCoupledGrid","pdeList" ,"acoustic" ) ) {
      writeGridFile_ = TRUE;
      Info->PrintF("acoustic","Writing grid def. file of coupled domain\n");
    }
    if( params->IsSet( "writeSrcFileperTS","pdeList" ,"acoustic" ) ) {
      writeSrcFileperTS_ = TRUE;
      Info->PrintF("acoustic","Writing coarse sources of coupled domain in time files (NrFiles=NrTimeSteps)\n");
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
    if( params->HasValue( "type", "nodalSrc", pdename_, "flowData" ) ) {
      nodalSrc_ = TRUE;
      Info->PrintF(pdename_, " Using FlowData as RHS nodal source\n" );
    }

    params->GetList( "name", regionNames, pdename_, "region" );
    params->GetList( "name", coupledRegionNames,  "MpCCI-flownoise",
                     "coupledregion" );
    ptgrid_->RegionNameToId( subdoms_, regionNames );
    ptgrid_->RegionNameToId( couplSubDomId_, coupledRegionNames );

    Boolean OnlyLinNodes=TRUE;
    Boolean Find=FALSE;
    //    std::cerr << "\nName of all regions:\n" << regionNames << std::endl;
    //    std::cerr << "Coupled regions:\n" << coupledRegionNames  << std::endl;
    

    for (UInt i=0; i<subdoms_.GetSize(); i++)
      {
        for (UInt j=0; j<couplSubDomId_.GetSize(); j++)
          {
            if (couplSubDomId_[j] == subdoms_[i])
              {

                ptgrid_->GetVolElems(elemssd,couplSubDomId_[j]);
		
                if (dim_ == 3)
                  { //TO SEND ONLY LINEAR TETRAS TO MPCCI SINCE NO QUAD TETRAS
                    //ARE ALLOWED
                    ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, FALSE);
                    ptgrid_->GetNodesOfElemList(mapSD_onlyLinNodes_, elemssd, OnlyLinNodes);
                    ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_onlyLinNodes_);
                  }
                else
                  {
                    ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, FALSE);
                    //ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_.GetSize() );
                    ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_allNodes_);
                  }
		
                Find=TRUE;

              }
          }
      }
    if (!Find) 
      {
        std::string msg="Subdom to be coupled with MpCCI ";
        msg+="is not in list of PDE subdoms. Please, check .xml-file";
        Error(msg.c_str(),__FILE__,__LINE__);
      }
    // On CFS side we always used the complete coupled region and later give a zero src to the
    // quadratic (or mid) nodes in case of 3D quadratric acoustic grid.
    MpCCInodes_=mapSD_allNodes_.GetSize();
    
    
    ptMpCCIexch_->PutExchangeGrid2MpCCI(couplSubDomId_);
#else
    if( params->HasValue( "type", "vortexSrc", pdename_, "flowData" ) ) {
      vortexSrc_ = TRUE;
      params->GetList( "name", regionNames, pdename_, "region" );
      params->GetList( "name", coupledRegionNames,  "MpCCI-flownoise",
                       "coupledregion" );
      ptgrid_->RegionNameToId( subdoms_, regionNames );
      ptgrid_->RegionNameToId( couplSubDomId_, coupledRegionNames );
      std::cerr << "\nSubdoms:\n" << regionNames << std::endl;
      std::cerr << "CoupledSubdoms:\n" << coupledRegionNames << std::endl;
      Info->PrintF(pdename_, " Using Vortex function as RHS nodal source\n" );
    }
    else 
      {
        if( params->HasValue( "type", "nodalSrc", pdename_, "flowData" ) ) {
          nodalSrc_ = TRUE;
          // Now verify that the type of analysis is HARMONIC
          // Determine type of analysis
          std::string analysis;
          AnalysisType analysisType;
          params->Get( "type", analysis, "analysis" );
          String2Enum( analysis, analysisType );
          if ( analysisType==HARMONIC )
            {
              isHarmonic_=TRUE;
              ComputeRHSforHarm_=TRUE;    
              Info->PrintF(pdename_, "Using FlowData from dataset as RHS nodal source\n" );
              Info->PrintF(pdename_, "Computing using nodal frequency files (No MpCCI used)\n" );



              // Now getting the correct size of coupled domain. It must have the same
              // numbering as the nodal frequency file, if not freqfunc will complain!!

              params->GetList( "name", regionNames, pdename_, "region" );
              params->GetList( "name", coupledRegionNames, pdename_, "coupledRegion" );
              ptgrid_->RegionNameToId( subdoms_, regionNames );
              ptgrid_->RegionNameToId( couplSubDomId_, coupledRegionNames );

              Boolean Find=FALSE;
              //              std::cerr << "Name of all regions:\n" << regionNames << std::endl;
              //              std::cerr << "Coupled regions:\n" << coupledRegionNames  << std::endl;
    

              for (UInt i=0; i<subdoms_.GetSize(); i++)
                {
                  for (UInt j=0; j<couplSubDomId_.GetSize(); j++)
                    {
                      if (couplSubDomId_[j] == subdoms_[i])
                        {
                          ptgrid_->GetVolElems(elemssd,couplSubDomId_[j]);
                          ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, FALSE);
                          //TRUE returns only the corner nodes (for 3D MpCCI linear coupling)
                          if (dim_==3)
                            ptgrid_->GetNodesOfElemList(mapSD_onlyLinNodes_, elemssd, TRUE);
                          Find=TRUE;
                        }
                    }
                }
              if (!Find) 
                {
                  std::string msg="Subdom to be coupled with NodalFreqFilesSet ";
                  msg+="is not in list of PDE subdoms. Please, check .xml-file";
                  Error(msg.c_str(),__FILE__,__LINE__);
                }
              Integer CoupledNodes;
              if (dim_==3)
                CoupledNodes=mapSD_onlyLinNodes_.GetSize();
              else 
                CoupledNodes=mapSD_allNodes_.GetSize();
              flowdata_.Resize(1, CoupledNodes);
            }
          else
            {
              Error( "Only Nodal Frequency Sources or Vortex Source are allowed if not using MpCCI"
                     , __FILE__, __LINE__ );
            }
        }
        
      }


    


    //     if (vortexSrc_)
    //       {
    // //         StdVector<Elem*> elemssd;
    // //         StdVector<std::string> regionNames, coupledRegionNames;
    // //         params->GetList( "name", regionNames, pdename_, "region" );
    // //         params->GetList( "name", coupledRegionNames,  "MpCCI-flownoise",
    // //                          "coupledregion" );
    // //         ptgrid_->RegionNameToId( subdoms_, regionNames );
    // //         ptgrid_->RegionNameToId( couplSubDomId_, coupledRegionNames );
    // //         std::cerr << "Subdoms:\n" << subdoms_ << std::endl;
    // //         std::cerr << "CoupledSubdoms:\n" << couplSubDomId_  << std::endl;
    //       }
    //     else
    //       Error( "No MpCCI coupling, need to define VortexSrc", __FILE__, __LINE__ );    


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

    solveStep_ = new SolveStepAcouFlowNoise(*this);
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
    Double valmult;

    static UInt timestep=1;

    std::cout<<"timestep counter in ComputeRHS: "<<timestep<<std::endl;

#ifdef MpCCI
    double starttime, endtime;
    //std::cout<<"MpCCInodes_: "<< MpCCInodes_ << " dimension: " 
    //             << dim_ << std::endl;
    if (nodalSrc_ == TRUE)
      //we get already the integrated acoustic source term
      flowdata_.Resize(1, MpCCInodes_);
    else
      flowdata_.Resize(1+dim_, MpCCInodes_);

    starttime = CCI_Wtime();
      
    Integer timestep = 0;
    ptMpCCIexch_->CouplCompPhase(flowdata_, atime);


    endtime = CCI_Wtime();
      

    //    std::cout<<"Transfer of Data CouplCompPhase() for 1 time step took: "
    //         <<(endtime-starttime)<<" seconds"<<std::endl;

    // Quadrupole computation
    std::cout << "Receiving "<<MpCCInodes_<<" acoustic nodal sources for time t = "<< atime<<"s"<< std::endl;


#else
    // If data from fluid file, call to get fluid flow data in flowdata_  
    std::string flowdata;
    //    Warning( "Using no MpCCI, vortexSrc is assumed.", __FILE__, __LINE__ );
      
//     conf->get("acousrc_file",flowdata);
//     // ReadFlowData(flowdata.c_str(), timestep, flowdata_); 

#endif 


// Variables for ramping
    Double xfmin, yfmin, zfmin, xfmax, yfmax, zfmax, facRampXmin, facRampYmin, 
      facRampZmin, facRampXmax, facRampYmax, facRampZmax;
    Double bndoffsetXmin, bndoffsetYmin, bndoffsetZmin, bndoffsetXmax, bndoffsetYmax, bndoffsetZmax ;
#ifdef MpCCI
    params->Get("xfmin",xfmin, "MpCCI-flownoise");
    params->Get("yfmin",yfmin, "MpCCI-flownoise");
    params->Get("xfmax",xfmax, "MpCCI-flownoise");
    params->Get("yfmax",yfmax, "MpCCI-flownoise");
    params->Get("facrampXmin",facRampXmin, "MpCCI-flownoise");
    params->Get("facrampYmin",facRampYmin, "MpCCI-flownoise");
    params->Get("facrampXmax",facRampXmax, "MpCCI-flownoise");
    params->Get("facrampYmax",facRampYmax, "MpCCI-flownoise");

    bndoffsetXmin=facRampXmin*xfmin;
    bndoffsetYmin=facRampYmin*yfmin; 
    bndoffsetXmax=facRampXmax*xfmax;
    bndoffsetYmax=facRampYmax*yfmax;

    if (dim_==3)
      {
    params->Get("zfmin",zfmin, "MpCCI-flownoise");
    params->Get("zfmax",zfmax, "MpCCI-flownoise");
    params->Get("facrampZmin",facRampZmin, "MpCCI-flownoise");
    params->Get("facrampZmax",facRampZmax, "MpCCI-flownoise");
    bndoffsetZmin=facRampZmin*zfmin;
    bndoffsetZmax=facRampZmax*zfmax;
      }
#endif 

    // Correct valmult value is -1.0, 
    // if plugging in source (ddTij/dxidxj) directly in weak form then 1.0
#ifdef MpCCI
    starttime = CCI_Wtime();
#endif 
    if (nodalSrc_ == FALSE) {
      
      if (vortexSrc_ == TRUE)
        {
          std::cout<<"Calling ComputeRHSwithVortexSource(atime)"<<std::endl;
          ComputeRHSwithVortexSource(atime);
        }
      else
        {
          valmult=-1.0;
      
          for (i=0; i<couplSubDomId_.GetSize(); i++)
            {
              ptgrid_->GetVolElems(elemssd,couplSubDomId_[i]);
              
              for (j=0; j< elemssd.GetSize(); j++)
                {
                  ptEl=elemssd[j]->ptElem;
                  BaseForm * linear_load = new LinearFlowNoiseInt(ptEl);
                  
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
                  //                   if (ptCoordNodes[0][ii]<bndoffsetXmin)
                  //                     {
                  //                       elemvec[ii]-=elemvec[ii]*
                  //                   (ptCoordNodes[0][ii]-bndoffsetXmin)/(xfmin-bndoffsetXmin);
                  //                     }
                  
                  //                   else
                  //                     if (ptCoordNodes[0][ii]>bndoffsetXmax)
                  //                       elemvec[ii]-=elemvec[ii]*
                  //                   (ptCoordNodes[0][ii]-bndoffsetXmax)/(xfmax-bndoffsetXmax);
                  //                   if (ptCoordNodes[1][ii]<bndoffsetYmin)
                  //                     elemvec[ii]-=elemvec[ii]*
                  //                   (ptCoordNodes[1][ii]-bndoffsetYmin)/(yfmin-bndoffsetYmin);
                  //                   else    
                  //                     if (ptCoordNodes[1][ii]>bndoffsetYmax)
                  //                       elemvec[ii]-=elemvec[ii]*
                  //                   (ptCoordNodes[1][ii]-bndoffsetYmax)/(yfmax-bndoffsetYmax);
                  //                 }
                  
                  //end ramping
                  
                  // CHANGE connecth
                  //Mesh2PDENode(connect_PDE,connecth,mesh2PDENode_);
                  eqnData_->Node2EQN(connecth, connect_PDE);
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
    
    else {     //for assigning nodal source
      UInt eqnDof, node, dof;
      Integer eqnNr;
      StdVector<UInt> connect(1);
 
      valmult=-1.0;    
     
      eqnDof = 1;
      dof    = 1;

      if (!isHarmonic_)
        {
          for (UInt idx=0; idx<flowdata_.GetSizeCol() ; idx++)
            {
              Double val = flowdata_[0][idx];
              node = idx + 1;

              ///////////////////////////////////////////////////////////////////////
                // Ramping before adding to RHS vector (NOW IT MAKES ZERO THE RHS ENTRY!)
                Matrix<Double> ptCoordNodes;
                connecth.Resize(1);
                connecth[0] = node;
                ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);         

                if (ptCoordNodes[0][0]<bndoffsetXmin)
                  //val -= val*(ptCoordNodes[0][0]-bndoffsetXmin)/(xfmin-bndoffsetXmin);
                  val = 0;
                else
                  if (ptCoordNodes[0][0]>bndoffsetXmax)
                    //val -= val*(ptCoordNodes[0][0]-bndoffsetXmax)/(xfmax-bndoffsetXmax);
                    val = 0;
                if (ptCoordNodes[1][0]<bndoffsetYmin)
                  //val -= val*(ptCoordNodes[1][0]-bndoffsetYmin)/(yfmin-bndoffsetYmin);
                  val = 0;
                else      
                  if (ptCoordNodes[1][0]>bndoffsetYmax)
                    //val -= val*(ptCoordNodes[1][0]-bndoffsetYmax)/(yfmax-bndoffsetYmax);
                    val = 0;
                if(dim_==3)
                  {
                    if (ptCoordNodes[2][0]<bndoffsetZmin)
                      //val -= val*(ptCoordNodes[2][0]-bndoffsetZmin)/(zfmin-bndoffsetZmin);
                      val = 0;
                    else      
                      if (ptCoordNodes[2][0]>bndoffsetZmax)
                        //val -= val*(ptCoordNodes[2][0]-bndoffsetZmax)/(zfmax-bndoffsetZmax);
                        val = 0;
                  }
                //end ramping
                /////////////////////////////////////////////////////////////////////////

                  val*=valmult;

                  //add to RHS
                  eqnData_->Node2EQN(node,dof,eqnNr,eqnDof);
                  algsys_->SetNodeRHS(val, pdeId_, eqnNr, eqnDof);      
            }
      
        } 
      else
        {
          std::cout << "Using frequency source files..."<< std::endl;
          FreqFunc * ptFreqFunc = new FreqFunc();
          StdVector<Double> Ampl_Phase;
          Ampl_Phase.Resize(2);

          //Getting current freq (received as atime)
          for (UInt idx=0; idx<flowdata_.GetSizeCol() ; idx++) {
          actFreq = atime;
            if (dim_==3)
              node = mapSD_onlyLinNodes_[idx];
            else
              node = mapSD_allNodes_[idx];    
          
              Ampl_Phase=ptFreqFunc->NodalFreqFuncAtFreq(actFreq,"freqsrcfile.node",
                                                         node);              
            Double valAmpl = Ampl_Phase[0];
            Double valPhase = Ampl_Phase[1];
#ifdef MpCCI            
              ///////////////////////////////////////////////////////////////////////
                // Ramping before adding to RHS vector (NOW IT MAKES ZERO THE RHS ENTRY!)
                Matrix<Double> ptCoordNodes;
                connecth.Resize(1);
                connecth[0] = node;
                ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);         

                if (ptCoordNodes[0][0]<bndoffsetXmin)
                  //valAmpl -= valAmpl*(ptCoordNodes[0][0]-bndoffsetXmin)/(xfmin-bndoffsetXmin);
                  valAmpl = 0;
                else
                  if (ptCoordNodes[0][0]>bndoffsetXmax)
                    //valAmpl -= valAmpl*(ptCoordNodes[0][0]-bndoffsetXmax)/(xfmax-bndoffsetXmax);
                    valAmpl = 0;
                if (ptCoordNodes[1][0]<bndoffsetYmin)
                  //valAmpl -= valAmpl*(ptCoordNodes[1][0]-bndoffsetYmin)/(yfmin-bndoffsetYmin);
                  valAmpl = 0;
                else      
                  if (ptCoordNodes[1][0]>bndoffsetYmax)
                    //valAmpl -= valAmpl*(ptCoordNodes[1][0]-bndoffsetYmax)/(yfmax-bndoffsetYmax);
                    valAmpl = 0;
                if(dim_==3)
                  {
                    if (ptCoordNodes[2][0]<bndoffsetZmin)
                      //valAmpl -= valAmpl*(ptCoordNodes[2][0]-bndoffsetZmin)/(zfmin-bndoffsetZmin);
                      valAmpl = 0;
                    else      
                      if (ptCoordNodes[2][0]>bndoffsetZmax)
                        //valAmpl -= valAmpl*(ptCoordNodes[2][0]-bndoffsetZmax)/(zfmax-bndoffsetZmax);
                        valAmpl = 0;
                  }
                //end ramping
                /////////////////////////////////////////////////////////////////////////
#endif
            valAmpl*=valmult;
            Complex complexValue( valAmpl * cos( valPhase / 180 * PI ),
                                  valAmpl * sin( valPhase / 180 * PI ) );
            //add to RHS
            //it was: eqnData_->Node2EQN(node,dof,eqnNr,eqnDof);
            //Now to try, since there is no sequential order of nodes:
            if (dim_ == 3)
              {
                eqnData_->Node2EQN(mapSD_onlyLinNodes_[idx],dof,eqnNr,eqnDof);
                //std::cout<<"EqnDOF= "<<eqnDof<<" EqnNr= "<<eqnNr<<" dof= "<<dof
                //         <<" NodeNumber= "<<mapSD_onlyLinNodes_[idx]<<" Ampl= "<<Ampl_Phase[0]<<std::endl;
              }
            else
              eqnData_->Node2EQN(mapSD_allNodes_[idx],dof,eqnNr,eqnDof);
            algsys_->SetNodeRHS(complexValue, pdeId_, eqnNr, eqnDof);  
          }
       }//end else for frequency analysis
    }//end else in case nodalSrc is TRUE
    
#ifdef MpCCI
    endtime = CCI_Wtime();
#endif    
  
    if (plotRHS_ && !isHarmonic_ && vortexFlag_!=6){
      ///////// For plotting the RHS as solution for analysing it
      std::cout<<"Init rhs_ in ComputeRHS(). vortexFlag="<<vortexFlag_<<std::endl;
      rhs_.SetNumSolutions(1);
      rhs_.SetNumNodes(numPDENodes_);
      rhs_.SetSolutionType(ACOU_RHSVAL);
      rhs_.SetNumDofs(1);
      rhs_.SetPtrEQNData(eqnData_, ptgrid_);
      rhs_.Init(0.0);
      
      Double *ptRHS;
      algsys_->GetRHSVal( ptRHS );
      rhs_.CopyFromAlgSysDataPointer(ptRHS);

      //For writing out the topology and RHSval (only for Vortex simulation)
      if (writeGridFile_ && writeSrcFileperTS_)
        {


            // static variable for suffix of output src file
            static Integer filenum=1;
            std::string filename;
            filename = "srcfile_ultrafine_";
            if ( filenum < 10 ) filename.append( "000" );
            else if ( filenum < 100 ) filename.append( "00" );
            else if ( filenum < 1000 ) filename.append( "0" );
            else if ( filenum > 10000 ) {
              Info->Error( "Number of src file exceeds 9999!",
                           __FILE__, __LINE__ );
            }
            filename.append( Info->GenStr( filenum ) );
            filename.append( ".dat" );
            filenum++;
            outsrcfile_ = new std::ofstream(filename.c_str());



          StdVector<Elem*> elemssd;  
          Point<2> ptNodalCoord2D;   
          for (UInt i=0; i<couplSubDomId_.GetSize(); i++)
            {
              UInt eqnDof, node, dof;
              Integer eqnNr;
              eqnDof = 1;
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
                      for (UInt ii=0; ii<elsize; ii++)
                        {
                          connecth[ii]=(elemssd[j]->connect)[ii];
                          if ((ii<(elsize-1)))
                            {
                              (*outelemfile_) << connecth[ii]<<" ";
                            }	      
                          else
                            {
                              if (writeGridFile_)
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


              ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, FALSE);
                               
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
                  eqnData_->Node2EQN(node,dof,eqnNr,eqnDof);
                  rhs_.Get(idx,1,RHSnodalVal);
                  (*outsrcfile_) <<  RHSnodalVal<<std::endl;
                  //std::cout<<"node: "<<node<<" eqnNr: "<<eqnNr<<" RHSnodalVal: "<<RHSnodalVal<<std::endl;
                }
                
              
            }
        }
    }
    
      
  
         timestep=timestep+1;

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
   params->Get("outType", vortexFlag_, "vortexSrc");

   for (UInt i=0; i<couplSubDomId_.GetSize(); i++)
     {
       std::cerr << "Current time: "<<atime;
       std::cerr <<std::endl;
       ptgrid_->GetVolElems(elemssd,couplSubDomId_[i]);
       if (vortexFlag_==1 || vortexFlag_==5 || vortexFlag_==6)
         //Getting Analytical solution (P_ak) or Tangential velocity (arg) as RHSval
         {
           std::cout<<"Getting Analytical solution (P_ak) or Tangential velocity (arg) as RHSval"<<std::endl;
           UInt eqnDof, node, dof;
           Integer eqnNr;
           StdVector<UInt> connect(1);
           eqnDof = 1;
           dof    = 1;
           if (vortexFlag_==6)
             {
               plotRHSVel_= TRUE;
               std::cout<<"Init rhs_ for putting vel vector field in RHSval for visualization"<<std::endl;
               rhs_.SetNumSolutions(1);
               rhs_.SetNumNodes(numPDENodes_);
               rhs_.SetSolutionType(ACOU_RHSVAL);
               rhs_.SetNumDofs(1);
               rhs_.SetPtrEQNData(eqnData_, ptgrid_);
               rhs_.Init(0.0);
               
               rhs2_.SetNumSolutions(1);
               rhs2_.SetNumNodes(numPDENodes_);
               rhs2_.SetSolutionType(ACOU_POT_NRBC);
               rhs2_.SetNumDofs(1);
               rhs2_.SetPtrEQNData(eqnData_, ptgrid_);
               rhs2_.Init(0.0);

             }

           ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, FALSE);
           for (UInt idx=0; idx<mapSD_allNodes_.GetSize() ; idx++)
             {
               Double val = mapSD_allNodes_[idx];
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
               tempVelX[0]=vectorVal[0];
               tempVelY[0]=vectorVal[1];
               
               //add to RHS
               eqnData_->Node2EQN(node,dof,eqnNr,eqnDof);
               algsys_->SetNodeRHS(val, pdeId_, eqnNr, eqnDof); 
                             
               if (vortexFlag_==6)
                 {
                   //std::cout<<"eqnNr: "<<eqnNr<<std::endl;
                   rhs_.SetNodalResult(eqnNr,tempVelX);
                   rhs2_.SetNodalResult(eqnNr,tempVelY);
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
               nodalval.Resize(elsize);
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
                  
                  
               for (ii=0; ii<elsize; ii++)
                 {
                   Double r_sqr=((ptCoordNodes[0][ii])*(ptCoordNodes[0][ii])+
                                 (ptCoordNodes[1][ii])*(ptCoordNodes[1][ii]));

                   //For square 100times100 domain                      
                   //                       if (r_sqr<=((1000./81.)*(1000./81.)))
                   //                         {
                   //                           nodalval[ii]=0;
                   //                         }
                      
                   //Original core (about 2.5m)  //Give also shifted results (above zero)
                   //if (r_sqr<=((200./81.)*(200./81.)))
                     //Smaller core source (about 1.5m core)  // Give better results
                   if (r_sqr<=((140./81.)*(140./81.))) 
                   //Bigger core (about 12m)
                   //if (r_sqr<=((1000./81.)*(1000./81.))) //Give shifted results (above zero)
                   //Full core   // Gives very high amplitudes
                   //      if (r_sqr<=-1)
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

                       //std::cout << "In else after calling VortexAnalytical, srcVal= "<<srcVal<<std::endl;
                       //                           //For debugging the quadratic middle nodes
                       //                           if ((j+1)==467)
                       //                           if (connecth[ii]==1708)
                       //                             {
                       //                               std::cerr <<"  "<< srcVal;
                       //                             }
                       //                           if ((j+1)==464)
                       //                           if (connecth[ii]==1356)
                       //                             {
                       //                               std::cerr <<"  "<< srcVal;
                       //                             }                          
                       //                           if ((j+1)==434)
                       //                           if (connecth[ii]==1354)
                       //                             {
                       //                               std::cerr <<"  "<< srcVal;
                       //                             }    
                       //                           if ((j+1)==437)
                       //                           if (connecth[ii]==1706)
                       //                             {
                       //                               std::cerr <<"  "<< srcVal;
                       //                             }    

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
                       //nodalval[ii]=1;
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
           

                   eqnData_->Node2EQN(connecth, connect_PDE);

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
                       ptEl->GetShFncAtIp(Sf, actIntPt);
                       valueAtIP[actIntPt-1]= nodalval*Sf;
                     }
                   Double jacDet;
                   for (UInt actIntPt=1; actIntPt<=nrIntPts;  actIntPt++) {
                        
                     jacDet = ptEl->CalcJacobianDetAtIp(actIntPt, ptCoordNodes);
                     ptEl -> GetShFncAtIp(Sf, actIntPt);
                        
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
                   eqnData_->Node2EQN(connecth, connect_PDE);
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
         
    //   if (timestep > 199)
    //     {
    //       timestep=timestep-200;
    //       strcat(anameloc,"8");
    //     }  
    //   else
    //     {
    //       if (timestep > 99)
    //      {
    //        timestep=timestep-100;
    //        strcat(anameloc,"7"); 
    //      }
    //       else
    //      {
    //        strcat(anameloc,"6"); 
    //      }      
    //     }
  
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
    UInt i,j;

    UInt buffernodenum = 0;
    for (i=0; i < maxnumnodes; i++)
      {
        flowdatafile >> buffernodenum;
        for (j=0; j < maxnumqtts; j++)
          flowdatafile >> flowdata_[j][buffernodenum-1];
      }

    flowdatafile.close();

    // With this we can printout the pressure data in 
    // a file with the format of a unverg file.
    /*  testflowf.open("test.flowf", std::ios::out);
        flowdatafile.seekg(pos, std::ios::beg);
        //testflowf << buffer;
        //testflowf << buffer2 <<  std::endl;
        //testflowf << "MAXNODES is: " << maxnumnodes << std::endl;

        for (i=0; i < maxnumnodes; i++)
        {
        testflowf << std::setw(10) << i+1;
        testflowf << std::endl;
        //for (j=0; j < maxnumqtts; j++)
        testflowf<< std::setiosflags(std::ios::uppercase | std::ios::scientific)
        << " " << flowdata_[0][i];
        testflowf << std::endl;
        }
        testflowf.close();*/
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
    

    params->Get("p0Phys", p0Phys, "vortexSrc");
    params->Get("rho0Phys", rho0Phys, "vortexSrc");
    params->Get("r_0Phys", r_0Phys, "vortexSrc");
    params->Get("GammaZirkPhys", GammaZirkPhys, "vortexSrc");
    params->Get("Equationgamma", Equationgamma, "vortexSrc");
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
    // Schritt  :   Berechne den Wirbel (dimensionsbehaftet)
    //--------------------------------------------------------------------------
    //
    r     = sqrt(XPhys*XPhys + YPhys*YPhys);             // Komplexe Koordinate
    std::complex<Double> z( XPhys       , YPhys      );  // des Wirbelpaars

    arg   = omegaPhys*tPhys;
    std::complex<Double> b( r_0Phys * cos(arg) ,           
              r_0Phys * sin(arg)  );        // Komplexe Bahn des Wirbelpaars

    std::complex<double> zq_bq(pow(z, 2) - pow(b, 2));
    
    // Komplexe Potentialfunktion
                                  
    std::complex<double> w(GammaZirkPhys/(2.*PI*Im) * log(zq_bq));
    
    // Ableitungen nach z, 1. Ordnung
      
    //Original:                            
    std::complex<double> w_z(GammaZirkPhys*z/(PI*Im*zq_bq));
    //Following Scully vortex core model:                            
    //std::complex<double> w_z(GammaZirkPhys*z/(PI*Im*(pow(z, 2) + pow(b, 2))));
    U_inc  =    real(w_z);                    
    V_inc  = - imag(w_z);

    // Ableitungen nach z, 2. Ordnung
    //Original:                                 
        std::complex<double> w_zz(-GammaZirkPhys*(pow(z,2) + pow(b,2)) / (PI*Im*pow(zq_bq,2)));
    //Following Scully vortex core model:  
    //std::complex<double> w_zz(GammaZirkPhys*(pow(b,2) - pow(z,2)) / (PI*Im*pow((pow(b,2) + pow(z,2)),2)));
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
                  if (abs(x)<1e-8) 
                    {
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
                      theta = atan(y / x);   // Polarkoordinaten r, theta
                    }
            
                  // Hankel Funktion 2.Ordg. und 2. Art H2_2
                  // mit Bessel Funkt. 2.Ordg., 1.Art JJ2 und 2.Art YY2
            
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
          srcVal=-P_t;
          //srcVal = 0.;
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
          //dTij_di*= density
          break;
        }
      case 4:
        {
          // Sending back only the velocity components for evaluation of
          // dTijdi in linearForm.cc
          dTij_di.Resize(2);   
          dTij_di[0] = U_inc; 
          dTij_di[1] = V_inc; 
          Double density=1.0;
          //dTij_di*= density
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
        {// Sending back the tangential velocity as srcVal
          srcVal=sqrt(U_inc*U_inc+V_inc*V_inc);
          P_ak = 0.;
          //srcVal=omegaPhys*r;
          dTij_di.Resize(2);   
          dTij_di[0] = U_inc; 
          dTij_di[1] = V_inc; 
          break;
        }
      default:
        (*error) << "Need to give an outType value between 1 and 4 "
                 << "Current outType value: " << outType;
        Error( __FILE__,__LINE__ );
      }
    

  }
  
  //___________________________________
  //
  // Besselfunktion 2. Ordg. und 1. Art
  //___________________________________
  
  Double AcouFlowNoise::JJ2(const Double x)
  {
    
    //--------------------------------------------------------------------------
    
    Double JJ2;
    Double RESULT;
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
  // Besselfunktion 2. Ordg. und 2. Art 
  //___________________________________

  Double AcouFlowNoise::YY2(const Double x)
  {
    //--------------------------------------------------------------------------
    Double YY2, RESULT;
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
