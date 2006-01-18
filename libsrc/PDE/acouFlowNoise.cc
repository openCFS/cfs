#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <complex>

#include "acouFlowNoise.hh"

#include "DataInOut/Unverg/outUnverg.hh"
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
  
    if( params->IsSet( "valRHS","pdeList" ,"acoustic" ) ) {
      plotRHS_ = TRUE;
      Info->PrintF(pdename_, "Writing RHS as solution of problem\n" );
    }

#ifdef MpCCI
    //check type of flow data
    if( params->HasValue( "type", "nodalSrc", pdename_, "flowData" ) ) {
      nodalSrc_ = TRUE;
      Info->PrintF(pdename_, " Using FlowData as RHS nodal source\n" );
    }
    StdVector<Elem*> elemssd;
    StdVector<std::string> regionNames, coupledRegionNames;
    params->GetList( "name", regionNames, pdename_, "region" );
    params->GetList( "name", coupledRegionNames,  "MpCCI-flownoise",
                     "coupledregion" );
    ptgrid_->RegionNameToId( subdoms_, regionNames );
    ptgrid_->RegionNameToId( couplSubDomId_, coupledRegionNames );

    Boolean OnlyLinNodes=TRUE;
    Boolean Find=FALSE;
    std::cerr << "Subdoms:\n" << subdoms_ << std::endl;
    std::cerr << "CoupledSubdoms:\n" << couplSubDomId_  << std::endl;
    

    for (int i=0; i<subdoms_.GetSize(); i++)
      {
        for (int j=0; j<couplSubDomId_.GetSize(); j++)
          {
            if (couplSubDomId_[j] == subdoms_[i])
              {

                ptgrid_->GetVolElems(elemssd,couplSubDomId_[j]);
		
                if (dim_ == 3)
                  { //TO SEND ONLY LINEAR TETRAS TO MPCCI SINCE NO QUAD TETRAS
                    //ARE ALLOWED
                    ptgrid_->GetNodesOfElemList(mapSD_allNodes_, elemssd, FALSE);
                    ptgrid_->GetNodesOfElemList(mapSD_, elemssd, OnlyLinNodes);
                    ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_);
                  }
                else
                  {
                    ptgrid_->GetNodesOfElemList(mapSD_, elemssd, FALSE);
                    //ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_.GetSize() );
                    ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_);
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
    // This is due to workaround from quad to lin elements for MpCCI
    if (dim_ == 3)
      MpCCInodes_=mapSD_allNodes_.GetSize();
    else
      MpCCInodes_=mapSD_.GetSize();
    
    ptMpCCIexch_->PutExchangeGrid2MpCCI(couplSubDomId_);
#else
    if( params->HasValue( "type", "vortexSrc", pdename_, "flowData" ) ) {
      vortexSrc_ = TRUE;
      StdVector<Elem*> elemssd;
      StdVector<std::string> regionNames, coupledRegionNames;
      params->GetList( "name", regionNames, pdename_, "region" );
      params->GetList( "name", coupledRegionNames,  "MpCCI-flownoise",
                       "coupledregion" );
      ptgrid_->RegionNameToId( subdoms_, regionNames );
      ptgrid_->RegionNameToId( couplSubDomId_, coupledRegionNames );
      std::cerr << "Subdoms:\n" << subdoms_ << std::endl;
      std::cerr << "CoupledSubdoms:\n" << couplSubDomId_  << std::endl;
      Info->PrintF(pdename_, " Using Vortex function as RHS nodal source\n" );
    }
    else
      Error( "No MpCCI coupling, need to define VortexSrc", __FILE__, __LINE__ );    

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

    Vector<Double> coeffMass, coeffDamp;
    Vector<Double> elemvec, nodalval;

    UInt i;

//     // get maximum number of elements from grid
//     UInt maxnumelem=ptgrid_->GetNumElems(subdoms_);

//     Double val;
    Matrix<Double> ptCoordNodes;
    // For ObstSurf
    Matrix<Double> ptCoordNodSurf;
    // For set of elements corresponding to surface elements
    Matrix<Double> ptCoordNodBelongSE;
    Matrix<Double> deriv;
    StdVector<UInt> connecth;
    // For changing connecth to PDE
    StdVector<Integer> connect_PDE; 
    Vector<UInt> connObstSurf; // For ObstSurf
    Vector<UInt> connBelongSE;
    // vector of 1D-elements (ObstSurf) from mesh-file
    StdVector<Elem*> ObstSurf;
    BaseFE * ptEl;
    Double actTime;

    UInt j;
    UInt elsize = 0;
    StdVector<Elem*> elemssd;     
 

    Double valmult;

    //Getting current time
    actTime = solveStep_->GetActTime();

    //  std::cout<<"timestep counter in ComputeRHS: "<<timestep<<std::endl;

#ifdef MpCCI
    double starttime, endtime;
    std::cout<<"MpCCInodes_: "<< MpCCInodes_ << " dimension: " 
             << dim_ << std::endl;
    if (nodalSrc_ == TRUE)
      //we get already the integrated acoustic source term
      flowdata_.Resize(1, MpCCInodes_);
    else
      flowdata_.Resize(1+dim_, MpCCInodes_);

    starttime = CCI_Wtime();
      

    Integer timestep = 0;
    ptMpCCIexch_->CouplCompPhase(flowdata_, actTime);


    endtime = CCI_Wtime();
      

    std::cout<<"Transfer of Data CouplCompPhase() for 1 time step took: "
             <<(endtime-starttime)<<" seconds"<<std::endl;


#else
    // If data from fluid file, call to get fluid flow data in flowdata_  
    std::string flowdata;
    //    Warning( "Using no MpCCI, vortexSrc is assumed.", __FILE__, __LINE__ );
      
//     conf->get("acousrc_file",flowdata);
//     // ReadFlowData(flowdata.c_str(), timestep, flowdata_); 

//     std::string analysis;
//     params->Get( "type", analysis, "analysis" );

//     if ( analysis=="harmonic" ) {
//       Double actFreq;   
//       actFreq = solveStep_->GetActFreq();
//     }



//     // Calling to get load (RHS acoustic source) either in time or frequency domain
//     // So step can be a time or a frequency value
//     ReadFlowData(flowdata.c_str(), actFreq, flowdata_); 

#endif 



  
    // Quadrupole computation
    std::cout << "Processing RHS volume elems for quadrupole... "<< std::endl;


    // Variables for ramping
    Double xfmin, yfmin, zfmin, xfmax, yfmax, zfmax, facRampXmin, facRampYmin, 
      facRampZmin, facRampXmax, facRampYmax, facRampZmax;
    Double bndoffsetXmin, bndoffsetYmin, bndoffsetZmin, bndoffsetXmax, bndoffsetYmax, bndoffsetZmax ;

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
    

    // Correct valmult value is -1.0, 
    // if plugging in source (ddTij/dxidxj) directly in weak form then 1.0
#ifdef MpCCI
    starttime = CCI_Wtime();
#endif 
    if (nodalSrc_ == FALSE) {
      
      if (vortexSrc_ == TRUE)
        {
          for (i=0; i<couplSubDomId_.GetSize(); i++)
            {
              std::cerr << atime;
              ptgrid_->GetVolElems(elemssd,couplSubDomId_[i]);
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
                  for (ii=0; ii<elsize; ii++)
                    {
                      Double r_sqr=((ptCoordNodes[0][ii])*(ptCoordNodes[0][ii])+
                                    (ptCoordNodes[1][ii])*(ptCoordNodes[1][ii]));

                      //For square 100times100 domain                      
//                       if (r_sqr<=((1000./81.)*(1000./81.)))
//                         {
//                           nodalval[ii]=0;
//                         }
                      
                      if (r_sqr<=((200./81.)*(200./81.)))
                        {
                          nodalval[ii]=0;
                        }
                      else
                        {
                          Double press=0;
                          //std::cout << "In else before calling VortexAnalytical"<<std::endl;
                          VortexAnalytical(press, ptCoordNodes[0][ii],
                                           ptCoordNodes[1][ii], actTime, 2);
                          //std::cout << "In else after calling VortexAnalytical, press= "<<press<<std::endl;
//                           //For debugging the quadratic middle nodes
//                           if ((j+1)==467)
//                           if (connecth[ii]==1708)
//                             {
//                               std::cerr <<"  "<< press;
//                             }
//                           if ((j+1)==464)
//                           if (connecth[ii]==1356)
//                             {
//                               std::cerr <<"  "<< press;
//                             }                          
//                           if ((j+1)==434)
//                           if (connecth[ii]==1354)
//                             {
//                               std::cerr <<"  "<< press;
//                             }    
//                           if ((j+1)==437)
//                           if (connecth[ii]==1706)
//                             {
//                               std::cerr <<"  "<< press;
//                             }    

                          nodalval[ii]=press;
                          //nodalval[ii]=1;
                        }
                      }


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
              //ALSO TO FINISH THE OUTPUT OF CONTROL DATA 
              std::cerr <<std::endl;
            }

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
    
    else {

      UInt eqnDof, node, dof;
      Integer eqnNr;
      StdVector<UInt> connect(1);
 
      valmult=-1.0;    
     
      eqnDof = 1;
      dof    = 1;
      for (UInt idx=0; idx<flowdata_.GetSizeCol() ; idx++) {
        Double val = flowdata_[0][idx];
        node = idx + 1;

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

        val*=valmult;

        //add to RHS
        eqnData_->Node2EQN(node,dof,eqnNr,eqnDof);
        algsys_->SetNodeRHS(val, pdeId_, eqnNr, eqnDof);      
      }

    
    
    }
#ifdef MpCCI
    endtime = CCI_Wtime();
#endif    
  
    if (plotRHS_){
      ///////// For plotting the RHS as solution for analysing it

      rhs_.SetNumSolutions(1);
      rhs_.SetNumNodes(numPDENodes_);
      rhs_.SetSolutionType(ACOU_RHSVAL);
      rhs_.SetNumDofs(1);
      rhs_.SetPtrEQNData(eqnData_, ptgrid_);
      rhs_.Init(0.0);
        
      Double *ptRHS;
      algsys_->GetRHSVal( ptRHS );
      rhs_.CopyFromAlgSysDataPointer(ptRHS);
    }
      
  
    //     timestep=timestep+1;

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

  void AcouFlowNoise::VortexAnalytical(Double &press, const Double x,
                                       const Double y, const Double t, 
                                       const UInt outType)
  {
    ENTER_FCN( "AcouFlowNoise::VortexAnalytical" );

    // From Vortex Implementation of HydSol by username iagschwa

    // local variable declaration

    Double      tPhys, XPhys, YPhys, arg;
    Double      p0Phys, rho0Phys, r_0Phys, GammaZirkPhys;
    Double      c0Phys, MaRot, iMaRot, qPhys, omegaPhys, akkPhys;
    Double Equationgamma; // before it was EQN%gamma
     
    Double      U_inc;                         // Variablen, Stroemungs-
    Double      V_inc;                         // groessen aus Potential-
    Double      P_inc;                      // theorie
    Double      P_ak;                       // Akust. Druck aus MAE
    Double      Phi_t, Phi_tt, Phi_ttt;     // und seine Ableitungen

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
    std::complex<double>    zq_bq;                         // z**2-b**2
    //std::complex<double>    Im;                            // i defined below
    std::complex<double>    w;                             // Komplexe Potential-
    // funktion
    std::complex<double>    w_z, w_t, w_zz, w_zt;          // und ihre Ableitungen
    std::complex<double>    w_zzz, w_tt, w_ttt, w_ztt;
    std::complex<double>    w_zzt;
    std::complex<double>    H2_2;
    // dimensionsbehaftete Var. die noch nicht in common definiert sind
    Double     KappaPhys, T_umlaufPhys, lambdaPhys;                    

 

    //--------------------------------------------------------------------------
    // ...Phys kennzeichnet die dimensionsbehafteten Groessen
    // ...     (ohne Endung) sind entdimensionalisierte Groessen
    //--------------------------------------------------------------------------
    //
    // Gegeben sind die folgenden Groessen :         
    //                                               
    p0Phys        = 0.714285714;                          // Umgeb.druck  [kg/(m*s^2)]
    rho0Phys      = 1.;                        // Umgeb.dichte     [kg/m^3]
    r_0Phys       = 1.;                           // halb. Wirbelabst. [  m  ]
    GammaZirkPhys = 1.003;                          // Zirkulation       [m^2/s]
    Equationgamma = 1.4;                        // Need to be defined here as well

    params->Get("p0Phys", p0Phys, "vortexSrc");
    params->Get("rho0Phys", rho0Phys, "vortexSrc");
    params->Get("r_0Phys", r_0Phys, "vortexSrc");
    params->Get("GammaZirkPhys", GammaZirkPhys, "vortexSrc");
    params->Get("Equationgamma", Equationgamma, "vortexSrc");
    
    c0Phys        = sqrt(Equationgamma*p0Phys/rho0Phys);  // Schallgeschw.[ m/s ]
    MaRot         = GammaZirkPhys/4./PI /r_0Phys/c0Phys;// Rotations-Mach-Zahl
    iMaRot        = 1. / MaRot;                      // 1/Ma_rot
    
    qPhys         = MaRot * c0Phys;                  // induz. Geschw.    [ m/s ]
    KappaPhys     = GammaZirkPhys / (2.*PI);     // Wirbelstaerke pro 
    // Wirbel            [m^2/s]
    T_umlaufPhys  = 2. * PI * r_0Phys / qPhys;   // Zeit fuer einen
    // Wirbelumlauf      [  s  ]
    omegaPhys     = 2. * PI / T_umlaufPhys;      // Winkelgeschw.des
    // Wirbel            [ 1/s ]
    akkPhys       = 2. * omegaPhys / c0Phys;         // Argument,Berechnung 
    // akust.Druck,MAE   [ 1/m ]
    lambdaPhys    = c0Phys * T_umlaufPhys;           // Wellenlaenge      [  m  ]

    std::complex<Double> Im(0., 1.);
      
    tPhys = t;
    XPhys = x;
    YPhys = y;

    //--------------------------------------------------------------------------
    // Schritt  :   Berechne den Wirbel (dimensionsbehaftet)
    //--------------------------------------------------------------------------
    //
    r     = sqrt(XPhys*XPhys + YPhys*YPhys);             // radius


    std::complex<Double> z( XPhys       , YPhys      );       // Wirbelpaars

    arg   = omegaPhys*tPhys;
    std::complex<Double> b( r_0Phys * cos(arg) ,           
              r_0Phys * sin(arg)  );        // Komplexe Bahn des Wirbelpaars

    zq_bq = pow(z, 2) - pow(b, 2);
    
    // Komplexe Potentialfunktion
                                  
    w = GammaZirkPhys/(2.*PI*Im) * log(zq_bq);
    
    // Ableitungen nach z, 1. Ordnung
                                  
    w_z    = GammaZirkPhys*z/(PI*Im*zq_bq);
    U_inc  =    real(w_z);                    
    V_inc  = - imag(w_z);  
    

    // Ableitungen nach z, 2. Ordnung
                                  
    w_zz  = -GammaZirkPhys*(pow(z,2) + pow(b,2)) / (PI*Im*pow(zq_bq,2));
    U_x   =   real(w_zz);
    V_x   = -imag(w_zz);
    U_y   =  V_x;
    V_y   = -U_x;
    
    // Ableitungen nach z, 3. Ordnung
                                  
    w_zzz = -GammaZirkPhys*(-2.*pow(z,3) - 6.*z*pow(b,2))/(PI*Im*pow(zq_bq,3));
    U_xx =  real( w_zzz);
    U_xy = -imag(w_zzz);
    
    V_yy =  imag(w_zzz);
    V_xy = -real( w_zzz);
    
    U_yy =  V_xy;
    V_xx =  U_xy;
    
    // Ableitungen nach z, nach t, 2. Ordnung
                                  
    w_zt  =  2.* omegaPhys* GammaZirkPhys* z* pow((b),2) / (PI*pow((zq_bq),2));
    U_t  =  real( w_zt);
    V_t  = -imag(w_zt);
    
    // Ableitungen nach t, 1. und 2. und 3. Ordnung
                                  
    w_t     = -omegaPhys*GammaZirkPhys * pow(b,2) / ( PI*zq_bq );
    Phi_t   = real(w_t);
    
    w_tt    = -(2.*Im*pow(omegaPhys,2)*GammaZirkPhys/PI)*pow(b,2)*pow(z,2)/pow(zq_bq,2);
    Phi_tt  = real(w_tt);
    
    w_ttt   = (4*pow(omegaPhys,3)*GammaZirkPhys/PI)*pow(b,2)*pow(z,2)*(pow(z,2)+pow(b,2))/pow(zq_bq,3);
    Phi_ttt = real(w_ttt);
    
    // Ableitungen nach z, nach t, Gesamtordnung 3
                                  
    w_ztt     = (4.*Im*pow(omegaPhys,2)*GammaZirkPhys/PI)*(pow(b,2)*pow(z,3) + z*pow(b,4))/pow(zq_bq,3);
    U_tt =  real( w_ztt);
    V_tt = -imag(w_ztt);
    
    w_zzt     = -(2.*omegaPhys*GammaZirkPhys/PI)*(3.*(pow(b,2))*(pow(z,2)) + pow(b,4))/pow(zq_bq,3);
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
    





    




    //-----------------Berechne den akkustischen Druck--------------------------

    if (outType == 1) 
      {              // Wenn outAkk gleich 1, dann soll der
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
            
            H2_2      = ( JJ2(akkPhys*r), -YY2(akkPhys*r) ) ;
            
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
        press=P_ak;
      }
    else
      {
        P_ak = 0.;
        press=-P_t;
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
