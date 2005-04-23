#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <stdio.h>

#include "acouFlowNoise.hh"
#include <Domain/bcs.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <MpCCIcpl/MpCCIexch.hh>
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Driver/solveStepAcouFlowNoise.hh"

#ifdef MpCCI
#include <cci.h>
#endif

namespace CoupledField
{

AcouFlowNoise::AcouFlowNoise(Grid *aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
			     FileType *aptFileType, WriteResults *aptOut)
    :AcousticPDE(aptgrid, aptbcs, aptTimeFunc, aptFileType, 
                 aptOut)
  {
    ENTER_FCN( "AcouFlowNoise::AcouFlowNoise" );

    // pdename_ is also acoustic for this case
    pdename_          = "acoustic";
    pdematerialclass_ = "fluid";
    nodalSrc_ = FALSE;
    plotRHS_ = FALSE;  

#ifdef MpCCI
    StdVector<Elem*> elemssd;

    params->GetList( "name", subdoms_, pdename_, "region" );
    params->GetList( "name", couplSubDomName_, "MpCCI-flownoise",
		     "coupledregion" );

    //check type of flow data
    if( params->HasValue( "type", "nodalSrc", pdename_, "flowData" ) ) {
      nodalSrc_ = TRUE;
      Info->PrintF(pdename_, " Using FlowData as RHS nodal source\n" );
    }
  
    if( params->IsSet( "valRHS","pdeList" ,"acoustic" ) ) {
      plotRHS_ = TRUE;
      Info->PrintF(pdename_, "Writing RHS as solution of problem\n" );
    }
    
    Boolean Find=FALSE;
    for (int i=0; i<subdoms_.GetSize(); i++)
      {
        for (int j=0; j<couplSubDomName_.GetSize(); j++)
          {
            if (couplSubDomName_[j] == subdoms_[i])
              {
                ptgrid_->GetElemSD(elemssd,couplSubDomName_[j],actlevel_);
                ptgrid_->CalcNumberOfNodesInPatch(elemssd,mapSD_);
                ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_.GetSize() );
                Find=TRUE;
              }
          }
      }
    if (!Find) 
      {
        std::string msg="Subdom to be coupled is not in list of PDE subdoms.";
	msg+="Please, check .xml-file";
        Error(msg.c_str(),__FILE__,__LINE__);
      }

    MpCCInodes_=mapSD_.GetSize();
    ptMpCCIexch_->PutExchangeGrid2MpCCI(couplSubDomName_);

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
    Vector<Double> elemvec;

    Integer i;

    Integer level=0;

    // get maximum number of elements from grid
    Integer maxnumelem=ptgrid_->GetMaxnumElem(level,subdoms_);


  double starttime, endtime;

    Double val;
    Matrix<Double> ptCoordNodes;
    // For ObstSurf
    Matrix<Double> ptCoordNodSurf;
    // For set of elements corresponding to surface elements
    Matrix<Double> ptCoordNodBelongSE;
    Matrix<Double> deriv;
    StdVector<Integer> connecth;
    // For changing connecth to PDE
    StdVector<Integer> connect_PDE; 
    Vector<Integer> connObstSurf; // For ObstSurf
    Vector<Integer> connBelongSE;
    // vector of 1D-elements (ObstSurf) from mesh-file
    StdVector<Elem*> ObstSurf;
    BaseFE * ptEl;
    BaseFE * ptElSurf;
    BaseFE * ptElBelongSE;
 
  
    Integer j,ii, elsize = -1;
    StdVector<Elem*> elemssd;     
    static Integer timestep=0;
    static int auxtime=0;
 

    Double valmult;


    //  std::cout<<"timestep counter in ComputeRHS: "<<timestep<<std::endl;

#ifdef MpCCI
    std::cout<<"MpCCInodes_: "<< MpCCInodes_ << " dimension: " 
	     << dim_ << std::endl;
    if (nodalSrc_ == TRUE)
      //we get already the integrated acoustic source term
      flowdata_.Resize(1, MpCCInodes_);
    else
      flowdata_.Resize(1+dim_, MpCCInodes_);

      starttime = CCI_Wtime();
      


       ptMpCCIexch_->CouplCompPhase(flowdata_, timestep);


      endtime = CCI_Wtime();
      

      std::cout<<"Transfer of Data CouplCompPhase() for 1 time step took: "
	       <<(endtime-starttime)<<" seconds"<<std::endl;


#else
    // If data from fluid file call to get fluid flow data in flowdata_  
    std::string flowdata;
    Error("Need change to XML-Standard");
      
    //      conf->get("acousrc_file",flowdata);
    //      ReadFlowData(flowdata.c_str(), timestep, flowdata_); 
#endif 



  
    // Quadrupole computation
    std::cout << "Processing RHS volume elems for quadrupole... "<< std::endl;


    // Variables for ramping
    Double xfmin, yfmin, xfmax, yfmax, facRampXmin, facRampYmin,
      facRampXmax, facRampYmax ;
    Double bndoffsetXmin, bndoffsetYmin, bndoffsetXmax, bndoffsetYmax ;

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

    // Correct valmult value is -1.0, 
    // if plugging in source (ddTij/dxidxj) directly in weak form then 1.0
#ifdef MpCCI
       starttime = CCI_Wtime();
#endif 
    if (nodalSrc_ == FALSE) {
      
      valmult=-1.0;
      
      for (i=0; i<couplSubDomName_.GetSize(); i++)
        {
          ptgrid_->GetElemSD(elemssd,couplSubDomName_[i],level);
        
          for (j=0; j< elemssd.GetSize(); j++)
            {
              ptEl=elemssd[j]->ptElem;
              BaseForm * linear_load = new LinearFlowNoiseInt(ptEl);
            
              Integer ii;
              elsize=(elemssd[j]->connect).GetSize();
              connecth.Resize(elsize);
              for (ii=0; ii<elsize; ii++)
                connecth[ii]=(elemssd[j]->connect)[ii];
            
              Matrix<Double> ptCoordNodes;
              ptgrid_->GetCoordNodesElemMat(connecth,  ptCoordNodes, level);        
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

    else {
      Integer eqnDof, eqnNr, node, dof;
      StdVector<Integer> connect(1);
 
      valmult=-1.0;    
     
      eqnDof = 1;
      dof    = 1;
      for (Integer idx=0; idx<flowdata_.GetSizeCol() ; idx++) {
        Double val = flowdata_[0][idx];
        node = idx + 1;

      // Ramping before adding to RHS vector (NOW IT MAKES ZERO THE RHS ENTRY!)
      Matrix<Double> ptCoordNodes;
      connecth.Resize(1);
      connecth[0] = node;
      ptgrid_->GetCoordNodesElemMat(connecth,  ptCoordNodes, level);	  

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
	rhs_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
	rhs_.Init(0.0);
	
	Double *ptRHS;
	algsys_->GetRHSVal( ptRHS );
	rhs_.CopyFromAlgSysDataPointer(ptRHS);
      }
      
  
      //     timestep=timestep+1;

  } 


  void AcouFlowNoise::ReadFlowData(const Char * aname, Integer timestep,
				   Matrix<Double> &flowdata_)
  {
    ENTER_FCN( "AcouFlowNoise::ReadFlowData" );

    std::ifstream flowdatafile;
    Integer maxnumnodes;
    Integer maxnumqtts;
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

    Integer ibuf;

    /* Set pointer to beginning of file: */
    std::string::size_type pos=0;
    flowdatafile.seekg(pos, std::ios::beg);

    flowdatafile >> ibuf >> maxnumqtts >> maxnumnodes >> ibuf;

    flowdata_.Resize(maxnumqtts,maxnumnodes);
    Integer i,j;

    Integer buffernodenum = 0;
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


  void AcouFlowNoise::WriteResultsInFile()
  {
    ENTER_FCN( "AcouFlowNoise::WriteResultsInFile" );

    NodeStoreSol<Double> sol_der1Array, sol_der2Array;
    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;
  
    if (analysistype_ == TRANSIENT)
      {
        sol_der1Array.SetNumSolutions(1);
        sol_der1Array.SetNumNodes(numPDENodes_);
        sol_der1Array.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
        sol_der1Array.SetNumDofs(dofspernode_);
        sol_der1Array.Init(0.0);
        sol_der1Array.SetAlgSysVector(getS1());
      
        sol_der2Array.SetNumSolutions(1);
        sol_der2Array.SetNumNodes(numPDENodes_);
        sol_der2Array.SetSolutionType(ACOU_POTENTIAL_DERIV_2);
        sol_der2Array.SetNumDofs(dofspernode_);
        sol_der2Array.Init(0.0);
        sol_der2Array.SetAlgSysVector(getS2());
 
	  std::cout<<"In AcouFlowNoise. Writing solution..."<<std::endl;
     

	if (plotRHS_)
	  {
	  std::cout<<"In AcouFlowNoise. Writing RHS as solution..."<<std::endl;
	         outFile_->WriteNodeSolutionTransient(rhs_,laststepcalc_,
	  					     lasttimecalc_);
		outFile_->WriteNodeSolutionTransient(rhs_, laststepcalc_, 
						     lasttimecalc_);
	  }
	//	else
	//  {
	  std::cout<<"In AcouFlowNoise. Writing solution..."<<std::endl;
	    solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
	    //  }
    
        outFile_->WriteNodeSolutionTransient(*solTransient,laststepcalc_,
					     lasttimecalc_);
        outFile_->WriteNodeSolutionTransient(sol_der1Array,
                                             laststepcalc_,lasttimecalc_);
        outFile_->WriteNodeSolutionTransient(sol_der2Array,
                                             laststepcalc_,lasttimecalc_);
      }
    else if (analysistype_ == HARMONIC)
      {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
        outFile_->WriteNodeSolutionHarmonic(*solHarmonic, actFreqStep_, 
                                            actFrequency_, complexFormat_);
      }
    else
      Error("AcouFlowNoisePDE: Only transient and harmonic results possible",
            __FILE__, __LINE__);
  }

} // end of namespace
