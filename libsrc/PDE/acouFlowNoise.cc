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

#ifdef MpCCI
   StdVector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,subdoms_[1],actlevel_);
       ptgrid_->CalcNumberOfNodesInPatch(elemssd,mapSD_);
 ptMpCCIexch_ = new MpCCIexch(ptgrid_,mapSD_.size() );
#endif
//  ReadBCs(pdename_);
 preComputeRHS();
}

AcouFlowNoise::~AcouFlowNoise()
{
  ENTER_FCN( "AcouFlowNoise::AcouFlowNoise" );
#ifdef MpCCI
  delete ptMpCCIexch_;
#endif
}

void AcouFlowNoise::preComputeRHS()
{
  ENTER_FCN( "AcouFlowNoise::preComputeRHS" );
  
  //first check if coupling per MpCCI or if a flow-result file is provided
  MpCCI_=FALSE;
  std::string strMpCCI="no";
  conf->ifget("MpCCI",strMpCCI,pdename_);
  if (strMpCCI == "yes")
    {
#ifdef MpCCI
      MpCCI_ = TRUE;
      std::cout << "DO COUPLING via MpCCI" << std::endl;
      //      MpCCIexch * ptMpCCIexch_ = new MpCCIexch(ptgrid_);
      //MpCCInodes_=ptgrid_->GetMaxnumnodes(0);

      MpCCInodes_=mapSD_.size();
      ptMpCCIexch_->PutExchangeGrid2MpCCI(subdoms_);
      flowdata_.Resize(3, MpCCInodes_);
      //  delete ptMpCCIexch_;
      SetRHSFnc = FALSE;
      SetRHSFlowSrc = FALSE;
      std::string isthererhs;
      conf->ifget("load_force",isthererhs,"acoustic");
  
      if (isthererhs=="FlowSrc" ) 
	{
	  std::cout<<"In flownoise precomputeRHS"<<std::endl;  
	  conf->getliststr("rhs_surfaces",rhs_surfaces_,"acoustic");
	  SetRHSFlowSrc=TRUE;
	}
#endif
    }
  else
    {
      SetRHSFnc = FALSE;
      SetRHSFlowSrc = FALSE;
      std::string isthererhs;
      conf->ifget("load_force",isthererhs,"acoustic");
  
      if (isthererhs=="FlowSrc" ) 
	{
	  std::cout<<"In flownoise precomputeRHS -no MpCCI used-"<<std::endl; 
	  conf->getliststr("rhs_surfaces",rhs_surfaces_,"acoustic");
	  SetRHSFlowSrc=TRUE;
	}
    }
}


void AcouFlowNoise::ComputeRHS(const Double atime)
{
  ENTER_FCN( "AcouFlowNoise::ComputeRHS" );

  Vector<Double> coeffMass, coeffDamp;
  Vector<Double> elemvec;
  Vector<Double> elemvecdip;
  Integer i;

  Integer level=0;

  // get maximum number of elements from grid
  Integer maxnumelem=ptgrid_->GetMaxnumElem(level,subdoms_);

  Double val;
  Matrix<Double> ptCoordNodes;
  Matrix<Double> ptCoordNodSurf; // For ObstSurf
  Matrix<Double> ptCoordNodBelongSE; // For set of elements corresponding to surface elements
  Matrix<Double> deriv;
  StdVector<Integer> connecth;
  StdVector<Integer> connect_PDE;	// For changing connecth to PDE
  Vector<Integer> connObstSurf;	// For ObstSurf
  Vector<Integer> connBelongSE;
  StdVector<Elem*> ObstSurf;  // vector of 1D-elements (ObstSurf) from mesh-file
  BaseFE * ptEl;
  BaseFE * ptElSurf;
  BaseFE * ptElBelongSE;
  Matrix<Double> flowdata_; // Where the nodal data is going to be stored
  
  Integer j,ii, elsize = -1;
  StdVector<Elem*> elemssd;     
  static Integer timestep=0;
  static int auxtime=0;
 
  std::string flowdata;
  conf->get("acousrc_file",flowdata);
    
  //For the subdomain next to boundary of obstacle (Next2Surf)
  StdVector<Elem*> Next2Surf;  // vector of Volume elements next to surface from mesh-file      
  BaseFE * ptElNext2Surf;      
  StdVector<Elem*> belongSE;
 
 
//   ObstSurf=ptBCs_->getEdgesBC(rhs_surfaces_[0],level);
//   std::cout<<"Number of surfelems: "<<ObstSurf.size()<<std::endl;
  
//   Next2Surf=ptBCs_->getNeighElemsForSurfaces(rhs_surfaces_[1],level);
//   ptgrid_->DefineBelonging4Elems(ObstSurf,Next2Surf,belongSE);
  
  Double valmult;

  valmult = 3.7; // 3.7 is to give units to the fluid data and coeffRHS is from inhom. wave eq.

  std::cout<<"timestep counter in ComputeRHS: "<<timestep<<std::endl;
  
  if(MpCCI_)
    {
#ifdef MpCCI
      std::cout<<"MpCCInodes_: "<< MpCCInodes_ << "dimension: " << dim_ << std::endl;
      //     MpCCIexch * ptMpCCIexch_ = new MpCCIexch(ptgrid_);
      flowdata_.Resize(1+dim_, MpCCInodes_);
      ptMpCCIexch_->CouplCompPhase(flowdata_, timestep);
      //  delete ptMpCCIexch_;
#endif
    }
  else
    // If data from fluid file call to get fluid flow data in flowdata_  
    ReadFlowData(flowdata.c_str(), timestep, flowdata_); 

  std::cout << "Processing RHS surface elems for dipole... "<< std::endl;

//   // This is for the loop over the surface elements
//   for (j=0; j< ObstSurf.size(); j++)
//     { 
//       ptElBelongSE=belongSE[j]->ptElem;
//       // This will be done inside the 2d element next to the 1d element to get gradP at the center
//       Integer n=ptElBelongSE->GetNumNodes(); // This returns number of integration points      

//       elsize=(belongSE[j]->connect).GetSize(); // Get element number of nodes 
//       connBelongSE.Resize(elsize);
//       for (ii=0; ii<elsize; ii++)
// 	connBelongSE[ii]=(belongSE[j]->connect)[ii];

//       ptgrid_->GetCoordNodesElemMat(connBelongSE,  ptCoordNodBelongSE, level);
  
    
//       StdVector<Double> gradN_x_P; // This is done due to the different parameter type Vector and StdVector
//       gradN_x_P.resize(dim_);
//       gradN_x_P*=0;
//       std::vector<Double> LCoord(dim_,0);
//       Double jacDet;
      
//       // TO COMMENT OUT ONLY WHILE USING FILES WITH GRADIENT VALUE!!!
//       // Gradient of P at center by average of value at four nodes of neighbour 2d element
//       //  for (ii=0; ii<n; ii++)
//       //     {  
// 	//ptElBelongSE->GetGradientShFncAtCenter(help[ii],ii+1);

//       jacDet=ptElBelongSE->CalcJacobianDet(LCoord, ptCoordNodBelongSE);
// 		std::cout<<"jacDet:"<<jacDet<<std::endl;      
	       
// 	ptElBelongSE->GetGlobDerivShFnc  (deriv, LCoord, ptCoordNodBelongSE);
// 	//deriv.Transpose(derivTrans);

// 	for (ii=0; ii<n; ii++)
// 	  {  
// 	    for(i=0;i<dim_;i++)
// 	      {
// 		gradN_x_P[i]+=jacDet*deriv[ii][i]*(flowdata_[0][connBelongSE[ii]-1]);
// 		std::cout<<"gradN_x_P["<<i<<"] :"<<gradN_x_P[i]<<std::endl;
// 		std::cout<<std::endl;
// 	      }
// 	  }

// 	ptElSurf=ObstSurf[j]->ptElem;
// 	std::cout<<"connect: "<<ObstSurf[j]->connect<<std::endl;
	
// 	BaseForm * linear_loaddipole = new LinearFlowNoiseInt(ptElSurf);
	  
// 	connObstSurf=ObstSurf[j]->connect;

// 	ptgrid_->GetCoordNodesElemMat(connObstSurf, ptCoordNodSurf, level);
//         std::cout<<"coordinates :"<<ptCoordNodSurf<<std::endl;
// 	linear_loaddipole->CalcElemVector4Dip(ptCoordNodSurf, connObstSurf,elemvecdip,gradN_x_P);
// 	elemvecdip*=valmult;



// 	// CHANGE connecth
// 	Mesh2PDENode(connect_PDE,connObstSurf,mesh2PDENode_);
// 	std::cout<<"connObstSurf :"<<connObstSurf<<std::endl;
//  	  std::cout<<"elemvect DIPOLE: "<<elemvecdip<<std::endl;
// 	// including the dipole contribution for testing!!!!!!!!
// 	algsys_->SetElementRHS(&elemvecdip[0], connect_PDE.GetPointer(), connect_PDE.GetSize());

// 	delete linear_loaddipole;

//       }

  
  // Quadrupole computation
  std::cout << "Processing RHS volume elems for quadrupole... "<< std::endl;


  // Variables for ramping
  Double xfmin, yfmin, xfmax, yfmax, facRampXmin, facRampYmin, facRampXmax, facRampYmax ;
  Double bndoffsetXmin, bndoffsetYmin, bndoffsetXmax, bndoffsetYmax ;
  conf->get("xfmin",xfmin,"acoustic"); // minimum x coord. of fluid domain
  conf->get("yfmin",yfmin,"acoustic"); // minimum y coord. of fluid domain	   
  conf->get("xfmax",xfmax,"acoustic"); // maximum x coord. of fluid domain
  conf->get("yfmax",yfmax,"acoustic"); // maximum y coord. of fluid domain
  conf->get("facrampXmin",facRampXmin,"acoustic"); // factor for starting ramping
  conf->get("facrampYmin",facRampYmin,"acoustic"); // factor for starting ramping
  conf->get("facrampXmax",facRampXmax,"acoustic"); // factor for starting ramping
  conf->get("facrampYmax",facRampYmax,"acoustic"); // factor for starting ramping
  bndoffsetXmin=facRampXmin*xfmin;
  bndoffsetYmin=facRampYmin*yfmin; 
  bndoffsetXmax=facRampXmax*xfmax;
  bndoffsetYmax=facRampYmax*yfmax;

  valmult=-1.0;
      
  for (i=(subdoms_.GetSize())-1; i<subdoms_.GetSize(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

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
	  linear_load->CalcElemVector4Quad(ptCoordNodes, connecth, flowdata_, elemvec);
	  elemvec*=valmult;

	  // Ramping before adding elemrhs to global vector to avoid spurious effect at bnd. of fluid dom.

	  for (ii=0; ii<elsize; ii++)
	    {
	      if (ptCoordNodes[0][ii]<bndoffsetXmin)
		{
		  
		  elemvec[ii]-=elemvec[ii]*(ptCoordNodes[0][ii]-bndoffsetXmin)/(xfmin-bndoffsetXmin);
		}
	       
	      else
		if (ptCoordNodes[0][ii]>bndoffsetXmax)
		  elemvec[ii]-=elemvec[ii]*(ptCoordNodes[0][ii]-bndoffsetXmax)/(xfmax-bndoffsetXmax);
	      if (ptCoordNodes[1][ii]<bndoffsetYmin)
		elemvec[ii]-=elemvec[ii]*(ptCoordNodes[1][ii]-bndoffsetYmin)/(yfmin-bndoffsetYmin);
	      else	
		if (ptCoordNodes[1][ii]>bndoffsetYmax)
		  elemvec[ii]-=elemvec[ii]*(ptCoordNodes[1][ii]-bndoffsetYmax)/(yfmax-bndoffsetYmax);
	    }
	   
	  //end ramping

	  // CHANGE connecth
	  //Mesh2PDENode(connect_PDE,connecth,mesh2PDENode_);
	  eqnData_->Node2EQN(connecth, connect_PDE);

	  //linear_load->CalcElemVector(ptCoordNodes, elemvec); // for setting with homogeneous rhs
	  
	  // if (j>100){ 
	  //  valmult=-10000.0;
	    //std::cout<<"elemvect quad: "<<elemvec<<std::endl;
	  
	  //}
	  
	
	  algsys_->SetElementRHS(&elemvec[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
	  
	  delete linear_load;
	}
    }


  
  //     timestep=timestep+1;

} 


void AcouFlowNoise::ReadFlowData(const Char * aname, Integer timestep, Matrix<Double> &flowdata_)
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
// 	{
// 	  timestep=timestep-100;
// 	  strcat(anameloc,"7");	
// 	}
//       else
// 	{
// 	  strcat(anameloc,"6");	
// 	}      
//     }
  
  sprintf(aux,"%i",timestep);

  //   if (timestep/10 < 1) strcat(anameloc,"0");
     
  strcat(anameloc,aux);
  strcat(anameloc,".dat");
  std::cout<<"name of currentfile: "<<anameloc<<std::endl;
  

  //flowdatafile.open(aname, std::ios::in | std::ios::binary); /* Open as a binary file */

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

  //With this we can printout the pressure data in a file with the format of a unverg file.
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
			testflowf << std::setiosflags(std::ios::uppercase | std::ios::scientific) << " " << flowdata_[0][i];
		testflowf << std::endl;
	}
	testflowf.close();*/
}


void AcouFlowNoise::SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
				     const Boolean reset)
{
  ENTER_FCN( "AcouFlowNoise::SolveStepTrans" );

  lasttimecalc_= asteptime;
  Boolean Recalc=FALSE;

  if (laststepcalc_==kstep && kstep!=0) Recalc=TRUE;
  else laststepcalc_= kstep;

  Double * ptsol;
  Integer job = 3;

  //perform predictor step
  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  
  TS_alg_->Predictor(solhelp->GetAlgSysVector());

  if (kstep==1)
    {
      job = 1;
      assemble_->AssembleMatrices(level);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      algsys_->InitRHS();
      assemble_->AssembleSrcRHS(level,lasttimecalc_);

      ComputeRHS(lasttimecalc_);
      TS_alg_->UpdateRHS();
    }
  else if (reset)
    {
      job    = 1;

      algsys_->InitMatrix(SYSTEM);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      algsys_->InitRHS();
      assemble_->AssembleSrcRHS(level,lasttimecalc_);
      ComputeRHS(lasttimecalc_);
      TS_alg_->UpdateRHS();
    }
  else
    {
      job    = 3;
      algsys_->InitRHS();
      assemble_->AssembleSrcRHS(level,lasttimecalc_);
      ComputeRHS(lasttimecalc_);
      TS_alg_->UpdateRHS();
    };

  SetBCs(level, lasttimecalc_);

#ifdef USE_OLAS
  algsys_->SetupPrecond(job);
#else
  algsys_->CalcPrecond(job);
#endif

  algsys_->Solve();

  // Save solution
  ptsol = algsys_->GetSolutionVal();
  sol_->SetAlgSysDataPointer(ptsol);
  
  //perform corrector step 
  TS_alg_->Corrector(solhelp->GetAlgSysVector());
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
      
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
      outFile_->WriteNodeSolutionTransient(*solTransient,laststepcalc_,lasttimecalc_);
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
