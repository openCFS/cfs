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
#ifdef TRACE
  (*trace) << "entering AcouFlowNoise::AcouFlowNoise " << std::endl;
#endif

 ReadBCs(pdename_);
 preComputeRHS();
}

AcouFlowNoise::~AcouFlowNoise()
{
}

void AcouFlowNoise::preComputeRHS()
{
#ifdef TRACE
  (*trace) << "entering AcouFlowNoise::preComputeRHS" << std::endl;
#endif

  //first check if coupling per MpCCI or if a flow-result file is provided
  MpCCI_=FALSE;
  std::string strMpCCI="no";
  conf->ifget("MpCCI",strMpCCI,pdename_);
  if (strMpCCI == "yes")
    {
#ifdef MpCCI
      MpCCI_ = TRUE;
      std::cout << "DO COUPLING via MpCCI" << std::endl;
      MpCCIexch * ptMpCCIexch_ = new MpCCIexch(ptgrid_);
      MpCCInodes_=ptgrid_->GetMaxnumnodes(0);

      ptMpCCIexch_->PutExchangeGrid2MpCCI(subdoms_);
      flowdata_.Resize(3, MpCCInodes_);

      SetRHSFnc = FALSE;
      SetRHSFlowSrc = FALSE;
      std::string isthererhs;
      conf->ifget("load_force",isthererhs,"acoustic");
  
      if (isthererhs=="FlowSrc" ) 
	{
	  std::cout<<"Aqui en flownoise"<<std::endl;  
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
	  conf->getliststr("rhs_surfaces",rhs_surfaces_,"acoustic");
	  SetRHSFlowSrc=TRUE;
	}
    }
}


void AcouFlowNoise::ComputeRHS(const Double atime)
{
#ifdef TRACE
  (*trace) << "entering AcouFlowNoise::ComputeRHS" << std::endl;
#endif

  Vector<Double> coeffMass, coeffDamp;
  Vector<Double> elemvec;
  Integer i;

  Integer level=0;

  // mass matrix part
  coeffMass = sol_old_*a0_+sol_der1_old_*a2_+sol_der2_old_*a3_;
  algsys_->UpdateRHS(MASS,coeffMass.get());

  // damping matrix part
  if (with_absBCs_) 
    {
      
      coeffDamp = -sol_der1_old_-sol_der2_old_*a6_;   
 
      algsys_->UpdateRHS(DAMPING,coeffDamp.get());

      coeffDamp = sol_old_*a1_+sol_der1_old_*a2_*a7_+sol_der2_old_*a7_*a3_;  
       
      algsys_->UpdateRHS(DAMPING,coeffDamp.get());
   }

   
  // get maximum number of elements from grid
  Integer maxnumelem=ptgrid_->GetMaxnumElem(level,subdoms_);

  Double integrShFnc,val,multiplier;
  Matrix<Double> ptCoordNodes;
  Matrix<Double> ptCoordNodSurf; // For ObstSurf
  Matrix<Double> ptCoordNodBelongSE; // For set of elements corresponding to surface elements
  Matrix<Double> deriv;
  Vector<Integer> connecth;
  Vector<Integer> connect_PDE;	// For changing connecth to PDE
  Vector<Integer> connObstSurf;	// For ObstSurf
  Vector<Integer> connBelongSE;
  std::vector<Elem*> ObstSurf;  // vector of 1D-elements (ObstSurf) from mesh-file
  BaseFE * ptEl;
  BaseFE * ptElSurf;
  BaseFE * ptElBelongSE;
  Matrix<Double> flowdata_; // Where the nodal data is going to be stored
  
  Integer j,ii, elsize = -1;
  std::vector<Elem*> elemssd;     
  static Integer timestep=0;
  static int auxtime=0;
 
  std::string flowdata;
  conf->get("acousrc_file",flowdata);
    
  //For the subdomain next to boundary of obstacle (Next2Surf)
  std::vector<Elem*> Next2Surf;  // vector of Volume elements next to surface from mesh-file      
  BaseFE * ptElNext2Surf;      
  std::vector<Elem*> belongSE;
 
 
  ObstSurf=ptBCs_->getEdgesBC(rhs_surfaces_[0],level);
  Next2Surf=ptBCs_->getNeighElemsForSurfaces(rhs_surfaces_[1],level);
  ptgrid_->DefineBelonging4Elems(ObstSurf,Next2Surf,belongSE);
  
  Double valmult;

  valmult = 3.7; // 3.7 is to give units to the fluid data and coeffRHS is from inhom. wave eq.

  std::cout<<"timestep counter in ComputeRHS: "<<timestep<<std::endl;
  
  if(MpCCI_)
    {
#ifdef MpCCI

      std::cout<<"MpCCInodes_: "<<MpCCInodes_<<std::endl;
      MpCCIexch * ptMpCCIexch_ = new MpCCIexch(ptgrid_);
      flowdata_.Resize(1+Dim_, MpCCInodes_);
      ptMpCCIexch_->CouplCompPhase(flowdata_, timestep);
#endif
    }
  else
    // If data from fluid file call to get fluid flow data in flowdata_  
    ReadFlowData(flowdata.c_str(), timestep, flowdata_); 






  std::cout << "Processing RHS 1D elems... "<< std::endl;

  // This is for the loop over the surface elements
  for (j=0; j< ObstSurf.size(); j++)
    { 
      ptElBelongSE=belongSE[j]->ptElem;
      // This will be done inside the 2d element next to the 1d element to get gradP at the center
      Integer n=ptElBelongSE->GetNumNodes(); // This returns number of integration points      

      elsize=(belongSE[j]->connect).size(); // Get element number of nodes 
      connBelongSE.Resize(elsize);
      for (ii=0; ii<elsize; ii++)
	connBelongSE[ii]=(belongSE[j]->connect)[ii];

      ptgrid_->GetCoordNodesElemMat(connBelongSE,  ptCoordNodBelongSE, level);
  
    
      std::vector<Double> gradN_x_P; // This is done due to the different parameter type Vector and std::vector
      gradN_x_P.resize(Dim_);// For 2D case only

      std::vector<Double> LCoord(Dim_,0);
      

      // TO COMMENT OUT ONLY WHILE USING FILES WITH GRADIENT VALUE!!!
      // Gradient of P at center by average of value at four nodes of neighbour 2d element
      //  for (ii=0; ii<n; ii++)
      //     {  
	//ptElBelongSE->GetGradientShFncAtCenter(help[ii],ii+1);

       
	ptElBelongSE->GetGlobDerivShFnc  (deriv, LCoord, ptCoordNodBelongSE /*, Double &    jacDet*/);
	//deriv.Transpose(derivTrans);

	for (ii=0; ii<n; ii++)
	  {  
	    for(i=0;i<Dim_;i++)
	      {
		gradN_x_P[i]+=deriv[ii][i]*(flowdata_[0][connBelongSE[ii]-1]);
	      }
	  }

	ptElSurf=ObstSurf[j]->ptElem;

	BaseForm * linear_loaddipole = new LinearFlowNoiseInt(ptElSurf);
	  
	connObstSurf=ObstSurf[j]->connect;

	ptgrid_->GetCoordNodesElemMat(connObstSurf, ptCoordNodSurf, level);
        std::cout<<"coordinates :"<<ptCoordNodSurf<<std::endl; 
	linear_loaddipole->CalcElemVector4Dip(ptCoordNodSurf, connObstSurf,elemvec,gradN_x_P);
	elemvec*=valmult;



	// CHANGE connecth
	Mesh2PDENode(connect_PDE,connObstSurf,Mesh2PDENode_);
	std::cout<<"connObstSurf :"<<connObstSurf<<std::endl; 

        algsys_->SetElementRHS(&elemvec[0], &connect_PDE[0], connect_PDE.size());

	delete linear_loaddipole;

      }

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



  // Quadrupole computation
  std::cout << "Processing RHS 2D elems... "<< std::endl;
  valmult=-1.0;
      
  for (i=(subdoms_.size())-1; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl=elemssd[j]->ptElem;
	  BaseForm * linear_load = new LinearFlowNoiseInt(ptEl);

	  Integer ii;
	  elsize=(elemssd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(elemssd[j]->connect)[ii];



      //  ptCoordNodes=new Point<2>[connecth.size()];
     Matrix<Double> ptCoordNodes;
	//   ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,level);
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
	Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);
	  //linear_load->CalcElemVector(ptCoordNodes, elemvec); // for setting with homogeneous rhs
	std::cout<<"elemvect quad: "<<elemvec<<std::endl;
	
	  algsys_->SetElementRHS(elemvec.get(), connect_PDE.get(), connect_PDE.size());

	  delete linear_load;
	  // delete [] ptCoordNodes;
	  
	}
    }


   timestep=timestep+1;

} 


void AcouFlowNoise::ReadFlowData(const Char * aname, Integer timestep, Matrix<Double> &flowdata_)
{
#ifdef TRACE
  (*trace) << "entering AcouFlowNoise ::ReadFlowData" << std::endl;
#endif

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
#ifdef TRACE
  (*trace) << "entering AcouFlowNoise::SolveStepTrans" << std::endl;
#endif

  lasttimecalc_= asteptime;
  Boolean Recalc=FALSE;

  if (laststepcalc_==kstep && kstep!=0) Recalc=TRUE;
  else laststepcalc_= kstep;

  Double * ptsol;
  Integer update,job;

  if (kstep==0)
    {
      update = 0;
      job = 1;
      SetupMatrices(level);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
    }
  else if (reset)
    {
      update = 1;
      job    = 1;

      algsys_->InitRHS();
      algsys_->InitMatrix(SYSTEM);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      ComputeRHS(lasttimecalc_);
    }
  else
    {
      update = 1;
      job    = 3;
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
    };

 
  SetBCs(level,update,lasttimecalc_); 
  algsys_->CalcPrecond(job);
  algsys_->Solve();
  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer i;
  for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
    sol_[0][i]=ptsol[i];


}

void AcouFlowNoise::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering AcouFlowNoise::WriteResultsInFile" << std::endl;
#endif
  Integer Dim = 2;

  Array<Double> arraysol,arraysol_der1,arraysol_der2;

  // transform solution vector for acoustic potential
//   arraysol=sol_;
//   arraysol_der1=sol_der1_;
//   arraysol_der2 = sol_der2_;
  
  
  TransformNodeSolution(arraysol,sol_,PDE2MeshNode_);
  TransformNodeSolution(arraysol_der1,sol_der1_,PDE2MeshNode_);
  TransformNodeSolution(arraysol_der2,sol_der2_,PDE2MeshNode_);

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteNodeSolution(arraysol,laststepcalc_,lasttimecalc_,"vp");
//       OutFile_->WriteNodeSolution(arraysol_der1,laststepcalc_,lasttimecalc_,"vp_1der");
//       OutFile_->WriteNodeSolution(arraysol_der2,laststepcalc_,lasttimecalc_,"vp_2der");
    }
  else
    {
      OutFile_->WriteNodeSolution(arraysol,laststepcalc_,lasttimecalc_,"fluid potential");
      //OutFile_->WriteNodeSolution(sol_der1,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ",1);
      //OutFile_->WriteNodeSolution(sol_der2,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ",1);
    }
}

} // end of namespace
