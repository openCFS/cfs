#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>

#include "acou2dFlowNoise.hh"
#include <Domain/bcs.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>

namespace CoupledField
{

Acou2dFlowNoise::Acou2dFlowNoise(Grid *aptgrid, BCs *aptbcs, Material *ptMaterial, TimeFunc *aptTimeFunc, 
				 FileType *aptFileType, WriteResults *aptOut)
:Acoustic2dPDE(aptgrid, aptbcs, ptMaterial, aptTimeFunc, aptFileType, aptOut)
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::Acou2dFlowNoise " << std::endl;
#endif

 ReadBCs(pdename_);
 preComputeRHS();
}

Acou2dFlowNoise::~Acou2dFlowNoise()
{
}

void Acou2dFlowNoise::preComputeRHS()
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::preComputeRHS" << std::endl;
#endif

  SetRHSFnc = FALSE;
  SetRHSFlowSrc = FALSE;
  std::string isthererhs;
  conf->ifget("load_force",isthererhs,"Acoustic2d");
  
  if (isthererhs=="FlowSrc" ) 
    {
      conf->getliststr("rhs_surfaces",rhs_surfaces_,"Acoustic2d");
      SetRHSFlowSrc=TRUE;
    }
}

void Acou2dFlowNoise::ComputeRHS(const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::ComputeRHS" << std::endl;
#endif

  Vector<Double> coeffMass, coeffDamp;
  Vector<Double> elemvec;
  Integer i;

  Integer level=0;

  // mass matrix part
  coeffMass = sol_old_*a0_+sol_der1_old_*a2_+sol_der2_old_*a3_;
  algsys_->UpdateRHS(MASS,coeffMass.get());

  // damping matrix part
  coeffDamp = -sol_der1_old_-sol_der2_old_*a6_;  
  algsys_->UpdateRHS(DAMPING,coeffDamp.get());

  coeffDamp = sol_old_*a1_+sol_der1_old_*a2_*a7_+sol_der2_old_*a7_*a3_;  
  algsys_->UpdateRHS(DAMPING,coeffDamp.get());

   
  // get maximum number of elements from grid
  Integer maxnumelem=ptgrid_->GetMaxnumElem(level,subdoms_);

  Double integrShFnc,val,multiplier;
  Point<2> * ptCoordNodes;
  Point<2> * ptCoordNodSurf; // For ObstSurf
  Point<2> * ptCoordNodBelongSE;
  Vector<Integer> connecth;
  Vector<Integer> connObstSurf;	// For ObstSurf
  Vector<Integer> connBelongSE;
  std::vector<Elem*> ObstSurf;  // vector of 1D-elements (ObstSurf) from mesh-file
  BaseElem * ptEl;
  BaseElem * ptElSurf;
  BaseElem * ptElBelongSE;
  Matrix<Double> nodedata; // Where the nodal data is going to be stored
  
  Integer j,ii, elsize = -1;
  std::vector<Elem*> elemssd;     
  static Integer timestep=0;
  static int auxtime=0;
  Double coeffRHS;

  coeffRHS  = 1.0;

  std::string flowdata;
  conf->get("acousrc_file",flowdata);
    
  //For the subdomain next to boundary of obstacle (Next2Surf)
  std::vector<Elem*> Next2Surf;  // vector of 2D-elements from mesh-file      
  BaseElem * ptElNext2Surf;      
  std::vector<Elem*> belongSE;

  
  ObstSurf=ptBCs_->getEdgesBC(rhs_surfaces_[0],level);
  Next2Surf=ptBCs_->getNeighElemsForSurfaces(rhs_surfaces_[1],level);
  ptgrid_->DefineBelonging4Elems(ObstSurf,Next2Surf,belongSE);

  Double valmult;

  valmult = 3.7; // 3.7 is to give units to the fluid data and coeffRHS is from inhom. wave eq.

  // Call to get fluid flow data in nodedata  
  ReadFlowData(flowdata.c_str(), timestep, nodedata); 

// PATCH TO NODE 80182 IN THE CHANNEL WITH SQR. OBST. WHICH HAD A NULL VALUE FROM LSTM.
  for (i=0;i<3;i++)
    nodedata[i][80182-1]=(nodedata[i][80183-1] + nodedata[i][80181-1])/2;
// END OF PATCH

  std::cout << "Processing RHS 1D elems... "<< std::endl;

  // This is for the loop over the surface elements
  for (j=0; j< ObstSurf.size(); j++)
    { 
      ptElBelongSE=belongSE[j]->ptElem;
      // This will be done inside the 2d element next to the 1d element to get gradP at the center
      Integer n=ptElBelongSE->GetNumNodes();      
      Vector<Double> * help=new Vector<Double>[n];
      
      elsize=(belongSE[j]->connect).size();
      connBelongSE.Resize(elsize);
      for (ii=0; ii<elsize; ii++)
	connBelongSE[ii]=(belongSE[j]->connect)[ii];
      
      ptCoordNodBelongSE=new Point<2>[connBelongSE.size()];
      ptgrid_->GetCoordNodesElem(connBelongSE,ptCoordNodBelongSE,level);
      
      std::vector<Double> gradN_x_P; // This is done due to the different parameter type Vector and std::vector
      gradN_x_P.resize(2);// For 2D case only
      
      Vector<Double> JinvX, JinvY;
      Jacobian<2> J;// For 2D case only
      ptElBelongSE->CalcJacobianAtCenter(J,ptCoordNodBelongSE);
      J.GetJinvX(JinvX);
      J.GetJinvY(JinvY);

      // TO COMMENT OUT ONLY WHILE USING FILES WITH GRADIENT VALUE!!!
      // Gradient of P at center by average of value at four nodes of neighbour 2d element
      for (ii=0; ii<n; ii++)
	{  
	  ptElBelongSE->GetGradientShFncAtCenter(help[ii],ii+1);
	  gradN_x_P[0]+=(help[ii]*JinvX)*(nodedata[0][connBelongSE[ii]-1]);	    
	  gradN_x_P[1]+=(help[ii]*JinvY)*(nodedata[0][connBelongSE[ii]-1]);
	}

      ptElSurf=ObstSurf[j]->ptElem;
      BaseForm<2> * linear_loaddipole = new LinearForm<2>(ptElSurf);
	  
      connObstSurf=ObstSurf[j]->connect;

      ptCoordNodSurf=new Point<2>[connObstSurf.size()];
      
      ptgrid_->GetCoordNodesElem(connObstSurf,ptCoordNodSurf,level);
      
      linear_loaddipole->CalcElemVector4FlowSrcDip(ptCoordNodSurf, connObstSurf,elemvec,gradN_x_P);
      elemvec*=valmult;
      algsys_->SetElementRHS(elemvec.get(), connObstSurf.get(), connObstSurf.size());

      delete linear_loaddipole;
      delete [] ptCoordNodSurf;
      delete [] ptCoordNodBelongSE;
      
    }


  // Variables for ramping
  Double xfmin, yfmin, xfmax, yfmax, facRampXmin, facRampYmin, facRampXmax, facRampYmax ;
  Double bndoffsetXmin, bndoffsetYmin, bndoffsetXmax, bndoffsetYmax ;
  conf->get("xfmin",xfmin,"Acoustic2d"); // minimum x coord. of fluid domain
  conf->get("yfmin",yfmin,"Acoustic2d"); // minimum y coord. of fluid domain	   
  conf->get("xfmax",xfmax,"Acoustic2d"); // maximum x coord. of fluid domain
  conf->get("yfmax",yfmax,"Acoustic2d"); // maximum y coord. of fluid domain
  conf->get("facrampXmin",facRampXmin,"Acoustic2d"); // factor for starting ramping
  conf->get("facrampYmin",facRampYmin,"Acoustic2d"); // factor for starting ramping
  conf->get("facrampXmax",facRampXmax,"Acoustic2d"); // factor for starting ramping
  conf->get("facrampYmax",facRampYmax,"Acoustic2d"); // factor for starting ramping
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
	  BaseForm<2> * linear_load = new LinearForm<2>(ptEl);

	  Integer ii;
	  elsize=(elemssd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(elemssd[j]->connect)[ii];

	  ptCoordNodes=new Point<2>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,level);
	  
	  linear_load->CalcElemVector4FlowSrcQuad(ptCoordNodes, connecth, nodedata, elemvec);
	  elemvec*=valmult;

	  // Ramping before adding elemrhs to global vector to avoid spurious effect at bnd. of fluid dom.

	  for (ii=0; ii<elsize; ii++)
	    {
	      if (ptCoordNodes[ii][0]<bndoffsetXmin)
		{
		  
		  elemvec[ii]-=elemvec[ii]*(ptCoordNodes[ii][0]-bndoffsetXmin)/(xfmin-bndoffsetXmin);
		 }
	       
	      else
		if (ptCoordNodes[ii][0]>bndoffsetXmax)
		  elemvec[ii]-=elemvec[ii]*(ptCoordNodes[ii][0]-bndoffsetXmax)/(xfmax-bndoffsetXmax);
	      if (ptCoordNodes[ii][1]<bndoffsetYmin)
		elemvec[ii]-=elemvec[ii]*(ptCoordNodes[ii][1]-bndoffsetYmin)/(yfmin-bndoffsetYmin);
	      else	
		if (ptCoordNodes[ii][1]>bndoffsetYmax)
		  elemvec[ii]-=elemvec[ii]*(ptCoordNodes[ii][1]-bndoffsetYmax)/(yfmax-bndoffsetYmax);
	    }
	   
	  //end ramping

	  //linear_load->CalcElemVector(ptCoordNodes, elemvec); // for setting with homogeneous rhs	      
	  algsys_->SetElementRHS(elemvec.get(), connecth.get(), connecth.size());

	  delete linear_load;
	  delete [] ptCoordNodes;
	  
	}
    }


  timestep=timestep+1;

} 


void Acou2dFlowNoise::ReadFlowData(const Char * aname, Integer timestep, Matrix<Double> &nodedata)
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise ::ReadFlowData" << std::endl;
#endif

  std::ifstream flowdatafile;
  Integer maxnumnodes;
  Integer maxnumqtts;
  std::ofstream testflowf;
  Char * aux=new Char[2];
  Char * anameloc=new Char[30];

  strcpy(anameloc,aname);
 	 
  if (timestep > 199)
    {
      timestep=timestep-200;
      strcat(anameloc,"8");
    }  
  else
    {
      if (timestep > 99)
	{
	  timestep=timestep-100;
	  strcat(anameloc,"7");	
	}
      else
	{
	  strcat(anameloc,"6");	
	}      
    }
  
  sprintf(aux,"%i",timestep);

  if (timestep/10 < 1) strcat(anameloc,"0");
     
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

  nodedata.Resize(maxnumqtts,maxnumnodes);
  Integer i,j;

  Integer buffernodenum = 0;
  for (i=0; i < maxnumnodes; i++)
    {
      flowdatafile >> buffernodenum;
      for (j=0; j < maxnumqtts; j++)
	flowdatafile >> nodedata[j][buffernodenum-1];
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
			testflowf << std::setiosflags(std::ios::uppercase | std::ios::scientific) << " " << nodedata[0][i];
		testflowf << std::endl;
	}
	testflowf.close();*/
}


void Acou2dFlowNoise::SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
				     const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::SolveStepTrans" << std::endl;
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
    sol_[i]=ptsol[i];

}

void Acou2dFlowNoise::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::WriteResultsInFile" << std::endl;
#endif

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"vp");
      OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"vp_1der");
      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"vp_2der");
    }
  else
    {
      OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"fluid potential");
      //OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ");
      //OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
    }
}


void Acou2dFlowNoise::PutExchangeGrid2MpCCI()
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::PutExchangeGrid2MpCCI" << std::endl;
#endif


  Integer maxnumelem=ptgrid_->GetMaxnumElem(actlevel_,subdoms_);

  Point<2> * ptCoordNodes;
  Vector<Integer> connecth;
  Integer i,j,ii;
  Integer elsize = -1;

  static float * NODEDATA=new float[3*size_];
  static int * TOPOLOGYDATA=new int[4*maxnumelem];
  std::vector<Elem*> elemssd;     

  for (i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],actlevel_);
      int k=0;
      for (j=0; j< elemssd.size(); j++)
	{
	  elsize=(elemssd[j]->connect).size();
	  connecth=elemssd[j]->connect;
	  ptCoordNodes=new Point<2>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,actlevel_);
	  for (ii=0; ii<elsize; ii++, k++)
	    {
	      TOPOLOGYDATA[k]=connecth[ii];
	      NODEDATA[3*(connecth[ii]-1)]=ptCoordNodes[ii][0];		      
	      NODEDATA[3*(connecth[ii]-1)+1]=ptCoordNodes[ii][1];
	      NODEDATA[3*(connecth[ii]-1)+2]= 0.0; // z-component for the 2d case is zero		      
	    }
	  delete [] ptCoordNodes;	
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


}
      
} // end of namespace
