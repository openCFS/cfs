#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>

#include "acou2dFlowNoise.hh"
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>

namespace CoupledField
{


Acou2dFlowNoise::Acou2dFlowNoise(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut):Acoustic2dPDE(ptalgsys, aptgrid, ptMaterial, aptTimeFunc, aptFileType, aptOut)   
  {
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::Acou2dFlowNoise " << std::endl;
#endif

 ReadBCs("Acoustic");
 preComputeRHS();
}

void Acou2dFlowNoise::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::SetBCs" << std::endl;
#endif

  Integer node;
  Double val, valTF;
  Double time=atime;

  //system matrix: id = 1
  Integer matrix_id = 1;

  // set dirichlet boundary conditions
  //val=ptTimeFunc_->TimeFuncAtTime(time,level);
  Integer i,j=0;

  // set dirichlet boundary conditions
  val=ptTimeFunc_->TimeFuncAtTime(time,level);

// if (time>1e-5)
//   val=0.0;
//val = 0.0;

  //Double val_hd = 0.0;
  //std::list<Integer> nodes;
  Integer sizebc=bcs_hd_.size();
  std::list<Integer> nodes_hd; 

  //! homogeneous dirichlet BCs
  for (i=0; i< bcs_hd_.size(); i++)
    {
#ifdef DEBUG
       (*debug) << " Dirichlet BC" << std::endl;
#endif
     nodes_hd=ptBCs->GetNodesLevel(bcs_hd_[i]);

     for (std::list<Integer>::const_iterator p=nodes_hd.begin(); p!=nodes_hd.end(); p++, j++)

	{
	  node=*p;
	  if (update==1)
	    {
#ifdef DEBUG
       (*debug) << " node: " << node << " val: " << val << std::endl;
#endif
            ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
	    }
	  else
          {
#ifdef DEBUG
      (*debug) << " node: " << node << " val: " << val << std::endl;
#endif
	   ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
	    }
	}
    }

  //! inhomogeneous dirichlet BCs
//  valTF=ptTimeFunc_->TimeFuncAtTime(atime,level);
//  for (i=0; i<bcs_id_.size(); i++)
//    {
//      nodes=ptBCs->GetNodesLevel(bcs_id_[i], level);

//      val=val_id_[i]*valTF;

//      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
//	{
//	  node=*p;
//          if (update==1)
//            {
//              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
//            }
//          else
//            {
//              ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,as_sysid_, matrix_id);
//            }
//	}
//    }


#ifdef TRACE
  (*trace) << "leaving Acou2dFlowNoise::SetBCs" << std::endl;
#endif
}


void Acou2dFlowNoise::preComputeRHS()
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::preComputeRHS" << std::endl;
#endif

  SetRHSFnc = FALSE;
  SetRHSFlowSrc = FALSE;
   //  if (conf->is_there("rhs_surfaces")) {
   std::string isthererhs;
   conf->ifget("load_force",isthererhs,"Acoustic");

   if (isthererhs=="FlowSrc" ) {
		conf->getliststr("rhs_surfaces",rhs_surfaces_,"Acoustic");
    SetRHSFlowSrc=TRUE;
  }
}

void Acou2dFlowNoise::ComputeRHS(const Double atime, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::ComputeRHS" << std::endl;
#endif

  Integer matrix_id;
  Vector<Double> coeffMass, coeffDamp;
  Vector<Double> elemvec;
  Integer i;

  Integer level=0;
  // mass matrix part
  matrix_id = 5; // Setting id for Mass matrix
  coeffMass = sol_old_*a0_+sol_der1_old_*a2_+sol_der2_old_*a3_;

  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffMass.get());

  // damping matrix part
  matrix_id = 3; // Setting id for Damping matrix	
  coeffDamp = -sol_der1_old_-sol_der2_old_*a6_;
  
  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffDamp.get());

  matrix_id = 3; // Setting id for Damping matrix	
  coeffDamp = sol_old_*a1_+sol_der1_old_*a2_*a7_+sol_der2_old_*a7_*a3_;
  
  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffDamp.get());

   
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
  Boolean reorg=TRUE;
  static float * NODEDATA=new float[3*size_];
  static int * TOPOLOGYDATA=new int[4*maxnumelem];
  std::vector<Elem*> elemssd;     
  static Integer timestep=1;
  static int auxtime=0;
  Double coeffRHS;
  Integer matnum;
  // Getting coeff. for RHS
  conf->get(subdoms_[(subdoms_.size())-1],matnum,"list_subdomains");
  // read density and compress with material number matnum
  Double density, compress;
  MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid");

  //      coeffRHS  = (density*density/compress);
  coeffRHS  = 1.0;

  //THIS IS A TEMPORAL OPTION WHEN GENERATING GRAD P FOR A 2 TIMES COARSER GRID
  //       Integer ctrelem=1;
  //       std::string::size_type pos=0;     
  //       Double gradAveX, gradAveY, gradsumX, gradsumY;
  //       Char * aux=new Char[2];
  //     Char * gradfileTS=new Char[30];
  //   sprintf(aux,"%i",timestep);
  //   sprintf(gradfileTS,"gradfile");
  //   strcat(gradfileTS,aux);
  //   strcat(gradfileTS,".dat");
  //   std::cout<<"name of currentfile: "<<gradfileTS<<std::endl;    

  //Openening file for saving average of gradient

  //std::ofstream gradfile;
  //gradfile.open(gradfileTS, std::ios::out);
  //      gradfile.seekg(pos, std::ios::beg);     
  // END OF OPENING GRADFILE


  // Reorganizing grid info so that exchange through MpCCi can be performed

  if (reorg)
    {
      for (i=0; i<subdoms_.size(); i++)
	{
	  ptgrid_->GetElemSD(elemssd,subdoms_[i],level);
	  int k=0;
	  for (j=0; j< elemssd.size(); j++)
	    {
	      elsize=(elemssd[j]->connect).size();
	      connecth=elemssd[j]->connect;
	      ptCoordNodes=new Point<2>[connecth.size()];
	      ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,level);
	      for (ii=0; ii<elsize; ii++, k++)
		{
		  TOPOLOGYDATA[k]=connecth[ii];
		  NODEDATA[3*(connecth[ii]-1)]=ptCoordNodes[ii][0];		      
		  NODEDATA[3*(connecth[ii]-1)+1]=ptCoordNodes[ii][1];
		  NODEDATA[3*(connecth[ii]-1)+2]= 0.0; // Z comp for the 2d case is zero		      
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
      reorg=FALSE;
    }
      
  std::string flowdata;
  conf->get("acousrc_file",flowdata);
    
  //For the subdomain next to boundary of obstacle (Next2Surf)
  std::vector<Elem*> Next2Surf;  // vector of 2D-elements from mesh-file      
  BaseElem * ptElNext2Surf;      
  std::vector<Elem*> belongSE;

  
  ObstSurf=ptBCs->getEdgesBC(rhs_surfaces_[0],level);
  Next2Surf=ptBCs->getNeighElemsForSurfaces(rhs_surfaces_[1],level);
  ptgrid_->DefineBelonging4Elems(ObstSurf,Next2Surf,belongSE);

  Double valmult;
  valmult = sin(atime*10000*2*3.1416); // value of 10kHz multiplier at the timestep
  //valmult = 1.2*coeffRHS; // 1.2 is to give units to the fluid data and coeffRHS is from inhom. wave eq.
  ///std::cout << "Time: "<< atime << std::endl; 


    ReadFlowData(flowdata.c_str(), timestep, nodedata); // Calls to get data in nodedata

    // PATCH TO NODE 80182 IN THE CHANNEL WITH SQR. OBST. WHICH HAD A NULL VALUE FROM LSTM.

//     for (i=0;i<3;i++)
//       nodedata[i][80182-1]=(nodedata[i][80183-1] + nodedata[i][80181-1])/2;

    // END OF PATCH

    std::cout << "Processing RHS 1D elems... "<< std::endl;
    // This is for the loop over the surface elements
 

    for (j=0; j< ObstSurf.size(); j++)
      { // loop over surface elements
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

	// TEMPORARY ALSO WHEN USING FILES!!!
	// 	      gradN_x_P[0]=nodedata[0][2*j];
	//  	      gradN_x_P[1]=nodedata[0][2*j+1];



	// THIS IS TO BE DONE ONLY ONCE EVERY SIMULATION!!
	// THIS IS JUST FOR SAVING THE AVERAGE OF THE GRADIENTS OF P FOR N ELEMENT SO THAT LATER WE COARSE THE GRID!!!
	// With this we can printout the average of grad(P) data in a file.
	  
	// 	  gradsumX+= gradN_x_P[0];
	// 	  gradsumY+= gradN_x_P[1];
	// 	  if (ctrelem==2)
	// 	    {
	// 	      ctrelem=0;
	// 	      gradAveX=gradsumX/2;
	// 	      gradsumX=0;
	// 	      gradAveY=gradsumY/2;
	// 	      gradsumY=0;
	// 	      gradfile << std::setiosflags(std::ios::uppercase | std::ios::scientific) << gradAveX;
	//               gradfile << std::endl;
	//               gradfile << std::setiosflags(std::ios::uppercase | std::ios::scientific) << gradAveY;
	//               gradfile << std::endl;
	// 	    }
	// 	  ctrelem++;


	// END OF SAVING GRADIENT (THIS IS A TEMPORAL OPTION)



	ptElSurf=ObstSurf[j]->ptElem;
	BaseForm<2> * linear_loaddipole = new LinearForm<2>(ptElSurf);
	  
	connObstSurf=ObstSurf[j]->connect;
	ptCoordNodSurf=new Point<2>[connObstSurf.size()];
 	  
	ptgrid_->GetCoordNodesElem(connObstSurf,ptCoordNodSurf,level);

	linear_loaddipole->CalcElemVector4FlowSrcDip(ptCoordNodSurf, connObstSurf,elemvec,gradN_x_P);
	elemvec*=valmult;
	ptalgsys_->PutElemRHS(elemvec.get(), connObstSurf.get(), connObstSurf.size(), as_sysid_);
	delete linear_loaddipole;
	delete [] ptCoordNodSurf;
	delete [] ptCoordNodBelongSE;

      }


    //THIS IS PART OF GENERATING GRADIENT FILE
    //	  gradfile.close();

    // Variables for ramping
	   Double xfmin, yfmin, xfmax, yfmax, facRampXmin, facRampYmin, facRampXmax, facRampYmax ;
	   Double bndoffsetXmin, bndoffsetYmin, bndoffsetXmax, bndoffsetYmax ;
	   conf->get("xfmin",xfmin,"Acoustic"); // minimum x coord. of fluid domain
	   conf->get("yfmin",yfmin,"Acoustic"); // minimum y coord. of fluid domain	   
	   conf->get("xfmax",xfmax,"Acoustic"); // maximum x coord. of fluid domain
	   conf->get("yfmax",yfmax,"Acoustic"); // maximum y coord. of fluid domain
	   conf->get("facrampXmin",facRampXmin,"Acoustic"); // factor for starting ramping
	   conf->get("facrampYmin",facRampYmin,"Acoustic"); // factor for starting ramping
	   conf->get("facrampXmax",facRampXmax,"Acoustic"); // factor for starting ramping
	   conf->get("facrampYmax",facRampYmax,"Acoustic"); // factor for starting ramping
	   bndoffsetXmin=facRampXmin*xfmin;
	   bndoffsetYmin=facRampYmin*yfmin; 
	   bndoffsetXmax=facRampXmax*xfmax;
	   bndoffsetYmax=facRampYmax*yfmax;
	  
      // Quadrupole computation
      std::cout << "Processing RHS 2D elems... "<< std::endl;
      valmult=coeffRHS;
      
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
	   std::cout<<"Aqui"<<std::endl;

	      linear_load->CalcElemVector4FlowSrcQuad(ptCoordNodes, connecth, nodedata, elemvec);
           elemvec*=valmult;

	   // Ramping before adding elemrhs to global vector to avoid spurious effect at bnd. of fluid dom.

	   for (ii=0; ii<elsize; ii++)
	     {
	       if (ptCoordNodes[ii][0]<bndoffsetXmin)
		 {
		   
		 elemvec[ii]-=elemvec[ii]*(ptCoordNodes[ii][0]-bndoffsetXmin)/(xfmin-bndoffsetXmin);
		 std::cout<<"nodenumber in ramping zone of Xmin: "<<connecth[ii]<<std::endl;
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
	      ptalgsys_->PutElemRHS(elemvec.get(), connecth.get(), connecth.size(), as_sysid_);
	      delete linear_load;
	      delete [] ptCoordNodes;
	    }
	}

    //        auxtime++;
    //        if (auxtime==9)
    //  	{
    //  	  auxtime=0;
    timestep=timestep+1;
    //  	}
      
     

} 



void Acou2dFlowNoise::ReadFlowData(const Char * aname,const Integer timestep, Matrix<Double> &nodedata)
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
  

  sprintf(aux,"%i",timestep);
  strcpy(anameloc,aname);
  strcat(anameloc,aux);
  strcat(anameloc,".dat");
  std::cout<<"name of currentfile: "<<anameloc<<std::endl;
  

  //flowdatafile.open(aname, std::ios::in | std::ios::binary); /* Open as a binary file */

  flowdatafile.open(anameloc,std::ios::in);
  if (!flowdatafile) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
			") Can't open Time-FlowSrc-File:" << anameloc << std::endl;
  exit(1);}

  std::string buffer, buffer2;
  //std::getline(flowdatafile,buffer);
  //std::getline(flowdatafile,buffer2);

  // for maxnumqtts, currently for 2D we are working with 3, p and vx, vy

  Integer ibuf;
  /* Set pointer to beginning of file: */
  std::string::size_type pos=0;
  flowdatafile.seekg(pos, std::ios::beg);
  //if (!pos) For later... when one datafile per timestep
  //{
  //std::getline(flowdatafile,buffer);
  flowdatafile >> ibuf >> maxnumqtts >> maxnumnodes >> ibuf;
  //}
  //std::cout << ibuf << ' ' << maxnumnodes << std::endl;
  nodedata.Resize(maxnumqtts,maxnumnodes);
  Integer i,j;

  Integer buffernodenum = 0;
  for (i=0; i < maxnumnodes; i++)
    {
      flowdatafile >> buffernodenum;
      for (j=0; j < maxnumqtts; j++)
	flowdatafile >> nodedata[j][buffernodenum-1];
    }
	//pos=flowdatafile.tellg(); // Getting the new position  for next time step after reading all data
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


// // THIS IS TEMPORARY FOR GENERATING THE GRAD OF P FOR A 2TIMES COARSER GRID!!
// void Acou2dFlowNoise :: ReadFlowData(const Char * aname,const Integer timestep, Matrix<Double> &nodedata)
// {
// #ifdef TRACE
//   (*trace) << "entering Acou2dFlowNoise ::ReadFlowData" << std::endl;
// #endif


//   std::ifstream flowdatafile;
//   Integer maxnumnodes;
//   Integer maxnumqtts;
//   std::ofstream testflowf;
//   Char * aux=new Char[2];
//   Char * anameloc=new Char[30];
  

//   sprintf(aux,"%i",timestep);
//   strcpy(anameloc,aname);
//   strcat(anameloc,aux);
//   strcat(anameloc,".dat");
//   std::cout<<"name of currentfile: "<<anameloc<<std::endl;
  

//   //flowdatafile.open(aname, std::ios::in | std::ios::binary); /* Open as a binary file */

//   flowdatafile.open(anameloc,std::ios::in);
//   if (!flowdatafile) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
// 			") Can't open Time-FlowSrc-File:" << aname << std::endl;
//   exit(1);}


//   /* Set pointer to beginning of file: */
//   std::string::size_type pos=0;
//   flowdatafile.seekg(pos, std::ios::beg);
//    nodedata.Resize(1,80);
//   Integer i;

//   for (i=0; i < 80; i++)
//     {
 
// 	flowdatafile >> nodedata[0][i];
//     }

//   flowdatafile.close();

// }


void Acou2dFlowNoise::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acou2dFlowNoise::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  SetupMatrices(level, ptBCs);
  SetBCs(ptBCs,level,update,0);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
}

void Acou2dFlowNoise::SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean reset)
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
      SetupMatrices(level,ptBCs);
      ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);
      ptalgsys_->ResetRHS(as_sysid_);
      ComputeRHS(lasttimecalc_,ptBCs);
    }
  else if (reset)
    {
      update = 1;
      job    = 1;

      ptalgsys_->ResetRHS(as_sysid_);
      ptalgsys_->ResetMatrix(0,0,1);
      ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);
      ComputeRHS(lasttimecalc_,ptBCs);
    }
  else
    {
      update = 1;
      job    = 3;
      ptalgsys_->ResetRHS(as_sysid_);
      ComputeRHS(lasttimecalc_,ptBCs);
    };

  SetBCs(ptBCs,level,update,lasttimecalc_);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
  ptsol = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Integer i;
  for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
    sol_[i]=ptsol[i];
  //std::cout << "maxnode:" <<  ptgrid_->GetMaxnumnodes(level) << " level:" << level << std::endl;
  // calculation of derivatives of solution
  CalcDerSol(); 

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
      //      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
    }
}

Acou2dFlowNoise::~Acou2dFlowNoise()
{
 
 
}

} // end of namespace
