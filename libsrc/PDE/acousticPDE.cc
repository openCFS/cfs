#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
#include "newmark.hh"
#include "newmarkdamp.hh"


namespace CoupledField
{

AcousticPDE::AcousticPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
			 WriteResults *aptOut)
:BasePDE(aptgrid,aptbcs,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::AcousticPDE " << std::endl;
#endif

  dofspernode_=1;

  pdename_    ="acoustic";
  pdematerialclass_ = "fluid";

  laststepcalc_=0;

  conf->getsubdompde(subdoms_,pdename_);

  AssignPDENodeNumbers(Mesh2PDENode_, PDE2MeshNode_, subdoms_);  
  NumPDENodes_ = PDE2MeshNode_.size();

  size_ = NumPDENodes_;

  sol_.reshape(dofspernode_, NumPDENodes_);
  sol_.init();

  //check, if problem is axisymmetric
  isaxi_ = FALSE;
  std::string subtype;
  conf->ifget("subtype",subtype,pdename_);
  if (subtype == "axi")
    isaxi_ = TRUE;

  with_absBCs_=FALSE;
  std::string absBCs="no";
  conf->ifget("absorbingBCs",absBCs,pdename_);
  if (absBCs == "yes")
      with_absBCs_ = TRUE;

  with_fracdamping_=FALSE;
  std::string frac_damping ="no";
  conf->ifget("frac_damping",frac_damping,pdename_);

  if (frac_damping == "yes")
    {
       with_fracdamping_ = TRUE;
       conf->get("frac_memory",frac_memory_,pdename_);
      (*infofile) << "\n Attenuation according to power law, number of memory is " << frac_memory_ 
		  << std::endl << std::endl;
    }

  ReadBCs(pdename_);

  if (analysistype_==HARMONIC)
    {
      conf->get("frequency", freq_, pdename_);
      solIm_.reshape(dofspernode_, NumPDENodes_);
      solIm_.init();
    }
}


void AcousticPDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE ::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_ = RSCALAR;
  GraphType_  = NODEGRAPH; 

  if (analysistype_==HARMONIC)
    SystemMatrix_     = TRUE;
  else
    {
      SystemMatrix_     = TRUE;
      StiffnessMatrix_  = TRUE;
      MassMatrix_       = TRUE;
      
      if (with_absBCs_ || with_fracdamping_)
	DampingMatrix_  = TRUE;
    }
}


void AcousticPDE::ComputeRHS(const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::ComputeRHS" << std::endl;
#endif
  
}

void AcousticPDE::SolveStepTrans(const Integer kstep, const Double asteptime, 
				   const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::SolveStepTrans" << std::endl;
#endif

  lasttimecalc_= asteptime;
  Boolean Recalc=FALSE;

  if (laststepcalc_==kstep && kstep!=0) Recalc=TRUE;
  else laststepcalc_= kstep;

  Double * ptsol;
  Integer update,job;

  //perform predictor step
  TS_alg_->Predictor(sol_);

  if (kstep==0)
    {
      update = 0;
      job = 1;
      SetupMatrices(level);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
      TS_alg_->UpdateRHS();
    }
  else if (reset)
    {
      update = 1;
      job    = 1;

      algsys_->InitRHS();
      algsys_->InitMatrix(SYSTEM);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      ComputeRHS(lasttimecalc_);
      TS_alg_->UpdateRHS();
    }
  else
    {
      update = 1;
      job    = 3;
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
      TS_alg_->UpdateRHS();
    };

  SetBCs(level,update,lasttimecalc_);
  algsys_->CalcPrecond(job);
  algsys_->Solve();
  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer k = 0;
  
  for (Integer i=0; i<NumPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      sol_[dim][i] = ptsol[k++];

  //perform corrector step  
  TS_alg_->Corrector(sol_);
}

void AcousticPDE::SolveStepHarmonic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::SolveStepHarmonic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 4;
  
  Double * ptsol;

  //compute and assemble element matrices
  SetupMatrices(level);
  
  // calculate source term
  //  SetupRHS(level);
  
  //account for bcs
  SetBCs(level,update,0);

  algsys_->CalcPrecond();    
  algsys_->Solve();
  
  ptsol = algsys_->GetSolutionVal();

  // save solution
  for (Integer i=0; i<NumPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      {
	sol_[dim][i] = ptsol[i];
	//sol_[dim][i] = ptsol[2*i];
	//	solIm_[dim][i] = ptsol[2*i+1];
      }
  
}


void AcousticPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::WriteResultsInFile" << std::endl;
#endif

  Array<Double> sol_mesh, solder1_mesh, solder2_mesh, solIm_mesh;
  Array<Double> sol_der1Array, sol_der2Array;

  if (analysistype_==HARMONIC)
    {
      TransformNodeSolution(sol_mesh,sol_,PDE2MeshNode_);
      TransformNodeSolution(solIm_mesh,solIm_,PDE2MeshNode_);      
      OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential, cw realpart,");
      OutFile_->WriteNodeSolution(solIm_mesh,laststepcalc_,lasttimecalc_,"fluid potential, cw imagpart, ");
    }
  else
    {  
      sol_der1Array = getS1();
      sol_der2Array = getS2();

      
      TransformNodeSolution(sol_mesh,sol_,PDE2MeshNode_);
      TransformNodeSolution(solder1_mesh,sol_der1Array,PDE2MeshNode_);
      TransformNodeSolution(solder2_mesh,sol_der2Array,PDE2MeshNode_);
      
      if (OutFile_->IsGMV())
	{
	  OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"vp");
	  //       OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"vp_1der");
	  //       OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"vp_2der");
	}
      else
	{
	  OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential");
	  OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ");
	  //      OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
	}
    }
}

void AcousticPDE :: InitTimeStepping(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::InitTimeStepping" << std::endl;
#endif

  if (with_fracdamping_)
    {
      //currently the parameter y is taken from the first subdomain
      //=> currently just one subdomain makes sense
      Double y = materialData_[0].GetDampingBeta();
      TS_alg_ = new NewmarkDamp(pdename_, algsys_, dofspernode_, NumPDENodes_, DampingMatrix_,
			      frac_memory_,y);
    }
  else
    TS_alg_ = new Newmark(pdename_, algsys_, dofspernode_, NumPDENodes_, DampingMatrix_);

  TS_alg_->Init(matrix_factor_, dt);

}


void AcousticPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
    (*trace) << "entering AcousticPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;
  Matrix<Double> ptCoord;


  BaseFE * ptEl;
  Double coeffstiff, coeffmass, coeffdamp;
  Vector<Integer> connecth, connect_PDE;
  std::vector<Elem*> elemssd;

  //for harmonic case
  std::vector<Double> harmVec;

  Integer i, j;
 
  for (i=0; i<subdoms_.size(); i++)
    {
      //compute material coefficient
      const Double density = materialData_[i].GetDensity();
      const Double compressibility = materialData_[i].GetCompressibility();
      const Double alpha0 =  materialData_[i].GetDampingAlfa();
      const Double y =  materialData_[i].GetDampingBeta();

      coeffmass  = density*density/compressibility;
      coeffstiff = density;
      coeffdamp  = 2*density*alpha0*sqrt(density/compressibility)/sin((y-1)*PI/2);

      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl = elemssd[j]->ptElem;
    
	  BaseForm * bilinear_mass  = new MassInt(ptEl, coeffmass, isaxi_);
	  BaseForm * bilinear_stiff = new LaplaceInt(ptEl, coeffstiff, isaxi_);

	  connecth=elemssd[j]->connect;
	  GetElemCoords(connecth, ptCoord, level); 

	  // CHANGE connecth
	  Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif     
	  // store stiffness matrix to special vector needed by algebraic system
	  // for harmonic analysis
	  Integer k=0;
	  if(analysistype_==HARMONIC)
	    algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), SYSTEM);
// 	    {
// 	      harmVec.resize(2*connecth.size()*connecth.size());
// 	      for(Integer iii=0; iii<elemmat.size_row(); iii++)
// 		for(Integer jjj=0; jjj < elemmat.size_row(); jjj++)
// 		  {    
// 		    harmVec[k] = elemmat[iii][jjj];
// 		    k++;
// 		  }
// 	    }
	  else
	    algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), STIFFNESS);

	  // mass part
	  bilinear_mass->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  if(analysistype_==HARMONIC)
	    {
	      elemmat *= -4*PI*PI*freq_*freq_;
	      *data << elemmat << std::endl;
	      algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), SYSTEM);
	    }
// 	    {
// 	      if (k!=elemmat.size_row()*elemmat.size_col())
// 		Error("k is wrong!!!!!!!!!!!!!!!!!!!!!",__FILE__,__LINE__);

// 	      for(Integer iii=0; iii<elemmat.size_row(); iii++)
// 		for(Integer jjj=0; jjj < elemmat.size_row(); jjj++)
// 		  {    
// 		    harmVec[k] =  -4*PI*PI*freq_*freq_*elemmat[iii][jjj];
// 		    // k initially  set in for loop above!!
// 		    k++;
// 		  }
// 	      algsys_->SetElementMatrix(&harmVec[0], connect_PDE.get(), connect_PDE.size(),SYSTEM);
// 	    }
	  else
	    algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), MASS);

	  //Damping part
	  if (with_fracdamping_ && analysistype_!=HARMONIC)
	    {
	      BaseForm * bilinear_damp  = new MassInt(ptEl, coeffdamp, isaxi_);
	      bilinear_damp->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	      (*debug) << "Damping matrix, ElementNumber  " << i << std::endl;
	      (*debug) << elemmat << std::endl;
#endif

	      algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), DAMPING);
	      delete bilinear_damp;
	    }

	  delete bilinear_stiff;
	  delete bilinear_mass;

	}
    }


  // BEGIN DAMPING MATRIX PART: Absorbing boundaries
  if (with_absBCs_ && analysistype_!=HARMONIC) {
    std::vector<Elem*>  DomainBnd;
    Double coeffdamp;
 
    conf->getliststr("bnd_for_absBCs",bnd_absBCs_,"acoustic");
    DomainBnd=ptBCs_->getFacesBC(bnd_absBCs_[0],level);

    for (j=0; j< DomainBnd.size(); j++)
      {
	const Double density = materialData_[0].GetDensity();
	const Double compressibility = materialData_[0].GetCompressibility();
	coeffdamp = density/((sqrt(compressibility/density)));
      
	ptEl=DomainBnd[j]->ptElem;
	// Here MassInt is used to calculate the damping matrix from the surface elements
	BaseForm * linear_damp = new MassInt(ptEl,coeffdamp,isaxi_);

	Integer ii;
	Integer elsize=(DomainBnd[j]->connect).size();
	connecth.Resize(elsize);
	for (ii=0; ii<elsize; ii++)
	  connecth[ii]=(DomainBnd[j]->connect)[ii];

	GetElemCoords(connecth, ptCoord, level);
	// CHANGE connecth
	Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

	linear_damp->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	(*debug) << "Connection array  " << std::endl;
	(*debug)  << connecth  << std::endl;

	(*debug) << "Dampingmatrix, ElementNumber  " <<   j << std::endl;
	(*debug) << elemmat << std::endl;
#endif


	algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), DAMPING);

	delete linear_damp;

      }
  }
  // END DAMPING MATRIX PART



#ifdef TRACE
    (*trace) << "Leaving AcousticPDE::SetupMatrices" << std::endl;
#endif
  }

#ifdef ADAPTGRID
Boolean AcousticPDE :: TestError(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::TestError" << std::endl;
#endif

  if (analysistype_!=HARMONIC)
    return TestError(level);

  //! this part is only for harmonic analysis
  if (!ptError_)
    ConstructorError();

  //Berechnung der Fehlerkarte
  Double           totalErr;
  Vector<Double>   solVecRe; // transform from Array to Vector format
  Vector<Double>   solVecIm; // - " -
  Integer          i,ssize;


  ssize = sol_.size();
  solVecRe.Resize(ssize);
  solVecIm.Resize(ssize);

  for (i=0; i<ssize; i++)
    {
      solVecRe[i] = sol_[0][i];
      solVecIm[i] = solIm_[0][i];
    }

  ptError_->CalcErrorMapHarmonic(solVecRe,solVecIm,subdoms_,
				 ptgrid_,errorMap_,totalErr,level);

  (*infofile) << " space error: " << totalErr <<
    " tolerance: " << tolSpaceErr_ << std::endl;

  if (totalErr > tolSpaceErr_) return TRUE;
  else return FALSE;
  
}
#endif //#ifdef ADAPTGRID

}


