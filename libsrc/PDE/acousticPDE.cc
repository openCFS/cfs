#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

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

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);

  laststepcalc_=0;

  AssignPDENodeNumbers();
  size_ = NumPDENodes_;

  sol_.Resize(size_);
  sol_.Init(0);

  sol_der1_.Resize(size_);
  sol_der1_.Init(0);

  sol_der1_old_.Resize(size_);
  sol_der1_old_.Init(0);

  sol_der1_.Resize(size_);
  sol_der1_.Init(0);

  sol_der2_.Resize(size_);
  sol_der2_.Init(0);

  sol_old_.Resize(size_);
  sol_old_.Init(0);

  sol_der2_old_.Resize(size_);
  sol_der2_old_.Init(0);

  alpha_ = 0.0;
  beta_  = 0.25;
  gamma_ = 0.5;

  //check if integration parameters are defined in conf-file
  conf->ifget("alpha_NM",alpha_,pdename_); 
  conf->ifget("beta_NM",beta_,pdename_); 
  conf->ifget("gamma_NM",gamma_,pdename_);

  with_absBCs_=FALSE;
  std::string absBCs="no";
  conf->ifget("absorbingBCs",absBCs,pdename_);
  if (absBCs == "yes")
    {
      with_absBCs_ = TRUE;
      Error("Currently no absorbing BCs supported",__FILE__,__LINE__);
    }

  SetMatrixFactors();

  //currently not available
  // preComputeRHS();
}


void AcousticPDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE ::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_ = RSCALAR;
  GraphType_  = NODEGRAPH; 

  SystemMatrix_     = TRUE;
  StiffnessMatrix_  = TRUE;
  MassMatrix_       = TRUE;

  if (with_absBCs_)
    DampingMatrix_  = TRUE;
}


void AcousticPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;       // factor for stiffness matrix

  matrix_factor_[1] = 0.0; 
  if (with_absBCs_)
    matrix_factor_[1] = 1.0*a1_; // factor for damping matrix

  matrix_factor_[2] = 0.0;       // factor for convection matrix
  matrix_factor_[3] = 1.0*a0_;   // factor for mass matrix
}


void AcousticPDE::ComputeRHS(const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::ComputeRHS" << std::endl;
#endif

  Integer n;
  Integer node;

  Vector<Double> coeffMass, coeffDamp;
  Vector<Double> elemvec;
  std::list<Integer> nodes_hd;
  Integer i;

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

  // save solution
  ptsol = algsys_->GetSolutionVal();
  Vector<Double> transsol(NumPDENodes_, ptsol);
  sol_=transsol;

  if (InfoPrint)
    (*infofile) << "maxnode:" <<  ptgrid_->GetMaxnumnodes(level) << std::endl;

  // calculation of derivatives of solution
  CalcDerSol(); 
}

void AcousticPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::WriteResultsInFile" << std::endl;
#endif

  Vector<Double> sol_mesh, solder1_mesh, solder2_mesh;  
  TransformNodeSolution(sol_mesh,sol_);
  TransformNodeSolution(solder1_mesh,sol_der1_);
  TransformNodeSolution(solder2_mesh,sol_der2_);

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteSolution(sol_mesh,laststepcalc_,lasttimecalc_,"vp");
      OutFile_->WriteSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"vp_1der");
      OutFile_->WriteSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"vp_2der");
    }
  else
    {
      OutFile_->WriteSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential");
      //      OutFile_->WriteSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ");
      //      OutFile_->WriteSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
    }

}

void AcousticPDE :: CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::CalcParameters" << std::endl;
#endif

  a2_=1.0/(beta_*dt);
  a0_=a2_*(1/dt);
  a1_=gamma_*a2_;
  a3_=0.5/(beta_)-1.0;
  a4_=gamma_/beta_-1.0;
  a5_=0.5*dt*(a4_-1.0);
  a6_=dt*(1-gamma_);
  a7_=gamma_*dt;
}

void AcousticPDE :: CalcDerSol()
{
#ifdef TRACE
  (*trace) << " entering  AcousticPDE :: CalcDerSol() " << std::endl;
#endif

  sol_der2_=(sol_ - sol_old_)*a0_ - (sol_der1_old_)*a2_ - sol_der2_old_*a3_;
  sol_der1_=sol_der1_old_+sol_der2_old_*a6_+sol_der2_*a7_;
}

void AcousticPDE::SaveSolAsPrevStep()
{
#ifdef TRACE
  (*trace) << " entering  AcousticPDE::SaveSolAsPrevStep() " << std::endl;
#endif

  sol_old_=sol_;
  sol_der1_old_=sol_der1_;
  sol_der2_old_=sol_der2_;  
}


void AcousticPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
    (*trace) << "entering AcousticPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;
  Matrix<Double> ptCoord;

  // This is a smaller matrix since it is just for linear 1D/2D elements.
  // Matrix<Double> elemmatbnd;

  BaseFE * ptEl;
  Double coeffstiff, coeffmass;
  Vector<Integer> connecth, connect_PDE;
  std::vector<Elem*> elemssd;

  Integer i, j;
 
  for (i=0; i<subdoms_.size(); i++)
    {
      //compute material coefficient
      const Double density = materialData_[i].GetDensity();
      const Double compressibility = materialData_[i].GetCompressibility();

      coeffmass  = density*density/compressibility;
      coeffstiff = density;

      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl = elemssd[j]->ptElem;
    
	  BaseForm * bilinear_mass  = new MassInt(ptEl, coeffmass);
	  BaseForm * bilinear_stiff = new LaplaceInt(ptEl, coeffstiff);

	  connecth=elemssd[j]->connect;
	  ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level); 

	  // CHANGE connecth
	  Mesh2PDENode(connect_PDE,connecth);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif     

	  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connecth.size(), STIFFNESS);

	  // mass part
	  bilinear_mass->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif
      
	  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connecth.size(), MASS);
  
	  delete bilinear_stiff;
	  delete bilinear_mass;

	}
    }

#ifdef TRACE
    (*trace) << "Leaving AcousticPDE::SetupMatrices" << std::endl;
#endif
  }


}


