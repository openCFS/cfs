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

AcousticPDE::AcousticPDE(Grid * aptgrid, BCs *aptbcs, Material *ptMaterial, TimeFunc *aptTimeFunc, 
		     FileType *aptFileType, WriteResults *aptOut)
:BasePDE(aptgrid,aptbcs,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::AcousticPDE " << std::endl;
#endif

  dofspernode_=1;

  if (!MatFile_) Error("You didn't specified material file. Check your config-file.");
 
  laststepcalc_=0;
  size_=ptgrid_->GetMaxnumnodes(0);

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
    with_absBCs_ = TRUE;

  SetMatrixFactors();
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

void AcousticPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, 
				    Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::SpecifyMatrices" << std::endl;
#endif

  matrixtype = RSCALAR; 
  graphtype  = NODEGRAPH; 

  matrixsystype[0] = SYSTEM;      // memory for the system matrix
  matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix

  matrixsystype[2] = 0;
  if (with_absBCs_)
    matrixsystype[2] = DAMPING;   // memory for the damping matrix

  matrixsystype[3] = 0;           // memory for the convection matrix
  matrixsystype[4] = MASS;        // memory for the mass matrix


  numdofpernode  = dofspernode_;
  numdirichlets  = GetNumRestraints(actlevel_);
  numconstraints = 0;
}


void AcousticPDE::SetAlgSys(const Integer as_sysid)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::SetAlgSys" << std::endl;
#endif

  as_sysid_ = as_sysid;

  //allocate according algebraic system
  algsys_ = new StandardSystem();

  //set solver parameters  
  SetSolverParameters();

  //set the graph type used for the system matrices
  Integer numnode = ptgrid_->GetMaxnumnodes(actlevel_);
  SetupMatrixGraph(numnode,NODEGRAPH);

  //allocate the necessary matrices as well as solver and preconditioner
  CreateMatrices_Solver();
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
  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer i;
  for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
    sol_[i]=ptsol[i];

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

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"vp");
      OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"vp_1der");
      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"vp_2der");
    }
  else
    {
      OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"fluid potential");
      OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ");
      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
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


void AcousticPDE ::CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff, Vector<Double> & coeffdamp)
{
  if (!MatFile_) Error("You didn't specialize material file. Check your config-file.");

  coeffmass.Resize(subdoms_.size());
  coeffstiff.Resize(subdoms_.size());
  if (with_absBCs_) coeffdamp.Resize(subdoms_.size());

  Integer i,matnum;
  for (i=0; i<subdoms_.size(); i++)
    {
      conf->get(subdoms_[i],matnum,"list_subdomains");

      // read density and compress with material number matnum
      Double density, compress;
      MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid");

      coeffmass[i]  = density*density/compress;
      coeffstiff[i] = density;
      if (with_absBCs_)
      coeffdamp[i]  = density/((sqrt(compress/density)));
    }
}

}


