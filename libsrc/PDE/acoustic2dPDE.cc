#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acoustic2dPDE.hh"
#include "actimeerror.hh"
#include "acspaceadapt.hh"
#include "outUnverg.hh"
#include "outGMV.hh"
#include "forms_header.hh"
 
namespace CoupledField
{

Acoustic2dPDE::Acoustic2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::Acoustic2dPDE " << std::endl;
#endif

  dofspernode_=1;
  ptgrid_=aptgrid; 

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

  s_oldold_.Resize(size_);
  s_oldold_.Init(0);

  conf->get("alpha_NM",alpha_,"Acoustic"); 
  conf->get("beta_NM",beta_,"Acoustic"); 
  conf->get("gamma_NM",gamma_,"Acoustic");

  conf->getsubdompde(subdoms_,"Acoustic");
  ReadBCs("Acoustic");
}

void Acoustic2dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, Integer &maxnumit, Integer &numeqcoarse)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SpecifySolver" << std::endl;
#endif

  conf->get("eps",eps,"Acoustic"); // relative accuracy in the precond. energy
  conf->get("dampiter",dampiter,"Acoustic"); // damping parameter for Jacobi, SSOR
  conf->get("maxnumit",maxnumit,"Acoustic"); // max number of iterations
  conf->get("solvertype",solvertype,"Acoustic"); // Richardson or CG
  conf->get("precondtype", precondtype, "Acoustic"); //ID or MG
  conf->get("numeqcoarse",numeqcoarse,"Acoustic"); // number of equation for coarsing

}


void Acoustic2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;
  matrix_factor_[1] = 0.0;
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 1.0*a0_;
}

void Acoustic2dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets,
Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SpecifyMatrices" << std::endl;
#endif

  matrixtype = 1;     // NOCLASS = 0
                      // RSPARSE = 1
                      // CSPARSE = 2
                      // RBLOCK  = 3
                      // CBLOCK  = 4
                      // RFULL   = 5
                      // CFULL   = 6
                      // MIXED   = 7

  /* matrixsystype: NOTYPE     = 0
                    SYSTEM     = 1
                    STIFFNESS  = 2
                    DAMPING    = 3
                    CONVECTION = 4
                    MASS       = 5
  */

  matrixsystype[0] = 1;   // memory for the system matrix
  matrixsystype[1] = 2;   // memory for the stiffness matrix
  matrixsystype[2] = 0;   // memory for the damping matrix
  matrixsystype[3] = 0;   // memory for the convection matrix
  matrixsystype[4] = 5;   // memory for the mass matrix

  graphtype = 1; // NOGRAPH = 0
                 // NODE   = 1
                 // EDGE   = 2
                 // FACE   = 3
                 // VOLUME = 4

  numdofpernode  = 1;
  numdirichlets = 1;
  numconstraints = 0;
}

void Acoustic2dPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;
  Point2D * ptCoord; 

  BaseElem * ptEl;

  Integer matrix_stiff=2;
  Integer matrix_mass=5;

  Vector<Double> coeffm, coeffst;
  CalcCoeff(coeffm, coeffst);

  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl=elemssd[j]->ptElem;

	  BaseForm<Point2D> * bilinear_stiff = new LaplaceInt<Point2D>(ptEl,1);
	  BaseForm<Point2D> * bilinear_mass  = new MassInt<Point2D>(ptEl,1);

	  Integer ii;
	  Integer elsize=(elemssd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(elemssd[j]->connect)[ii];
  
	  ptCoord=new Point2D[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level); 

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);

	  elemmat*=coeffst[i];

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;

	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;


#endif     

	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_stiff);

	  // mass part
	  bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

	  elemmat*=coeffm[i];

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;

	  (*debug) << elemmat << std::endl;
#endif
      
	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_,matrix_mass);

  
	  delete bilinear_stiff;
	  delete bilinear_mass;
	  delete [] ptCoord;
	}
    }
#ifdef TRACE
  (*trace) << "Leaving Acoustic2dPDE::SetupMatrices" << std::endl;
#endif
}

void Acoustic2dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val, valueTF;

  //system matrix: id = 1
  Integer matrix_id = 1;

  val=ptTimeFunc_->TimeFuncAtTime(atime,level);

  // set dirichlet boundary conditions
  Integer i,j=0;
  std::list<Integer> nodes_hd;
  Integer sizebc=bcs_hd_.size();
  for (i=0; i< bcs_hd_.size(); i++)
    {
      nodes_hd=ptBCs->GetNodesLevel(bcs_hd_[i]);
   
      for (std::list<Integer>::const_iterator p=nodes_hd.begin(); p!=nodes_hd.end(); p++, j++)
	{
	  node=*p;
	  if (update==1)
	    {
	      ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
	    }
	  else
	    {
    
	      ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
	    }
	}  
    }

#ifdef TRACE
  (*trace) << "leaving Acoustic2dPDE::SetBCs" << std::endl;
#endif
}

void Acoustic2dPDE::ComputeRHS(const Double atime, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::ComputeRHS" << std::endl;
#endif
  Integer n;
  Integer node;
  Integer matrix_id;
  Vector<Double> coeffMass;
  std::list<Integer> nodes_hd;

  matrix_id = 5;
  coeffMass=sol_old_*a0_+sol_der1_old_*a2_+sol_der2_old_*a3_;

  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffMass.get());

  //  Double  val=ptTimeFunc_->TimeFuncWaveSt(atime);
  // Double val=ptTimeFunc_->TimeFuncAtTime(atime,0);
  // if (!ptBCs) Error("You didn't provide pointer to BCs in function ComputeRHS");

  //      nodes_hd=ptBCs->GetNodesLevel("vp-load");

  //   for (std::list<Integer>::const_iterator p=nodes_hd.begin(); p!=nodes_hd.end(); p++)
  //   {
  //     node=*p;
  //    val = 1;
  //     std::cerr << " node: " << node << " val: " << val << " time: " << atime << std::endl;
  //     ptalgsys_->AddRHS(val,0,as_sysid_);
  //    ptalgsys_->AddRHS(val,167,as_sysid_);
  //   }        
   
}

void Acoustic2dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  SetupMatrices(level);
  SetBCs(ptBCs,level,update,0);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
}

void Acoustic2dPDE::SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepTrans" << std::endl;
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

  std::cout << "maxnode:" <<  ptgrid_->GetMaxnumnodes(level) << " level:" << level << std::endl;

  // calculation of derivatives of solution
  CalcDerSol(); 

}

void Acoustic2dPDE::SolveStepTransNewMesh(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepTransNewMesh" << std::endl;
#endif

  lasttimecalc_= asteptime;
  Boolean Recalc=TRUE;

  //   Char * name="testref";
  //   WriteResults * ptOut=new WriteResultsGMV(name);
  //   ptOut->Init(ptgrid_);
  //   ptOut->WriteGrid(0);
  //   if (ptOut) delete ptOut; 

  Double * ptsol;

  Integer update,job;

  update = 0;
  job = 1;
  SetupMatrices(level);
  ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);

  if (kstep) {
    ptalgsys_->ResetRHS(as_sysid_);
    ComputeRHS(lasttimecalc_);
  }
  std::cerr << " we compute RHS " << std::endl;

  SetBCs(ptBCs,level,update,lasttimecalc_);
  std::cerr << " we set BC " << std::endl;
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
  ptsol = ptalgsys_->GetSolution(as_sysid_);
  std::cerr << " we get solution " << std::endl;

  // save solution
  Integer i;
  for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
    sol_[i]=ptsol[i];

  std::cerr << sol_ << std::endl;
  // calculation of derivatives of solution
  CalcDerSol();
}

void Acoustic2dPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::WriteResultsInFile" << std::endl;
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

void Acoustic2dPDE :: CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::CalcParameters" << std::endl;
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

void Acoustic2dPDE :: CalcDerSol()
{
#ifdef TRACE
  (*trace) << " entering  Acoustic2dPDE :: CalcDerSol() " << std::endl;
#endif

  sol_der2_=(sol_ - sol_old_)*a0_ - (sol_der1_old_)*a2_ - sol_der2_old_*a3_;
  sol_der1_=sol_der1_old_+sol_der2_old_*a6_+sol_der2_*a7_;
}

void Acoustic2dPDE::SaveSolAsPrevStep()
{
#ifdef TRACE
  (*trace) << " entering  Acoustic2dPDE::SaveSolAsPrevStep() " << std::endl;
#endif

  sol_old_=sol_;
  sol_der1_old_=sol_der1_;
  sol_der2_old_=sol_der2_;  
}

Double Acoustic2dPDE::CalcEnergyNorm()
{
 Double help1, help2;
 help1=ptalgsys_->CalcEnergyNorm(0,0,2,sol_.get());
 help2=ptalgsys_->CalcEnergyNorm(0,0,5,sol_.get());
 
 return sqrt(help1+help2);
}

void Acoustic2dPDE::CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff)
{
  if (!MatFile_) Error("You didn't specialize material file. Check your config-file.");

  coeffmass.Resize(subdoms_.size());
  coeffstiff.Resize(subdoms_.size());
 
  Integer i,matnum;
  for (i=0; i<subdoms_.size(); i++)
    {
      conf->get(subdoms_[i],matnum,"list_subdomains");

      // read density and compress with material number matnum
      Double density, compress;
      MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid");

      coeffmass[i]  = density*density/compress;
      coeffstiff[i] = density;
    }
}

void Acoustic2dPDE::RestoreSol()
{
#ifdef TRACE
  (*trace) << " entering Acoustic2dPDE::RestoreSol() " << std::endl;
#endif

  Vector<Double> help;
 
  help=sol_old_; 
  ptgrid_->ProlongSol(sol_old_,help,0);
  sol_old_=help;

  help=sol_der1_old_;
  ptgrid_->ProlongSol(sol_der1_old_,help,0);
  sol_der1_old_=help;
  help=sol_der2_old_;
  ptgrid_->ProlongSol(sol_der2_old_,help,0);
  sol_der2_old_=help;

  Integer sizesol=sol_old_.size();
  sol_.Resize(sizesol);
  sol_der1_.Resize(sizesol);
  sol_der2_.Resize(sizesol);
}

TimeErrorEstimator * Acoustic2dPDE::CreatePtTimeError()
{
  return ptTimeError_=new AcousticTimeErrorEstimator(this);
}

SpaceErrorEstimator * Acoustic2dPDE::CreatePtSpaceError()
{
  return ptSpaceError_=new AcousticSpaceErrorEstimator(this,ptgrid_);
}

Acoustic2dPDE::~Acoustic2dPDE()
{
  if (ptTimeError_) delete ptTimeError_;
  if (ptSpaceError_) delete ptSpaceError_;
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !! old stuff, which is not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void Acoustic2dPDE :: CalculationDerivativesSol(const Boolean Recalc)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::CalculationDerivativesSol" << std::endl;
#endif
  
  if (Recalc) std::cout << laststepcalc_ << std::endl;

  if (!Recalc)
  { sol_der2_old_=sol_der2_;
    sol_der1_old_=sol_der1_; 
  }

  if (!Recalc)
  sol_der2_=(sol_-sol_old_)*a0_-sol_der1_old_*a2_-sol_der2_old_*a3_;
  else
  sol_der2_=(sol_-s_oldold_)*a0_-sol_der1_old_*a2_-sol_der2_old_*a3_;

  sol_der1_=sol_der1_old_+sol_der2_old_*a6_+sol_der2_*a7_;


  if (!Recalc) { s_oldold_=sol_old_;}
  sol_old_=sol_;
}

void Acoustic2dPDE::CalcThirdDerivateFromEquation(Vector<Double> & result)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::CalcThirdDerivateFromEquation" << std::endl;
#endif

  Double * ptres;

  Double mat_factor[4];
  Integer update,job;

  mat_factor[0] = 0.0;
  mat_factor[1] = 0.0;
  mat_factor[2] = 0.0;
  mat_factor[3] = 1.0;

  job    = 1;

  ptalgsys_->ResetRHS(as_sysid_);
  ptalgsys_->ResetMatrix(0,0,1);
  ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,mat_factor);

  Integer matrix_id;
  Vector<Double> coeffMass;

  matrix_id = 2;
  Vector<Double> help=-sol_der1_;

  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id, help.get());

  update=1; // use boundary condition at this time step

  //  SetBCs(ptBCs,level,update,lasttimecalc_);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
  ptres = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Integer level=0; //   !!!!!!!!!!!!!!

  Vector<Double> transres(ptgrid_->GetMaxnumnodes(level), ptres);
  result=transres;

  (*infofile) << " result " << transres << std::endl;
}

} // end of namespace

// void Acoustic2dPDE::SetupMatrices(const Integer level)
// {
// #ifdef TRACE
//   (*trace) << "entering Acoustic2dPDE::SetupMatrices" << std::endl;
// #endif
 
//   Vector<Double> coeffm, coeffst;
//   CalcCoeff(coeffm, coeffst);

//  Integer i;
//  for (i=0; i<subdoms_.size(); i++)
// {
//   PutElemMatInAlgSys putelmatalgsys(ptalgsys_,ptgrid_,coeffm[i],coeffst[i],as_sysid_,level);

//   ptgrid_->forEachElemSd(putelmatalgsys,subdoms_[i]);
// }
// }

/*
void PutElemMatInAlgSys::operator()(Elem t)
{
  Matrix<Double> elemmat;

  BaseForm<Point2D> * bilinear_stiff = new LaplaceInt<Point2D>(t.ptElem,1);
  BaseForm<Point2D> * bilinear_mass  = new MassInt<Point2D>(t.ptElem,1);

  Point2D * ptCoord=new Point2D[t.connect.size()];
  ptgrid_->GetCoordNodesElem(t.connect,ptCoord,level_);

  // stiffness part
  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);

  elemmat*=coeffs_;

#ifdef DEBUG
      (*debug) << "Stiffnessmatrix, ElementNumber  "  << std::endl;

      (*debug) << elemmat << std::endl;
#endif

  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), t.connect.get(), t.connect.size(), sysid_, sysid_, matrix_stiff_);

      // mass part
  bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

  elemmat*=coeffm_;

#ifdef DEBUG
      (*debug) << "Massmatrix, ElementNumber  "  << std::endl;

      (*debug) << elemmat << std::endl;
#endif

  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), t.connect.get(), t.connect.size(), sysid_, sysid_,matrix_mass_);

  delete bilinear_stiff;
  delete bilinear_mass;
  delete [] ptCoord;
}
*/
