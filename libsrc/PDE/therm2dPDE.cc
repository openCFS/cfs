#include <fstream>
#include <iostream>
#include <string>

#include "therm2dPDE.hh"
#include "forms_header.hh"
 
namespace CoupledField
{

Therm2dPDE::Therm2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * agrid, Material * aMatFile, TimeFunc * aptTimeFunc, FileType * aInFile, WriteResults * aptOut)
:BasePDE(ptalgsys, aMatFile, aInFile, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::Therm2dPDE " << std::endl;
#endif

  dofspernode_ = 1;
  ptgrid_=agrid;

  std::string formulation;
  conf->get("formulation", formulation, "Thermal");
  
  if (formulation=="effmass") formulation_=0;
  else { 
    if (formulation=="effstiff") formulation_=1;
    else Error("This type of formulation is wrong",__FILE__,__LINE__);
  }

  size_ = ptgrid_->GetMaxnumnodes(0);
  sol_.Resize(size_);
  sol_.Init(0);
  sol_der1_.Resize(size_);
  sol_der1_.Init(0);
  sol_old_.Resize(size_);
  sol_old_.Init(0);
  sol_der1_old_.Resize(size_);
  sol_der1_old_.Init(0);
  sol_pred_.Resize(size_);
  sol_pred_.Init(0);

  conf->get("gamma_parabolic",gamma_parab_,"Thermal");

  conf->getsubdompde(subdoms_,"Thermal");
  ReadBCs("Thermal"); 
}

Therm2dPDE::~Therm2dPDE()
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::~Therm2dPDE" << std::endl;
#endif

  ;
}


void Therm2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::SetMatrixFactors" << std::endl;
#endif

  switch(formulation_) {
  case 0:
    matrix_factor_[0] = 1.0*a0_;
    matrix_factor_[1] = 0.0;
    matrix_factor_[2] = 0.0;
    matrix_factor_[3] = 1.0;
    break;
  case 1:
     matrix_factor_[0] = 1.0;
     matrix_factor_[1] = 0.0;
     matrix_factor_[2] = 0.0;
     matrix_factor_[3] = 1.0*factor1_;
     break;
  }
}

void Therm2dPDE::CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::CalcParameters" << std::endl;
#endif

 a0_ = gamma_parab_*dt;
 factor1_ = 1.0/(gamma_parab_*dt);
 factor2_ = (1-gamma_parab_)/gamma_parab_;
 factor_pred_ = (1-gamma_parab_)*dt;
 
}

void Therm2dPDE::Predictor()
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::Predictor" << std::endl;
#endif

  sol_pred_ = sol_old_+ sol_der1_old_*factor_pred_;

}

void Therm2dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter,  Integer &maxnumit, Integer &numeqcoarse, Double &coarsealpha)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::SpecifySolver" << std::endl;
#endif

  conf->get("eps",eps,"Thermal"); // relative accuracy in the precond. energy
  conf->get("dampiter",dampiter,"Thermal"); // damping parameter for Jacobi, SSOR
  conf->get("maxnumit",maxnumit,"Thermal"); // max number of iterations
  conf->get("solvertype",solvertype,"Thermal"); // Richardson or CG
  conf->get("precondtype", precondtype, "Thermal"); //ID or MG
  conf->get("numeqcoarse",numeqcoarse,"Thermal"); // number of equation for coarsing
  conf->get("coarsealpha",coarsealpha,"Thermal"); // coarsing parameter for AMG
}

void Therm2dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints) 
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::SpecifyMatrices" << std::endl;
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
  numdirichlets  = 1;
  numconstraints = 0;

}

void Therm2dPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;
  Point2D * ptCoord; 

  BaseElem * ptElem;

  Integer matrix_stiff=2;
  Integer matrix_mass=5;

  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Integer i,j; 
  for (i=0; i<subdoms_.size(); i++)
  {
    ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

    for (j=0; j<elemssd.size(); j++)
      {
	ptElem=elemssd[j]->ptElem;
	
	BaseForm<Point2D> * bilinear_stiff = new LaplaceInt<Point2D>(ptElem,1);
	BaseForm<Point2D> * bilinear_mass  = new MassInt<Point2D>(ptElem,1);

	if (!bilinear_stiff) Error("Problems with allocation of object Laplace");
	if (!bilinear_mass)  Error("Problems with allocation of object Mass");

	 Integer ii;
	  Integer elsize=(elemssd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(elemssd[j]->connect)[ii];
  
	  ptCoord=new Point2D[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level); 

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
 
	  
	  std::cout << " stiffness: " << elemmat << std::endl;
#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;

	  (*debug) << elemmat << std::endl;
	  (*debug) << "Connection array " << std::endl;
	  (*debug)  << connecth  << std::endl;

#endif
  
      ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_stiff);

      // mass part
      bilinear_mass->CalcElemMatrix(ptCoord, elemmat);
 
      std::cout << " mass: " << elemmat << std::endl;

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
  (*trace) << "Leaving Therm2dPDE::SetupMatrices" << std::endl;
#endif
}

void Therm2dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val, valueTF;

  //system matrix: id = 1
  Integer matrix_id = 1;

  val=ptTimeFunc_->TimeFuncAtTime(atime,level);
  //  val=20;

  // set homogeneous boundary conditions
  Integer i,j=0;
  std::list<Integer> nodes_hd;
  Integer sizebc=bcs_hd_.size();
  for (i=0; i< bcs_hd_.size(); i++)
    {
      val = 0;
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

  // set dirichlet boundary conditions
  val=1.0;
  j=0;
  std::list<Integer> nodes_id;
  sizebc=bcs_id_.size();
  for (i=0; i< bcs_id_.size(); i++)
    {
      nodes_id=ptBCs->GetNodesLevel(bcs_id_[i]);
   
      for (std::list<Integer>::const_iterator p=nodes_id.begin(); p!=nodes_id.end(); p++, j++)
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

void Therm2dPDE::ComputeRHS(const Double atime, BCs * ptBCs=NULL)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::ComputeRHS" << std::endl;
#endif

  Integer matrix_id;
  Vector<Double> help;

  switch(formulation_) {
  case 0:
    matrix_id = 2;
    help=sol_pred_*(-1);
    ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,help.get());
    break;

  case 1:
    matrix_id = 5;
    help=sol_pred_*factor1_;
    ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,help.get());
    break;
  }
}


void Therm2dPDE::SolveStepStatic(BCs * ptBCs, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::SolveStepStatic" << std::endl;
#endif

  Integer type = 0;
  Integer update = 0;
  Integer job = 1;

  SetupMatrices(type);
  SetBCs(ptBCs,level,update,0);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);

}


void Therm2dPDE::SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::SolveStepTrans" << std::endl;
#endif

  lasttimecalc_= asteptime;
  laststepcalc_= kstep;

  Double * ptsol;

  Integer update,job;

  if (kstep==0) 
    {
      update = 0;
      job = 1;
      SetupMatrices(level);
      ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);
    }
  else if (reset) {   
      update = 1;
      job    = 1;

      Predictor();

      ptalgsys_->ResetRHS(as_sysid_);
      ptalgsys_->ResetMatrix(0,0,1);
      ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);
      ComputeRHS(lasttimecalc_,ptBCs);
    }
  else {
    update = 1;
    job = 3;

    Predictor();
    ptalgsys_->ResetRHS(as_sysid_);
    ComputeRHS(lasttimecalc_,ptBCs);
    }

  SetBCs(ptBCs,level,update,lasttimecalc_);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
  ptsol = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Integer i;
  switch(formulation_) {

  case 0:      // effective mass formulation
    for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
      sol_der1_[i]=ptsol[i];
    break;

  case 1:      // effective stiffness formulation
    for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
      sol_[i]=ptsol[i];
    break;
  }

  CalcDerSol();
}

void Therm2dPDE::CalcDerSol()
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::CalcDerSol" << std::endl;
#endif

  switch(formulation_) {
  case 0:        // effective mass formulation
    sol_ = sol_pred_ + sol_der1_*a0_;
    break;

  case 1:        // effective stiffness formulation
    sol_der1_ = (sol_-sol_pred_)*factor1_;
    break;
  }
}

void Therm2dPDE::SaveSolAsPrevStep()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SaveSolAsPrevStep" << std::endl;
#endif

  sol_old_=sol_;
  sol_der1_old_=sol_der1_;
}

void Therm2dPDE:: WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::WriteResultsInFile" << std::endl;
#endif

  OutFile_->WriteSolution(sol_, laststepcalc_, lasttimecalc_," temperature"); 
}

} // end of namespace
