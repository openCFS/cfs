#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "interface_linalg.hh"
#include "therm2dPDE.hh"

 
namespace CoupledField
{

Therm2dPDE::Therm2dPDE(AbstractAlgebraicSys * ptalgsys, Grid<Point2D> * agrid, Material * aMatFile, TimeFunc * aptTimeFunc, FileType * aInFile, WriteResults<Point2D> * aptOut)
:BasePDE(ptalgsys, aMatFile, aInFile, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::Therm2dPDE " << std::endl;
#endif

  dofspernode_ = 1;
  doftype_     = 14;
  ptgrid_=agrid;

  conf->get("formulation", effsysmat_, "Thermal");

  size_=ptgrid_->GetMaxnumnodes(0);
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

  InFile_->ReadParabolicParam(gamma_parab_);
  std::cout << "gamma_parab_" << std::endl;
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

if (effsysmat_==1)
{
 matrix_factor_[0] = 1.0*a0_;
  matrix_factor_[1] = 0.0;
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 1.0;
}
else 
{
  matrix_factor_[0] = 1.0;
  matrix_factor_[1] = 0.0;
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 1.0*factor1_;
}

}

void Therm2dPDE::CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::CalcParameters" << std::endl;
#endif

 a0_=gamma_parab_*dt;
 factor1_ = 1.0/(gamma_parab_*dt);
 factor2_ = (1-gamma_parab_)/gamma_parab_;
 factor_pred_ = (1-gamma_parab_)*dt;
 
}

void Therm2dPDE::Predictor()
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::Predictor" << std::endl;
#endif

  sol_pred_=sol_+ sol_der1_*factor_pred_;

}

void Therm2dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter,  Integer &maxnumit, Integer &numeqcoarse)
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

//  if (statickey==0)
//    {
//      matrixsystype[1] = 2;   
//      matrixsystype[4] = 5;   
//    }

  graphtype = 1; // NOGRAPH = 0
                 // NODE   = 1
                 // EDGE   = 2
                 // FACE   = 3 
                 // VOLUME = 4 

  numdofpernode  = 1;
  numdirichlets  = 1;
  numconstraints = 0;

}

void Therm2dPDE::SetupMatrices(const Integer type)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::SetupMatrices" << std::endl;
#endif

  Integer numnodeelem;
  numnodeelem=ptgrid_->GetNumNodesPerElem(0,0);
  Integer * help=new Integer[numnodeelem];
  Matrix<Double> elemmat;

  Point2D * ptCoord=new Point2D[numnodeelem];

  Integer numelem=ptgrid_->GetMaxnumElem(0);
 
  BaseElem * ptElem;

  BaseElem ** ptArrayElem=ptgrid_->getptArrayElem();

  Integer matrix_stiff=2;
  Integer matrix_mass=5;

  Integer i; 
  for (i=0; i<numelem; i++)
  {
    ptElem=ptArrayElem[i];

    ptElem->test();

   BaseForm<Point2D> * bilinear_stiff = new LaplaceInt<Point2D>(ptElem,1);
   BaseForm<Point2D> * bilinear_mass  = new MassInt<Point2D>(ptElem,1);

   if (!bilinear_stiff) Error("Problems with allocation of object Laplace");
   if (!bilinear_mass)  Error("Problems with allocation of object Mass");

   ptgrid_->GetConnection(help,0,i,numnodeelem);
   ptgrid_->GetCoordOfNodesElem(i,0,numnodeelem,ptCoord);

       // stiffness part
      std::cout << "before" << std::endl;
      bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);

      std::cout << "intermediate" << std::endl;
      ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), help, numnodeelem, AS_sysid_, AS_sysid_, matrix_stiff);

       std::cout << "intermediate1" << std::endl;
      // mass part
      bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

        std::cout << "intermediate2" << std::endl;
      ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), help, numnodeelem, AS_sysid_, AS_sysid_,matrix_mass);

       std::cout << "end" << std::endl;

      delete bilinear_stiff;
      delete bilinear_mass;
     }

   delete [] ptCoord;
   delete [] help;

#ifdef TRACE
  (*trace) << "Leaving Therm2dPDE::SetupMatrices" << std::endl;
#endif
}

void Therm2dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::SetBCs" << std::endl;
#endif
 
  Integer num, node, type, tfunc;
  Double val, valueTF;

  //system matrix: id = 1
  Integer matrix_id = 1;

  std::list<NodeRestraint> restr;
  ptBCs->GetRestraints(restr,level);

  val=ptTimeFunc_->TimeFuncAtTime(atime,level);

  Integer i=0;
  for (std::list<NodeRestraint>::const_iterator p=restr.begin(); p!=restr.end(); p++, i++)
  {
        Integer node=p->nodalnum;
        if (p->dof==doftype_)
	  {
        if (update==1)
           {
        ptalgsys_->SetBCDirichletUpdate(i+1, node, val, dofspernode_, AS_sysid_, AS_sysid_, matrix_id);
           }
          else
           {
              ptalgsys_->SetBCDirichlet(i+1, node, val, dofspernode_, AS_sysid_,AS_sysid_, matrix_id);
           }
          }
   }  
}

void Therm2dPDE::ComputeRHS()
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::ComputeRHS" << std::endl;
#endif
  Integer n;
  Integer matrix_id;

  if (effsysmat_ == 1)
    {
      matrix_id = 2;
      Vector<Double> help=sol_pred_*(-1); 
      ptalgsys_->UpdateRHS(AS_sysid_,AS_sysid_,matrix_id,help.get());
      
    }
  else
    {
      matrix_id = 5;
      Vector<Double> help=sol_pred_*factor1_;
      ptalgsys_->UpdateRHS(AS_sysid_,AS_sysid_,matrix_id,help.get());
    }
  std::cout << " end " << std::endl;
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
  ptalgsys_->ComputePrecond(job,AS_sysid_);
  ptalgsys_->SolveAlgSys(AS_sysid_);

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
      Integer type=0;
      update = 0;
      job = 1;
      SetupMatrices(type);
      ptalgsys_->ComputeSysMatrix(AS_sysid_,AS_sysid_,matrix_factor_);
    }
  else if (reset)
    {
      update = 1;
      job    = 0;

   //   if (effsysmat_ == 2)
   	Predictor();

      ptalgsys_->ResetRHS(AS_sysid_);
      ptalgsys_->ResetMatrix(0,0,1);
      ptalgsys_->ComputeSysMatrix(AS_sysid_,AS_sysid_,matrix_factor_);
      ComputeRHS();
    }
      else
    { update = 1;
      job = 3;

    //  if (effsysmat_ == 2)
        Predictor();

      ptalgsys_->ResetRHS(AS_sysid_);
      ComputeRHS();
    }

  SetBCs(ptBCs,level,update,lasttimecalc_);
  ptalgsys_->ComputePrecond(job,AS_sysid_);
  ptalgsys_->SolveAlgSys(AS_sysid_);
  ptsol = ptalgsys_->GetSolution(AS_sysid_);

 Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);

  //save solution
  if (effsysmat_==1)
    {
      sol_der1_old_=sol_der1_;
      sol_der1_=transsol;
      sol_old_=sol_;

      CalculationDerivativesSol();
    }
  else
    {
     sol_old_=sol_;
     sol_=transsol;
     sol_der1_old_=sol_der1_;
     CalculationDerivativesSol();
    }
}

void Therm2dPDE::CalculationDerivativesSol()
{
#ifdef TRACE
  (*trace) << "entering Therm2dPDE::CalculationDerivativesSol" << std::endl;
#endif

  if (effsysmat_==1)
   {
    sol_+=sol_der1_*a0_+sol_der1_old_*factor_pred_;
   }
  else
   {
    sol_der1_=(sol_-sol_pred_)*factor1_;
   }
}

void Therm2dPDE:: WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::WriteResultsInFile" << std::endl;
#endif

  OutFile_->WriteSolution(sol_, laststepcalc_, lasttimecalc_," temperature"); 
}

} // end of namespace
