#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acoustic2dPDE.hh"
#include "actimeerror.hh"
#include "outUnverg.hh"
#include "forms_header.hh"

 
namespace CoupledField
{

Acoustic2dPDE::Acoustic2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::Acoustic2dPDE " << std::endl;
#endif

//  doftype_=5;
  doftype_=vp_restraint;
  dofspernode_=1;
  ptgrid_=aptgrid;

  Char * name="result";
  ptOutput=new WriteResultsUnverg(name);
  ptOutput->Init(ptgrid_);
  ptOutput->WriteGrid(0);  

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

  conf->get("alpha_NM",alpha_,"Acoustic"); 
  conf->get("beta_NM",beta_,"Acoustic"); 
  conf->get("gamma_NM",gamma_,"Acoustic");

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

void Acoustic2dPDE::SetupMatrices(const Integer type)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetupMatrices" << std::endl;
#endif

  Integer level=0;

  Matrix<Double> elemmat;
  Point2D * ptCoord; 

  Integer numelems=ptgrid_->GetMaxnumElem(0);

  BaseElem * ptElem;

  Integer matrix_stiff=2;
  Integer matrix_mass=5;

  Double coeffmass,coeffstiff;

/// !!!!!!!!
  CalcCoeff(coeffmass, coeffstiff,0);

 Vector<Integer> connecth;
 Integer i;
 for (i=0; i< numelems; i++)
{
  ptElem=ptgrid_->GetptElem(i);

  BaseForm<Point2D> * bilinear_stiff = new LaplaceInt<Point2D>(ptElem,1);
  BaseForm<Point2D> * bilinear_mass  = new MassInt<Point2D>(ptElem,1);

  ptgrid_->GetConnection(connecth,i,level);
  
  ptCoord=new Point2D[connecth.size()];
  ptgrid_->GetCoordOfNodesElem(i,0,connecth.size(),ptCoord);

  // stiffness part
  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
  elemmat*=coeffstiff;

#ifdef DEBUG
      (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;

      (*debug) << elemmat << std::endl;
#endif     

  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), AS_sysid_, AS_sysid_, matrix_stiff);

      // mass part
  bilinear_mass->CalcElemMatrix(ptCoord, elemmat);
  elemmat*=coeffmass;

#ifdef DEBUG
      (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;

      (*debug) << elemmat << std::endl;
#endif
      
  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), AS_sysid_, AS_sysid_,matrix_mass);

  
  delete bilinear_stiff;
  delete bilinear_mass;
  delete [] ptCoord;
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
              ptalgsys_->SetBCDirichlet(i+1, node, val, dofspernode_, AS_sysid_, AS_sysid_, matrix_id);
            }
	    }
    }
#ifdef TRACE
  (*trace) << "leaving Acoustic2dPDE::SetBCs" << std::endl;
#endif
}

void Acoustic2dPDE::ComputeRHS()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::ComputeRHS" << std::endl;
#endif
  Integer n;
  Integer matrix_id;
  Vector<Double> coeffMass;

      matrix_id = 5;
      coeffMass=sol_*a0_+sol_der1_*a2_+sol_der2_*a3_;

      ptalgsys_->UpdateRHS(AS_sysid_,AS_sysid_,matrix_id,coeffMass.get());
}

void Acoustic2dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepStatic" << std::endl;
#endif

  Integer type = 0;
  Integer update = 0;
  Integer job = 1;

  SetupMatrices(type);
  SetBCs(ptBCs,level,update,0);
  ptalgsys_->ComputePrecond(job,AS_sysid_);
  ptalgsys_->SolveAlgSys(AS_sysid_);
}

void Acoustic2dPDE::SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepTrans" << std::endl;
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
        job    = 1;

        ptalgsys_->ResetRHS(AS_sysid_);
        ptalgsys_->ResetMatrix(0,0,1);
        ptalgsys_->ComputeSysMatrix(AS_sysid_,AS_sysid_,matrix_factor_);
        ComputeRHS();
       }
       else
        {
         update = 1;
         job    = 3;
         ptalgsys_->ResetRHS(AS_sysid_);
         ComputeRHS();
        };

  SetBCs(ptBCs,level,update,lasttimecalc_);
  ptalgsys_->ComputePrecond(job,AS_sysid_);
  ptalgsys_->SolveAlgSys(AS_sysid_);
  ptsol = ptalgsys_->GetSolution(AS_sysid_);

    // save solution
//  sol_.TransformInVector(ptgrid_->GetMaxnumnodes(level),ptsol);
   Integer i;
   for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
    sol_[i]=ptsol[i];

  std::cout << "maxnode:" <<  ptgrid_->GetMaxnumnodes(level) << " level:" << level << std::endl;

   // calculation of derivatives of solution
  CalculationDerivativesSol();

//  ptOutput->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"fluid potential");
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

  ptalgsys_->ResetRHS(AS_sysid_);
  ptalgsys_->ResetMatrix(0,0,1);
  ptalgsys_->ComputeSysMatrix(AS_sysid_,AS_sysid_,mat_factor);

  Integer matrix_id;
  Vector<Double> coeffMass;

  matrix_id = 2;
  Vector<Double> help=-sol_der1_;
  std::cout << help << std::endl;

  ptalgsys_->UpdateRHS(AS_sysid_,AS_sysid_,matrix_id, help.get());

  update=1; // use boundary condition at this time step
//  SetBCs(ptBCs,level,update,lasttimecalc_);
  ptalgsys_->ComputePrecond(job,AS_sysid_);
  ptalgsys_->SolveAlgSys(AS_sysid_);
  ptres = ptalgsys_->GetSolution(AS_sysid_);

 // save solution
  Integer level=0; //   !!!!!!!!!!!!!!

  Vector<Double> transres(ptgrid_->GetMaxnumnodes(level), ptres);
  result=transres;

  (*infofile) << " result " << transres << std::endl;
}

void Acoustic2dPDE:: WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::WriteResultsInFile" << std::endl;
#endif

  if (OutFile_->IsGMV())
{
  OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"fluid_potential");
  OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"fluid_potential, 1st deriv., ");
  OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"fluid_potential, 2nd deriv., ");
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

void Acoustic2dPDE :: CalculationDerivativesSol()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::CalculationDerivativesSol" << std::endl;
#endif

  sol_der2_old_=sol_der2_;

  sol_der2_=(sol_-sol_old_)*a0_-sol_der1_*a2_-sol_der2_old_*a3_;

  sol_der1_old_=sol_der1_;
  sol_der1_+=sol_der2_old_*a6_+sol_der2_*a7_;

  sol_old_=sol_;
}

Double Acoustic2dPDE::CalcEnergyNorm()
{
 Double help1, help2;
 help1=ptalgsys_->CalcEnergyNorm(0,0,2,sol_.get());
 help2=ptalgsys_->CalcEnergyNorm(0,0,5,sol_.get());
 
 return sqrt(help1+help2);
}

void Acoustic2dPDE::CalcCoeff(Double & coeffmass, Double & coeffstiff, const Integer numsubdom)
{
 // get material number for subdomain with number numsubdom from config-file
 Integer matnum;
// conf->getmatnum(matnum,numsubdom);
 conf->get("acoustic2d",matnum);

 // read density and compress with material number matnum
 Double density, compress;
 if (!MatFile_) Error("You didn't specialize material file. Use option -m");
 MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid"); 

 coeffmass  = density*density/compress;
 coeffstiff = density;

 std::cout << " calc coeff is over " << std::endl;
}

TimeErrorEstimator * Acoustic2dPDE::CreatePtTimeError()
{
 return ptTimeError_=new AcousticTimeErrorEstimator(this);
}

Acoustic2dPDE::~Acoustic2dPDE()
{
 if (ptTimeError_) delete ptTimeError_;
}

} // end of namespace
