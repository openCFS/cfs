#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "interface_linalg.hh"
#include "acoustic2dPDE.hh"
 
namespace CoupledField
{

Acoustic2dPDE::Acoustic2dPDE(Grid<Point2D> * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults<Point2D> * aptOut)
:BasePDE(ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::Acoustic2dPDE " << std::endl;
#endif

  dofspernode_=1;
  ptgrid_=aptgrid;

  Double density, compress;
  if (MatFile_)
  {
    MatFile_->ReadDensityAndCompress(density,compress);
    coeff_=sqrt(compress/density);
  }
  else coeff_=1;

   sol_.Resize(ptgrid_->GetMaxnumnodes(0));
   sol_.Init(0);
   sol_der1_.Resize(ptgrid_->GetMaxnumnodes(0));
   sol_der1_.Init(0);
   sol_der2_.Resize(ptgrid_->GetMaxnumnodes(0));
   sol_der2_.Init(0);
    sol_old_.Resize(ptgrid_->GetMaxnumnodes(0));
   sol_old_.Init(0);
  InFile_->ReadIntegrationParam(alpha_, beta_, gamma_);

}

void Acoustic2dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, Integer &maxnumit)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SpecifySolver" << std::endl;
#endif

  conf->get("eps",eps,"Acoustic"); // relative accuracy in the precond. energy
  conf->get("dampiter",dampiter,"Acoustic"); // damping parameter for Jacobi, SSOR
  conf->get("maxnumit",maxnumit,"Acoustic"); // max number of iterations
  conf->get("solvertype",solvertype,"Acoustic"); // Richardson or CG
  conf->get("precondtype", precondtype, "Acoustic"); //ID or MG
}

void Acoustic2dPDE::SetAlgSys_id(Integer as_sysid)
{
 AS_sysid_ = as_sysid;
}

void Acoustic2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;
  matrix_factor_[1] = 0.0;
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 1.0*coeff_*a0_;
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
  numconstraints = 0;
}

void Acoustic2dPDE::SetupMatrices(AbstractAlgebraicSys * algsys, Integer type)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetupMatrices" << std::endl;
#endif

  Integer k,l;
  Integer i,ii,iii;
  Integer irow,icln;

  Integer numnodeelem=ptgrid_->GetNumNodesPerElem(0,0);
  Integer * help=new Integer[numnodeelem];
  Matrix<Double> elemmat;

  std::cout << "numnodeleme" << numnodeelem << std::endl;
  Double * pilesmat = new Double[numnodeelem*numnodeelem];

  Point2D * ptCoord=new Point2D[numnodeelem];

  BaseElem * ptElem;
  switch(numnodeelem)
  {
    case 3:
       ptElem=new Triangle1(GaussOrder3);
       break;

    case 4:
       ptElem=new Quad1(GaussOrder2);
       break;

    default:
       Error("Number of nodes per element is strange",__FILE__,__LINE__);
  }

   Integer numelem=ptgrid_->GetMaxnumElem(0);

   BaseForm<Point2D> * bilinear_stiff = new LaplaceInt<Point2D>(ptElem,1);
   BaseForm<Point2D> * bilinear_mass  = new MassInt<Point2D>(ptElem,1);

   Integer matrix_stiff=2;
   Integer matrix_mass=5;

   for (i=0; i<numelem; i++)
     {
       ptgrid_->GetConnection(help,0,i,numnodeelem);
       ptgrid_->GetCoordOfNodesElem(i,0,numnodeelem,ptCoord);

       // stiffness part
       bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);

       std::cout << "heh" << std::endl;
       ii = 0;
       for (k=0;k<numnodeelem;k++)
         {
           for (l=0;l<numnodeelem;l++)
             {
               pilesmat[ii]=elemmat[k][l];
               ii = ii+1;
             }
         }
       algsys->PutElemMatAlgSys(pilesmat, help, numnodeelem, AS_sysid_, AS_sysid_, matrix_stiff);

    std::cout << "stiff" << std::endl;
           // mass part
           bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

           ii = 0;
           for (k=0;k<numnodeelem;k++)
             {
               for (l=0;l<numnodeelem;l++)
                 {
                   pilesmat[ii]=elemmat[k][l];
                   ii = ii+1;
                 }
             }
           algsys->PutElemMatAlgSys(pilesmat, help, numnodeelem, AS_sysid_, AS_sysid_,matrix_mass);
     }

    std::cout << " mass " <<std::endl;
   delete [] ptCoord;
   delete [] help;
   delete [] pilesmat;
}

void Acoustic2dPDE::SetBCs(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, const Integer level, const Integer update, const Double atime)
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
    for (list<NodeRestraint>::const_iterator p=restr.begin(); p!=restr.end(); p++, i++)
   {
          Integer node=p->nodalnum;
          if (update==1)
            {
              ptalgsys->SetBCDirichletUpdate(i+1, node, val, dofspernode_, AS_sysid_, AS_sysid_, matrix_id);
            }
          else
            {
              ptalgsys->SetBCDirichlet(i+1, node, val, dofspernode_, AS_sysid_,
AS_sysid_, matrix_id);
            }
    }
}

void Acoustic2dPDE::ComputeRHS(AbstractAlgebraicSys *ptalgsys)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::ComputeRHS" << std::endl;
#endif
  Integer n;
  Integer matrix_id;
  Vector<Double> coeffMass;

      matrix_id = 5;
      coeffMass=sol_*a0_+sol_der1_*a2_+sol_der2_*a3_;

      ptalgsys->UpdateRHS(AS_sysid_,AS_sysid_,matrix_id,coeffMass.get());
}

void Acoustic2dPDE::SolveStepStatic(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepStatic" << std::endl;
#endif

  Integer type = 0;
  Integer update = 0;
  Integer job = 1;

  SetupMatrices(ptalgsys, type);
  SetBCs(ptalgsys, ptBCs,level,update,0);
  ptalgsys->ComputePrecond(job,AS_sysid_);
  ptalgsys->SolveAlgSys(AS_sysid_);
}

void Acoustic2dPDE::SolveStepTrans(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean reset)
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
  SetupMatrices(ptalgsys, type);
  ptalgsys->ComputeSysMatrix(AS_sysid_,AS_sysid_,matrix_factor_);
     }
  else if (reset)
       {
        update = 1;
        job    = 0;

        ptalgsys->ResetRHS(AS_sysid_);
        ptalgsys->ResetMatrix(0,0,1);
        ptalgsys->ComputeSysMatrix(AS_sysid_,AS_sysid_,matrix_factor_);
        ComputeRHS(ptalgsys);
       }
       else
        {
      update = 1;
      job    = 3;
      ptalgsys->ResetRHS(AS_sysid_);
      ComputeRHS(ptalgsys);
        };

  SetBCs(ptalgsys, ptBCs,level,update,lasttimecalc_);
  ptalgsys->ComputePrecond(job,AS_sysid_);
  ptalgsys->SolveAlgSys(AS_sysid_);
  ptsol = ptalgsys->GetSolution(AS_sysid_);

    // save solution
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);
  sol_=transsol;
    
   // calculation of derivatives of solution
   CalculationDerivativesSol();
}

void Acoustic2dPDE::SolveStepTransNewAssemble(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, const Integer level)
{
#ifdef TRACE
 (*trace) << "entering Acoustic2dPDE::SolveStepTransNewAssemble" << std::endl;
#endif

  Integer type = 0;
  Double * ptsol;

  Integer update = 0;
  Integer job = 1;
  SetupMatrices(ptalgsys, type);
  ptalgsys->ComputeSysMatrix(AS_sysid_,AS_sysid_,matrix_factor_);
  SetBCs(ptalgsys, ptBCs,level,update,lasttimecalc_);
  ptalgsys->ComputePrecond(job,AS_sysid_);
  ptalgsys->SolveAlgSys(AS_sysid_);
  ptsol = ptalgsys->GetSolution(AS_sysid_);

   // save solution
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);

  sol_=transsol;
}

void Acoustic2dPDE::SolveStepTransContinue(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepTransContinue" << std::endl;
#endif

  Integer type = 0;
  Double * ptsol;

      Integer update = 1;
      Integer job    = 3;

      ptalgsys->ResetRHS(AS_sysid_);
      ComputeRHS(ptalgsys);
      SetBCs(ptalgsys, ptBCs,level,update,lasttimecalc_);
      ptalgsys->ComputePrecond(job,AS_sysid_);
      ptalgsys->SolveAlgSys(AS_sysid_);
      ptsol = ptalgsys->GetSolution(AS_sysid_);

   // save solution
   Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);

   sol_=transsol;
}

void Acoustic2dPDE::SolveStepTransReset(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepTransReset" << std::endl;
#endif

  Integer type = 0;
  Double * ptsol;

  Integer update = 1;
  Integer job    = 0;

  ptalgsys->ResetRHS(AS_sysid_);
  ptalgsys->ResetMatrix(0,0,1);
  ptalgsys->ComputeSysMatrix(AS_sysid_,AS_sysid_,matrix_factor_);
  ComputeRHS(ptalgsys);
  SetBCs(ptalgsys, ptBCs,level,update,lasttimecalc_);
  ptalgsys->ComputePrecond(job,AS_sysid_);
  ptalgsys->SolveAlgSys(AS_sysid_);
  ptsol = ptalgsys->GetSolution(AS_sysid_);

   // save solution
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);

  sol_=transsol;
}

void Acoustic2dPDE:: WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::WriteResultsInFile" << std::endl;
#endif

  OutFile_->WriteSolution(sol_, laststepcalc_, lasttimecalc_); /// !!!!!!!!!!
  OutFile_->WriteFirstDerSolution(sol_der1_, laststepcalc_,lasttimecalc_);
  OutFile_->WriteSecondDerSolution(sol_der2_,laststepcalc_,lasttimecalc_);

}

void Acoustic2dPDE :: CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::CalcParamForNewmarkMethod" << std::endl;
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
  Vector<Double> sol_der2_old=sol_der2_;

  sol_der2_=(sol_-sol_old_)*a0_-sol_der1_*a2_-sol_der2_old*a3_;
  sol_der1_+=sol_der2_old*a6_+sol_der2_*a7_;

  sol_old_=sol_;
}

} // end of namespace
