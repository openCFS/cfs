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

  dofspernode=1;
  ptgrid=aptgrid;
  doftype=14;

  Double density, compress,c;
  if (ptMaterial)
  {
    ptMaterial->ReadDensityAndCompress(density,compress);
    coeff=sqrt(compress/density);
  }
  else coeff=1;

   sol.Resize(ptgrid->GetMaxnumnodes(0));
   sol.Init(0);
   sol_der1.Resize(ptgrid->GetMaxnumnodes(0));
   sol_der1.Init(0);
   sol_der2.Resize(ptgrid->GetMaxnumnodes(0));
   sol_der2.Init(0);
    sol_old.Resize(ptgrid->GetMaxnumnodes(0));
   sol_old.Init(0);
//  Double dt0=0.1e-7;                          /// !!!!!!!!!!!!
//  CalcParamForNewmarkMethod(dt0);
  InFile->ReadIntegrationParam(alpha, beta, gamma);

  // read a tolerance for the algebraic system from config-file
//  Double epsilon;
//  conf->get("epsilon",epsilon);
//  ptWork=new InterfaceAlgSys<Point2D>(ptgrid,0,epsilon);

//  ptWork->AssembleSysMatrix(CoefLaplace,CoefMass);
//  ptWork->SetRHS();
  
//  ptWork->SetNodesBoundaryCondition(aptFileType);
//  ptWork->SetDirichletBoundaryCondSysMat_PenaltyMethod();

//  size=ptWork->getSize();

//  sol.Resize(size);
//  sol_der1.Resize(size);
//  sol_der2.Resize(size);
}

void Acoustic2dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter,
                     Integer &maxnumit)
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
 AS_sysid = as_sysid;
}

void Acoustic2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetMatrixFactors" << std::endl;
#endif
 
  Double firstdt=0.5e-7; 

  CalcParamForNewmarkMethod(firstdt);

  matrix_factor[0] = 1.0;
  matrix_factor[1] = 0.0;
  matrix_factor[2] = 0.0;
  matrix_factor[3] = 1.0*coeff*a0;
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
  Integer irow,icln,matdim;

  Integer numnodeelem;
  numnodeelem=ptgrid->GetNumNodesPerElem(0,0);
  Integer * help=new Integer[numnodeelem];
  Matrix<Double> elemmat;
  matdim=numnodeelem;

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

   Integer numelem=ptgrid->GetMaxnumElem(0);

   bilinear_stiff = new LaplaceInt<Point2D>(ptElem,1);
   bilinear_mass  = new MassInt<Point2D>(ptElem,1);

   Integer matrix_stiff=2;
   Integer matrix_mass=5;

   for (i=0; i<numelem; i++)
     {
       ptgrid->GetConnection(help,0,i,numnodeelem);
       ptgrid->GetCoordOfNodesElem(i,0,numnodeelem,ptCoord);

       // stiffness part
       bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);

       ii = 0;
       for (k=0;k<matdim;k++)
         {
           for (l=0;l<matdim;l++)
             {
               pilesmat[ii]=elemmat[k][l];
               ii = ii+1;
             }
         }
       algsys->PutElemMatAlgSys(pilesmat, help, numnodeelem, AS_sysid, AS_sysid, matrix_stiff);

    std::cout << "stiff" << std::endl;
           // mass part
           bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

           ii = 0;
           for (k=0;k<matdim;k++)
             {
               for (l=0;l<matdim;l++)
                 {
                   pilesmat[ii]=elemmat[k][l];
                   ii = ii+1;
                 }
             }
           algsys->PutElemMatAlgSys(pilesmat, help, numnodeelem, AS_sysid, AS_sysid,matrix_mass);
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

  Integer i, num, node, type, tfunc;
  Double val, valueTF;


  //system matrix: id = 1
  Integer matrix_id = 1;

  num = ptBCs->GetNumDirichlet(level);

  for (i=0;i<num;i++)
    {
          node = ptBCs->GetDirichletNode(i,level);
          std::cout << "before time"  << node << std::endl;
          val=ptTimeFunc->TimeFuncAtTime(atime,level);

	  std::cout << val << " val" << std::endl;

          if (update==1)
            {
              ptalgsys->SetBCDirichletUpdate(i+1, node, val, dofspernode, AS_sysid, AS_sysid, matrix_id);
            }
          else
            {
              ptalgsys->SetBCDirichlet(i+1, node, val, dofspernode, AS_sysid, AS_sysid, matrix_id);
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
      coeffMass=sol*a0+sol_der1*a2+sol_der2*a3;

      ptalgsys->UpdateRHS(AS_sysid,AS_sysid,matrix_id,coeffMass.get());
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
  ptalgsys->ComputePrecond(job,AS_sysid);
  ptalgsys->SolveAlgSys(AS_sysid);
}

void Acoustic2dPDE::SolveStepTrans(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, Integer kstep, Double asteptime, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepTrans" << std::endl;
#endif

  Integer n,update;
  Integer type = 0;
  Integer job;
  Double * ptsol;

  lasttimecalc= asteptime;
  StepTime=asteptime;
  laststepcalc= kstep;

  std::cout << "kstep" << kstep << std::endl;
  if (kstep==0)
    {
      update = 0;
      job = 1;
      std::cout << StepTime << std::endl;
      SetupMatrices(ptalgsys, type);
      ptalgsys->ComputeSysMatrix(AS_sysid,AS_sysid,matrix_factor);
      SetBCs(ptalgsys, ptBCs,level,update,StepTime);
      ptalgsys->ComputePrecond(job,AS_sysid);
      ptalgsys->SolveAlgSys(AS_sysid);
      ptsol = ptalgsys->GetSolution(AS_sysid);
    }
  else
    {
      update = 1;
      job    = 3;

      ptalgsys->ResetRHS(AS_sysid);
      ComputeRHS(ptalgsys);
      SetBCs(ptalgsys, ptBCs,level,update,StepTime);
      ptalgsys->ComputePrecond(job,AS_sysid);
      ptalgsys->SolveAlgSys(AS_sysid);
      ptsol = ptalgsys->GetSolution(AS_sysid);
    }

   // save solution
   Integer ii;

   Vector<Double> transsol(ptgrid->GetMaxnumnodes(level), ptsol);

   sol=transsol;

   // calculation of derivatives of solution
   CalculationDerivativesSol();
}

void Acoustic2dPDE:: WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::WriteResultsInFile" << std::endl;
#endif

  std::cout << " Attention " << sol << std::endl;
  std::cout << "Step&Time" << laststepcalc << " " << lasttimecalc << std::endl;
  OutFile->WriteSolution(sol, laststepcalc, lasttimecalc); /// !!!!!!!!!!!!!!
  OutFile->WriteFirstDerSolution(sol_der1, laststepcalc,lasttimecalc);
  OutFile->WriteSecondDerSolution(sol_der2,laststepcalc,lasttimecalc);

}

void Acoustic2dPDE :: CalcParamForNewmarkMethod(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::CalcParamForNewmarkMethod" << std::endl;
#endif

 a2=1.0/(beta*dt);
 a0=a2*(1/dt);
 a1=gamma*a2;
 a3=0.5/(beta)-1.0;
 a4=gamma/beta-1.0;
 a5=0.5*dt*(a4-1.0);
 a6=dt*(1-gamma);
 a7=gamma*dt;
}

void Acoustic2dPDE :: CalculationDerivativesSol()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::CalculationDerivativesSol" << std::endl;
#endif
  Vector<Double> sol_der2_old=sol_der2;

  sol_der2=(sol-sol_old)*a0-sol_der1*a2-sol_der2_old*a3;
  sol_der1+=sol_der2_old*a6+sol_der2*a7;

  sol_old=sol;
}
} // end of namespace
