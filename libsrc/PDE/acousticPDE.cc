#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "interface_linalg.hh"
#include "acousticPDE.hh"
 
namespace CoupledField
{

template<class Dim>
AcousticPDE<Dim>::AcousticPDE(const Double dt0, Grid<Dim> * ptgrid, const Integer level,Material * aptMaterial, FileType * aptFileType)
:PDE(aptFileType,aptMaterial)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::AcousticPDE " << std::endl;
#endif

  ptTimeFunc=new TimeFunc(aptFileType);

  Double density, compress,c;
  if (ptMaterial)
  {
    ptMaterial->ReadDensityAndCompress(density,compress);
    std::cout << "density " << density << "compri " << compress << std::endl;
    c=sqrt(compress/density);
  }
  else c=1;

  CalcParamForNewmarkMethod(dt0);

  CoefLaplace=1.0;
  CoefMass=a0*c;

  // read a tolerance for the algebraic system from config-file
  Double epsilon;
  conf->get("epsilon",epsilon);
  ptWork=new InterfaceAlgSys<Dim>(ptgrid,level,epsilon);

  ptWork->AssembleSysMatrix(CoefLaplace,CoefMass);
  ptWork->SetRHS();
  
  ptWork->SetNodesBoundaryCondition(aptFileType);
  ptWork->SetDirichletBoundaryCondSysMat_PenaltyMethod();

  size=ptWork->getSize();

  sol.Resize(size);
  sol_der1.Resize(size);
  sol_der2.Resize(size);
}

template<class Dim>
void AcousticPDE<Dim>::SolveNewmarkMethodStatic(const Double atime)
{
#ifdef TRACE
  (*trace) << "entering PDE::SolveNewmarkMethod" << std::endl;
#endif
 
   Double valueTF=ptTimeFunc->TimeFuncAtTime(atime,0);

   ptWork->SetDirichletBoundaryCondRHS_PenaltyMethod(valueTF);
 
   ptWork->CG(100, Jacobi);

#ifdef DEBUG
   std::string title=" System matrix after applying boundary condition ";
   ptWork->printAb(debug, title);
#endif
 
  /// Calculation of derivatives of solution
  Vector<Double> sol_der2_old=sol_der2;
  sol_der2=(ptWork->getSolution()-sol)*a0-sol_der1*a2-sol_der2_old*a3;
  sol_der1+=sol_der2_old*a6+sol_der2*a7;
 
#ifdef DEBUG
   ptWork->printx(debug,atime);
#endif
 
   Vector<Double> VectM;
   VectM=sol*a0+sol_der1*a2+sol_der2*a3;

   ptWork->SetRHS(VectM);
   sol=ptWork->getSolution();
}

} // end of namespace
