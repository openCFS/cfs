#include <fstream.h>
#include <iostream.h>
#include <string>
#include <math.h>

#include <general_head.hh> 
#include <utils_head.hh>
#include <datainout_head.hh>
#include <domain_head.hh>
#include <linalg_head.hh>
#include "pde.hh"
#include "acousticPDE.hh"
 
namespace CoupledField
{

AcousticPDE::AcousticPDE(const Double epsilon, const Double dt0, Grid<Point2D> * ptgrid, Material * aptMaterial, FileType * aptFileType)
:PDE(aptFileType,aptMaterial)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::AcousticPDE " << endl;
#endif

  ptTimeFunc=new TimeFunc(aptFileType);

  Double density, compress,c;
  if (ptMaterial)
  {
    ptMaterial->ReadDensityAndCompress(density,compress);
    c=sqrt(compress/density);
  }
  else c=1;

  CalcParamForNewmarkMethod(dt0);

  CoefLaplace=1.0;
  CoefMass=a0*c;

  ptWork=new WorkWithSysMat<Point2D, Matrix<Double> >(ptgrid, epsilon);  

  ptWork->AssembleSysMatrix(CoefLaplace,CoefMass);
  ptWork->SetRHS();
  
  ptWork->SetNodesBoundaryCondition(aptFileType);
  ptWork->SetDirichletBoundaryCondSysMat_PenaltyMethod();
  
  size=ptWork->getSize();
 
  sol.Resize(size);
  sol_der1.Resize(size);
  sol_der2.Resize(size);
}

void AcousticPDE::SolveNewmarkMethodStatic(const Double atime)
{
#ifdef TRACE
  (*trace) << "entering PDE::SolveNewmarkMethod" << endl;
#endif
 
   Boolean NeedRestore=FALSE;
   Double valueTF=ptTimeFunc->TimeFuncAtTime(atime,0);
 
   if (valueTF==0)
     { NeedRestore=TRUE; ptWork->SetDirichletBoundaryCondZero_Cut();}
   else
     ptWork->SetDirichletBoundaryCondRHS_PenaltyMethod(valueTF);
 
   ptWork->CG(100, Jacobi);
 
#ifdef DEBUG
   string title=" System matrix after applying boundary condition ";
   ptWork->printAb(debug, title);
#endif
 
  if (NeedRestore) ptWork->Restore();
 
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
