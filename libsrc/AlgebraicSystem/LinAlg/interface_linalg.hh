#ifndef FILE_INTERFACE_LINALG_2002
#define FILE_INTERFACE_LINALG_2002

#include "work.hh"
#include "abstractAS.hh"
 
namespace CoupledField
{

template<class Dim>
class InterfaceAlgSys: public AbstractAlgSys<Dim>
{
public:
  /// Constructor with parameter - pointer to Grid and tolerance
  InterfaceAlgSys(Grid<Dim> * aptgrid, const Integer level, const Double aeps); 

 //! Restore solution after applying penalty method for zero boundary condition
 virtual void Restore(){ ptWork->Restore() ;} 

  //! Assemble system matrix from element matrix and boundary condition
virtual void AssembleSysMatrix(const Double CoefLaplace, const Double CoefMass)
 { ptWork-> AssembleSysMatrix(CoefLaplace, CoefMass);}

   //! Assemble right hand side
  virtual void SetRHS(const Vector<Double> & CoefM=Vector<Double>(), const Vector<Double> & CoefS=Vector<Double>(), const Vector<Double> & f=Vector<Double>()){ ptWork->SetRHS(CoefM, CoefS, f);}

  //! Set Dirichlet condition to right hand side, penalty method
  virtual void SetDirichletBoundaryCondRHS_PenaltyMethod(const Double val_tf)
  {ptWork->SetDirichletBoundaryCondRHS_PenaltyMethod(val_tf);}

   //! Set Dirichlet condition to system matrix, penalty method
 virtual void SetDirichletBoundaryCondSysMat_PenaltyMethod()
 {ptWork->SetDirichletBoundaryCondSysMat_PenaltyMethod();}

    //! Set zero Dirichlet boundary condition, cut row and colomn
  virtual void SetDirichletBoundaryCondZero_Cut()
  {ptWork->SetDirichletBoundaryCondZero_Cut();}

  //! Define array of nodes at which we have Dirichlet boundary condition
  virtual void SetNodesBoundaryCondition(FileType * aptFileType)
  {ptWork->SetNodesBoundaryCondition(aptFileType);}

    //! Conjugate Gradient for positive symmetric matrix
  virtual Boolean CG(const Integer maxIter, enum precond typePrecond)
 {return ptWork->CG(maxIter,typePrecond);}

  //! Generalized minimum residual method
  virtual Boolean GMRes_m(const Integer maxIter, enum precond typePrecond, const Integer m){ return ptWork->GMRes_m(maxIter,typePrecond,m); }

  //! BiCGSTAB method, can be used for all types of matrix
//  virtual Boolean BiCGSTAB(const Integer maxIter, enum precond typePrecond)=0;

  //! Return size of system matrix
  Integer getSize(){ return ptWork->getSize();}
 
  // Print system matrix and right hand side
  void printAb(std::ostream * out, const std::string & title) const 
  { ptWork->printAb(out,title);}
 
  //Print solution
  void printx(std::ostream * out, const Double time) const 
  { ptWork->printx(out,time);}
 
  // Return solution
   Vector<Double> getSolution()
  { return ptWork->getSolution(); }

private:
  WorkWithSysMat<Dim, Matrix<Double> > * ptWork;
} ;

template<class Dim>
inline InterfaceAlgSys<Dim>::InterfaceAlgSys(Grid<Dim> * aptgrid, const Integer level, const Double aeps)
: AbstractAlgSys<Dim>(aptgrid,level,aeps)
{
#ifdef TRACE
 (*trace) << "Entering InterfaceAlgSys::InterfaceAlgSys" << std::endl;
#endif
   ptWork=new WorkWithSysMat<Dim, Matrix<Double> >(ptGrid,level,eps);
}

} // end of namespace

#endif
