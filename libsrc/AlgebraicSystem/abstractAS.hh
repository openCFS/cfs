#ifndef FILE_ABSTRACTAS_2001
#define FILE_ABSTRACTAS_2001
 
namespace CoupledField
{

template<class T_leaftype>
class AbstractAlgSys
{
public:

  T_leaftype & asLeaf()
  { return static_cast<T_leaftype&>(*this); }

 //! Restore solution after applying penalty method for zero boundary condition
 virtual void Restore(){ asLeaf.Restore(); }

  //! Assemble system matrix from element matrix and boundary condition
 virtual void AssembleSysMatrix(const Double CoefLaplace, const Double CoefMass){ asLeaf.AssembleSysMatrix(CoefLaplace, CoefMass); }

   //! Assemble right hand side
  virtual void SetRHS(const Vector<Double> & CoefM=Vector<Double>(), const Vector<Double> & CoefS=Vector<Double>(), const Vector<Double> & f=Vector<Double>())
{ asLeaf.SetRHS(CoefM=Vector<Double>(), CoefS=Vector<Double>(), f=Vector<Double>() ); }

  //! Set Dirichlet condition to right hand side, penalty method
  virtual void SetDirichletBoundaryCondRHS_PenaltyMethod(const Double val_tf)
 { asLeaf.SetDirichletBoundaryCondRHS_PenaltyMethod(val_tf); }

   //! Set Dirichlet condition to system matrix, penalty method
 virtual void SetDirichletBoundaryCondSysMat_PenaltyMethod() { asLeaf.SetDirichletBoundaryCondSysMat_PenaltyMethod(); }

    //! Set zero Dirichlet boundary condition, cut row and colomn
  virtual void SetDirichletBoundaryCondZero_Cut(){ asLeaf.SetDirichletBoundaryCondZero_Cut(); }

  //! Define array of nodes at which we have Dirichlet boundary condition
  virtual void SetNodesBoundaryCondition(FileType * aptFileType) { asLeaf.SetNodesBoundaryCondition(aptFileType); }

    //! Conjugate Gradient for positive symmetric matrix
  virtual Boolean CG(const Integer maxIter, enum precond typePrecond){ asLeaf.CG( maxIter, typePrecond); }

  //! Generalized minimum residual method
  virtual Boolean GMRes_m(const Integer maxIter, enum precond typePrecond, const Integer m){ asLeaf.GMRes_m( maxIter, typePrecond, m) ;}

  //! BiCGSTAB method, can be used for all types of matrix
//  virtual Boolean BiCGSTAB(const Integer maxIter, enum precond typePrecond)=0;

} ;

} // end of namespace

#endif
