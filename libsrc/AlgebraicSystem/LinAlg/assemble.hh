#ifndef FILE_ASSEMBLE_2001
#define FILE_ASSEMBLE_2001

#include "clock.hh"
#include "systemmatrix.hh"

#include "elements_header.hh"
#include "matrix_header.hh"

#include "laplaceint.hh"
#include "massint.hh"
#include "linearform.hh"

#include "domain.hh"
 
namespace CoupledField
{

  //! Class for assembling global system matrix
  /*! This class provides methods for assembling global system matrix according to the weak form of differential equations.
  */

template< Integer Dim, class T_Matrix> 
class Assemble: virtual public SystemMatrix<T_Matrix>  
{
public:

  //! stiffness matrix
  T_Matrix S;

  //! mass matrix
  T_Matrix M;

  //! indicator for mass matrix
  Boolean IsCalcM;

  //! indicator for stiffnes matrix
  Boolean IsCalcS;

  //! constructor 
  Assemble(Grid * aptgrid, const Integer level);

  //! deconstructor
  ~Assemble();
    
  //! Restore solution after applying penalty method for zero boundary condition
  void Restore();

  //! Assemble system matrix from element matrix and boundary condition
  void AssembleSysMatrix(const Double CoefLaplace, const Double CoefMass);

  //! Assemble right hand side
  void SetRHS(const Vector<Double> & CoefM=Vector<Double>(), const Vector<Double> & CoefS=Vector<Double>(), const Vector<Double> & f=Vector<Double>());

  //! Assemble global matrix from element matrix
  template<class typeBaseForm>
  void AssembleGlobal(T_Matrix & Mat) const;

  //! Set Dirichlet condition to right hand side, penalty method
  void SetDirichletBoundaryCondRHS_PenaltyMethod(const Double val_tf); 

  //! Set Dirichlet condition to system matrix, penalty method
  void SetDirichletBoundaryCondSysMat_PenaltyMethod();

  //! Set zero Dirichlet boundary condition, cut row and colomn 
  void SetDirichletBoundaryCondZero_Cut();

  //! Define array of nodes at which we have Dirichlet boundary condition
  void SetNodesBoundaryCondition(FileType * aptFileType);

  //! Set matrix A: only for test LinSystem
  void SetAb();

  //! Print, only for testing
  void Print();

private:

  //! 
  Grid * ptgrid;

  //!
  Integer level;

  //!
  Integer size;

  //!
  Double BigConst;

  //!
  Vector<Integer> nodesDirBC; 
};

#ifdef __GNUC__
template class Assemble<2, Matrix<Double> >;
template class Assemble<2, SymMatrix<Double> >;
template class Assemble<3, Matrix<Double> >;
template class Assemble<3, SymMatrix<Double> >;
template class Assemble<2, SparseMatrix<Double> >;
template class Assemble<3, SparseMatrix<Double> >;
#endif

}

#endif  
  
