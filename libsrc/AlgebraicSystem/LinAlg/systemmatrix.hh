#ifndef FILE_SYSTEMMATRIX_2001
#define FILE_SYSTEMMATRIX_2001

#include "matrix_header.hh"

namespace CoupledField
{

  //! Class for storing matrix and right hand side
  /*! This class is created only for storing matrix and right hand side, which are arised from finite element discretization of solving problem. It is base class for classes LinSystem and Assemble. It is done in order to store these data once in program and provides explicit access to them from methods of classes LinSystem and Assemble. Owing to implementation of this class as template, we can define type of system matrix (ex. sparse, symmetric and so on), when we determine object of the class.
   */

template <class T_Matrix>
class SystemMatrix
{
protected:
  //! system matrix
  T_Matrix A;

  //! RHS
  Vector<Double> b;

  //! size of system
  Integer size;

  //! solution
  Vector<Double> x;
   
  //! deconstructor
  virtual ~SystemMatrix(){};

public:
  //! print on screen: only for check
  void print(){ std::cout  << A; std::cout << b; std::cout << x;}
  //! print solution on screen
  void printxscreen() { std::cout << " solution " << x << std::endl;}
  //! print in file matrix A and RHS
  void printAb(std::ostream * out, const std::string & title) const;
  //! print solution 
  void printx(std::ostream * out, const Double time) const;
  //! Return size of matrix A
  Integer getSize(){return A.getSize();}
};

#ifdef __GNUC__
  template class SystemMatrix< Matrix<Double> >;
  template class SystemMatrix< SymMatrix<Double> >;
  template class SystemMatrix< SparseMatrix<Double> >;
//  template class SystemMatrix< SymBandMatrix<Double> >; 
#endif


}
#endif
