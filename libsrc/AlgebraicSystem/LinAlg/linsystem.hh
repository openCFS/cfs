#ifndef FILE_LINSYSTEM_2001
#define FILE_LINSYSTEM_2001
 
namespace CoupledField
{

  //! This class provides methods for solving linear algebraic system
  /*! This class provides methods for solving linear algebraic system. This class is derived from SystemMatrix, so it works with matrix and right hand side, which are stored in it.
  */

template< class T, class T_Matrix> 
class LinSystem: virtual public SystemMatrix<T_Matrix> 
{
public:
  //! constructor : tolerance
  LinSystem(const Double);

  //! deconstructor
  ~LinSystem();

  //! set A and b only for check function
  void Set();
  void check()
  { std::cout << " check " << A*x << std::endl << " b " << std::endl << b << std::endl;}

  //! Conjugate Gradient for positive symmetric matrix
  Boolean CG(const Integer maxIter, enum precond typePrecond);

  //! Generalized minimum residual method
  Boolean GMRes_m(const Integer maxIter, enum precond typePrecond, const Integer m);

  //! BiCGSTAB method, can be used for all types of matrix
  Boolean BiCGSTAB(const Integer maxIter, enum precond typePrecond);

  //!
  void printsol(){ std::out << x << std::endl;}
private:

  //! Tolerance 
  Double eps; 
};

#ifdef __GNUC__
template class LinSystem<Double, Matrix<Double> >;
template class LinSystem<Double, SymMatrix<Double> >; 
template class LinSystem<Double, SymBandMatrix<Double> >;
#endif

}

#endif  
  
