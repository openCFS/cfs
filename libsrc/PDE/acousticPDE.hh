#ifndef FILE_ACOUSTICPDE_2001
#define FILE_ACOUSTICPDE_2001

#include "timefunc.hh"
#include "pde.hh"
#include "domain.hh"
#include "abstractAS.hh"
 
namespace CoupledField
{

  //! Class for acoustic equation
  /*! 
    This class is derived from class PDE. It is used for solving acoustic equation on one time step. In this class we work with complex WorkWithSystemMatrix. When we create object of this class, we determine what kind of matrix should be system matrix, for example, sparse matrix, band matrix or so on. We set rules for assembling global system matrix according to weak form of PDE, define right hand side and set boundary conditions. Then we cause one of methods of LinSystem for solving linear system. On the last step we calculate first and second derivatives of the solution.
  */

template<class Dim>
class AcousticPDE: virtual public PDE
{
public:

  //!
  AcousticPDE(const Double , const Double , Grid<Dim> * , const Integer level, Material * , FileType * ptFileType);

  //!
  ~AcousticPDE();

  //!
  Vector<Double> & getS() { return sol;}

  //!
  Vector<Double> & getS1() { return sol_der1;}

  //!
  Vector<Double> & getS2() { return sol_der2;}

  //!
  void SolveNewmarkMethodStatic(const Double atime);

private:

  //!
  Double CoefLaplace;

  //!
  Double CoefMass;   

  //!
  TimeFunc * ptTimeFunc;

  //!
   AbstractAlgSys<Dim> * ptWork;
  //!
  Vector<Double> sol, sol_der2, sol_der1;

  //!
  Integer size;  

};

template<class Dim>
inline AcousticPDE<Dim>::~AcousticPDE()
{
 if (ptWork) delete ptWork;
 if (ptTimeFunc) delete ptTimeFunc;
}

#ifdef __GNUC__
template class AcousticPDE<Point2D>;
template class AcousticPDE<Point3D>;
#endif

} // end of namespace
#endif
