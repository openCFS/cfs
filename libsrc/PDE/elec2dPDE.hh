#ifndef FILE_ELEC2dPDE_2001
#define FILE_ELEC2dPDE_2001

#include "timefunc.hh"
#include "basepde.hh"
#include "domain.hh"
#include "material.hh"
#include "laplaceint.hh"

namespace CoupledField
{

  //! Class for acoustic equation
  /*! 
    This class is derived from class PDE. It is used for solving electrostatic equation 
  */

class Elec2dPDE: virtual public BasePDE
{
public:

  //!
  Elec2dPDE(Grid<Point2D> * , Material * , FileType * ptFileType);

  //!
  ~Elec2dPDE();

  //!
  void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps, Double &adampiter, 
                     Integer &amaxnumit);
  //!
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode,  
                       Integer &numdirichlets, Integer &numconstraints);

  //!
  virtual void CalcElemStiff();

  //!
  virtual void AssembleGlobal(AbstractAlgebraicSys * algsys);

  //!
  Vector<Double> & getS() { return sol;}

  //!
  void SolveStatic();

private:

  //!
  Grid<Point2D> * ptgrid;

  //!
  BaseForm<Point2D> * bilinear;

  //!
  TimeFunc * ptTimeFunc;

  //!
  Vector<Double> sol;

  //!
  Integer size;  

};

inline Elec2dPDE::~Elec2dPDE()
{
  // if (ptWork) delete ptWork;
  // if (ptTimeFunc) delete ptTimeFunc;
}

} // end of namespace
#endif
