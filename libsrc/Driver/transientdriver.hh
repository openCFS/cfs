#ifndef FILE_TRANSIENTDRIVER_2001
#define FILE_TRANSIENTDRIVER_2001

#include "basedriver.hh"

namespace CoupledField
{

//! driver for transient problems.it is derived from BaseDriver;
  class TransientDriver : virtual public BaseDriver
{
public:
  //!  constructor
  /*!
    \param adomain pointer to class Domain
  */
  TransientDriver(Domain * adomain);

   //! deconstructor 
  virtual ~TransientDriver();

  //!  main method, where time-stepping is implemented. it is for transient and static problem
  virtual void SolveProblem();

  //! method with adaptivity in time
  virtual void SolveProblemAdapt();

  //! method with adaptivity in space 
  virtual void SolveProblemAdaptSpace();

  //!  to setup matrices of PDE. we call according method of class PDE for setup matrices of PDE in assembling procedure.
  /*!
    \param pdenumber (input) number of PDE
    \param matrixtype (input) type of system matrix in algebraic system
  */
  void SetupMatricesPDE(const Integer pdenumber, const Integer matrixtype);

protected:

private:
  //!
  Integer numstep_,isavebegin_,isaveincr_,isaveend_;

  //!
  Double firstdt_;

};

}

#endif // FILE_TRANSIENTDRIVER
