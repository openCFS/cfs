#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "basedriver.hh"

namespace CoupledField
{

//! driver for static problems. it is derived from BaseDriver
class StaticDriver : virtual public BaseDriver
{
public:
  //! constructor
  /*!
    \param adomain pointer to class Domain
  */
  StaticDriver(Domain * adomain);

   //! deconstructor 
  virtual ~StaticDriver();
  
  //!  main method, where time-stepping is implemented. it is for transient and static problem
  virtual void SolveProblem();

  //! method with adaptivity in space
  void SolveProblemAdaptSpace();

  //!  to setup matrices of PDE. we call according method of class PDE for setup matrices of PDE in assembling procedure.
  /*!
    \param pdenumber number of PDE
    \param matrixtype type of matrix
  */
  void SetupMatricesPDE(const Integer pdenumber, const Integer matrixtype);

protected:

private:

};

}

#endif // FILE_STATICDRIVER
