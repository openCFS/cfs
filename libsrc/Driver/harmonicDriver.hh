#ifndef FILE_HARMONICDRIVER_2001
#define FILE_HARMONICDRIVER_2001

#include "basedriver.hh"

namespace CoupledField
{

//! driver for static problems. it is derived from BaseDriver
class HarmonicDriver : virtual public BaseDriver
{
public:
  //! constructor
  /*!
    \param adomain pointer to class Domain
  */
  HarmonicDriver(Domain * adomain);

   //! deconstructor 
  virtual ~HarmonicDriver();
  
  //!  main method, where time-stepping is implemented. it is for transient and static problem
  virtual void SolveProblem();

  //!  to setup matrices of PDE. we call according method of class PDE for setup matrices of PDE in assembling procedure.
  /*!
    \param pdenumber number of PDE
    \param matrixtype type of matrix
  */
  virtual void SetupMatricesPDE(const Integer pdenumber, const Integer matrixtype);
};

}

#endif // FILE_HARMONICDRIVER
