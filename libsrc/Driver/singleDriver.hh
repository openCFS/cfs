#ifndef FILE_SINGLEDRIVER_2004
#define FILE_SINGLEDRIVER_2004

#include "basedriver.hh"

namespace CoupledField
{

//! Abstract base class for sinlge driver (static, transient, harmonic)
  class SingleDriver : public BaseDriver
{
public:
  
  //! constructor
  /*!
    \param adomain pointer to class Domain
  */
  SingleDriver(Domain * adomain) :
    BaseDriver(adomain)
    {};

   //! deconstructor
  virtual ~SingleDriver(){};
  
  //! main method, where time-stepping is implemented. it is for transient and static problem
  virtual void SolveProblem()=0;

  //! to setup matrices of PDE. we call according method of class PDE for setup matrices of PDE in assembling procedure.
  /*!
    \param pdenumber number of PDE
    \param matrixtype type of matrix
  */
  virtual void SetupMatricesPDE(Integer pdenumber, const Integer matrixtype)
  { Error("SetupMatricesPDE not implemented in base class!",__FILE__,__LINE__); };
  
private:
   
  };

}

#endif 
