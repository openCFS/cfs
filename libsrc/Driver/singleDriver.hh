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
  //! \param adomain pointer to class Domain
  //! \param stepOffset offset for starting (time)step
  //! \param timeOffset offset for starting time
  //! \param driverTag tag for current driver section
  //! \param isPartOfSequence true, if driver is part of  multiSequence
  SingleDriver(Domain * adomain, 
	       Integer stepOffset, 
	       Double timeOffset, 
	       std::string driverTag,
	       Boolean isPartOfSequence) :
    BaseDriver(adomain), 
    driverTag_(driverTag),
    stepOffset_(stepOffset), 
    timeOffset_(timeOffset),
    isPartOfSequence_(isPartOfSequence)
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
  
protected:
  
  //! true, if driver is part of  multiSequence
  Boolean isPartOfSequence_;

  //! tag for driver section
  std::string driverTag_;

  //! offset for first timestep
  Integer stepOffset_;

  //! offset for first time
  Double timeOffset_;

private:
   
  };

}

#endif 
