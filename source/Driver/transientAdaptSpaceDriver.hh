#ifndef FILE_TRANSIENTADAPTIVESPACEDRIVER_2003
#define FILE_TRANSIENTADAPTIVESPACEDRIVER_2003

#include "basedriver.hh"

namespace CoupledField
{

  //! driver for transient problems.it is derived from BaseDriver;
  class TransientAdaptSpaceDriver : virtual public TransientDriver
  {
  public:
    //!  constructor
    /*!
      \param adomain pointer to class Domain
    */
    TransientAdaptSpaceDriver(Domain * adomain);

    //! deconstructor 
    virtual ~TransientAdaptSpaceDriver();

    //!  main method, where time-stepping is implemented. it is for transient and static problem
    virtual void SolveProblem();

    //!  to setup matrices of PDE. we call according method of class PDE for setup matrices of PDE in assembling procedure.
    /*!
      \param pdenumber (input) number of PDE
      \param matrixtype (input) type of system matrix in algebraic system
    */
    virtual void SetupMatricesPDE(const Integer pdenumber, const Integer matrixtype);

  protected:

  private:
    //!
    Integer numstep_,isavebegin_,isaveincr_,isaveend_;

    //!
    Double firstdt_;

  };

#endif
