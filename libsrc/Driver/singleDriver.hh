#ifndef FILE_SINGLEDRIVER_2004
#define FILE_SINGLEDRIVER_2004

#include "basedriver.hh"

#include "PDE/basePDE.hh"
#include "stdSolveStep.hh"

namespace CoupledField {


  //! Abstract base class for sinlge driver (static, transient, harmonic)
  class SingleDriver : public BaseDriver {

  public:
  
  //! Constructor
  //! \param adomain pointer to class Domain
  //! \param stepOffset offset for starting (time)step
  //! \param timeOffset offset for starting time
  //! \param driverTag tag for current driver section
  //! \param isPartOfSequence true, if driver is part of  multiSequence
  SingleDriver(Domain * adomain, 
	       Integer stepOffset = 0, 
	       Double timeOffset = 0.0, 
	       std::string driverTag ="anyTag",
	       Boolean isPartOfSequence = FALSE);

    //! Default destructor
    virtual ~SingleDriver();
  
    //! set the pde, which gets to be solved
    void SetPDE(BasePDE * pde) {
      ptPDE_ = pde;
    }

    //! main method, where time-stepping is implemented. it is for transient and static problem
  virtual void SolveProblem()=0;

    //! to setup matrices of PDE. we call according method of class PDE for setup matrices of PDE in assembling procedure.
    //! \param pdenumber number of PDE
    //! \param matrixtype type of matrix
    virtual void SetupMatricesPDE(Integer pdenumber, const Integer matrixtype){
      Error( "SetupMatricesPDE must be implemented by derived class!",
	     __FILE__, __LINE__ );
    };
  
    //! return the actua frequency within a harmonic analysis
    virtual Double GetActFrequency() {return actFreq_; };

  protected:
  
    //! intialize all PDEs

    //! intialize all PDEs: A more detailed description what and how the PDEs
    //! are initialized would be much appreciated!
    void GetMyPDEs();

    //! pointer to basePDE 
    BasePDE * ptPDE_;

    //! true, if driver is part of a multiSequence
    Boolean isPartOfSequence_;

    //! tag for driver section
    std::string driverTag_;

    //! offset for first timestep
    Integer stepOffset_;

    //! offset for first time
    Double timeOffset_;

    //!
    Double actFreq_;

  };

}

#endif 
