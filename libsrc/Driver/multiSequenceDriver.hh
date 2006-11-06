#ifndef FILE_MULTISEQUENCEDRIVER_2004
#define FILE_MULTISEQUENCEDRIVER_2004

#include "basedriver.hh"
#include "singleDriver.hh"
#include "Utils/StdVector.hh"
#include "PDE/pdememento.hh"

namespace CoupledField
{

  // forward class declaration
  class BasePDE;

  //! Driver class for multi sequence simulations

  //! Driver class for multi sequence simulations.
  //! Multo sequence here means, that either different types
  //! of simulation (static, transient) can be computed one after
  //! another (e.g. prestressed materials with addtitional source
  //! terms) or that mutliple instances of the same anlysis-type can be
  //! computed (e.g. pule-echo simulation in transient case)

  class MultiSequenceDriver : public BaseDriver
  {
  public:
  
    //! constructor
    MultiSequenceDriver();

    //! destructir
    virtual ~MultiSequenceDriver();

    //! Initialization method
    void Init();
  
    //! main method, where time-stepping is implemented. 
    //! it is for transient and static problem
    void SolveProblem();

    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename );

    //! Return current singleDriver 
    SingleDriver* GetSingleDriver() { return actDriver_; }

  private:

    //! number of sequence steps
    UInt numSteps_;

    //! current sequence step
    UInt curSequenceStep_;

    //! current time/frequency step
    UInt actStep_;
  
    //! current time
    Double actTime_;

    //! current singleDriver object
    SingleDriver * actDriver_;
  
    //! stores for each step the participating pdes as name
    StdVector<StdVector<std::string> > pdesPerStep_;

    //! stores for each step the usage of the values for each pde
    StdVector<StdVector<PDEMemento::ValueUsageType> > valueUsagePerStep_;

    //! stores for each step the participating pdes as pointer

    //! stores for each step the tags of each pde
    StdVector<StdVector<std::string> > tagsPerStep_;

    //! stores for each step the analisystype of each pde
    StdVector<StdVector<AnalysisType> > analysisPerStep_;

    //! stores the current set of drivers
    StdVector<SingleDriver *> drivers;
  };

}

#endif 
