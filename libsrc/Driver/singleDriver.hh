#ifndef FILE_SINGLEDRIVER_2004
#define FILE_SINGLEDRIVER_2004

#include "basedriver.hh"

namespace CoupledField {

  // forward class declaration
  class BasePDE;

  //! Abstract base class for sinlge driver (static, transient, harmonic)
  class SingleDriver : public BaseDriver {
    
  public:
    
    //! Constructor
    //! \param driverTag tag for current driver section
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    SingleDriver( std::string driverTag ="anyTag",
                  bool isPartOfSequence = false );
    
    //! Default destructor
    virtual ~SingleDriver();
    
    //! set the pdes, which gets to be solved
    void SetPDE(BasePDE * pde);
    
    //! main method, where time-stepping is implemented. it is for transient and static problem
    virtual void SolveProblem() = 0;
  
  protected:
  
    //! Trigger reading of restart
    virtual void ReadRestart( ) {};
    
    //! Initialize PDEs
    void InitializePDEs();

    //! pointer to basePDE 
    BasePDE * ptPDE_;

    //! true, if driver is part of a multiSequence
    bool isPartOfSequence_;

    //! tag for driver section
    std::string driverTag_;

  };

}

#endif 
